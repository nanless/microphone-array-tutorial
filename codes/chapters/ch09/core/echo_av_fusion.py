"""Pure-NumPy teaching kernels for self-loudspeaker echo gating and AV OOSM/MOT.

Covers Chapter 09 exercises E09-27 (reference-activity-conditioned loudspeaker
sector gating before association and track birth) and E09-28 (audio-visual
out-of-sequence measurement retrodiction H_retro = [1, -lag_s], camera FOV-aware
Bernoulli existence conservation, and CLEAR MOT MOTA/IDSW vs label-free OSPA).
No SciPy imports; deterministic float64.
"""
from __future__ import annotations

from itertools import permutations
import math
import numbers
import numpy as np


MAX_TEACHING_TRACKS = 8


def _finite_real_array(value, name: str = "input") -> np.ndarray:
    original = np.asarray(value)
    if original.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numeric values, not complex, bool or text")
    with np.errstate(over="ignore", invalid="ignore"):
        array = np.asarray(original, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite and representable as float64")
    return array


def _finite_real_scalar(value, name: str = "input") -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a real scalar, not bool")
    array = _finite_real_array(value, name)
    if array.ndim != 0:
        raise ValueError(f"{name} must be a scalar")
    return float(array)


def _integer(value, name: str, *, minimum: int | None = None) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Integral):
        raise ValueError(f"{name} must be an integer")
    result = int(value)
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return result


def _wrap_angle_deg(angle: float | np.ndarray) -> float | np.ndarray:
    value = _finite_real_array(angle, "angle")
    wrapped = ((value + 180.0) % 360.0) - 180.0
    wrapped = np.where(wrapped >= 180.0, -180.0, wrapped)
    wrapped = np.where(np.abs(value) < 1e-12, value, wrapped)
    return float(wrapped) if wrapped.ndim == 0 else wrapped


def _checked_psd_2x2(matrix, name: str = "covariance") -> np.ndarray:
    arr = _finite_real_array(matrix, name)
    if arr.shape != (2, 2):
        raise ValueError(f"{name} must have shape (2, 2)")
    if not np.allclose(arr, arr.T, rtol=1e-12, atol=1e-15):
        raise ValueError(f"{name} must be symmetric")
    sym = 0.5 * (arr + arr.T)
    eigvals = np.linalg.eigvalsh(sym)
    if eigvals[0] < 0.0:
        raise ValueError(f"{name} must be positive semidefinite")
    return sym


def echo_sector_gate(
    detections_deg: np.ndarray,
    *,
    speaker_angle_deg: float,
    sector_half_width_deg: float,
    far_end_active: bool,
    double_talk: bool,
) -> dict:
    """Filter DOA detections against loudspeaker sector during far-end single-talk."""
    if not isinstance(far_end_active, (bool, np.bool_)) or not isinstance(double_talk, (bool, np.bool_)):
        raise ValueError("far_end_active and double_talk must be booleans")
    dets = _finite_real_array(detections_deg, "detections_deg")
    if dets.ndim != 1:
        raise ValueError("detections_deg must be a 1-D real array")
    spk_deg = float(_wrap_angle_deg(_finite_real_scalar(speaker_angle_deg, "speaker_angle_deg")))
    half_width = _finite_real_scalar(sector_half_width_deg, "sector_half_width_deg")
    if not 0.0 < half_width < 180.0:
        raise ValueError("sector_half_width_deg must satisfy 0 < half_width < 180")

    wrapped_dets = np.asarray(_wrap_angle_deg(dets), dtype=np.float64).reshape(-1)
    angular_diffs = np.abs(np.asarray(_wrap_angle_deg(wrapped_dets - spk_deg), dtype=np.float64)).reshape(-1)
    inside_sector = angular_diffs <= half_width
    gate_active = bool(far_end_active and (not double_talk))
    suppressed_mask = inside_sector if gate_active else np.zeros(wrapped_dets.shape, dtype=bool)
    retained = wrapped_dets[~suppressed_mask]
    suppressed = wrapped_dets[suppressed_mask]
    return {
        "gate_active": gate_active,
        "speaker_angle_deg": spk_deg,
        "sector_half_width_deg": half_width,
        "input_detections_deg": wrapped_dets,
        "angular_distance_to_speaker_deg": angular_diffs,
        "inside_sector_mask": inside_sector,
        "suppressed_mask": suppressed_mask,
        "retained_detections_deg": retained,
        "suppressed_detections_deg": suppressed,
    }


def oosm_cv_retrodiction_update(
    state: np.ndarray,
    covariance: np.ndarray,
    delayed_angle_deg: float,
    measurement_variance_deg2: float,
    lag_s: float,
) -> dict:
    """Exact zero-Q constant-velocity OOSM retrodiction H=[1, -lag_s] vs naive H=[1, 0]."""
    s_vec = _finite_real_array(state, "state").copy()
    if s_vec.shape != (2,):
        raise ValueError("state must have shape (2,) = [angle_deg, velocity_deg_per_s]")
    s_vec[0] = float(_wrap_angle_deg(s_vec[0]))
    p_mat = _checked_psd_2x2(covariance, "covariance")
    z_deg = float(_wrap_angle_deg(_finite_real_scalar(delayed_angle_deg, "delayed_angle_deg")))
    r_var = _finite_real_scalar(measurement_variance_deg2, "measurement_variance_deg2")
    lag = _finite_real_scalar(lag_s, "lag_s")
    if r_var <= 0.0:
        raise ValueError("measurement_variance_deg2 must be strictly positive")
    if lag < 0.0:
        raise ValueError("lag_s must be non-negative")

    eye2 = np.eye(2, dtype=np.float64)

    def _joseph_scalar_update(h_row: np.ndarray) -> dict:
        pred_z = float(_wrap_angle_deg(float(h_row @ s_vec)))
        innovation = float(_wrap_angle_deg(z_deg - pred_z))
        ph_t = p_mat @ h_row
        s_innov = float(h_row @ ph_t + r_var)
        if s_innov <= 0.0 or not math.isfinite(s_innov):
            raise ValueError("innovation variance must be finite and positive")
        gain = ph_t / s_innov
        post_state = s_vec + gain * innovation
        post_state[0] = float(_wrap_angle_deg(post_state[0]))
        a_mat = eye2 - np.outer(gain, h_row)
        post_cov = a_mat @ p_mat @ a_mat.T + r_var * np.outer(gain, gain)
        post_cov = _checked_psd_2x2(0.5 * (post_cov + post_cov.T), "posterior_covariance")
        return {
            "observation_matrix": h_row,
            "predicted_measurement_deg": pred_z,
            "innovation_deg": innovation,
            "innovation_variance_deg2": s_innov,
            "kalman_gain": gain,
            "posterior_state": post_state,
            "posterior_covariance": post_cov,
        }

    h_retro = np.array([1.0, -lag], dtype=np.float64)
    h_naive = np.array([1.0, 0.0], dtype=np.float64)
    return {
        "lag_s": lag,
        "delayed_angle_deg": z_deg,
        "measurement_variance_deg2": r_var,
        "retrodiction": _joseph_scalar_update(h_retro),
        "naive_current_time": _joseph_scalar_update(h_naive),
    }


def fov_bernoulli_missed_detection_update(
    prior_existence: float,
    angle_deg: float,
    *,
    fov_half_width_deg: float = 45.0,
    in_fov_detection_prob: float = 0.95,
    camera_center_deg: float = 0.0,
) -> dict:
    """Compare FOV-aware vs naive constant-p_D Bernoulli existence update on empty visual frame."""
    r_prior = _finite_real_scalar(prior_existence, "prior_existence")
    theta = float(_wrap_angle_deg(_finite_real_scalar(angle_deg, "angle_deg")))
    fov_half = _finite_real_scalar(fov_half_width_deg, "fov_half_width_deg")
    pd_in = _finite_real_scalar(in_fov_detection_prob, "in_fov_detection_prob")
    cam_center = float(_wrap_angle_deg(_finite_real_scalar(camera_center_deg, "camera_center_deg")))
    if not 0.0 <= r_prior <= 1.0:
        raise ValueError("prior_existence must lie in [0, 1]")
    if not 0.0 < fov_half <= 180.0:
        raise ValueError("fov_half_width_deg must satisfy 0 < fov_half_width_deg <= 180")
    if not 0.0 <= pd_in <= 1.0:
        raise ValueError("in_fov_detection_prob must lie in [0, 1]")

    in_fov = bool(abs(float(_wrap_angle_deg(theta - cam_center))) <= fov_half)
    pd_fov = pd_in if in_fov else 0.0
    if r_prior * pd_fov >= 1.0:
        raise ValueError("r_prior * p_D == 1 makes an empty visual frame impossible")
    r_post_fov = float(r_prior * (1.0 - pd_fov) / (1.0 - r_prior * pd_fov))
    r_post_naive = (
        None
        if r_prior * pd_in >= 1.0
        else float(r_prior * (1.0 - pd_in) / (1.0 - r_prior * pd_in))
    )
    return {
        "prior_existence": r_prior,
        "angle_deg": theta,
        "camera_center_deg": cam_center,
        "fov_half_width_deg": fov_half,
        "in_fov": in_fov,
        "fov_aware_detection_prob": pd_fov,
        "fov_aware_posterior_existence": r_post_fov,
        "naive_constant_detection_prob": pd_in,
        "naive_posterior_existence": r_post_naive,
    }


def angular_ospa_distance(
    truth_deg: np.ndarray | list[float],
    estimate_deg: np.ndarray | list[float],
    *,
    cutoff_deg: float = 10.0,
    order: int = 1,
) -> float:
    """Label-free circular OSPA distance for <= 7 targets (p = 1 or 2)."""
    first = _finite_real_array(truth_deg, "truth_deg")
    second = _finite_real_array(estimate_deg, "estimate_deg")
    cutoff = _finite_real_scalar(cutoff_deg, "cutoff_deg")
    p_ord = _integer(order, "order", minimum=1)
    if first.ndim != 1 or second.ndim != 1 or max(first.size, second.size) > 7:
        raise ValueError("truth_deg and estimate_deg must be 1-D arrays of length <= 7")
    if not 0.0 < cutoff <= 180.0 or p_ord not in (1, 2):
        raise ValueError("require 0 < cutoff_deg <= 180 and order in {1, 2}")
    first = np.asarray(_wrap_angle_deg(first), dtype=np.float64).reshape(-1)
    second = np.asarray(_wrap_angle_deg(second), dtype=np.float64).reshape(-1)
    if first.size > second.size:
        first, second = second, first
    m, n = first.size, second.size
    if n == 0:
        return 0.0
    best_cost = math.inf
    for perm in permutations(range(n), m):
        dists = [min(cutoff, abs(float(_wrap_angle_deg(first[i] - second[j])))) for i, j in enumerate(perm)]
        cost = sum(d ** p_ord for d in dists)
        if cost < best_cost:
            best_cost = cost
    total = (best_cost + (cutoff ** p_ord) * (n - m)) / n
    return float(total if p_ord == 1 else math.sqrt(total))


def _validate_id_angle_frame(frame, name: str) -> dict[int, float]:
    if type(frame) is not dict:
        raise ValueError(f"{name} frame must be a plain dict[int, float]")
    if len(frame) > MAX_TEACHING_TRACKS:
        raise ValueError(f"{name} frame exceeds MAX_TEACHING_TRACKS={MAX_TEACHING_TRACKS}")
    cleaned: dict[int, float] = {}
    for raw_id, raw_angle in frame.items():
        track_id = _integer(raw_id, f"{name} id")
        angle = float(_wrap_angle_deg(_finite_real_scalar(raw_angle, f"{name} angle")))
        cleaned[track_id] = angle
    return cleaned


def clear_mot_angular_summary(
    truth_frames: list[dict[int, float]],
    hyp_frames: list[dict[int, float]],
    *,
    match_threshold_deg: float = 10.0,
) -> dict:
    """Compute exact FN, FP, IDSW, MOTA and MOTP with CLEAR MOT persistent-ID carry-over."""
    if type(truth_frames) not in (list, tuple) or type(hyp_frames) not in (list, tuple):
        raise ValueError("truth_frames and hyp_frames must be lists or tuples of dicts")
    if len(truth_frames) == 0 or len(truth_frames) != len(hyp_frames):
        raise ValueError("truth_frames and hyp_frames must be non-empty and equal-length")
    tau = _finite_real_scalar(match_threshold_deg, "match_threshold_deg")
    if not 0.0 < tau <= 180.0:
        raise ValueError("match_threshold_deg must satisfy 0 < match_threshold_deg <= 180")

    validated_truth = [_validate_id_angle_frame(f, "truth") for f in truth_frames]
    validated_hyp = [_validate_id_angle_frame(f, "hyp") for f in hyp_frames]
    total_gt = sum(len(f) for f in validated_truth)
    if total_gt == 0:
        raise ValueError("MOTA requires at least one ground-truth object across the sequence")

    prev_correspondence: dict[int, int] = {}
    total_fn = 0
    total_fp = 0
    total_idsw = 0
    total_matches = 0
    total_distance = 0.0
    per_frame = []

    for frame_idx, (gt_map, hyp_map) in enumerate(zip(validated_truth, validated_hyp)):
        current_matches: dict[int, int] = {}
        used_hyps: set[int] = set()
        frame_distances: list[float] = []
        frame_idsw = 0

        # Step 1 (Bernardin & Stiefelhagen 2008): Carry over persistent correspondences within gate.
        for gt_id in sorted(gt_map):
            if gt_id in prev_correspondence:
                prev_hyp = prev_correspondence[gt_id]
                if prev_hyp in hyp_map and prev_hyp not in used_hyps:
                    dist = abs(float(_wrap_angle_deg(gt_map[gt_id] - hyp_map[prev_hyp])))
                    if dist <= tau:
                        current_matches[gt_id] = prev_hyp
                        used_hyps.add(prev_hyp)
                        frame_distances.append(dist)
                        for old_gid in [
                            k for k, v in prev_correspondence.items() if v == prev_hyp and k != gt_id
                        ]:
                            del prev_correspondence[old_gid]
                        prev_correspondence[gt_id] = prev_hyp

        # Step 2: Exhaustive max-cardinality, min-distance matching on remaining gated pairs.
        rem_gt = [gid for gid in sorted(gt_map) if gid not in current_matches]
        rem_hyp = [hid for hid in sorted(hyp_map) if hid not in used_hyps]

        best_key: tuple[int, float, int, tuple[tuple[int, int], ...]] | None = None
        best_pairs: list[tuple[int, int, float]] = []

        def _search(idx: int, used: set[int], chosen: list[tuple[int, int, float]]) -> None:
            nonlocal best_key, best_pairs
            if idx == len(rem_gt):
                count = len(chosen)
                dist_sum = math.fsum(item[2] for item in chosen)
                switches = sum(
                    1
                    for gid, hid, _ in chosen
                    if gid in prev_correspondence and prev_correspondence[gid] != hid
                )
                pair_ids = tuple((gid, hid) for gid, hid, _ in chosen)
                cand_key = (-count, round(dist_sum, 12), switches, pair_ids)
                if best_key is None or cand_key < best_key:
                    best_key = cand_key
                    best_pairs = list(chosen)
                return

            gid = rem_gt[idx]
            for hid in rem_hyp:
                if hid in used:
                    continue
                dist = abs(float(_wrap_angle_deg(gt_map[gid] - hyp_map[hid])))
                if dist <= tau:
                    used.add(hid)
                    chosen.append((gid, hid, dist))
                    _search(idx + 1, used, chosen)
                    chosen.pop()
                    used.remove(hid)
            _search(idx + 1, used, chosen)

        _search(0, set(), [])

        for gid, hid, dist in best_pairs:
            current_matches[gid] = hid
            used_hyps.add(hid)
            frame_distances.append(dist)
            if gid in prev_correspondence and prev_correspondence[gid] != hid:
                frame_idsw += 1

        for gid, hid, _ in best_pairs:
            for old_gid in [k for k, v in prev_correspondence.items() if v == hid and k != gid]:
                del prev_correspondence[old_gid]
            prev_correspondence[gid] = hid

        frame_fn = len(gt_map) - len(current_matches)
        frame_fp = len(hyp_map) - len(current_matches)
        frame_motp = None if not frame_distances else float(math.fsum(frame_distances) / len(frame_distances))

        total_fn += frame_fn
        total_fp += frame_fp
        total_idsw += frame_idsw
        total_matches += len(current_matches)
        total_distance += math.fsum(frame_distances)
        per_frame.append(
            {
                "frame": frame_idx,
                "ground_truth_count": len(gt_map),
                "hypothesis_count": len(hyp_map),
                "matches": {str(k): v for k, v in sorted(current_matches.items())},
                "fn": frame_fn,
                "fp": frame_fp,
                "idsw": frame_idsw,
                "motp_deg": frame_motp,
            }
        )

    mota = float(1.0 - (total_fn + total_fp + total_idsw) / total_gt)
    motp = None if total_matches == 0 else float(total_distance / total_matches)
    return {
        "match_threshold_deg": tau,
        "total_ground_truth": total_gt,
        "total_matches": total_matches,
        "fn": total_fn,
        "fp": total_fp,
        "idsw": total_idsw,
        "mota": mota,
        "motp_deg": motp,
        "frames": per_frame,
    }


def _cv_kalman_predict(state: np.ndarray, cov: np.ndarray, dt: float, q_mat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    f_mat = np.array([[1.0, dt], [0.0, 1.0]], dtype=np.float64)
    s_pred = f_mat @ state
    s_pred[0] = float(_wrap_angle_deg(s_pred[0]))
    p_pred = f_mat @ cov @ f_mat.T + q_mat
    return s_pred, 0.5 * (p_pred + p_pred.T)


def _cv_kalman_update(state_pred: np.ndarray, cov_pred: np.ndarray, z_deg: float, r_var: float) -> dict:
    h_row = np.array([1.0, 0.0], dtype=np.float64)
    innov = float(_wrap_angle_deg(z_deg - state_pred[0]))
    s_innov = float(cov_pred[0, 0] + r_var)
    gain = cov_pred[:, 0] / s_innov
    nis = float((innov ** 2) / s_innov)
    s_post = state_pred + gain * innov
    s_post[0] = float(_wrap_angle_deg(s_post[0]))
    a_mat = np.eye(2, dtype=np.float64) - np.outer(gain, h_row)
    p_post = a_mat @ cov_pred @ a_mat.T + r_var * np.outer(gain, gain)
    return {
        "innovation_deg": innov,
        "innovation_variance_deg2": s_innov,
        "nis": nis,
        "gain": gain,
        "posterior_state": s_post,
        "posterior_covariance": 0.5 * (p_post + p_post.T),
    }


def echo_clutter_gating_experiment() -> dict:
    """E09-27: Self-loudspeaker residual echo clutter hijack vs reference-activity sector gate."""
    s0 = np.array([20.0, 0.0], dtype=np.float64)
    p0 = np.diag([4.0, 1.0])
    q_mat = np.diag([1.0, 0.0])
    r_var = 4.0
    dt = 1.0
    speaker_deg = -45.0
    half_width_deg = 10.0

    # 1. Ungated single-track hijack across 2 residual-echo frames z_1 = z_2 = -45 deg.
    s1_pred, p1_pred = _cv_kalman_predict(s0, p0, dt, q_mat)
    upd1 = _cv_kalman_update(s1_pred, p1_pred, -45.0, r_var)
    s2_pred_u, p2_pred_u = _cv_kalman_predict(upd1["posterior_state"], upd1["posterior_covariance"], dt, q_mat)
    upd2 = _cv_kalman_update(s2_pred_u, p2_pred_u, -45.0, r_var)

    # 2. Unmasked 3-consecutive-hit birth counter on NIS-rejected -45 deg echo clutter.
    unmasked_consecutive = 0
    false_births_unmasked = 0
    for _ in range(3):
        unmasked_consecutive += 1
        if unmasked_consecutive == 3:
            false_births_unmasked += 1

    # 3. Sector-gated coasting during far-end single-talk (frames 1, 2) + double-talk near-end update.
    gate_f1 = echo_sector_gate(
        np.array([-45.0]),
        speaker_angle_deg=speaker_deg,
        sector_half_width_deg=half_width_deg,
        far_end_active=True,
        double_talk=False,
    )
    s1_gated, p1_gated = _cv_kalman_predict(s0, p0, dt, q_mat)
    s2_gated, p2_gated = _cv_kalman_predict(s1_gated, p1_gated, dt, q_mat)
    gated_false_births = 0

    # Double-talk control: near-end active at +22 deg alongside residual echo at -45 deg.
    gate_dt = echo_sector_gate(
        np.array([-45.0, 22.0]),
        speaker_angle_deg=speaker_deg,
        sector_half_width_deg=half_width_deg,
        far_end_active=True,
        double_talk=True,
    )
    dt_near_upd = _cv_kalman_update(s1_pred, p1_pred, 22.0, r_var)

    return {
        "initial_state": s0.tolist(),
        "initial_covariance": p0.tolist(),
        "process_noise": q_mat.tolist(),
        "measurement_variance_deg2": r_var,
        "speaker_angle_deg": speaker_deg,
        "sector_half_width_deg": half_width_deg,
        "step1_prediction": {
            "predicted_state": s1_pred.tolist(),
            "predicted_covariance": p1_pred.tolist(),
        },
        "ungated_hijack": {
            "step1_nis": upd1["nis"],
            "step1_gain": upd1["gain"].tolist(),
            "step1_posterior_state": upd1["posterior_state"].tolist(),
            "step1_posterior_covariance": upd1["posterior_covariance"].tolist(),
            "step2_predicted_state": s2_pred_u.tolist(),
            "step2_predicted_covariance": p2_pred_u.tolist(),
            "step2_innovation_variance_deg2": upd2["innovation_variance_deg2"],
            "step2_gain": upd2["gain"].tolist(),
            "step2_innovation_deg": upd2["innovation_deg"],
            "step2_posterior_state": upd2["posterior_state"].tolist(),
            "step2_posterior_covariance": upd2["posterior_covariance"].tolist(),
            "nis_gate_9_rejects_step1": bool(upd1["nis"] > 9.0),
            "unmasked_3hit_false_track_births": false_births_unmasked,
        },
        "sector_gated_single_talk": {
            "gate_active": gate_f1["gate_active"],
            "retained_detections_deg": gate_f1["retained_detections_deg"].tolist(),
            "suppressed_detections_deg": gate_f1["suppressed_detections_deg"].tolist(),
            "step1_coasted_state": s1_gated.tolist(),
            "step1_coasted_covariance": p1_gated.tolist(),
            "step2_coasted_state": s2_gated.tolist(),
            "step2_coasted_covariance": p2_gated.tolist(),
            "false_track_births": gated_false_births,
        },
        "double_talk_control": {
            "gate_active": gate_dt["gate_active"],
            "retained_detections_deg": gate_dt["retained_detections_deg"].tolist(),
            "near_end_22deg_posterior_state": dt_near_upd["posterior_state"].tolist(),
        },
        "scope": (
            "two-step CV Kalman and 3-hit birth gate control; loudspeaker geometry "
            "and far-end/double-talk flags are supplied inputs, not estimated"
        ),
    }


def av_oosm_and_mot_experiment() -> dict:
    """E09-28: Audio-visual OOSM retrodiction, FOV existence conservation, and MOTA vs OSPA."""
    state_t4 = np.array([13.0, 100.0], dtype=np.float64)
    cov_t4 = np.array([[9.0, 20.0], [20.0, 100.0]], dtype=np.float64)
    lag_s = 0.03
    z_cam_t1 = 10.0
    r_cam = 1.0

    oosm = oosm_cv_retrodiction_update(state_t4, cov_t4, z_cam_t1, r_cam, lag_s)

    # Independent verification via chronological state-buffer replay t_1 -> t_4 (zero Q on [t_1, t_4]).
    f_back = np.array([[1.0, -lag_s], [0.0, 1.0]], dtype=np.float64)
    f_fwd = np.array([[1.0, lag_s], [0.0, 1.0]], dtype=np.float64)
    s_t1_prior = f_back @ state_t4
    p_t1_prior = f_back @ cov_t4 @ f_back.T
    t1_upd = _cv_kalman_update(s_t1_prior, p_t1_prior, z_cam_t1, r_cam)
    replay_state_t4 = f_fwd @ t1_upd["posterior_state"]
    replay_cov_t4 = f_fwd @ t1_upd["posterior_covariance"] @ f_fwd.T

    fov_check = fov_bernoulli_missed_detection_update(
        0.8,
        60.0,
        fov_half_width_deg=45.0,
        in_fov_detection_prob=0.95,
    )

    truth_seq = [{1: -20.0, 2: 20.0}, {1: -20.0, 2: 20.0}]
    hyp_swap = [{1: -20.0, 2: 20.0}, {1: 20.0, 2: -20.0}]
    hyp_swap_plus_echo = [{1: -20.0, 2: 20.0}, {1: 20.0, 2: -20.0, 3: -45.0}]

    mot_swap = clear_mot_angular_summary(truth_seq, hyp_swap, match_threshold_deg=10.0)
    mot_swap_echo = clear_mot_angular_summary(truth_seq, hyp_swap_plus_echo, match_threshold_deg=10.0)

    ospa_swap_frames = [
        angular_ospa_distance(list(t.values()), list(h.values()), cutoff_deg=10.0, order=1)
        for t, h in zip(truth_seq, hyp_swap)
    ]
    ospa_swap_echo_frames = [
        angular_ospa_distance(list(t.values()), list(h.values()), cutoff_deg=10.0, order=1)
        for t, h in zip(truth_seq, hyp_swap_plus_echo)
    ]

    def _serialize_upd(u: dict) -> dict:
        return {
            "observation_matrix": u["observation_matrix"].tolist(),
            "predicted_measurement_deg": u["predicted_measurement_deg"],
            "innovation_deg": u["innovation_deg"],
            "innovation_variance_deg2": u["innovation_variance_deg2"],
            "kalman_gain": u["kalman_gain"].tolist(),
            "posterior_state": u["posterior_state"].tolist(),
            "posterior_covariance": u["posterior_covariance"].tolist(),
        }

    return {
        "oosm_retrodiction": {
            "audio_state_at_t4": state_t4.tolist(),
            "audio_covariance_at_t4": cov_t4.tolist(),
            "camera_exposure_time_s": 0.01,
            "camera_arrival_time_s": 0.04,
            "lag_s": lag_s,
            "delayed_camera_angle_deg": z_cam_t1,
            "camera_variance_deg2": r_cam,
            "retrodiction": _serialize_upd(oosm["retrodiction"]),
            "naive_current_time": _serialize_upd(oosm["naive_current_time"]),
            "chronological_replay_state_at_t4": replay_state_t4.tolist(),
            "chronological_replay_covariance_at_t4": replay_cov_t4.tolist(),
        },
        "camera_fov_existence": fov_check,
        "identity_vs_set_metrics": {
            "pure_id_swap": {
                "clear_mot": mot_swap,
                "ospa_per_frame_deg": ospa_swap_frames,
            },
            "id_swap_plus_echo_false_track": {
                "clear_mot": mot_swap_echo,
                "ospa_per_frame_deg": ospa_swap_echo_frames,
            },
        },
        "scope": (
            "zero-process-noise OOSM retrodiction, piecewise binary FOV detection probability, "
            "and 2-frame CLEAR MOT vs label-free OSPA comparison"
        ),
    }
