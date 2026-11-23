"""Multi-speaker diarization error rate (DER), overlap bounds, TS-VAD fusion, and GSS-MVDR.

Conventions:
- Binary speaker activity matrices have shape ``(K, T)`` where axis 0 indexes
  speakers and axis 1 indexes uniform frames of duration ``frame_dur_s`` seconds.
- Exact integer frame counts are matched before multiplying by ``frame_dur_s`` so
  floating-point rounding never alters optimal speaker permutations.
"""

from __future__ import annotations

import itertools
import numpy as np

from codes.chapters.ch02.core.conventions import (
    finite_real_array,
    finite_real_scalar,
    hermitian_part,
)


def _validate_binary_activity(
    activity: np.ndarray,
    name: str,
    *,
    allow_empty_speakers: bool = False,
) -> np.ndarray:
    """Validate a binary 0/1 speaker-by-frame matrix and return int64 ``(K, T)``."""
    raw = np.asarray(activity)
    if raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain numeric 0/1 values, not complex, bool, or text")
    arr = finite_real_array(raw, name)
    if arr.ndim != 2 or arr.shape[1] < 1:
        raise ValueError(f"{name} must have 2-D shape (K, T) with T >= 1")
    if not allow_empty_speakers and arr.shape[0] < 1:
        raise ValueError(f"{name} must contain at least one speaker row")
    if arr.size > 0 and not np.all((arr == 0.0) | (arr == 1.0)):
        raise ValueError(f"{name} entries must be binary 0 or 1")
    return arr.astype(np.int64)


def _hungarian_max_int(score: list[list[int]]) -> tuple[tuple[int, ...], int]:
    """Exact O(K^3) integer Kuhn-Munkres solver maximizing ``sum_i score[i][perm[i]]``."""
    k_dim = len(score)
    if k_dim == 0:
        return (), 0
    max_entry = max(max(row) for row in score)
    inf_cost = (max_entry + 1) * (k_dim + 2) + 1
    u = [0] * (k_dim + 1)
    v = [0] * (k_dim + 1)
    p = [0] * (k_dim + 1)
    way = [0] * (k_dim + 1)

    for i in range(1, k_dim + 1):
        p[0] = i
        j0 = 0
        minv = [inf_cost] * (k_dim + 1)
        used = [False] * (k_dim + 1)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = inf_cost
            j1 = 0
            for j in range(1, k_dim + 1):
                if not used[j]:
                    cur = (max_entry - score[i0 - 1][j - 1]) - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            for j in range(k_dim + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break

    perm = [0] * k_dim
    for j in range(1, k_dim + 1):
        if p[j] > 0:
            perm[p[j] - 1] = j - 1
    total = sum(score[i][perm[i]] for i in range(k_dim))
    return tuple(perm), total


def _optimal_speaker_assignment(score_matrix: np.ndarray) -> tuple[tuple[int, ...], int]:
    """Find optimal 1-to-1 speaker permutation in exact integer arithmetic."""
    k_dim = int(score_matrix.shape[0])
    score_list = [[int(score_matrix[i, j]) for j in range(k_dim)] for i in range(k_dim)]
    if k_dim <= 7:
        best_total = -1
        best_perm: tuple[int, ...] = tuple(range(k_dim))
        for perm in itertools.permutations(range(k_dim)):
            total = sum(score_list[i][perm[i]] for i in range(k_dim))
            if total > best_total:
                best_total = total
                best_perm = perm
        return best_perm, best_total
    return _hungarian_max_int(score_list)


def compute_der_components(
    ref_activity: np.ndarray,
    hyp_activity: np.ndarray,
    frame_dur_s: float = 0.5,
) -> dict:
    """Compute exact NIST-style framewise DER components under optimal global speaker bijection.

    Note that ``DER = (miss_s + false_alarm_s + speaker_confusion_s) / ref_speaker_time_s``
    can legitimately exceed ``1.0`` (100%) when false alarms plus confusions exceed
    total reference speaker-time; it is never artificially clipped to 1.0.
    """
    ref = _validate_binary_activity(ref_activity, "ref_activity", allow_empty_speakers=False)
    hyp = _validate_binary_activity(hyp_activity, "hyp_activity", allow_empty_speakers=True)
    dt = finite_real_scalar(frame_dur_s, "frame_dur_s")
    if dt <= 0.0:
        raise ValueError("frame_dur_s must be strictly positive")
    if ref.shape[1] != hyp.shape[1]:
        raise ValueError("ref_activity and hyp_activity must have the same frame count T")

    n_frames = ref.shape[1]
    ref_speaker_frames = int(np.sum(ref))
    if ref_speaker_frames <= 0:
        raise ValueError("ref_activity must have positive total active reference speaker-time")

    k_ref, k_hyp = ref.shape[0], hyp.shape[0]
    k_pad = max(k_ref, k_hyp)
    ref_pad = np.zeros((k_pad, n_frames), dtype=np.int64)
    hyp_pad = np.zeros((k_pad, n_frames), dtype=np.int64)
    ref_pad[:k_ref, :] = ref
    if k_hyp > 0:
        hyp_pad[:k_hyp, :] = hyp

    overlap_counts = ref_pad @ hyp_pad.T
    best_perm, total_correct_frames = _optimal_speaker_assignment(overlap_counts)
    hyp_aligned = hyp_pad[list(best_perm), :]

    n_ref_t = np.sum(ref_pad, axis=0)
    n_hyp_t = np.sum(hyp_aligned, axis=0)
    n_cor_t = np.sum((ref_pad == 1) & (hyp_aligned == 1), axis=0)

    miss_frames = int(np.sum(np.maximum(0, n_ref_t - n_hyp_t)))
    fa_frames = int(np.sum(np.maximum(0, n_hyp_t - n_ref_t)))
    conf_frames = int(np.sum(np.minimum(n_ref_t, n_hyp_t) - n_cor_t))
    total_error_frames = miss_frames + fa_frames + conf_frames

    miss_s = float(miss_frames * dt)
    fa_s = float(fa_frames * dt)
    conf_s = float(conf_frames * dt)
    ref_time_s = float(ref_speaker_frames * dt)
    rec_time_s = float(n_frames * dt)
    if not all(np.isfinite(v) for v in (miss_s, fa_s, conf_s, ref_time_s, rec_time_s)):
        raise ValueError("DER duration metrics exceed finite float64 range")

    return {
        "miss_frames": miss_frames,
        "false_alarm_frames": fa_frames,
        "speaker_confusion_frames": conf_frames,
        "correct_matched_frames": int(total_correct_frames),
        "ref_speaker_frames": ref_speaker_frames,
        "miss_s": miss_s,
        "false_alarm_s": fa_s,
        "speaker_confusion_s": conf_s,
        "ref_speaker_time_s": ref_time_s,
        "recording_duration_s": rec_time_s,
        "der": float(total_error_frames / ref_speaker_frames),
        "exceeds_unit_der": bool(total_error_frames > ref_speaker_frames),
        "optimal_permutation": list(best_perm),
    }


def single_label_overlap_der_bound(
    ref_activity: np.ndarray,
    frame_dur_s: float = 0.5,
) -> dict:
    """Compute the exact DER lower bound of any hard single-label clustering system.

    Because a single-label assigner predicts ``sum_k hyp[k, t] <= 1`` on every frame,
    any frame with ``m_t >= 2`` active speakers necessarily incurs at least ``m_t - 1``
    missed speaker-frames even under oracle VAD and zero speaker confusion.
    """
    ref = _validate_binary_activity(ref_activity, "ref_activity", allow_empty_speakers=False)
    dt = finite_real_scalar(frame_dur_s, "frame_dur_s")
    if dt <= 0.0:
        raise ValueError("frame_dur_s must be strictly positive")

    active_per_frame = np.sum(ref, axis=0)
    ref_speaker_frames = int(np.sum(active_per_frame))
    if ref_speaker_frames <= 0:
        raise ValueError("ref_activity must have positive total active reference speaker-time")

    unavoidable_miss_per_frame = np.maximum(0, active_per_frame - 1)
    min_miss_frames = int(np.sum(unavoidable_miss_per_frame))
    overlap_frame_count = int(np.sum(active_per_frame >= 2))

    min_miss_s = float(min_miss_frames * dt)
    ref_time_s = float(ref_speaker_frames * dt)
    if not (np.isfinite(min_miss_s) and np.isfinite(ref_time_s)):
        raise ValueError("DER bound durations exceed finite float64 range")

    return {
        "overlap_frames": overlap_frame_count,
        "min_miss_frames": min_miss_frames,
        "ref_speaker_frames": ref_speaker_frames,
        "min_miss_s": min_miss_s,
        "ref_speaker_time_s": ref_time_s,
        "der_lower_bound": float(min_miss_frames / ref_speaker_frames),
    }


def _stable_sigmoid(logits: np.ndarray) -> np.ndarray:
    """Numerically stable logistic sigmoid avoiding exp overflow for large negative logits."""
    out = np.empty_like(logits, dtype=np.float64)
    pos = logits >= 0.0
    out[pos] = 1.0 / (1.0 + np.exp(-logits[pos]))
    exp_neg = np.exp(logits[~pos])
    out[~pos] = exp_neg / (1.0 + exp_neg)
    return out


def tsvad_spatial_spectral_posterior(
    obs_ipd_cos: np.ndarray | float,
    emb_sim: np.ndarray | float,
    prior: float = 0.5,
    kappa_sp: float = 3.0,
    tau_sp: float = 0.25,
    kappa_emb: float = 4.0,
    tau_emb: float = 0.45,
) -> np.ndarray:
    """Fuse spatial IPD alignment cosine and spectral speaker-embedding cosine in log-odds space.

    Each speaker receives an independent Bernoulli posterior in ``(0, 1)``, allowing
    multiple speakers to be simultaneously active during overlapped speech.
    """
    sp = finite_real_array(obs_ipd_cos, "obs_ipd_cos")
    emb = finite_real_array(emb_sim, "emb_sim")
    if sp.shape != emb.shape:
        raise ValueError("obs_ipd_cos and emb_sim must have identical shapes")
    if np.any(np.abs(sp) > 1.0 + 1e-12) or np.any(np.abs(emb) > 1.0 + 1e-12):
        raise ValueError("cosine similarities obs_ipd_cos and emb_sim must lie within [-1, 1]")
    sp = np.clip(sp, -1.0, 1.0)
    emb = np.clip(emb, -1.0, 1.0)

    pi = finite_real_scalar(prior, "prior")
    k_sp = finite_real_scalar(kappa_sp, "kappa_sp")
    t_sp = finite_real_scalar(tau_sp, "tau_sp")
    k_emb = finite_real_scalar(kappa_emb, "kappa_emb")
    t_emb = finite_real_scalar(tau_emb, "tau_emb")
    if not (0.0 < pi < 1.0):
        raise ValueError("prior must lie strictly inside (0, 1)")
    if k_sp < 0.0 or k_emb < 0.0:
        raise ValueError("kappa_sp and kappa_emb must be non-negative")

    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        prior_logit = float(np.log(pi / (1.0 - pi)))
        logits = np.asarray(prior_logit + k_sp * (sp - t_sp) + k_emb * (emb - t_emb), dtype=np.float64)
    if not (np.isfinite(prior_logit) and np.all(np.isfinite(logits))):
        raise ValueError("TS-VAD log-odds exceed finite float64 range")
    return _stable_sigmoid(logits)


def gss_soft_gated_mvdr(
    observations: np.ndarray,
    posteriors: np.ndarray,
    target_index: int = 0,
    diag_loading: float = 1e-3,
    ref_mic: int = 0,
) -> dict:
    """Estimate target-exclusive and interference SCMs from TS-VAD posteriors and solve Souden MVDR.

    ``target_scm`` is a weighted solo-observation SCM, still containing ambient
    noise and soft-mask leakage. It approximates a rank-one target SCM only in
    the high-SNR control (E15-10); no general distortionless guarantee follows.

    Parameters
    ----------
    observations:
        Narrowband multichannel snapshots ``X`` of shape ``(M, T)``.
    posteriors:
        Multi-speaker activity posteriors ``Gamma`` of shape ``(K, T)`` in ``[0, 1]``.
    target_index:
        Index ``k`` of the target speaker to extract.
    diag_loading:
        Relative trace loading added to the interference-plus-noise covariance.
    ref_mic:
        Reference microphone index for the Souden target-image convention.
    """
    raw_x = np.asarray(observations)
    if raw_x.dtype.kind not in "iufc" or raw_x.ndim != 2 or min(raw_x.shape) < 1:
        raise ValueError("observations must be a non-empty 2-D numeric array (M, T)")
    x_mat = np.asarray(raw_x, dtype=np.complex128)
    if not np.all(np.isfinite(x_mat)):
        raise ValueError("observations must be finite")

    gamma = finite_real_array(posteriors, "posteriors")
    if gamma.ndim != 2 or gamma.shape[1] != x_mat.shape[1] or gamma.shape[0] < 1:
        raise ValueError("posteriors must have shape (K, T) matching observations frame count T")
    if np.any(gamma < 0.0) or np.any(gamma > 1.0):
        raise ValueError("posteriors entries must lie in [0, 1]")

    if isinstance(target_index, (bool, np.bool_)) or not isinstance(target_index, (int, np.integer)):
        raise ValueError("target_index must be an integer")
    if isinstance(ref_mic, (bool, np.bool_)) or not isinstance(ref_mic, (int, np.integer)):
        raise ValueError("ref_mic must be an integer")
    n_mics = x_mat.shape[0]
    n_spk = gamma.shape[0]
    if not (0 <= target_index < n_spk):
        raise ValueError("target_index out of range")
    if not (0 <= ref_mic < n_mics):
        raise ValueError("ref_mic out of range")
    load = finite_real_scalar(diag_loading, "diag_loading")
    if load < 0.0:
        raise ValueError("diag_loading must be non-negative")

    other_inactive = np.ones(gamma.shape[1], dtype=np.float64)
    for spk in range(n_spk):
        if spk != target_index:
            other_inactive *= (1.0 - gamma[spk, :])
    target_mask = gamma[target_index, :] * other_inactive
    interf_mask = 1.0 - gamma[target_index, :]

    tar_denom = float(np.sum(target_mask))
    int_denom = float(np.sum(interf_mask))
    if tar_denom <= 1e-12 or int_denom <= 1e-12:
        raise ValueError("insufficient target-exclusive or target-inactive posterior support")

    x_peak = float(np.max(np.abs(x_mat)))
    if x_peak <= 0.0:
        raise np.linalg.LinAlgError("observations are identically zero")
    x_scaled = np.empty_like(x_mat)
    x_scaled.real = x_mat.real / x_peak
    x_scaled.imag = x_mat.imag / x_peak

    r_ss_scaled = hermitian_part((x_scaled * target_mask[None, :]) @ x_scaled.conj().T / tar_denom)
    r_vv_raw_scaled = hermitian_part((x_scaled * interf_mask[None, :]) @ x_scaled.conj().T / int_denom)
    avg_diag_scaled = float(np.trace(r_vv_raw_scaled).real / n_mics)
    if avg_diag_scaled <= 0.0 or not np.isfinite(avg_diag_scaled):
        raise np.linalg.LinAlgError("interference covariance has non-positive trace")
    r_vv_scaled = hermitian_part(r_vv_raw_scaled + max(load, 0.0) * avg_diag_scaled * np.eye(n_mics))
    if np.linalg.eigvalsh(r_vv_scaled)[0] <= 0.0:
        raise np.linalg.LinAlgError("loaded interference covariance is not strictly positive definite")

    rinv_rss = np.linalg.solve(r_vv_scaled, r_ss_scaled)
    denom_trace = float(np.trace(rinv_rss).real)
    if denom_trace <= 1e-15 or not np.isfinite(denom_trace):
        raise np.linalg.LinAlgError("degenerate target spatial covariance trace in Souden MVDR")
    weights = rinv_rss[:, ref_mic] / denom_trace
    if not np.all(np.isfinite(weights)):
        raise np.linalg.LinAlgError("Souden MVDR weights contain non-finite values")

    with np.errstate(over="ignore", invalid="ignore"):
        enhanced = weights.conj() @ x_mat
        scale_sq = x_peak * x_peak
    if not np.all(np.isfinite(enhanced)):
        raise ValueError("enhanced_frames exceed finite float64 range")
    if not np.isfinite(scale_sq):
        raise ValueError("spatial covariance matrix exceeds finite float64 range")
    r_ss = r_ss_scaled * scale_sq
    r_vv = r_vv_scaled * scale_sq

    return {
        "target_exclusive_mask": target_mask,
        "interference_mask": interf_mask,
        "target_scm": r_ss,
        "interference_scm": r_vv,
        "weights": weights,
        "enhanced_frames": enhanced,
    }
