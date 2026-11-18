"""EEND-EDA variable-speaker attractor inference, fixed-slot DER bounds, and enrollment robustness.

Conventions:
- Frame embeddings ``E`` have shape ``(D, T)`` and candidate attractors ``A`` have
  shape ``(D, S_max)``. Sequential EEND-EDA retains attractors ``s = 0 .. S_hat - 1``
  up to the first index whose existence probability ``p_s = sigma(w^T a_s + b)`` drops
  below ``exist_threshold`` (default ``0.5``).
- Clean speaker enrollment vectors are unit-L2-norm vectors ``||e_k||_2 = 1``.
  Contamination by an interferer at linear leakage ratio ``alpha >= 0`` yields
  ``e_cont = (e_1 + alpha * e_2) / ||e_1 + alpha * e_2||_2``.
"""
from __future__ import annotations

import numpy as np

from codes.chapters.ch02.core.conventions import (
    finite_real_array,
    finite_real_scalar,
)
from codes.chapters.ch17.core.diarization_tsvad import (
    _stable_sigmoid,
    _validate_binary_activity,
)


def eend_eda_inference(
    frame_embeddings: np.ndarray,
    candidate_attractors: np.ndarray,
    w_exist: np.ndarray,
    b_exist: float = -1.0,
    *,
    exist_threshold: float = 0.5,
    act_threshold: float = 0.5,
) -> dict:
    """Run sequential EEND-EDA attractor existence stopping and multi-label diarization.

    Parameters
    ----------
    frame_embeddings : array_like of shape ``(D, T)``
        D-dimensional encoder frame representations across ``T >= 1`` frames.
    candidate_attractors : array_like of shape ``(D, S_max)``
        Sequentially decoded candidate speaker attractors ``[a_1, ..., a_{S_max}]``.
    w_exist : array_like of shape ``(D,)``
        Linear existence-head weight vector ``w_exist``.
    b_exist : float, default -1.0
        Existence-head scalar bias ``b``.
    exist_threshold : float, default 0.5
        Stopping threshold in ``(0, 1)``; retention halts at the first ``p_s < exist_threshold``.
    act_threshold : float, default 0.5
        Framewise posterior threshold in ``(0, 1)`` for binary speaker activity decisions.
    """
    emb = finite_real_array(frame_embeddings, "frame_embeddings")
    attr = finite_real_array(candidate_attractors, "candidate_attractors")
    w_vec = finite_real_array(w_exist, "w_exist")
    bias = finite_real_scalar(b_exist, "b_exist")
    tau_exist = finite_real_scalar(exist_threshold, "exist_threshold")
    tau_act = finite_real_scalar(act_threshold, "act_threshold")

    if emb.ndim != 2 or emb.shape[0] < 1 or emb.shape[1] < 1:
        raise ValueError("frame_embeddings must have shape (D, T) with D >= 1 and T >= 1")
    d_dim, n_frames = emb.shape
    if attr.ndim != 2 or attr.shape[0] != d_dim or attr.shape[1] < 1:
        raise ValueError("candidate_attractors must have shape (D, S_max) matching frame_embeddings D")
    if w_vec.ndim != 1 or w_vec.shape[0] != d_dim:
        raise ValueError("w_exist must have shape (D,) matching frame_embeddings D")
    if not (0.0 < tau_exist < 1.0):
        raise ValueError("exist_threshold must lie strictly inside (0, 1)")
    if not (0.0 < tau_act < 1.0):
        raise ValueError("act_threshold must lie strictly inside (0, 1)")

    with np.errstate(over="raise", invalid="raise"):
        try:
            existence_logits = w_vec @ attr + bias
        except FloatingPointError as exc:
            raise ValueError("EEND-EDA existence logits overflowed float64") from exc
    if not np.all(np.isfinite(existence_logits)):
        raise ValueError("EEND-EDA existence logits must be finite")

    existence_probs = _stable_sigmoid(existence_logits)
    stop_indices = np.flatnonzero(existence_probs < tau_exist)
    s_hat = int(stop_indices[0]) if stop_indices.size > 0 else int(attr.shape[1])

    active_attractors = attr[:, :s_hat].copy()
    if s_hat > 0:
        with np.errstate(over="raise", invalid="raise"):
            try:
                activity_logits = active_attractors.T @ emb
            except FloatingPointError as exc:
                raise ValueError("EEND-EDA activity logits overflowed float64") from exc
        if not np.all(np.isfinite(activity_logits)):
            raise ValueError("EEND-EDA activity logits must be finite")
        posteriors = _stable_sigmoid(activity_logits)
        binary_activity = (posteriors >= tau_act).astype(np.int64)
    else:
        activity_logits = np.zeros((0, n_frames), dtype=np.float64)
        posteriors = np.zeros((0, n_frames), dtype=np.float64)
        binary_activity = np.zeros((0, n_frames), dtype=np.int64)

    return {
        "existence_logits": existence_logits,
        "existence_probs": existence_probs,
        "estimated_num_speakers": s_hat,
        "active_attractors": active_attractors,
        "activity_logits": activity_logits,
        "posteriors": posteriors,
        "binary_activity": binary_activity,
    }


def fixed_slot_der_lower_bound(
    ref_activity: np.ndarray,
    max_slots: int,
    frame_dur_s: float = 0.5,
) -> dict:
    """Compute exact DER lower bound of any fixed-``max_slots`` PIT diarization model.

    Because NIST framewise error satisfies
    ``miss_t + fa_t + conf_t = max(n_ref_t, n_hyp_t) - n_cor_t >= n_ref_t - n_cor_t``
    on every frame ``t``, any model restricted to at most ``max_slots`` output tracks
    can match at most ``M = min(K_ref, max_slots)`` reference speakers under optimal
    1-to-1 permutation, leaving at least the ``max(0, K_ref - max_slots)`` least-active
    reference speakers completely unmatched.
    """
    ref = _validate_binary_activity(ref_activity, "ref_activity", allow_empty_speakers=False)
    if isinstance(max_slots, (bool, np.bool_)) or not isinstance(max_slots, (int, np.integer)):
        raise ValueError("max_slots must be an integer")
    slots = int(max_slots)
    if slots < 0:
        raise ValueError("max_slots must be non-negative")
    dt = finite_real_scalar(frame_dur_s, "frame_dur_s")
    if dt <= 0.0:
        raise ValueError("frame_dur_s must be strictly positive")

    speaker_frames = np.sum(ref, axis=1)
    ref_speaker_frames = int(np.sum(speaker_frames))
    if ref_speaker_frames <= 0:
        raise ValueError("ref_activity must have positive total active reference speaker-time")

    k_ref = int(ref.shape[0])
    order_desc = np.argsort(-speaker_frames, kind="mergesort")
    matched_count = min(k_ref, slots)
    matched_sorted = np.sort(order_desc[:matched_count])
    unmatched = order_desc[matched_count:]

    min_error_frames = int(np.sum(speaker_frames[unmatched]))
    min_error_s = float(min_error_frames * dt)
    ref_time_s = float(ref_speaker_frames * dt)
    if not (np.isfinite(min_error_s) and np.isfinite(ref_time_s)):
        raise ValueError("DER bound durations exceed finite float64 range")

    oracle_hyp = ref[matched_sorted, :].copy() if matched_count > 0 else np.zeros((0, ref.shape[1]), dtype=np.int64)
    return {
        "num_ref_speakers": k_ref,
        "max_slots": slots,
        "unmatched_speaker_count": int(unmatched.size),
        "matched_speaker_indices": [int(i) for i in matched_sorted],
        "unmatched_speaker_indices": [int(i) for i in unmatched],
        "speaker_active_frames": [int(v) for v in speaker_frames],
        "min_unmatched_frames": min_error_frames,
        "ref_speaker_frames": ref_speaker_frames,
        "min_unmatched_s": min_error_s,
        "ref_speaker_time_s": ref_time_s,
        "der_lower_bound": float(min_error_frames / ref_speaker_frames),
        "oracle_slot_hypothesis": oracle_hyp,
    }


def contaminated_enrollment_metrics(
    e_target: np.ndarray,
    e_interferer: np.ndarray,
    leakage_alpha: float = 0.5,
    *,
    gate_kappa: float = 8.0,
    gate_tau: float = 0.40,
) -> dict:
    """Quantify cosine drift and mask-gate leakage under contaminated speaker enrollment.

    For unit-norm embeddings ``||e_1||_2 = ||e_2||_2 = 1`` with inner product
    ``rho_12 = e_1^T e_2`` and interferer leakage ratio ``alpha >= 0``, evaluates
    ``e_cont = (e_1 + alpha * e_2) / ||e_1 + alpha * e_2||_2`` (Eq. 17-34) and the
    induced target/interferer mask gates ``g_k = sigma(kappa * (e^T e_k - tau))``.
    """
    e1 = finite_real_array(e_target, "e_target")
    e2 = finite_real_array(e_interferer, "e_interferer")
    if e1.ndim != 1 or e2.ndim != 1 or e1.shape != e2.shape or e1.size < 2:
        raise ValueError("e_target and e_interferer must be 1-D vectors of identical length D >= 2")

    n1 = float(np.linalg.norm(e1))
    n2 = float(np.linalg.norm(e2))
    if not (np.isfinite(n1) and np.isfinite(n2)) or abs(n1 - 1.0) > 1e-6 or abs(n2 - 1.0) > 1e-6:
        raise ValueError("e_target and e_interferer must have unit L2 norm")

    alpha = finite_real_scalar(leakage_alpha, "leakage_alpha")
    kappa = finite_real_scalar(gate_kappa, "gate_kappa")
    tau = finite_real_scalar(gate_tau, "gate_tau")
    if alpha < 0.0:
        raise ValueError("leakage_alpha must be non-negative")
    if kappa < 0.0:
        raise ValueError("gate_kappa must be non-negative")

    with np.errstate(over="raise", invalid="raise"):
        try:
            e_mix = e1 + alpha * e2
            mix_norm = float(np.linalg.norm(e_mix))
        except FloatingPointError as exc:
            raise ValueError("contaminated enrollment vector overflowed float64") from exc
    if not np.isfinite(mix_norm) or mix_norm <= 1e-12:
        raise ValueError("contaminated enrollment vector cancelled to near-zero norm")

    e_cont = e_mix / mix_norm
    rho_12 = float(np.clip(np.dot(e1, e2), -1.0, 1.0))
    tar_cos = float(np.clip(np.dot(e_cont, e1), -1.0, 1.0))
    int_cos = float(np.clip(np.dot(e_cont, e2), -1.0, 1.0))

    radicand = 1.0 + 2.0 * alpha * rho_12 + alpha * alpha
    denom = float(np.sqrt(max(radicand, 1e-30)))
    analytic_tar_cos = float(np.clip((1.0 + alpha * rho_12) / denom, -1.0, 1.0))
    analytic_int_cos = float(np.clip((rho_12 + alpha) / denom, -1.0, 1.0))
    drift_deg = float(np.rad2deg(np.arccos(tar_cos)))

    with np.errstate(over="raise", invalid="raise"):
        try:
            logits = np.array(
                [
                    kappa * (1.0 - tau),
                    kappa * (rho_12 - tau),
                    kappa * (tar_cos - tau),
                    kappa * (int_cos - tau),
                ],
                dtype=np.float64,
            )
        except FloatingPointError as exc:
            raise ValueError("enrollment gate logits overflowed float64") from exc
    if not np.all(np.isfinite(logits)):
        raise ValueError("enrollment gate logits must be finite")
    gates = _stable_sigmoid(logits)
    clean_tar_gate, clean_int_gate, cont_tar_gate, cont_int_gate = (float(v) for v in gates)
    inflation = float(cont_int_gate / clean_int_gate) if clean_int_gate > 0.0 else None

    return {
        "leakage_alpha": alpha,
        "leakage_db": float(20.0 * np.log10(alpha)) if alpha > 0.0 else None,
        "clean_speaker_cosine": rho_12,
        "unnormalized_norm": mix_norm,
        "e_contaminated": e_cont,
        "target_cosine": tar_cos,
        "interferer_cosine": int_cos,
        "analytic_target_cosine": analytic_tar_cos,
        "analytic_interferer_cosine": analytic_int_cos,
        "angular_drift_deg": drift_deg,
        "clean_target_gate": clean_tar_gate,
        "clean_interferer_gate": clean_int_gate,
        "contaminated_target_gate": cont_tar_gate,
        "contaminated_interferer_gate": cont_int_gate,
        "interferer_gate_inflation_ratio": inflation,
    }


def gated_enrollment_pooling(
    frame_embeddings: np.ndarray,
    frame_snr_db: np.ndarray,
    speaker_posteriors: np.ndarray,
    *,
    target_index: int = 0,
    min_snr_db: float = 10.0,
    min_exclusive_prob: float = 0.6,
) -> dict:
    """Compare energy-VAD-only enrollment pooling against joint SNR + diarization-exclusive gating.

    Target-exclusive frame probability follows Eq. 17-25:
    ``m_{k,t}^{tar} = gamma_{k,t} * prod_{j != k} (1 - gamma_{j,t})``,
    rejecting high-SNR pauses where a competing talker speaks into the enrollment buffer.
    """
    emb = finite_real_array(frame_embeddings, "frame_embeddings")
    snr = finite_real_array(frame_snr_db, "frame_snr_db")
    gamma = finite_real_array(speaker_posteriors, "speaker_posteriors")

    if emb.ndim != 2 or emb.shape[0] < 2 or emb.shape[1] < 1:
        raise ValueError("frame_embeddings must have shape (D, T_enr) with D >= 2 and T_enr >= 1")
    n_frames = emb.shape[1]
    if snr.ndim != 1 or snr.shape[0] != n_frames:
        raise ValueError("frame_snr_db must have 1-D shape (T_enr,) matching frame_embeddings")
    if gamma.ndim != 2 or gamma.shape[1] != n_frames or gamma.shape[0] < 1:
        raise ValueError("speaker_posteriors must have shape (K, T_enr) matching frame_embeddings")
    if np.any(gamma < 0.0) or np.any(gamma > 1.0):
        raise ValueError("speaker_posteriors entries must lie in [0, 1]")

    if isinstance(target_index, (bool, np.bool_)) or not isinstance(target_index, (int, np.integer)):
        raise ValueError("target_index must be an integer")
    t_idx = int(target_index)
    if not (0 <= t_idx < gamma.shape[0]):
        raise ValueError("target_index out of range")

    snr_thr = finite_real_scalar(min_snr_db, "min_snr_db")
    excl_thr = finite_real_scalar(min_exclusive_prob, "min_exclusive_prob")
    if not (0.0 < excl_thr <= 1.0):
        raise ValueError("min_exclusive_prob must lie in (0, 1]")

    other_inactive = np.ones(n_frames, dtype=np.float64)
    for spk in range(gamma.shape[0]):
        if spk != t_idx:
            other_inactive *= (1.0 - gamma[spk, :])
    target_excl = gamma[t_idx, :] * other_inactive

    vad_mask = snr >= snr_thr
    if not np.any(vad_mask):
        raise ValueError("no enrollment frames satisfy frame_snr_db >= min_snr_db")

    joint_mask = vad_mask & (target_excl >= excl_thr)
    if not np.any(joint_mask):
        raise ValueError("no enrollment frames survive joint SNR and target-exclusive gating")

    vad_mean = np.mean(emb[:, vad_mask], axis=1)
    vad_norm = float(np.linalg.norm(vad_mean))
    if not np.isfinite(vad_norm) or vad_norm <= 1e-12:
        raise ValueError("energy-VAD pooled enrollment vector has vanishing norm")

    joint_mean = np.mean(emb[:, joint_mask], axis=1)
    joint_norm = float(np.linalg.norm(joint_mean))
    if not np.isfinite(joint_norm) or joint_norm <= 1e-12:
        raise ValueError("joint-gated pooled enrollment vector has vanishing norm")

    return {
        "target_exclusive_prob": target_excl,
        "vad_selected_indices": [int(i) for i in np.flatnonzero(vad_mask)],
        "joint_selected_indices": [int(i) for i in np.flatnonzero(joint_mask)],
        "rejected_crosstalk_indices": [int(i) for i in np.flatnonzero(vad_mask & (~joint_mask))],
        "vad_selected_count": int(np.sum(vad_mask)),
        "joint_selected_count": int(np.sum(joint_mask)),
        "vad_only_embedding": vad_mean / vad_norm,
        "gated_embedding": joint_mean / joint_norm,
    }
