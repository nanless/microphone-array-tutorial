"""Run E15-01..15 with pure NumPy and in-memory PCM16; never write files or access network.

Usage: .venv/bin/python -m codes.chapters.ch15.chapter15_exercises [--exercise E15-02]
"""

from __future__ import annotations

if __name__ == "__main__" and not __package__:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import argparse
import json
import math
import numpy as np

from codes.chapters.ch00.core.audio_samples import pcm16_bytes, read_pcm16
from codes.chapters.ch15.core.av_alignment import map_video_to_audio_frames
from codes.chapters.ch15.core.conditional_tse import (
    film_conditioning,
    inspect_standard_sisdr,
    tse_dual_state_loss,
    zero_overlap_bypass_decision,
)
from codes.chapters.ch15.core.diarization_tsvad import (
    compute_der_components,
    gss_soft_gated_mvdr,
    single_label_overlap_der_bound,
    tsvad_spatial_spectral_posterior,
)
from codes.chapters.ch15.core.distance_bubble import (
    broadside_distance_bubble_weights,
    spherical_ula_steering,
)
from codes.chapters.ch15.core.eend_enrollment_robustness import (
    contaminated_enrollment_metrics,
    eend_eda_inference,
    fixed_slot_der_lower_bound,
    gated_enrollment_pooling,
)
from codes.chapters.ch15.core.tse_audio import (
    OUTPUT_KEYS,
    SAMPLE_RATE,
    measure_signal,
    run_experiment,
)

EXERCISE_IDS = tuple(f"E15-{i:02d}" for i in range(1, 16))


def _json(value):
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return _json(np.stack([value.real, value.imag], axis=-1).tolist())
        return [_json(v) for v in value.tolist()]
    if isinstance(value, (complex, np.complexfloating)):
        return [_json(float(value.real)), _json(float(value.imag))]
    if isinstance(value, np.generic):
        return _json(value.item())
    if isinstance(value, float):
        if math.isnan(value):
            raise ValueError("NaN is not permitted in Chapter 15 exercise JSON output")
        if math.isinf(value):
            return "-Infinity" if value < 0 else "Infinity"
        return value
    if isinstance(value, dict):
        return {key: _json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(item) for item in value]
    return value


def _meeting_reference() -> np.ndarray:
    """Six-frame (0.5 s/frame = 3.0 s) 2-speaker meeting with 2 overlap frames and 1 silent frame."""
    return np.array(
        [
            [1, 1, 1, 1, 0, 0],  # Speaker 1: 2.0 s active
            [0, 0, 1, 1, 1, 0],  # Speaker 2: 1.5 s active -> total ref speaker-time = 3.5 s
        ],
        dtype=np.int64,
    )


def _tsvad_meeting_cues() -> tuple[np.ndarray, np.ndarray]:
    """Deterministic 6-frame spatial IPD cosines and speaker-embedding similarities."""
    obs_ipd_cos = np.array(
        [
            [0.95, 0.92, 0.82, 0.80, -0.40, 0.00],  # Speaker 1 spatial alignment
            [-0.45, -0.40, 0.78, 0.81, 0.94, -0.05],  # Speaker 2 spatial alignment
        ],
        dtype=np.float64,
    )
    emb_sim = np.array(
        [
            [0.88, 0.86, 0.48, 0.47, 0.12, 0.08],  # Overlap degrades single-vector embedding sim!
            [0.10, 0.14, 0.46, 0.49, 0.89, 0.06],
        ],
        dtype=np.float64,
    )
    return obs_ipd_cos, emb_sim


def _ex15_film_and_visual_confidence_gating() -> dict:
    """E15-15: FiLM affine conditioning contrast swing & lip-occlusion visual confidence gating."""
    h_vec = np.array([2.0, -1.0], dtype=np.float64)
    w_gamma = np.array([[1.5, 0.2], [0.2, 1.5]], dtype=np.float64)
    w_beta = np.array([[0.1, 0.0], [0.0, 0.1]], dtype=np.float64)
    e_a = np.array([1.0, 0.0], dtype=np.float64)
    e_b = np.array([0.0, 1.0], dtype=np.float64)

    gamma_a = w_gamma @ e_a
    beta_a = w_beta @ e_a
    gamma_b = w_gamma @ e_b
    beta_b = w_beta @ e_b

    def _eval_scale(h_in: np.ndarray) -> dict:
        h_film_a = film_conditioning(h_in[:, None], gamma_a, beta_a)[:, 0]
        h_film_b = film_conditioning(h_in[:, None], gamma_b, beta_b)[:, 0]
        ratio_a = float((h_film_a[0] ** 2) / (h_film_a[1] ** 2))
        ratio_b = float((h_film_b[0] ** 2) / (h_film_b[1] ** 2))
        ratio_a_db = float(10.0 * math.log10(ratio_a))
        ratio_b_db = float(10.0 * math.log10(ratio_b))

        h_add_a = h_in + beta_a
        h_add_b = h_in + beta_b
        add_ratio_a = float((h_add_a[0] ** 2) / (h_add_a[1] ** 2))
        add_ratio_b = float((h_add_b[0] ** 2) / (h_add_b[1] ** 2))
        add_ratio_a_db = float(10.0 * math.log10(add_ratio_a))
        add_ratio_b_db = float(10.0 * math.log10(add_ratio_b))

        return {
            "input_h": h_in,
            "input_channel_energy_ratio": float((h_in[0] ** 2) / (h_in[1] ** 2)),
            "input_channel_energy_ratio_db": float(10.0 * math.log10((h_in[0] ** 2) / (h_in[1] ** 2))),
            "film_target_A": {
                "gamma": gamma_a,
                "beta": beta_a,
                "output_h": h_film_a,
                "channel_energy_ratio": ratio_a,
                "channel_energy_ratio_db": ratio_a_db,
            },
            "film_target_B": {
                "gamma": gamma_b,
                "beta": beta_b,
                "output_h": h_film_b,
                "channel_energy_ratio": ratio_b,
                "channel_energy_ratio_db": ratio_b_db,
            },
            "film_contrast_swing_db": ratio_a_db - ratio_b_db,
            "additive_bias_only_target_A": {
                "output_h": h_add_a,
                "channel_energy_ratio": add_ratio_a,
                "channel_energy_ratio_db": add_ratio_a_db,
            },
            "additive_bias_only_target_B": {
                "output_h": h_add_b,
                "channel_energy_ratio": add_ratio_b,
                "channel_energy_ratio_db": add_ratio_b_db,
            },
            "additive_bias_only_swing_db": add_ratio_a_db - add_ratio_b_db,
        }

    # Part 2: Log-odds AV speech-presence fusion with visual confidence gating fallback
    def _sigmoid(x: float) -> float:
        return 1.0 / (1.0 + math.exp(-x))

    ell_a = 2.0
    ell_v = -4.0
    c_min = 0.3
    c_max = 0.8

    def _gate_weight(c_v: float) -> float:
        return float(np.clip((c_v - c_min) / (c_max - c_min), 0.0, 1.0))

    def _fuse(w_v: float, c_v: float | None = None) -> dict:
        ell_av = ell_a + w_v * ell_v
        p_av = _sigmoid(ell_av)
        out = {
            "visual_weight": float(w_v),
            "fused_log_odds": float(ell_av),
            "fused_speech_prob": float(p_av),
            "target_detected": bool(p_av >= 0.5),
        }
        if c_v is not None:
            out["visual_confidence"] = float(c_v)
        return out

    w_crit = -ell_a / ell_v
    c_crit = c_min + w_crit * (c_max - c_min)

    return {
        "film_conditioning": {
            "nominal_scale": _eval_scale(h_vec),
            "scaled_plus20db_10x_amplitude": _eval_scale(10.0 * h_vec),
            "pure_multiplicative_asymptotic_swing_db": float(40.0 * math.log10(1.5 / 0.2)),
        },
        "visual_confidence_gating": {
            "audio_log_odds": ell_a,
            "audio_speech_prob": _sigmoid(ell_a),
            "visual_log_odds": ell_v,
            "visual_speech_prob": _sigmoid(ell_v),
            "c_min": c_min,
            "c_max": c_max,
            "critical_flip_weight": float(w_crit),
            "critical_flip_confidence": float(c_crit),
            "ungated_wv_1p0": _fuse(1.0),
            "severe_occlusion_cv_0p25": _fuse(_gate_weight(0.25), c_v=0.25),
            "partial_occlusion_cv_0p45": _fuse(_gate_weight(0.45), c_v=0.45),
        },
    }


def run_experiments() -> dict:
    # E15-01: Dual-state TSE loss & target-absent SI-SDR singularity
    s_ref_act = np.array([1.0, 1.0, -1.0, -1.0], dtype=np.float64)
    v_int = np.array([1.0, -1.0, 1.0, -1.0], dtype=np.float64)
    x_mix_act = s_ref_act + v_int
    s_hat_act = 0.5 * s_ref_act + 0.05 * v_int
    present_res = tse_dual_state_loss(s_hat_act, s_ref_act, x_mix_act)

    s_ref_zero = np.zeros(4, dtype=np.float64)
    sisdr_singularity = inspect_standard_sisdr(0.1 * v_int, s_ref_zero)
    absent_quiet = tse_dual_state_loss(np.zeros(4), s_ref_zero, v_int, tau_floor=1e-3)
    absent_leak = tse_dual_state_loss(np.sqrt(0.099) * v_int, s_ref_zero, v_int, tau_floor=1e-3)
    absent_pass = tse_dual_state_loss(v_int, s_ref_zero, v_int, tau_floor=1e-3)

    # E15-02: Zero-overlap bypass protection
    y_pri = np.array([1.0, 0.0, -1.0, 0.0], dtype=np.float64)
    y_sec_fricative = np.array([0.0, 0.1, 0.0, -0.1], dtype=np.float64)
    y_sec_overlap = np.array([0.0, 0.5, 0.0, -0.5], dtype=np.float64)
    bypass_single = zero_overlap_bypass_decision(y_pri, y_sec_fricative, threshold_db=-18.0)
    keep_overlap = zero_overlap_bypass_decision(y_pri, y_sec_overlap, threshold_db=-18.0)

    # E15-03: Broadside distance bubble (f = 3430 Hz, lambda = 0.10 m)
    radii_eval = np.array([0.20, 0.40, 0.80, 1.60, 4.00], dtype=np.float64)
    bub_10cm = broadside_distance_bubble_weights(0.10, 0.40, 3430.0, c=343.0, radii_m=radii_eval)
    bub_2cm = broadside_distance_bubble_weights(0.02, 0.40, 3430.0, c=343.0, radii_m=radii_eval)
    rej_10cm_db = float(-10.0 * np.log10(bub_10cm["radial_power_response"][-1]))
    rej_2cm_db = float(-10.0 * np.log10(bub_2cm["radial_power_response"][-1]))

    # E15-04: Phase-only vs 1/r spherical spreading model mismatch
    bub_10cm_rel = broadside_distance_bubble_weights(
        0.10, 0.40, 3430.0, c=343.0, radii_m=radii_eval, attenuation="relative_1_over_r"
    )
    a_true_rel = spherical_ula_steering(
        [-0.10, 0.0, 0.10], 0.40, 0.0, 3430.0, c=343.0, attenuation="relative_1_over_r"
    )
    mismatch_gain = float(abs(np.vdot(bub_10cm["weights"], a_true_rel)) ** 2)

    # E15-05: Tikhonov regularization sweep on tiny d = 0.02 m bubble
    reg_sweep = []
    for reg_val in (0.0, 1e-4, 5e-4, 2e-3):
        res = broadside_distance_bubble_weights(
            0.02, 0.40, 3430.0, c=343.0, reg=reg_val, radii_m=np.array([0.40, 4.00])
        )
        reg_sweep.append(
            {
                "reg": reg_val,
                "weight_norm_l2": float(np.sqrt(res["weight_norm_squared"])),
                "target_power_gain": res["target_power_gain"],
                "interferer_4m_power_gain": float(res["radial_power_response"][1]),
                "wng_db": res["wng_db"],
            }
        )

    # E15-06, E15-07, E15-08: Overlapped meeting DER, single-label lower bound, DER > 100%
    ref_act = _meeting_reference()
    hyp_imperfect = np.array([[1, 1, 1, 0, 1, 0], [0, 0, 0, 1, 0, 1]], dtype=np.int64)
    hyp_swapped = hyp_imperfect[::-1, :]
    der_imperfect = compute_der_components(ref_act, hyp_imperfect, frame_dur_s=0.5)
    der_swapped = compute_der_components(ref_act, hyp_swapped, frame_dur_s=0.5)

    overlap_bound = single_label_overlap_der_bound(ref_act, frame_dur_s=0.5)
    hyp_oracle_single = np.array([[1, 1, 1, 1, 0, 0], [0, 0, 0, 0, 1, 0]], dtype=np.int64)
    der_oracle_single = compute_der_components(ref_act, hyp_oracle_single, frame_dur_s=0.5)

    hyp_over_unit = np.array(
        [
            [0, 0, 0, 0, 1, 1],
            [0, 0, 1, 0, 1, 1],
            [1, 1, 0, 0, 1, 1],
        ],
        dtype=np.int64,
    )
    der_over_unit = compute_der_components(ref_act, hyp_over_unit, frame_dur_s=0.5)

    # E15-09 & E15-10: TS-VAD Bayesian fusion and GSS soft-gated MVDR
    obs_ipd_cos, emb_sim = _tsvad_meeting_cues()
    post_spec_only = tsvad_spatial_spectral_posterior(
        obs_ipd_cos, emb_sim, prior=0.5, kappa_sp=0.0, kappa_emb=4.0, tau_emb=0.45
    )
    post_joint = tsvad_spatial_spectral_posterior(
        obs_ipd_cos, emb_sim, prior=0.5, kappa_sp=3.0, tau_sp=0.25, kappa_emb=4.0, tau_emb=0.45
    )
    der_spec_only = compute_der_components(ref_act, (post_spec_only >= 0.5).astype(np.int64), 0.5)
    der_joint = compute_der_components(ref_act, (post_joint >= 0.5).astype(np.int64), 0.5)

    # 2-mic narrowband snapshots over 60 subframes (10 per 0.5 s block) for GSS-MVDR check
    a_spk1 = np.array([1.0 + 0.0j, np.exp(-1j * np.pi / 3.0)], dtype=np.complex128)
    a_spk2 = np.array([1.0 + 0.0j, np.exp(+1j * np.pi / 2.5)], dtype=np.complex128)
    sub_idx = np.arange(60, dtype=np.float64)
    blk_idx = (sub_idx // 10).astype(np.int64)
    s1_sub = ref_act[0, blk_idx] * np.exp(1j * 2.0 * np.pi * 0.13 * sub_idx)
    s2_sub = ref_act[1, blk_idx] * np.exp(1j * 2.0 * np.pi * 0.29 * sub_idx)
    n_sub = 0.05 * np.vstack(
        [np.exp(1j * 2.0 * np.pi * 0.41 * sub_idx), np.exp(1j * 2.0 * np.pi * 0.53 * sub_idx)]
    )
    x_sub = np.outer(a_spk1, s1_sub) + np.outer(a_spk2, s2_sub) + n_sub
    gamma_sub = post_joint[:, blk_idx]
    gss_res = gss_soft_gated_mvdr(x_sub, gamma_sub, target_index=0, diag_loading=1e-3, ref_mic=0)
    # Compare with naive un-excluded target mask (setting speaker-2 posterior to 0 in target mask)
    naive_gamma = np.vstack([gamma_sub[0], np.zeros_like(gamma_sub[1])])
    naive_res = gss_soft_gated_mvdr(x_sub, naive_gamma, target_index=0, diag_loading=1e-3, ref_mic=0)
    gss_int_leak = float(abs(np.vdot(gss_res["weights"], a_spk2)) ** 2)
    naive_int_leak = float(abs(np.vdot(naive_res["weights"], a_spk2)) ** 2)

    # E15-11: Camera lag onset consonant clipping sweep
    lag_sweep = []
    for lag_s in (0.0, 0.040, 0.080, 0.120):
        res_lag = map_video_to_audio_frames(
            30, audio_hop_s=0.016, video_fps=25.0, av_offset_s=lag_s, onset_frame=10, consonant_frames=8
        )
        lag_sweep.append(
            {
                "av_offset_ms": float(lag_s * 1000.0),
                "clipped_consonant_frames": res_lag["onset_audit"]["clipped_consonant_frames"],
                "clipped_duration_ms": float(res_lag["onset_audit"]["clipped_duration_s"] * 1000.0),
                "consonant_loss_fraction": res_lag["onset_audit"]["consonant_loss_fraction"],
            }
        )

    # E15-12: 6-file synthetic distance-bubble audio three-layer verification in memory
    audio_report, audio_arrays = run_experiment()
    decoded_pcm = {
        key: read_pcm16(pcm16_bytes(arr, SAMPLE_RATE))[1] for key, arr in audio_arrays.items()
    }
    pcm_scores = {
        key: measure_signal(decoded_pcm[key], key, reference=decoded_pcm["near_target_ref"][0], pcm=True)
        for key in OUTPUT_KEYS
    }

    # E15-13: EEND-EDA variable-speaker attractors vs fixed 2-slot PIT failure
    w_exist = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    attr_2spk = np.array([[2.0, 0.0, 0.1], [0.0, 2.0, 0.1], [0.0, 0.0, 0.0]], dtype=np.float64)
    emb_3spk = np.array(
        [
            [2.0, 2.0, 2.0, -2.0],
            [-2.0, 2.0, 2.0, -2.0],
            [-2.0, -2.0, 2.0, 2.0],
        ],
        dtype=np.float64,
    )
    eda_2spk_stop = eend_eda_inference(emb_3spk, attr_2spk, w_exist, b_exist=-1.0)

    ref_3spk = np.array(
        [
            [1, 1, 1, 0],
            [0, 1, 1, 0],
            [0, 0, 1, 1],
        ],
        dtype=np.int64,
    )
    attr_3spk = np.array(
        [[2.0, 0.0, 0.0, 0.1], [0.0, 2.0, 0.0, 0.1], [0.0, 0.0, 2.0, 0.1]],
        dtype=np.float64,
    )
    eda_3spk = eend_eda_inference(emb_3spk, attr_3spk, w_exist, b_exist=-1.0)
    der_eda_3spk = compute_der_components(ref_3spk, eda_3spk["binary_activity"], frame_dur_s=0.5)
    pit_2slot_bound = fixed_slot_der_lower_bound(ref_3spk, max_slots=2, frame_dur_s=0.5)
    der_pit_oracle = compute_der_components(
        ref_3spk, pit_2slot_bound["oracle_slot_hypothesis"], frame_dur_s=0.5
    )

    # E15-14: Enrollment crosstalk contamination & diarization-gated pooling remedy
    e_tar = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    e_int_orth = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    e_int_corr = np.array([0.3, math.sqrt(0.91), 0.0], dtype=np.float64)
    cont_orth = contaminated_enrollment_metrics(e_tar, e_int_orth, leakage_alpha=0.5)
    cont_corr = contaminated_enrollment_metrics(e_tar, e_int_corr, leakage_alpha=0.5)

    enr_frames = np.column_stack([e_tar, e_tar, e_tar, e_tar, e_int_orth, e_int_orth])
    enr_snr_db = np.array([20.0, 20.0, 20.0, 20.0, 15.0, 15.0], dtype=np.float64)
    enr_gamma = np.array(
        [
            [0.95, 0.95, 0.95, 0.95, 0.05, 0.05],
            [0.04, 0.04, 0.04, 0.04, 0.92, 0.92],
        ],
        dtype=np.float64,
    )
    enr_pool = gated_enrollment_pooling(
        enr_frames, enr_snr_db, enr_gamma, target_index=0, min_snr_db=10.0, min_exclusive_prob=0.6
    )

    exercises = {
        "E15-01": {
            "target_present": present_res,
            "standard_sisdr_singularity": sisdr_singularity,
            "target_absent_quiet": absent_quiet,
            "target_absent_moderate_leak": absent_leak,
            "target_absent_full_passthrough": absent_pass,
        },
        "E15-02": {
            "single_speaker_fricative_bypass": bypass_single,
            "two_speaker_overlap_retained": keep_overlap,
        },
        "E15-03": {
            "d_10cm": {
                "outer_phase_deg": bub_10cm["outer_phase_deg"],
                "weights": bub_10cm["weights"],
                "wng_db": bub_10cm["wng_db"],
                "rejection_4m_db": rej_10cm_db,
                "farfield_plane_power_gain": bub_10cm["farfield_plane_power_gain"],
            },
            "d_2cm": {
                "outer_phase_deg": bub_2cm["outer_phase_deg"],
                "weights": bub_2cm["weights"],
                "wng_db": bub_2cm["wng_db"],
                "rejection_4m_db": rej_2cm_db,
                "farfield_plane_power_gain": bub_2cm["farfield_plane_power_gain"],
            },
        },
        "E15-04": {
            "outer_amplitude_ratio_r0_over_r1": float(0.40 / np.hypot(0.40, 0.10)),
            "matched_relative_1_over_r_wng_db": bub_10cm_rel["wng_db"],
            "matched_relative_1_over_r_rej_4m_db": float(
                -10.0 * np.log10(bub_10cm_rel["radial_power_response"][-1])
            ),
            "phase_weights_on_attenuated_field_power_gain": mismatch_gain,
            "phase_weights_on_attenuated_field_db": float(10.0 * np.log10(mismatch_gain)),
        },
        "E15-05": {"d_2cm_regularization_sweep": reg_sweep},
        "E15-06": {
            "imperfect_2cluster": der_imperfect,
            "label_swapped_2cluster": der_swapped,
            "permutation_invariant_der": der_imperfect["der"] == der_swapped["der"],
        },
        "E15-07": {
            "theoretical_bound": overlap_bound,
            "oracle_single_label_check": der_oracle_single,
        },
        "E15-08": {"over_triggering_3track": der_over_unit},
        "E15-09": {
            "spectral_only_posteriors": post_spec_only,
            "joint_spatial_spectral_posteriors": post_joint,
            "spectral_only_der": der_spec_only["der"],
            "joint_der": der_joint["der"],
        },
        "E15-10": {
            "gss_weights": gss_res["weights"],
            "gss_interferer_power_gain_db": float(10.0 * np.log10(max(gss_int_leak, 1e-30))),
            "naive_interferer_power_gain_db": float(10.0 * np.log10(max(naive_int_leak, 1e-30))),
        },
        "E15-11": {"camera_lag_consonant_clipping_sweep": lag_sweep},
        "E15-12": {
            "analytic": {k: audio_report["analytic_expectations"][k] for k in OUTPUT_KEYS},
            "float": {k: audio_report["float_measurements"][k] for k in OUTPUT_KEYS},
            "pcm": {k: pcm_scores[k]["pcm_integer_measurements"] for k in OUTPUT_KEYS},
        },
        "E15-13": {
            "two_speaker_stopping_check": {
                "existence_logits": eda_2spk_stop["existence_logits"],
                "existence_probs": eda_2spk_stop["existence_probs"],
                "estimated_num_speakers": eda_2spk_stop["estimated_num_speakers"],
            },
            "three_speaker_eend_eda": {
                **eda_3spk,
                "der": der_eda_3spk["der"],
            },
            "fixed_2slot_pit_bound": pit_2slot_bound,
            "fixed_2slot_oracle_der": der_pit_oracle["der"],
        },
        "E15-14": {
            "orthogonal_speakers_alpha0p5": cont_orth,
            "correlated_rho0p3_alpha0p5": cont_corr,
            "six_frame_enrollment_remedy": {
                **enr_pool,
                "vad_only_target_cosine": float(np.dot(enr_pool["vad_only_embedding"], e_tar)),
                "vad_only_interferer_cosine": float(np.dot(enr_pool["vad_only_embedding"], e_int_orth)),
                "gated_target_cosine": float(np.dot(enr_pool["gated_embedding"], e_tar)),
                "gated_interferer_cosine": float(np.dot(enr_pool["gated_embedding"], e_int_orth)),
            },
        },
        "E15-15": _ex15_film_and_visual_confidence_gating(),
    }
    if set(exercises) != set(EXERCISE_IDS):
        raise AssertionError("Chapter 15 exercise ID set mismatch")
    return _json(
        {
            "schema_version": 1,
            "scope": "Chapter 15 deterministic TSE, near-field bubble, TS-VAD/GSS, and AV controls",
            "complex_encoding": "complex scalar -> [real,imag], complex arrays -> trailing [real,imag] axis",
            "exercises": exercises,
        }
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exercise", choices=EXERCISE_IDS)
    parser.add_argument("--json", action="store_true", help="emit JSON output (default)")
    parser.add_argument("--threshold-db", type=float, default=-18.0, help="bypass threshold in dB for E15-02 check")
    args = parser.parse_args(argv)
    if isinstance(args.threshold_db, bool) or not math.isfinite(args.threshold_db):
        raise ValueError("threshold_db must be a finite real number")
    result = run_experiments()
    if args.threshold_db != -18.0:
        y_pri = np.array([1.0, 0.0, -1.0, 0.0], dtype=np.float64)
        y_sec_fricative = np.array([0.0, 0.1, 0.0, -0.1], dtype=np.float64)
        y_sec_overlap = np.array([0.0, 0.5, 0.0, -0.5], dtype=np.float64)
        result["exercises"]["E15-02"] = _json(
            {
                "single_speaker_fricative_bypass": zero_overlap_bypass_decision(
                    y_pri, y_sec_fricative, threshold_db=args.threshold_db
                ),
                "two_speaker_overlap_retained": zero_overlap_bypass_decision(
                    y_pri, y_sec_overlap, threshold_db=args.threshold_db
                ),
            }
        )
    if args.exercise:
        result["exercises"] = {args.exercise: result["exercises"][args.exercise]}
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
