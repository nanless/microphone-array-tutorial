"""Run E15-01..18 with pure NumPy and in-memory PCM16; never write files or access network.

Usage: .venv/bin/python -m codes.chapters.ch15.chapter15_exercises [--exercise E15-05]
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
from codes.chapters.ch15.core.av_alignment import (
    av_confidence_gated_mask,
    map_video_to_audio_frames,
)
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
    fresnel_distance,
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

EXERCISE_IDS = tuple(f"E15-{i:02d}" for i in range(1, 19))


def _json(value):
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return np.stack([value.real, value.imag], axis=-1).tolist()
        return [_json(v) for v in value.tolist()]
    if isinstance(value, (complex, np.complexfloating)):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.generic):
        return _json(value.item())
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
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


def run_experiments() -> dict:
    # E15-01: FiLM conditioning
    h_feat = np.array([[1.0, -0.5, 2.0], [0.5, 1.5, -1.0]], dtype=np.float64)
    z_enrolled = film_conditioning(h_feat, [2.0, 0.5], [-0.5, 0.25])
    z_identity = film_conditioning(h_feat, [1.0, 1.0], [0.0, 0.0])
    z_zero = film_conditioning(h_feat, [0.0, 0.0], [0.0, 0.0])

    # E15-02: Dual-state TSE loss & target-absent SI-SDR singularity
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

    # E15-03: Zero-overlap bypass protection
    y_pri = np.array([1.0, 0.0, -1.0, 0.0], dtype=np.float64)
    y_sec_fricative = np.array([0.0, 0.1, 0.0, -0.1], dtype=np.float64)
    y_sec_overlap = np.array([0.0, 0.5, 0.0, -0.5], dtype=np.float64)
    bypass_single = zero_overlap_bypass_decision(y_pri, y_sec_fricative, threshold_db=-18.0)
    keep_overlap = zero_overlap_bypass_decision(y_pri, y_sec_overlap, threshold_db=-18.0)

    # E15-04 & E15-05: Fresnel distance & broadside distance bubble (f = 3430 Hz, lambda = 0.10 m)
    radii_eval = np.array([0.20, 0.40, 0.80, 1.60, 4.00], dtype=np.float64)
    bub_10cm = broadside_distance_bubble_weights(0.10, 0.40, 3430.0, c=343.0, radii_m=radii_eval)
    bub_2cm = broadside_distance_bubble_weights(0.02, 0.40, 3430.0, c=343.0, radii_m=radii_eval)
    rej_10cm_db = float(-10.0 * np.log10(bub_10cm["radial_power_response"][-1]))
    rej_2cm_db = float(-10.0 * np.log10(bub_2cm["radial_power_response"][-1]))

    # E15-06: Phase-only vs 1/r spherical spreading model mismatch
    bub_10cm_rel = broadside_distance_bubble_weights(
        0.10, 0.40, 3430.0, c=343.0, radii_m=radii_eval, attenuation="relative_1_over_r"
    )
    a_true_rel = spherical_ula_steering(
        [-0.10, 0.0, 0.10], 0.40, 0.0, 3430.0, c=343.0, attenuation="relative_1_over_r"
    )
    mismatch_gain = float(abs(np.vdot(bub_10cm["weights"], a_true_rel)) ** 2)

    # E15-07: Tikhonov regularization sweep on tiny d = 0.02 m bubble
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

    # E15-08, E15-09, E15-10: Overlapped meeting DER, single-label lower bound, DER > 100%
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

    # E15-11 & E15-12: TS-VAD Bayesian fusion and GSS soft-gated MVDR
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

    # E15-13, E15-14, E15-15: AV frame cadence, camera lag onset loss, occlusion fallback
    av_sync = map_video_to_audio_frames(15, audio_hop_s=0.016, video_fps=25.0, av_offset_s=0.0)
    counts_per_vid = [int(np.sum(av_sync["video_frame_indices"] == m)) for m in range(6)]
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

    mask_audio = np.array([0.82, 0.82, 0.82], dtype=np.float64)
    mask_av_clean_turn_occ = np.array([0.96, 0.65, 0.05], dtype=np.float64)
    alpha_schedule = np.array([0.90, 0.40, 0.05], dtype=np.float64)
    fused_masks = av_confidence_gated_mask(mask_audio, mask_av_clean_turn_occ, alpha_schedule)

    # E15-16: 6-file synthetic distance-bubble audio three-layer verification in memory
    audio_report, audio_arrays = run_experiment()
    decoded_pcm = {
        key: read_pcm16(pcm16_bytes(arr, SAMPLE_RATE))[1] for key, arr in audio_arrays.items()
    }
    pcm_scores = {
        key: measure_signal(decoded_pcm[key], key, reference=decoded_pcm["near_target_ref"][0], pcm=True)
        for key in OUTPUT_KEYS
    }

    # E15-17: EEND-EDA variable-speaker attractors vs fixed 2-slot PIT failure
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

    # E15-18: Enrollment crosstalk contamination & diarization-gated pooling remedy
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
            "features": h_feat,
            "enrolled_modulated": z_enrolled,
            "identity_modulated": z_identity,
            "zero_gated": z_zero,
        },
        "E15-02": {
            "target_present": present_res,
            "standard_sisdr_singularity": sisdr_singularity,
            "target_absent_quiet": absent_quiet,
            "target_absent_moderate_leak": absent_leak,
            "target_absent_full_passthrough": absent_pass,
        },
        "E15-03": {
            "single_speaker_fricative_bypass": bypass_single,
            "two_speaker_overlap_retained": keep_overlap,
        },
        "E15-04": {
            "d_10cm": {
                "aperture_m": 0.20,
                "fresnel_distance_m": fresnel_distance(0.20, 3430.0, c=343.0),
                "exact_extra_path_mm": float(bub_10cm["extra_path_m"] * 1000.0),
                "fresnel_quadratic_approx_mm": float((0.10 ** 2) / (2.0 * 0.40) * 1000.0),
            },
            "d_2cm": {
                "aperture_m": 0.04,
                "fresnel_distance_m": fresnel_distance(0.04, 3430.0, c=343.0),
                "exact_extra_path_mm": float(bub_2cm["extra_path_m"] * 1000.0),
                "fresnel_quadratic_approx_mm": float((0.02 ** 2) / (2.0 * 0.40) * 1000.0),
            },
        },
        "E15-05": {
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
        "E15-06": {
            "outer_amplitude_ratio_r0_over_r1": float(0.40 / np.hypot(0.40, 0.10)),
            "matched_relative_1_over_r_wng_db": bub_10cm_rel["wng_db"],
            "matched_relative_1_over_r_rej_4m_db": float(
                -10.0 * np.log10(bub_10cm_rel["radial_power_response"][-1])
            ),
            "phase_weights_on_attenuated_field_power_gain": mismatch_gain,
            "phase_weights_on_attenuated_field_db": float(10.0 * np.log10(mismatch_gain)),
        },
        "E15-07": {"d_2cm_regularization_sweep": reg_sweep},
        "E15-08": {
            "imperfect_2cluster": der_imperfect,
            "label_swapped_2cluster": der_swapped,
            "permutation_invariant_der": der_imperfect["der"] == der_swapped["der"],
        },
        "E15-09": {
            "theoretical_bound": overlap_bound,
            "oracle_single_label_check": der_oracle_single,
        },
        "E15-10": {"over_triggering_3track": der_over_unit},
        "E15-11": {
            "spectral_only_posteriors": post_spec_only,
            "joint_spatial_spectral_posteriors": post_joint,
            "spectral_only_der": der_spec_only["der"],
            "joint_der": der_joint["der"],
        },
        "E15-12": {
            "gss_weights": gss_res["weights"],
            "gss_interferer_power_gain_db": float(10.0 * np.log10(max(gss_int_leak, 1e-30))),
            "naive_interferer_power_gain_db": float(10.0 * np.log10(max(naive_int_leak, 1e-30))),
        },
        "E15-13": {
            "hops_per_video_frame": av_sync["hops_per_video_frame"],
            "video_frame_indices_first_15": av_sync["video_frame_indices"],
            "audio_frames_per_video_exposure": counts_per_vid,
        },
        "E15-14": {"camera_lag_consonant_clipping_sweep": lag_sweep},
        "E15-15": {
            "stages": ["clean_frontal", "profile_turn", "hand_occlusion"],
            "alpha_v": alpha_schedule,
            "mask_audio": mask_audio,
            "mask_av": mask_av_clean_turn_occ,
            "mask_fused": fused_masks,
            "ungated_av_error_at_occlusion": float(abs(1.0 - mask_av_clean_turn_occ[2])),
            "gated_fused_error_at_occlusion": float(abs(1.0 - fused_masks[2])),
        },
        "E15-16": {
            "analytic": {k: audio_report["analytic_expectations"][k] for k in OUTPUT_KEYS},
            "float": {k: audio_report["float_measurements"][k] for k in OUTPUT_KEYS},
            "pcm": {k: pcm_scores[k]["pcm_integer_measurements"] for k in OUTPUT_KEYS},
        },
        "E15-17": {
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
        "E15-18": {
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
    args = parser.parse_args(argv)
    result = run_experiments()
    if args.exercise:
        result["exercises"] = {args.exercise: result["exercises"][args.exercise]}
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
