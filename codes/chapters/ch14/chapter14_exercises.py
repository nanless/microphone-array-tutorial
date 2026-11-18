"""Run E14-01..18 with pure NumPy and in-memory PCM; never write files or use network.

Usage: python -m codes.chapters.ch14.chapter14_exercises [--exercise E14-07]
"""
from __future__ import annotations

if __name__ == "__main__" and not __package__:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import argparse
import json
import numpy as np

from codes.chapters.ch00.core.audio_samples import pcm16_bytes, read_pcm16
from codes.chapters.ch02.core.conventions import SPEED_OF_SOUND
from codes.chapters.ch14.core.ambisonic_sofa_rendering import (
    REGULAR_HEXAGON_2D_DIRS,
    REGULAR_OCTAHEDRON_3D_DIRS,
    audit_sofa_simple_freefield_geometry,
    foa_loudspeaker_decode_metrics,
    foa_max_re_order_weight,
    head_tracking_latency_audit,
    sofa_cartesian_to_spherical,
    sofa_spherical_to_cartesian,
)
from codes.chapters.ch14.core.binaural_matching import (
    binaural_coherence_with_mixing,
    bsm_directional_metrics,
    bsm_ls_weights,
    bsm_magls_weights,
)
from codes.chapters.ch14.core.spatial_parametric import (
    TETRA_DIRECTIONS,
    active_intensity_and_diffuseness,
    encode_tetrahedral_a_to_b,
    foa_rotation_matrix,
    open_sphere_foa_radial_weights,
    reorder_foa,
    rotate_foa,
    synthesize_tetrahedral_a_format,
    tetrahedral_a_to_b_matrix,
)
from codes.chapters.ch14.core.wearable_audio import (
    SAMPLE_RATE,
    measure_signal,
    run_experiment,
)
from codes.chapters.ch14.core.wearable_fusion import (
    corcos_wind_coherence,
    cross_modal_mvdr_mmse,
    hear_through_response,
)

EXERCISE_IDS = tuple(f"E14-{i:02d}" for i in range(1, 19))


def _json(value):
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return np.stack([value.real, value.imag], axis=-1).tolist()
        return value.tolist()
    if isinstance(value, (complex, np.complexfloating)):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(item) for item in value]
    return value


def _ex01_tetra_encoding() -> dict:
    enc_wxyz = tetrahedral_a_to_b_matrix(cardioid_alpha=0.5, order="WXYZ")
    enc_acn = tetrahedral_a_to_b_matrix(cardioid_alpha=0.5, order="ACN")
    gram = TETRA_DIRECTIONS.T @ TETRA_DIRECTIONS
    col_sum = np.sum(TETRA_DIRECTIONS, axis=0)

    cases = {}
    for label, az, el in (
        ("plus_x", 0.0, 0.0),
        ("plus_y", np.pi / 2.0, 0.0),
        ("plus_z", 0.0, np.pi / 2.0),
        ("az45_el0", np.pi / 4.0, 0.0),
    ):
        a_sig = synthesize_tetrahedral_a_format(az, el, amplitude=1.0, cardioid_alpha=0.5)
        b_wxyz = encode_tetrahedral_a_to_b(a_sig, cardioid_alpha=0.5, order="WXYZ")
        b_acn = reorder_foa(b_wxyz, src_order="WXYZ", dst_order="ACN")
        cases[label] = {"a_format": a_sig, "b_wxyz": b_wxyz, "b_acn": b_acn}

    return {
        "capsule_gram_matrix": gram,
        "capsule_direction_sum": col_sum,
        "encoder_wxyz": enc_wxyz,
        "encoder_acn": enc_acn,
        "plane_wave_checks": cases,
    }


def _ex02_open_sphere_radial() -> dict:
    radius_m = 0.015
    c = SPEED_OF_SOUND
    f_null_hz = c / (2.0 * radius_m)
    f_alias_hz = c / (2.0 * np.pi * radius_m)
    freqs = np.array([100.0, 1000.0, f_alias_hz, f_null_hz], dtype=np.float64)
    kr = 2.0 * np.pi * freqs * radius_m / c
    b0, b1 = open_sphere_foa_radial_weights(kr)
    dipole_eq_db_100hz = float(-20.0 * np.log10(np.abs(b1[0])))
    return {
        "radius_m": radius_m,
        "sound_speed_mps": c,
        "first_bessel_null_hz": float(f_null_hz),
        "foa_kr1_limit_hz": float(f_alias_hz),
        "test_frequencies_hz": freqs,
        "kr_values": kr,
        "b0_weights": b0,
        "b1_weights": b1,
        "dipole_boost_required_at_100hz_db": dipole_eq_db_100hz,
    }


def _ex03_foa_rotation() -> dict:
    src_px = np.array([1.0, 1.0, 0.0, 0.0], dtype=np.float64)
    src_py = np.array([1.0, 0.0, 1.0, 0.0], dtype=np.float64)
    scene_yaw90 = rotate_foa(src_px, yaw_rad=np.pi / 2.0, listener_head_rotation=False)
    head_yaw90_on_px = rotate_foa(src_px, yaw_rad=np.pi / 2.0, listener_head_rotation=True)
    head_yaw90_on_py = rotate_foa(src_py, yaw_rad=np.pi / 2.0, listener_head_rotation=True)

    r_zy = foa_rotation_matrix(yaw_rad=np.pi / 2.0, pitch_rad=np.pi / 2.0)
    r_yz = foa_rotation_matrix(pitch_rad=np.pi / 2.0) @ foa_rotation_matrix(yaw_rad=np.pi / 2.0)
    return {
        "scene_yaw90_from_plus_x": scene_yaw90,
        "listener_yaw90_world_plus_x": head_yaw90_on_px,
        "listener_yaw90_world_plus_y": head_yaw90_on_py,
        "commutator_frobenius_norm": float(np.linalg.norm(r_zy - r_yz, ord="fro")),
    }


def _ex04_intensity_bounds() -> dict:
    traveling = active_intensity_and_diffuseness(1.0, np.array([1.0, 0.0, 0.0]))
    standing = active_intensity_and_diffuseness(2.0, np.array([0.0, 0.0, 0.0]))
    reactive = active_intensity_and_diffuseness(1.0, np.array([1.0j, 0.0, 0.0]))
    return {
        "traveling_plus_x": traveling,
        "standing_wave_antinode": standing,
        "pure_reactive_quadrature": reactive,
    }


def _ex05_two_wave_interference() -> dict:
    cases = {}
    for deg in (0, 60, 90, 180):
        phi = np.deg2rad(deg)
        p = 1.0 + np.exp(1j * phi)
        u = np.array([1.0 + 0.0j, np.exp(1j * phi), 0.0 + 0.0j])
        res = active_intensity_and_diffuseness(p, u, rho0=1.0, c=1.0)
        az_deg = float(np.rad2deg(np.arctan2(res["intensity_unit"][1], res["intensity_unit"][0])))
        cases[f"delta_phi_{deg}deg"] = {
            **res,
            "intensity_norm": float(np.linalg.norm(res["active_intensity"])),
            "phantom_intensity_azimuth_deg": az_deg if res["valid_direction"] else None,
        }
    return cases


def _ex06_dirac_stream_split() -> dict:
    # Single wave vs two uncorrelated equal-power opposite waves (+x and -x) over 4 orthogonal phases
    phases = np.array([0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi])
    p_opp = 1.0 + np.exp(1j * phases)
    u_opp = np.vstack([1.0 - np.exp(1j * phases), np.zeros(4), np.zeros(4)])
    opp = active_intensity_and_diffuseness(p_opp, u_opp, rho0=1.0, c=1.0)

    p_quad = 1.0 + 1.0j
    u_quad = np.array([1.0 + 0.0j, 1.0j, 0.0 + 0.0j])
    quad = active_intensity_and_diffuseness(p_quad, u_quad, rho0=1.0, c=1.0)
    e_tot = quad["energy_density"]
    psi = quad["diffuseness"]
    return {
        "quadrature_two_wave": {
            "total_energy": e_tot,
            "diffuseness": psi,
            "directional_energy": (1.0 - psi) * e_tot,
            "diffuse_energy": psi * e_tot,
        },
        "opposite_uncorrelated_waves": opp,
    }


def _ex07_bsm_ls_notch() -> dict:
    h_mat = np.array([[1.0, 1.0, 1.0], [1.0j, 1.0, -1.0j]], dtype=np.complex128)
    p_ear = np.array([1.0, -1.0, 1.0], dtype=np.complex128)
    w_ls = bsm_ls_weights(h_mat, p_ear)
    metrics = bsm_directional_metrics(h_mat, p_ear, w_ls)
    return {
        "gram_inverse": np.linalg.inv(h_mat @ h_mat.conj().T),
        "w_ls": w_ls,
        **metrics,
    }


def _ex08_bsm_magls_monotonic() -> dict:
    h_mat = np.array([[1.0, 1.0, 1.0], [1.0j, 1.0, -1.0j]], dtype=np.complex128)
    p_ear = np.array([1.0, -1.0, 1.0], dtype=np.complex128)
    w_mag1, hist1 = bsm_magls_weights(h_mat, np.abs(p_ear), n_iters=1)
    metrics1 = bsm_directional_metrics(h_mat, p_ear, w_mag1)

    # 4-mic circular array (r=0.025 m) at 4000 Hz across 12 azimuths
    az = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    mic_angles = np.array([0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi])
    kr = 2.0 * np.pi * 4000.0 * 0.025 / SPEED_OF_SOUND
    h_12 = np.exp(1j * kr * np.cos(mic_angles[:, None] - az[None, :]))
    # Synthetic +90 deg left-ear target with head-shadow magnitude + rapid phase advance
    p_12 = (1.0 + 0.6 * np.sin(az)) * np.exp(1j * 2.0 * kr * np.sin(az))
    w_ls_12 = bsm_ls_weights(h_12, p_12, reg=0.05)
    w_mag_12, hist_12 = bsm_magls_weights(
        h_12, np.abs(p_12), reg=0.05, n_iters=25, init_phases=np.angle(p_12)
    )
    return {
        "anchor_1step_w": w_mag1,
        "anchor_1step_metrics": metrics1,
        "multi_dir_ls_mag_mse": bsm_directional_metrics(h_12, p_12, w_ls_12)["weighted_magnitude_mse"],
        "multi_dir_magls_mag_mse": bsm_directional_metrics(h_12, p_12, w_mag_12)["weighted_magnitude_mse"],
        "objective_history": hist_12,
        "monotonically_non_increasing": bool(np.all(np.diff(hist_12) <= 0.0)),
    }


def _ex09_bsm_weighting_regularization() -> dict:
    h_mat = np.array([[1.0, 1.0, 1.0], [1.0j, 1.0, -1.0j]], dtype=np.complex128)
    p_ear = np.array([1.0, -1.0, 1.0], dtype=np.complex128)
    frontal_w = np.array([1.0, 4.0, 1.0], dtype=np.float64)
    results = {}
    for label, w_om, reg in (
        ("uniform_reg0", None, 0.0),
        ("frontal4_reg0", frontal_w, 0.0),
        ("frontal4_reg0.5", frontal_w, 0.5),
        ("frontal4_reg2.0", frontal_w, 2.0),
    ):
        w = bsm_ls_weights(h_mat, p_ear, weights=w_om, reg=reg)
        m = bsm_directional_metrics(h_mat, p_ear, w, weights=w_om)
        results[label] = {
            "weights": w,
            "weight_norm_squared": float(np.real(np.vdot(w, w))),
            "broadside_magnitude": float(m["synthesized_magnitude"][1]),
            "weighted_complex_mse": m["weighted_complex_mse"],
        }
    return results


def _ex10_bmvdr_coherence_collapse() -> dict:
    gamma_in = 0.2
    r_nn = np.array([[1.0, gamma_in], [gamma_in, 1.0]], dtype=np.complex128)
    a = np.array([1.0, 1.0], dtype=np.complex128)
    r_inv_a = np.linalg.solve(r_nn, a)
    w_bmvdr = r_inv_a / np.real(np.vdot(a, r_inv_a))
    res = binaural_coherence_with_mixing(r_nn, w_bmvdr, w_bmvdr, beta=0.0, ref_idx=(0, 1))
    snr_gain_db = float(-10.0 * np.log10(res["var_L"]))
    return {
        "gamma_in": gamma_in,
        "w_bmvdr_L": w_bmvdr,
        "w_bmvdr_R": w_bmvdr,
        "output_noise_variance": res["var_L"],
        "snr_gain_db": snr_gain_db,
        "ic_out": res["ic_out"],
    }


def _ex11_reference_mixing_tradeoff() -> dict:
    gamma_in = 0.2
    r_nn = np.array([[1.0, gamma_in], [gamma_in, 1.0]], dtype=np.complex128)
    w_bmvdr = np.array([0.5, 0.5], dtype=np.complex128)
    sweep = {}
    for beta in (0.0, 0.2, 0.5, float(1.0 / np.sqrt(2.0)), 1.0):
        res = binaural_coherence_with_mixing(r_nn, w_bmvdr, w_bmvdr, beta=beta, ref_idx=(0, 1))
        key = f"beta_{beta:.5f}"
        sweep[key] = {
            "beta": beta,
            "var_out": res["var_L"],
            "cross_LR": float(np.real(res["cross_covariance"])),
            "ic_out": res["ic_out"],
            "snr_gain_db": float(-10.0 * np.log10(res["var_L"])),
        }
    return {"beta_for_ic_0p5": float(1.0 / np.sqrt(2.0)), "sweep": sweep}


def _ex12_wind_coherence() -> dict:
    freqs = np.array([100.0, 250.0, 500.0, 1000.0, 2000.0], dtype=np.float64)
    return corcos_wind_coherence(freqs, d_m=0.015, u_conv_mps=3.0, alpha_long=0.14)


def _ex13_cross_modal_fusion() -> dict:
    low = cross_modal_mvdr_mmse(2.0, np.diag([4.0, 0.04]), phi_s=1.0)
    high = cross_modal_mvdr_mmse(0.1, np.diag([0.05, 0.04]), phi_s=1.0)
    return {
        "low_250hz": low,
        "high_4000hz": high,
        "naive_eq_high_noise_boost_db": float(10.0 * np.log10(high["equalized_inear_noise_variance"] / 0.04)),
    }


def _ex14_hearthrough_comb() -> dict:
    freqs = np.array([125.0, 250.0, 500.0, 750.0, 1000.0, 2000.0, 4000.0], dtype=np.float64)
    resp = hear_through_response(freqs, tau_proc_s=0.002, g_leak=1.0, g_ht=1.0)
    tau_max_phase30_2khz_us = 1e6 / (12.0 * 2000.0)
    tau_max_notch_above_4khz_ms = 1e3 / (2.0 * 4000.0)
    return {
        **resp,
        "tau_max_phase30deg_at_2khz_us": float(tau_max_phase30_2khz_us),
        "tau_max_first_notch_above_4khz_ms": float(tau_max_notch_above_4khz_ms),
    }


def _ex15_wearable_resource_budget() -> dict:
    fs_comm = 16000
    n_fft = 256
    hop = 128
    bins = n_fft // 2 + 1
    m_chan = 3  # 2 outer air mics + 1 occluded in-ear mic
    cov_complex64_bytes = bins * (m_chan**2) * 8
    cov_complex128_bytes = bins * (m_chan**2) * 16
    return {
        "uplink_stft_pipeline": {
            "sample_rate_hz": fs_comm,
            "n_fft": n_fft,
            "hop_samples": hop,
            "window_duration_ms": 1000.0 * n_fft / fs_comm,
            "frame_advance_period_ms": 1000.0 * hop / fs_comm,
            "frequency_bins": bins,
            "channels": m_chan,
            "covariance_complex64_bytes": cov_complex64_bytes,
            "covariance_complex128_bytes": cov_complex128_bytes,
        },
        "hearthrough_timedomain_pipeline": {
            "sample_rate_hz": 96000,
            "fir_taps": 16,
            "sample_period_us": 1e6 / 96000.0,
            "macs_per_second_per_ear": 96000 * 16,
        },
    }


def _ex16_audio_three_tier_audit() -> dict:
    report, signals = run_experiment()
    decoded = {k: read_pcm16(pcm16_bytes(v, SAMPLE_RATE))[1] for k, v in signals.items()}
    pcm_ref = decoded["wearable_clean_ref"]
    pcm_scores = {
        k: measure_signal(decoded[k], k, reference=pcm_ref, pcm=True) for k in signals
    }
    return {
        "analytic_expectations": report["analytic_expectations"],
        "float_measurements": report["float_measurements"],
        "pcm_measurements": pcm_scores,
    }


def _ex17_ambisonics_max_re_decoding() -> dict:
    src_px = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    modes_2d = {}
    for mode in ("basic", "max_re", "in_phase"):
        g1 = foa_max_re_order_weight(spatial_dim=2, mode=mode)
        modes_2d[mode] = {
            "g1": g1,
            **foa_loudspeaker_decode_metrics(REGULAR_HEXAGON_2D_DIRS, src_px, spatial_dim=2, g1=g1),
        }

    modes_3d = {}
    for mode in ("basic", "max_re", "in_phase"):
        g1 = foa_max_re_order_weight(spatial_dim=3, mode=mode)
        modes_3d[mode] = {
            "g1": g1,
            **foa_loudspeaker_decode_metrics(REGULAR_OCTAHEDRON_3D_DIRS, src_px, spatial_dim=3, g1=g1),
        }

    rear_power_ratio_2d = float((modes_2d["max_re"]["gains"][3] / modes_2d["basic"]["gains"][3]) ** 2)
    rear_power_ratio_3d = float((modes_3d["max_re"]["gains"][1] / modes_3d["basic"]["gains"][1]) ** 2)
    return {
        "hexagon_2d": {
            "modes": modes_2d,
            "max_re_rear_power_suppression_db": float(10.0 * np.log10(rear_power_ratio_2d)),
        },
        "octahedron_3d": {
            "modes": modes_3d,
            "max_re_rear_power_suppression_db": float(10.0 * np.log10(rear_power_ratio_3d)),
        },
    }


def _ex18_sofa_and_head_tracking_latency() -> dict:
    src_sph = np.array([30.0, 45.0, 1.5], dtype=np.float64)
    src_cart = sofa_spherical_to_cartesian(src_sph)
    src_roundtrip = sofa_cartesian_to_spherical(src_cart)

    valid_sph_ears = np.array([[90.0, 0.0, 0.09], [-90.0, 0.0, 0.09]], dtype=np.float64)
    swapped_cart_ears = np.array([[0.0, -0.09, 0.0], [0.0, 0.09, 0.0]], dtype=np.float64)

    sofa_valid = audit_sofa_simple_freefield_geometry(
        [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], valid_sph_ears, coord_type="spherical"
    )
    sofa_ear_swap = audit_sofa_simple_freefield_geometry(
        [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], swapped_cart_ears, coord_type="cartesian"
    )
    sofa_missing_type = audit_sofa_simple_freefield_geometry(
        [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], valid_sph_ears, coord_type="cartesian"
    )

    low_lat = head_tracking_latency_audit(180.0, 0.020, maa_deg=2.0)
    high_lat_lookahead = head_tracking_latency_audit(
        180.0, 0.080, maa_deg=2.0, prediction_horizon_s=0.065
    )
    return {
        "source_position_spherical": src_sph,
        "source_position_cartesian_m": src_cart,
        "source_position_roundtrip": src_roundtrip,
        "sofa_valid_audit": sofa_valid,
        "sofa_ear_swap_audit": sofa_ear_swap,
        "sofa_missing_spherical_type_audit": sofa_missing_type,
        "latency_20ms": low_lat,
        "latency_80ms_with_65ms_lookahead": high_lat_lookahead,
    }


def run_experiments() -> dict:
    exercises = {
        "E14-01": _ex01_tetra_encoding(),
        "E14-02": _ex02_open_sphere_radial(),
        "E14-03": _ex03_foa_rotation(),
        "E14-04": _ex04_intensity_bounds(),
        "E14-05": _ex05_two_wave_interference(),
        "E14-06": _ex06_dirac_stream_split(),
        "E14-07": _ex07_bsm_ls_notch(),
        "E14-08": _ex08_bsm_magls_monotonic(),
        "E14-09": _ex09_bsm_weighting_regularization(),
        "E14-10": _ex10_bmvdr_coherence_collapse(),
        "E14-11": _ex11_reference_mixing_tradeoff(),
        "E14-12": _ex12_wind_coherence(),
        "E14-13": _ex13_cross_modal_fusion(),
        "E14-14": _ex14_hearthrough_comb(),
        "E14-15": _ex15_wearable_resource_budget(),
        "E14-16": _ex16_audio_three_tier_audit(),
        "E14-17": _ex17_ambisonics_max_re_decoding(),
        "E14-18": _ex18_sofa_and_head_tracking_latency(),
    }
    if set(exercises) != set(EXERCISE_IDS):
        raise AssertionError("exercise ID set differs")
    return _json(
        {
            "schema_version": 1,
            "scope": "original finite known-model spatial capture and wearable array teaching controls",
            "complex_encoding": "complex scalar -> [real,imag], complex arrays -> trailing [real,imag] axis",
            "exercises": exercises,
        }
    )


def main(argv=None) -> int:
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
