"""Independent hand-calculated unit tests for Chapter 16 spatial & wearable modules."""
from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np

from codes.chapters.ch16.chapter16_exercises import EXERCISE_IDS, main as cli_main, run_experiments
from codes.chapters.ch16.core.ambisonic_sofa_rendering import (
    REGULAR_HEXAGON_2D_DIRS,
    REGULAR_OCTAHEDRON_3D_DIRS,
    audit_sofa_simple_freefield_geometry,
    foa_loudspeaker_decode_metrics,
    foa_max_re_order_weight,
    head_tracking_latency_audit,
    sofa_cartesian_to_spherical,
    sofa_spherical_to_cartesian,
)
from codes.chapters.ch16.core.binaural_matching import (
    binaural_coherence_with_mixing,
    bsm_directional_metrics,
    bsm_ls_weights,
    bsm_magls_weights,
)
from codes.chapters.ch16.core.spatial_parametric import (
    TETRA_DIRECTIONS,
    active_intensity_and_diffuseness,
    encode_tetrahedral_a_to_b,
    foa_rotation_matrix,
    foa_to_pressure_velocity,
    open_sphere_foa_radial_weights,
    reorder_foa,
    rotate_foa,
    synthesize_tetrahedral_a_format,
    tetrahedral_a_to_b_matrix,
)
from codes.chapters.ch16.core.wearable_audio import (
    CHANNEL_COUNTS,
    FILE_NAMES,
    SAMPLE_RATE,
    SAMPLES,
    measure_signal,
    raised_cosine_envelope,
    run_experiment as run_audio_experiment,
)
from codes.chapters.ch16.core.wearable_fusion import (
    corcos_wind_coherence,
    cross_modal_mvdr_mmse,
    hear_through_response,
)
from codes.chapters.ch16.examples.generate_wearable_audio import (
    check_assets,
    generate_assets,
)


class TestSpatialParametric(unittest.TestCase):
    def test_tetrahedral_encoding_and_channel_reordering(self):
        # Regular tetrahedron geometry: U^T U = (4/3) I_3, sum_q u_q = 0.
        np.testing.assert_allclose(
            TETRA_DIRECTIONS.T @ TETRA_DIRECTIONS, (4.0 / 3.0) * np.eye(3), atol=1e-15
        )
        np.testing.assert_allclose(np.sum(TETRA_DIRECTIONS, axis=0), np.zeros(3), atol=1e-15)

        # Plane wave at azimuth=30 deg, elevation=45 deg with complex phasor 2 - 1j:
        az = np.pi / 6.0
        el = np.pi / 4.0
        amp = 2.0 - 1.0j
        nx = np.cos(el) * np.cos(az)  # sqrt(6)/4
        ny = np.cos(el) * np.sin(az)  # sqrt(2)/4
        nz = np.sin(el)  # sqrt(2)/2
        expected_wxyz = amp * np.array([1.0, nx, ny, nz])
        expected_acn = amp * np.array([1.0, ny, nz, nx])

        a_sig = synthesize_tetrahedral_a_format(az, el, amplitude=amp, cardioid_alpha=0.5)
        b_wxyz = encode_tetrahedral_a_to_b(a_sig, cardioid_alpha=0.5, order="WXYZ")
        b_acn = encode_tetrahedral_a_to_b(a_sig, cardioid_alpha=0.5, order="ACN")
        np.testing.assert_allclose(b_wxyz, expected_wxyz, atol=1e-14)
        np.testing.assert_allclose(b_acn, expected_acn, atol=1e-14)
        np.testing.assert_allclose(
            reorder_foa(b_wxyz, src_order="WXYZ", dst_order="ACN"), expected_acn, atol=1e-15
        )
        np.testing.assert_allclose(
            reorder_foa(b_acn, src_order="ACN", dst_order="WXYZ"), expected_wxyz, atol=1e-15
        )

        # Non-0.5 hypercardioid alpha = 0.75: g_W = 1/(4*0.25) = 1.0, g_1 = sqrt(3)/3 = 1/sqrt(3)
        enc_hyper = tetrahedral_a_to_b_matrix(cardioid_alpha=0.75, order="WXYZ")
        np.testing.assert_allclose(enc_hyper[0], np.ones(4), atol=1e-15)
        np.testing.assert_allclose(np.abs(enc_hyper[1:]), np.full((3, 4), 1.0 / np.sqrt(3.0)), atol=1e-15)
        a_hyper = synthesize_tetrahedral_a_format(az, el, amplitude=amp, cardioid_alpha=0.75)
        np.testing.assert_allclose(
            encode_tetrahedral_a_to_b(a_hyper, cardioid_alpha=0.75, order="WXYZ"),
            expected_wxyz,
            atol=1e-14,
        )

    def test_foa_rotation_scene_vs_listener_head(self):
        # Source at +x: [1, 1, 0, 0]. Scene yaw +90 deg moves it to +y: [1, 0, 1, 0].
        src_px = np.array([1.0, 1.0, 0.0, 0.0])
        np.testing.assert_allclose(
            rotate_foa(src_px, yaw_rad=np.pi / 2.0, listener_head_rotation=False),
            np.array([1.0, 0.0, 1.0, 0.0]),
            atol=1e-15,
        )
        # Listener turns head +90 deg left (+y): world +x source is now on listener's right (-y).
        np.testing.assert_allclose(
            rotate_foa(src_px, yaw_rad=np.pi / 2.0, listener_head_rotation=True),
            np.array([1.0, 0.0, -1.0, 0.0]),
            atol=1e-15,
        )
        # Orthogonality in ACN order as well: R_acn @ R_acn.T == I_4 and det == +1.
        r_acn = foa_rotation_matrix(0.4, -0.3, 0.7, order="ACN")
        np.testing.assert_allclose(r_acn @ r_acn.T, np.eye(4), atol=1e-15)
        self.assertAlmostEqual(float(np.linalg.det(r_acn)), 1.0, places=14)

    def test_ex16_1_active_intensity_and_diffuseness(self):
        # Case 1: single +x propagating plane wave (p=1, u=[1,0,0]) => psi = 0.
        single = active_intensity_and_diffuseness(1.0, np.array([1.0, 0.0, 0.0]), rho0=1.0, c=1.0)
        np.testing.assert_allclose(single["active_intensity"], [0.5, 0.0, 0.0], atol=1e-15)
        self.assertAlmostEqual(single["energy_density"], 0.5, places=15)
        self.assertAlmostEqual(single["diffuseness"], 0.0, places=15)
        np.testing.assert_allclose(single["intensity_unit"], [1.0, 0.0, 0.0], atol=1e-15)
        np.testing.assert_allclose(single["doa_unit"], [-1.0, 0.0, 0.0], atol=1e-15)

        # Case 2: two orthogonal plane waves in-phase (p=2, u=[1,1,0])
        # ||I_a|| = sqrt(2), E = 1.5, psi = 1 - 2*sqrt(2)/3
        in_phase = active_intensity_and_diffuseness(2.0, np.array([1.0, 1.0, 0.0]), rho0=1.0, c=1.0)
        np.testing.assert_allclose(in_phase["active_intensity"], [1.0, 1.0, 0.0], atol=1e-15)
        self.assertAlmostEqual(in_phase["energy_density"], 1.5, places=15)
        self.assertAlmostEqual(in_phase["diffuseness"], 1.0 - 2.0 * np.sqrt(2.0) / 3.0, places=14)
        np.testing.assert_allclose(
            in_phase["intensity_unit"], [1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0), 0.0], atol=1e-15
        )

        # Case 3: quadrature (Delta_phi = 90 deg: p = 1+j, u = [1, j, 0])
        # ||I_a|| = 1/sqrt(2), E = 1.0, psi = 1 - 1/sqrt(2)
        quad = active_intensity_and_diffuseness(
            1.0 + 1.0j, np.array([1.0 + 0.0j, 1.0j, 0.0 + 0.0j]), rho0=1.0, c=1.0
        )
        np.testing.assert_allclose(quad["active_intensity"], [0.5, 0.5, 0.0], atol=1e-15)
        np.testing.assert_allclose(quad["reactive_intensity"], [-0.5, 0.5, 0.0], atol=1e-15)
        self.assertAlmostEqual(quad["energy_density"], 1.0, places=15)
        self.assertAlmostEqual(quad["diffuseness"], 1.0 - 1.0 / np.sqrt(2.0), places=14)

    def test_foa_sign_and_scale_invariance_and_rejections(self):
        # Source at +x => FOA [W,X,Y,Z] = [2, 2, 0, 0] => physical u = [-2/(rho0*c), 0, 0]
        p, u = foa_to_pressure_velocity(np.array([2.0, 2.0, 0.0, 0.0]), rho0=1.2, c=343.0)
        res = active_intensity_and_diffuseness(p, u, rho0=1.2, c=343.0)
        np.testing.assert_allclose(res["doa_unit"], [1.0, 0.0, 0.0], atol=1e-15)
        self.assertAlmostEqual(res["diffuseness"], 0.0, places=15)

        # Zero-energy degenerate input returns valid_direction=False and finite zeros.
        zero_res = active_intensity_and_diffuseness(0.0, np.zeros(3))
        self.assertFalse(zero_res["valid_direction"])
        self.assertEqual(zero_res["diffuseness"], 1.0)

        # Scale invariance of diffuseness down to tiny underflow-prone amplitudes 1e-150:
        tiny = active_intensity_and_diffuseness(
            (1.0 + 1.0j) * 1e-150, np.array([1.0, 1.0j, 0.0]) * 1e-150
        )
        self.assertAlmostEqual(tiny["diffuseness"], 1.0 - 1.0 / np.sqrt(2.0), places=14)

        # Open-sphere Bessel null at kr = pi: b0(pi) = 0, b1(pi) = j / pi.
        b0, b1 = open_sphere_foa_radial_weights(np.array([0.0, np.pi]))
        np.testing.assert_allclose(b0, [1.0, 0.0], atol=1e-15)
        np.testing.assert_allclose(b1, [0.0j, 1.0j / np.pi], atol=1e-15)

        with self.assertRaises(ValueError):
            tetrahedral_a_to_b_matrix(cardioid_alpha=1.0)
        with self.assertRaises(ValueError):
            tetrahedral_a_to_b_matrix(cardioid_alpha=np.nan)
        with self.assertRaises(ValueError):
            reorder_foa(np.ones(3), src_order="WXYZ", dst_order="ACN")
        with self.assertRaises(ValueError):
            reorder_foa(np.ones(4), src_order="FUMA", dst_order="ACN")
        with self.assertRaises(ValueError):
            open_sphere_foa_radial_weights(np.array([-0.1, 1.0]))
        with self.assertRaises(ValueError):
            open_sphere_foa_radial_weights(np.array([np.inf]))
        with self.assertRaises(ValueError):
            active_intensity_and_diffuseness(1.0, np.array([1.0, np.nan, 0.0]))
        with self.assertRaises(ValueError):
            active_intensity_and_diffuseness(np.array([1.0, 2.0]), np.zeros((3, 3)))
        with self.assertRaises(ValueError):
            rotate_foa(np.array([1.0, 0.0, 0.0, 0.0]), yaw_rad=1.0 + 1.0j)  # type: ignore[arg-type]

        # Extreme kr = 1e150 must not overflow radial weights (b0 -> 0, b1 -> 0)
        b0_huge, b1_huge = open_sphere_foa_radial_weights(np.array([1e150]))
        self.assertTrue(np.all(np.isfinite(b0_huge)) and np.all(np.isfinite(b1_huge)))
        self.assertLess(float(np.abs(b1_huge[0])), 1e-149)


class TestBinauralMatching(unittest.TestCase):
    def test_ex16_2_bsm_ls_vs_magls_exact_anchor(self):
        h_mat = np.array([[1.0, 1.0, 1.0], [1.0j, 1.0, -1.0j]], dtype=np.complex128)
        p_ear = np.array([1.0, -1.0, 1.0], dtype=np.complex128)

        # Hand-derived: w_LS = [0.5, -0.5]^T => y_hat = [(1-j)/2, 0, (1+j)/2]^T
        w_ls = bsm_ls_weights(h_mat, p_ear)
        np.testing.assert_allclose(w_ls, [0.5, -0.5], atol=1e-15)
        m_ls = bsm_directional_metrics(h_mat, p_ear, w_ls)
        np.testing.assert_allclose(
            m_ls["synthesized_response"], [0.5 - 0.5j, 0.0 + 0.0j, 0.5 + 0.5j], atol=1e-15
        )
        np.testing.assert_allclose(
            m_ls["synthesized_magnitude"], [1.0 / np.sqrt(2.0), 0.0, 1.0 / np.sqrt(2.0)], atol=1e-15
        )
        self.assertAlmostEqual(m_ls["weighted_complex_mse"], 2.0 / 3.0, places=15)
        self.assertAlmostEqual(
            m_ls["weighted_magnitude_mse"], (4.0 - 2.0 * np.sqrt(2.0)) / 3.0, places=15
        )

        # Hand-derived 1-step MagLS from zero initial phases: w^(1) = [1, 0]^T => |y_hat| = [1, 1, 1]^T
        w_mag, hist = bsm_magls_weights(h_mat, np.abs(p_ear), n_iters=5)
        np.testing.assert_allclose(w_mag, [1.0, 0.0], atol=1e-15)
        np.testing.assert_allclose(hist, np.zeros(5), atol=1e-15)

    def test_magls_monotonic_non_increasing_and_singular_rejection(self):
        rng = np.random.default_rng(1602)
        h_mat = rng.normal(size=(4, 15)) + 1j * rng.normal(size=(4, 15))
        mag = rng.uniform(0.2, 2.0, size=15)
        q_weights = rng.uniform(0.5, 2.5, size=15)
        _, hist = bsm_magls_weights(h_mat, mag, weights=q_weights, reg=0.1, n_iters=30)
        self.assertEqual(hist.shape, (30,))
        self.assertTrue(np.all(np.diff(hist) <= 0.0))

        # Rank-deficient H with reg=0 must raise ValueError rather than silently failing.
        h_singular = np.ones((2, 3), dtype=np.complex128)
        with self.assertRaises(ValueError):
            bsm_ls_weights(h_singular, np.ones(3), reg=0.0)

    def test_ex16_3_bmvdr_coherence_collapse_and_mixing(self):
        r_nn = np.array([[1.0, 0.2], [0.2, 1.0]], dtype=np.complex128)
        w_bmvdr = np.array([0.5, 0.5], dtype=np.complex128)

        # Closed form: var(beta) = 0.6 + 0.4*beta^2, cross(beta) = 0.6 - 0.4*beta^2,
        # IC_out(beta) = (3 - 2*beta^2) / (3 + 2*beta^2)
        for beta, exp_var, exp_ic in (
            (0.0, 0.6, 1.0),
            (0.2, 0.616, 73.0 / 77.0),
            (0.5, 0.7, 5.0 / 7.0),
            (1.0, 1.0, 0.2),
        ):
            res = binaural_coherence_with_mixing(r_nn, w_bmvdr, w_bmvdr, beta=beta, ref_idx=(0, 1))
            self.assertAlmostEqual(res["var_L"], exp_var, places=15)
            self.assertAlmostEqual(res["var_R"], exp_var, places=15)
            self.assertAlmostEqual(res["ic_out"], exp_ic, places=14)

        # Scale invariance at extreme scales (1e200 and 1e-200) must never overflow/underflow IC!
        for scale in (1e200, 1e-200):
            res0 = binaural_coherence_with_mixing(
                scale * r_nn, w_bmvdr, w_bmvdr, beta=0.0, ref_idx=(0, 1)
            )
            self.assertAlmostEqual(res0["ic_out"], 1.0, places=14)
            res_half = binaural_coherence_with_mixing(
                scale * r_nn, w_bmvdr, w_bmvdr, beta=0.5, ref_idx=(0, 1)
            )
            self.assertAlmostEqual(res_half["ic_out"], 5.0 / 7.0, places=14)

        # Reject indefinite/non-Hermitian/NaN R_nn, invalid beta, invalid ref_idx, and MagLS overflow
        with self.assertRaises(ValueError):
            binaural_coherence_with_mixing(np.array([[1.0, 2.0], [2.0, 1.0]]), w_bmvdr, w_bmvdr)
        with self.assertRaises(ValueError):
            binaural_coherence_with_mixing(np.array([[1.0, 1.0j], [1.0j, 1.0]]), w_bmvdr, w_bmvdr)
        with self.assertRaises(ValueError):
            binaural_coherence_with_mixing(np.array([[np.nan, 0.0], [0.0, 1.0]]), w_bmvdr, w_bmvdr)
        with self.assertRaises(ValueError):
            binaural_coherence_with_mixing(r_nn, w_bmvdr, w_bmvdr, beta=1.2)
        with self.assertRaises(ValueError):
            binaural_coherence_with_mixing(r_nn, w_bmvdr, w_bmvdr, ref_idx=(True, 1))  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            bsm_magls_weights(np.eye(2), np.array([1e200, 1e200]), reg=0.0, n_iters=2)
        with self.assertRaises(ValueError):
            bsm_magls_weights(np.eye(2), np.array([1.0, -0.1]), reg=0.1, n_iters=2)
        with self.assertRaises(ValueError):
            bsm_magls_weights(np.eye(2), np.array([1.0, 1.0]), reg=0.1, n_iters=True)  # type: ignore[arg-type]


class TestWearableFusion(unittest.TestCase):
    def test_ex16_4_cross_modal_mvdr_mmse_exact_rationals(self):
        # Low band (250 Hz): c = [1, 2]^T, R_nn = diag(4.0, 0.04) => w_MVDR = [1/401, 200/401]^T, SINR = 100.25
        low = cross_modal_mvdr_mmse(2.0, np.diag([4.0, 0.04]), phi_s=1.0)
        np.testing.assert_allclose(low["w_mvdr"], [1.0 / 401.0, 200.0 / 401.0], atol=1e-15)
        self.assertAlmostEqual(low["output_sinr_linear"], 100.25, places=13)
        self.assertAlmostEqual(low["output_sinr_db"], 10.0 * np.log10(100.25), places=13)
        np.testing.assert_allclose(
            low["w_mmse"], (401.0 / 405.0) * np.array([1.0 / 401.0, 200.0 / 401.0]), atol=1e-15
        )
        self.assertAlmostEqual(low["mmse_total_mse"], 4.0 / 405.0, places=15)

        # High band (4000 Hz): c = [1, 0.1]^T, R_nn = diag(0.05, 0.04) => w_MVDR = [80/81, 10/81]^T, SINR = 20.25
        high = cross_modal_mvdr_mmse(0.1, np.diag([0.05, 0.04]), phi_s=1.0)
        np.testing.assert_allclose(high["w_mvdr"], [80.0 / 81.0, 10.0 / 81.0], atol=1e-15)
        self.assertAlmostEqual(high["output_sinr_linear"], 20.25, places=13)
        self.assertAlmostEqual(high["output_sinr_db"], 10.0 * np.log10(20.25), places=13)
        # Naive 10x inverse equalization amplifies 0.04 sensor floor 100x to 4.0 (+20 dB)
        self.assertAlmostEqual(high["equalized_inear_noise_variance"], 4.0, places=13)

        # Scale invariance at 1e200 when R_nn and phi_s scale together, and overflow rejection
        scaled = cross_modal_mvdr_mmse(2.0, 1e200 * np.diag([4.0, 0.04]), phi_s=1e200)
        np.testing.assert_allclose(scaled["w_mvdr"], [1.0 / 401.0, 200.0 / 401.0], atol=1e-15)
        self.assertAlmostEqual(scaled["output_sinr_linear"], 100.25, places=12)
        with self.assertRaises(ValueError):
            cross_modal_mvdr_mmse(2.0, 1e-200 * np.eye(2), phi_s=1e200)
        with self.assertRaises(ValueError):
            cross_modal_mvdr_mmse(np.nan, np.eye(2), phi_s=1.0)
        with self.assertRaises(ValueError):
            hear_through_response(np.array([250.0]), tau_proc_s=-1e-3)
        with self.assertRaises(ValueError):
            hear_through_response(np.array([250.0]), tau_proc_s=0.002, g_leak=1e200, g_ht=1e200)
        with self.assertRaises(ValueError):
            corcos_wind_coherence(np.array([250.0]), d_m=0.015, u_conv_mps=0.0)

    def test_hear_through_comb_and_corcos_wind(self):
        freqs = np.array([0.0, 250.0, 500.0, 750.0, 4000.0])
        ht = hear_through_response(freqs, tau_proc_s=0.002, g_leak=1.0, g_ht=1.0)
        np.testing.assert_allclose(ht["power_gain"], [4.0, 0.0, 4.0, 0.0, 4.0], atol=1e-14)
        expected_notches = np.array([250.0, 750.0, 1250.0, 1750.0, 2250.0, 2750.0, 3250.0, 3750.0])
        np.testing.assert_allclose(ht["comb_notch_frequencies_hz"], expected_notches, atol=1e-14)

        coh = corcos_wind_coherence(
            np.array([0.0, 250.0, 1000.0]),
            d_m=0.015,
            u_conv_mps=3.0,
            alpha_long=0.14,
            sound_speed=343.0,
        )
        np.testing.assert_allclose(
            coh["wind_coherence"],
            [1.0, np.exp(-0.14 * 2.5 * np.pi), np.exp(-0.14 * 10.0 * np.pi)],
            atol=1e-15,
        )
        expected_sinc_250 = np.sin(2.0 * np.pi * 250.0 * 0.015 / 343.0) / (
            2.0 * np.pi * 250.0 * 0.015 / 343.0
        )
        self.assertAlmostEqual(float(coh["acoustic_diffuse_coherence"][1]), expected_sinc_250, places=15)


class TestAmbisonicAndSofaRendering(unittest.TestCase):
    def test_e16_17_foa_max_re_vs_basic_and_in_phase_decoding(self):
        # Closed-form modal weights: Chebyshev T_2(g1) = 0 -> 1/sqrt(2), Legendre P_2(g1) = 0 -> 1/sqrt(3)
        self.assertEqual(foa_max_re_order_weight(2, "basic"), 1.0)
        self.assertAlmostEqual(foa_max_re_order_weight(2, "max_re"), 1.0 / np.sqrt(2.0), places=15)
        self.assertEqual(foa_max_re_order_weight(2, "in_phase"), 0.5)
        self.assertEqual(foa_max_re_order_weight(3, "basic"), 1.0)
        self.assertAlmostEqual(foa_max_re_order_weight(3, "max_re"), 1.0 / np.sqrt(3.0), places=15)
        self.assertAlmostEqual(foa_max_re_order_weight(3, "in_phase"), 1.0 / 3.0, places=15)

        src_px = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        # 2-D regular hexagon (L = 6): Basic r_E = 2/3, max-r_E r_E = 1/sqrt(2), in-phase r_E = 2/3
        hex_basic = foa_loudspeaker_decode_metrics(REGULAR_HEXAGON_2D_DIRS, src_px, spatial_dim=2, g1=1.0)
        hex_maxre = foa_loudspeaker_decode_metrics(
            REGULAR_HEXAGON_2D_DIRS, src_px, spatial_dim=2, g1=1.0 / np.sqrt(2.0)
        )
        hex_inph = foa_loudspeaker_decode_metrics(REGULAR_HEXAGON_2D_DIRS, src_px, spatial_dim=2, g1=0.5)
        self.assertTrue(hex_basic["is_isotropic_layout"])
        self.assertAlmostEqual(hex_basic["r_V_norm"], 1.0, places=15)
        self.assertAlmostEqual(hex_basic["r_E_norm"], 2.0 / 3.0, places=15)
        self.assertAlmostEqual(hex_basic["min_gain"], -1.0 / 6.0, places=15)
        self.assertTrue(hex_basic["has_negative_lobe"])

        self.assertAlmostEqual(hex_maxre["r_V_norm"], 1.0 / np.sqrt(2.0), places=15)
        self.assertAlmostEqual(hex_maxre["r_E_norm"], 1.0 / np.sqrt(2.0), places=15)
        self.assertAlmostEqual(hex_maxre["r_E_acos_deg"], 45.0, places=13)
        self.assertAlmostEqual(hex_maxre["min_gain"], (1.0 - np.sqrt(2.0)) / 6.0, places=15)
        # Rear speaker l = 3 (-x) power suppression relative to Basic (-1/6): (sqrt(2) - 1)^2
        rear_supp_2d_db = 10.0 * np.log10((np.sqrt(2.0) - 1.0) ** 2)
        self.assertAlmostEqual(
            10.0 * np.log10((hex_maxre["gains"][3] / hex_basic["gains"][3]) ** 2),
            rear_supp_2d_db,
            places=13,
        )
        self.assertFalse(hex_inph["has_negative_lobe"])
        self.assertAlmostEqual(hex_inph["min_gain"], 0.0, places=15)

        # 3-D regular octahedron (L = 6): Basic r_E = 1/2, max-r_E r_E = 1/sqrt(3), in-phase min_gain = 0
        oct_basic = foa_loudspeaker_decode_metrics(REGULAR_OCTAHEDRON_3D_DIRS, src_px, spatial_dim=3, g1=1.0)
        oct_maxre = foa_loudspeaker_decode_metrics(
            REGULAR_OCTAHEDRON_3D_DIRS, src_px, spatial_dim=3, g1=1.0 / np.sqrt(3.0)
        )
        oct_inph = foa_loudspeaker_decode_metrics(
            REGULAR_OCTAHEDRON_3D_DIRS, src_px, spatial_dim=3, g1=1.0 / 3.0
        )
        self.assertTrue(oct_basic["is_isotropic_layout"])
        self.assertAlmostEqual(oct_basic["r_V_norm"], 1.0, places=15)
        self.assertAlmostEqual(oct_basic["r_E_norm"], 0.5, places=15)
        self.assertAlmostEqual(oct_basic["min_gain"], -1.0 / 3.0, places=15)
        self.assertAlmostEqual(oct_maxre["r_E_norm"], 1.0 / np.sqrt(3.0), places=15)
        # Rear speaker l = 1 (-x): gain (1 - sqrt(3))/6 vs Basic -2/6 => ratio (sqrt(3) - 1)/2
        rear_supp_3d_db = 10.0 * np.log10(((np.sqrt(3.0) - 1.0) / 2.0) ** 2)
        self.assertAlmostEqual(
            10.0 * np.log10((oct_maxre["gains"][1] / oct_basic["gains"][1]) ** 2),
            rear_supp_3d_db,
            places=13,
        )
        self.assertFalse(oct_inph["has_negative_lobe"])
        self.assertAlmostEqual(oct_inph["min_gain"], 0.0, places=15)

    def test_e16_18_sofa_geometry_and_head_tracking_latency(self):
        # Exact SOFA spherical [30 deg, 45 deg, 1.5 m] -> [3*sqrt(6)/8, 3*sqrt(2)/8, 3*sqrt(2)/4]
        expected_cart = np.array(
            [3.0 * np.sqrt(6.0) / 8.0, 3.0 * np.sqrt(2.0) / 8.0, 3.0 * np.sqrt(2.0) / 4.0],
            dtype=np.float64,
        )
        actual_cart = sofa_spherical_to_cartesian(np.array([30.0, 45.0, 1.5]))
        np.testing.assert_allclose(actual_cart, expected_cart, atol=1e-15)
        np.testing.assert_allclose(sofa_cartesian_to_spherical(actual_cart), [30.0, 45.0, 1.5], atol=1e-13)

        valid_sph = np.array([[90.0, 0.0, 0.09], [-90.0, 0.0, 0.09]], dtype=np.float64)
        ok = audit_sofa_simple_freefield_geometry([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], valid_sph, coord_type="spherical")
        self.assertTrue(ok["passed"])
        np.testing.assert_allclose(ok["ear_left_projection_m"], [0.09, -0.09], atol=1e-15)

        # Bug 1: swapped left/right ears
        swapped = audit_sofa_simple_freefield_geometry(
            [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [[0.0, -0.09, 0.0], [0.0, 0.09, 0.0]], coord_type="cartesian"
        )
        self.assertFalse(swapped["passed"])
        self.assertIn("ear_sign_swap_or_nonpositive_left_ear", swapped["issues"])
        with self.assertRaises(ValueError):
            audit_sofa_simple_freefield_geometry(
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [[0.0, -0.09, 0.0], [0.0, 0.09, 0.0]],
                coord_type="cartesian",
                strict=True,
            )

        # Bug 2: spherical degrees parsed as Cartesian metres
        unconverted = audit_sofa_simple_freefield_geometry(
            [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], valid_sph, coord_type="cartesian"
        )
        self.assertFalse(unconverted["passed"])
        self.assertTrue(unconverted["likely_unconverted_spherical_degrees"])
        self.assertIn("implausible_head_radius_m", unconverted["issues"])

        # Head-tracking motion-to-sound latency at 180 deg/s: 20 ms -> 3.6 deg, 80 ms -> 14.4 deg
        lat20 = head_tracking_latency_audit(180.0, 0.020, maa_deg=2.0)
        self.assertAlmostEqual(lat20["angular_lag_deg"], 3.6, places=14)
        self.assertAlmostEqual(lat20["foa_chord_error"], 2.0 * np.sin(np.deg2rad(1.8)), places=15)
        self.assertTrue(lat20["exceeds_maa"])

        lat80 = head_tracking_latency_audit(180.0, 0.080, maa_deg=2.0, prediction_horizon_s=0.065)
        self.assertAlmostEqual(lat80["angular_lag_deg"], 14.4, places=13)
        self.assertAlmostEqual(lat80["compensated_lag_deg"], 2.7, places=13)
        self.assertAlmostEqual(lat80["abrupt_stop_overshoot_deg"], 11.7, places=13)

        # Defensive boundary rejections
        with self.assertRaises(ValueError):
            foa_max_re_order_weight(True, "max_re")  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            foa_loudspeaker_decode_metrics(REGULAR_OCTAHEDRON_3D_DIRS, [1.0, 0.0, 0.0], spatial_dim=2)
        with self.assertRaises(ValueError):
            sofa_spherical_to_cartesian([0.0, 95.0, 1.0])
        with self.assertRaises(ValueError):
            head_tracking_latency_audit(180.0, -0.01)


class TestWearableAudioAndExercises(unittest.TestCase):
    def test_audio_generation_and_strict_readonly_check(self):
        report, signals = run_audio_experiment()
        self.assertEqual(SAMPLE_RATE, 16000)
        self.assertEqual(SAMPLES, 32000)
        self.assertEqual(set(signals), set(FILE_NAMES))
        for k, arr in signals.items():
            self.assertEqual(arr.shape, (CHANNEL_COUNTS[k], 32000))
            self.assertLess(float(np.max(np.abs(arr))), 0.90)

        # Verify fused output achieves exact MVDR noise reduction on [1600, 30400)
        expected_fused_mse = 0.005 / 100.25 + 0.005 / 20.25
        self.assertAlmostEqual(
            report["float_measurements"]["wearable_fused_output"]["reference_mse_ch0"],
            expected_fused_mse,
            places=15,
        )

        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "wearable_audio"
            # Missing directory must fail --check without creating anything
            with self.assertRaises(ValueError):
                check_assets(out_dir)
            self.assertFalse(out_dir.exists())

            manifest = generate_assets(out_dir, check=False)
            members = {p.name for p in out_dir.iterdir()}
            self.assertEqual(len(members), 7)
            self.assertEqual(members, set(FILE_NAMES.values()) | {"MANIFEST.json"})

            checked = generate_assets(out_dir, check=True)
            self.assertEqual(checked, manifest)

            # Extra member must cause check_assets to fail
            extra = out_dir / "extra.txt"
            extra.write_text("bad", encoding="utf-8")
            with self.assertRaises(ValueError):
                check_assets(out_dir)
            extra.unlink()

            # Symlink member must cause check_assets to fail
            wav_path = out_dir / "wearable_clean_ref.wav"
            backup = out_dir.parent / "backup.wav"
            wav_path.replace(backup)
            os.symlink(backup, wav_path)
            with self.assertRaises(ValueError):
                check_assets(out_dir)

    def test_chapter16_exercises_cli_and_catalog(self):
        with self.assertRaises(ValueError):
            raised_cosine_envelope(n_samples=100, fade_samples=60)
        with self.assertRaises(ValueError):
            measure_signal(np.zeros((1, SAMPLES)), "wearable_clean_ref", reference=np.zeros((1, SAMPLES)), pcm=1)  # type: ignore[arg-type]

        payload = run_experiments()
        self.assertEqual(tuple(payload["exercises"].keys()), EXERCISE_IDS)
        self.assertEqual(len(EXERCISE_IDS), 18)
        ex = payload["exercises"]

        # E16-03: exact Frobenius norm of Rz(pi/2)Ry(pi/2) - Ry(pi/2)Rz(pi/2) is sqrt(6)
        self.assertAlmostEqual(ex["E16-03"]["commutator_frobenius_norm"], float(np.sqrt(6.0)), places=14)

        # E16-09: hand-derived frontal weighting W_Omega = diag(1, 4, 1)
        f0 = ex["E16-09"]["frontal4_reg0"]
        np.testing.assert_allclose(f0["weights"], [[0.2, 0.0], [-0.8, 0.0]], atol=1e-14)
        self.assertAlmostEqual(f0["weight_norm_squared"], 0.68, places=14)
        self.assertAlmostEqual(f0["broadside_magnitude"], 0.6, places=14)
        f2 = ex["E16-09"]["frontal4_reg2.0"]
        np.testing.assert_allclose(f2["weights"], [[0.0, 0.0], [-0.5, 0.0]], atol=1e-14)
        self.assertAlmostEqual(f2["weight_norm_squared"], 0.25, places=14)
        self.assertAlmostEqual(f2["broadside_magnitude"], 0.5, places=14)

        # E16-12: hand-derived 250 Hz wind exp(-0.35*pi) and 2000 Hz diffuse sinc(60/343)
        self.assertAlmostEqual(ex["E16-12"]["wind_coherence"][1], float(np.exp(-0.35 * np.pi)), places=14)
        arg_2k = 60.0 * np.pi / 343.0
        self.assertAlmostEqual(
            ex["E16-12"]["acoustic_diffuse_coherence"][4], float(np.sin(arg_2k) / arg_2k), places=14
        )

        # E16-16: causal fade-in tail difference on wearable_hearthrough_comb
        ht_ana = ex["E16-16"]["analytic_expectations"]["wearable_hearthrough_comb"]["reference_mse_ch0"]
        ht_flt = ex["E16-16"]["float_measurements"]["wearable_hearthrough_comb"]["reference_mse_ch0"]
        ht_pcm = ex["E16-16"]["pcm_measurements"]["wearable_hearthrough_comb"]["reference_mse_ch0"]
        self.assertAlmostEqual(ht_ana, 0.01, places=15)
        self.assertGreater(ht_ana - ht_flt, 7.0e-9)
        self.assertLess(ht_ana - ht_flt, 8.5e-9)
        self.assertAlmostEqual(ht_pcm, 0.010001495537419235, places=12)

        # E16-17 & E16-18 exact anchors
        self.assertAlmostEqual(
            ex["E16-17"]["octahedron_3d"]["modes"]["max_re"]["r_E_norm"], 1.0 / np.sqrt(3.0), places=14
        )
        self.assertTrue(ex["E16-18"]["sofa_valid_audit"]["passed"])
        self.assertFalse(ex["E16-18"]["sofa_ear_swap_audit"]["passed"])
        self.assertAlmostEqual(ex["E16-18"]["latency_20ms"]["angular_lag_deg"], 3.6, places=14)

        # Ensure JSON serializability without NaN/Inf
        encoded = json.dumps(payload, allow_nan=False)
        self.assertIn("E16-18", encoded)
        self.assertEqual(cli_main(["--exercise", "E16-17"]), 0)


if __name__ == "__main__":
    unittest.main()
