"""Independent hand-derived unit tests for Chapter 17 kernels, exercises, and audio assets."""

from __future__ import annotations

import io
import json
import math
from pathlib import Path
import tempfile
import unittest
import wave

import numpy as np

from codes.chapters.ch17 import chapter17_exercises as ex17
from codes.chapters.ch17.core import (
    av_alignment,
    conditional_tse,
    diarization_tsvad,
    distance_bubble,
    eend_enrollment_robustness,
    tse_audio,
)
from codes.chapters.ch17.examples import generate_tse_audio as tse_gen


class ConditionalTSETests(unittest.TestCase):
    def test_film_conditioning_hand_example_and_3d_broadcast(self):
        features = np.array([[1.0, -0.5, 2.0], [0.5, 1.5, -1.0]], dtype=np.float64)
        expected = np.array([[1.5, -1.5, 3.5], [0.5, 1.0, -0.25]], dtype=np.float64)
        actual = conditional_tse.film_conditioning(features, [2.0, 0.5], [-0.5, 0.25])
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-15)

        # Neutral (1, 0) identity and (0, 0) zero gate
        np.testing.assert_array_equal(
            conditional_tse.film_conditioning(features, [1.0, 1.0], [0.0, 0.0]), features
        )
        np.testing.assert_array_equal(
            conditional_tse.film_conditioning(features, [0.0, 0.0], [0.0, 0.0]),
            np.zeros_like(features),
        )

        # 3-D (C, F, T) feature map broadcast along channel axis 0
        feat_3d = np.stack([features, 2.0 * features], axis=1)
        out_3d = conditional_tse.film_conditioning(feat_3d, [2.0, 0.5], [-0.5, 0.25])
        np.testing.assert_allclose(out_3d[:, 0, :], expected, rtol=0.0, atol=1e-15)
        self.assertEqual(out_3d.shape, (2, 2, 3))

    def test_film_rejects_complex_and_shape_mismatch(self):
        h = np.ones((2, 3), dtype=np.float64)
        with self.assertRaises(ValueError):
            conditional_tse.film_conditioning(h + 0j, [1.0, 1.0], [0.0, 0.0])
        with self.assertRaises(ValueError):
            conditional_tse.film_conditioning(h, [1.0], [0.0, 0.0])

    def test_dual_state_loss_present_scale_invariance_and_absent_singularity(self):
        s_ref = np.array([1.0, 1.0, -1.0, -1.0], dtype=np.float64)
        v_int = np.array([1.0, -1.0, 1.0, -1.0], dtype=np.float64)
        x_mix = s_ref + v_int

        for scale in (0.25, 0.5, 1.0, 3.0):
            s_hat = scale * (s_ref + 0.1 * v_int)
            res = conditional_tse.tse_dual_state_loss(s_hat, s_ref, x_mix)
            self.assertTrue(res["target_present"])
            self.assertFalse(res["standard_sisdr_undefined"])
            self.assertAlmostEqual(res["si_sdr_db"], 20.0, places=12)
            self.assertAlmostEqual(res["loss_db"], -20.0, places=12)
            self.assertAlmostEqual(res["scale_alpha"], scale, places=12)

        # Target-absent: standard SI-SDR is 0/0 undefined, dual-state uses energy floor
        s_zero = np.zeros(4, dtype=np.float64)
        diag = conditional_tse.inspect_standard_sisdr(0.2 * v_int, s_zero)
        self.assertTrue(diag["projection_undefined"])
        self.assertIsNone(diag["standard_sisdr_db"])

        quiet = conditional_tse.tse_dual_state_loss(np.zeros(4), s_zero, v_int, tau_floor=1e-3)
        self.assertFalse(quiet["target_present"])
        self.assertTrue(quiet["standard_sisdr_undefined"])
        self.assertAlmostEqual(quiet["floor_energy"], 0.004, places=15)
        self.assertAlmostEqual(quiet["loss_db"], 0.0, places=15)

        # ||s_hat||^2 = 0.396 = 99 * 0.004 -> 10 * log10(1 + 99) = +20.0 dB exact
        leak = conditional_tse.tse_dual_state_loss(
            math.sqrt(0.099) * v_int, s_zero, v_int, tau_floor=1e-3
        )
        self.assertAlmostEqual(leak["loss_db"], 20.0, places=12)

        # Full interferer passthrough: 10 * log10(1 + 1000) = 10 * log10(1001)
        full_pass = conditional_tse.tse_dual_state_loss(v_int, s_zero, v_int, tau_floor=1e-3)
        self.assertAlmostEqual(full_pass["loss_db"], 10.0 * math.log10(1001.0), places=12)

        with self.assertRaises(ValueError):
            conditional_tse.tse_dual_state_loss(np.zeros(4), s_zero, np.zeros(4))

    def test_zero_overlap_bypass_fricative_recovery_and_scale_invariance(self):
        y_pri = np.array([1.0, 0.0, -1.0, 0.0], dtype=np.float64)
        y_fric = np.array([0.0, 0.1, 0.0, -0.1], dtype=np.float64)
        y_ovl = np.array([0.0, 0.5, 0.0, -0.5], dtype=np.float64)

        for gain in (1e-3, 1.0, 1e3):
            b1 = conditional_tse.zero_overlap_bypass_decision(gain * y_pri, gain * y_fric, -18.0)
            self.assertTrue(b1["bypass_separation"])
            self.assertFalse(b1["secondary_active"])
            self.assertAlmostEqual(b1["energy_ratio_db"], -20.0, places=12)
            np.testing.assert_allclose(b1["output"], gain * (y_pri + y_fric), rtol=0.0, atol=1e-15)

            b2 = conditional_tse.zero_overlap_bypass_decision(gain * y_pri, gain * y_ovl, -18.0)
            self.assertFalse(b2["bypass_separation"])
            self.assertTrue(b2["secondary_active"])
            self.assertAlmostEqual(b2["energy_ratio_db"], 10.0 * math.log10(0.25), places=12)
            np.testing.assert_allclose(b2["output"], gain * y_pri, rtol=0.0, atol=1e-15)


class DistanceBubbleTests(unittest.TestCase):
    def test_fresnel_distance_and_closed_form_bubble_metrics(self):
        # c = 343 m/s, f = 3430 Hz -> lambda = 0.10 m exact
        self.assertAlmostEqual(distance_bubble.fresnel_distance(0.20, 3430.0, 343.0), 0.80, places=14)
        self.assertAlmostEqual(distance_bubble.fresnel_distance(0.04, 3430.0, 343.0), 0.032, places=14)

        # Hand-derived closed form for d = 0.10 m, r1 = 0.40 m, r2 = 4.00 m
        delta_r1 = math.hypot(0.40, 0.10) - 0.40
        phi1 = -2.0 * math.pi * delta_r1 / 0.10
        expected_wng_10 = (8.0 / 3.0) * (math.sin(0.5 * phi1) ** 2)
        delta_r2 = math.hypot(4.00, 0.10) - 4.00
        phi2 = -2.0 * math.pi * delta_r2 / 0.10
        expected_b2_r2 = (math.sin(0.5 * phi2) / math.sin(0.5 * phi1)) ** 2
        w_out_expected = 1.0 / (2.0 * (complex(math.cos(-phi1), math.sin(-phi1)) - 1.0))
        expected_w10 = np.array([w_out_expected, -2.0 * w_out_expected, w_out_expected])

        res10 = distance_bubble.broadside_distance_bubble_weights(
            0.10, 0.40, 3430.0, c=343.0, radii_m=np.array([0.40, 4.00])
        )
        np.testing.assert_allclose(res10["weights"], expected_w10, rtol=0.0, atol=1e-13)
        self.assertAlmostEqual(res10["extra_path_m"] * 1000.0, 12.310562561766059, places=10)
        self.assertAlmostEqual(res10["outer_phase_deg"], math.degrees(phi1), places=11)
        self.assertAlmostEqual(res10["wng_db"], 10.0 * math.log10(expected_wng_10), places=11)
        self.assertAlmostEqual(res10["wng_db"], -4.2094, places=3)
        self.assertAlmostEqual(-10.0 * math.log10(res10["radial_power_response"][1]), -10.0 * math.log10(expected_b2_r2), places=10)
        self.assertAlmostEqual(-10.0 * math.log10(res10["radial_power_response"][1]), 19.6533, places=3)
        self.assertLess(res10["farfield_plane_power_gain"], 1e-28)

        # Tiny d = 0.02 m collapse check
        delta_r_tiny = math.hypot(0.40, 0.02) - 0.40
        phi_tiny = -2.0 * math.pi * delta_r_tiny / 0.10
        expected_wng_02 = 10.0 * math.log10((8.0 / 3.0) * (math.sin(0.5 * phi_tiny) ** 2))
        res02 = distance_bubble.broadside_distance_bubble_weights(0.02, 0.40, 3430.0, c=343.0)
        self.assertAlmostEqual(res02["wng_db"], expected_wng_02, places=10)
        self.assertAlmostEqual(res02["wng_db"], -31.8237, places=3)

    def test_spherical_steering_modes_and_two_tap_fir(self):
        pos = np.array([-0.10, 0.0, 0.10], dtype=np.float64)
        a_ph = distance_bubble.spherical_ula_steering(pos, 0.40, 0.0, 3430.0, attenuation="phase_only")
        a_rel = distance_bubble.spherical_ula_steering(pos, 0.40, 0.0, 3430.0, attenuation="relative_1_over_r")
        a_abs = distance_bubble.spherical_ula_steering(pos, 0.40, 0.0, 3430.0, attenuation="absolute_1_over_r")
        np.testing.assert_allclose(np.abs(a_ph), [1.0, 1.0, 1.0], rtol=0.0, atol=1e-15)
        ratio = 0.40 / math.hypot(0.40, 0.10)
        np.testing.assert_allclose(np.abs(a_rel), [ratio, 1.0, ratio], rtol=0.0, atol=1e-15)
        np.testing.assert_allclose(np.abs(a_abs), [1.0 / math.hypot(0.40, 0.10), 2.5, 1.0 / math.hypot(0.40, 0.10)], rtol=0.0, atol=1e-15)

        res10 = distance_bubble.broadside_distance_bubble_weights(0.10, 0.40, 3200.0, c=343.0)
        taps = distance_bubble.broadside_two_tap_real_fir(res10["weights"], 3200.0, 16000.0)
        omega = 2.0 * math.pi * 3200.0 / 16000.0
        freq_resp = taps[:, 0] + taps[:, 1] * np.exp(-1j * omega)
        np.testing.assert_allclose(freq_resp, np.conj(res10["weights"]), rtol=0.0, atol=1e-14)
        # Broadband broadside plane-wave null: column sums vanish
        np.testing.assert_allclose(np.sum(taps, axis=0), [0.0, 0.0], rtol=0.0, atol=1e-14)


class DiarizationTSVADTests(unittest.TestCase):
    def setUp(self):
        self.ref = np.array([[1, 1, 1, 1, 0, 0], [0, 0, 1, 1, 1, 0]], dtype=np.int64)

    def test_der_decomposition_overlap_bound_and_over_100_percent(self):
        bound = diarization_tsvad.single_label_overlap_der_bound(self.ref, frame_dur_s=0.5)
        self.assertEqual(bound["overlap_frames"], 2)
        self.assertAlmostEqual(bound["min_miss_s"], 1.0, places=15)
        self.assertAlmostEqual(bound["ref_speaker_time_s"], 3.5, places=15)
        self.assertAlmostEqual(bound["der_lower_bound"], 2.0 / 7.0, places=15)

        hyp_imp = np.array([[1, 1, 1, 0, 1, 0], [0, 0, 0, 1, 0, 1]], dtype=np.int64)
        res = diarization_tsvad.compute_der_components(self.ref, hyp_imp, frame_dur_s=0.5)
        self.assertAlmostEqual(res["miss_s"], 1.0, places=15)
        self.assertAlmostEqual(res["false_alarm_s"], 0.5, places=15)
        self.assertAlmostEqual(res["speaker_confusion_s"], 0.5, places=15)
        self.assertAlmostEqual(res["der"], 4.0 / 7.0, places=15)
        self.assertFalse(res["exceeds_unit_der"])

        # Label-swapped permutation invariance
        res_swap = diarization_tsvad.compute_der_components(self.ref, hyp_imp[::-1, :], frame_dur_s=0.5)
        self.assertAlmostEqual(res_swap["der"], 4.0 / 7.0, places=15)

        # Legitimate DER > 100%
        hyp_over = np.array([[0, 0, 0, 0, 1, 1], [0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 1]], dtype=np.int64)
        res_over = diarization_tsvad.compute_der_components(self.ref, hyp_over, frame_dur_s=0.5)
        self.assertAlmostEqual(res_over["miss_s"], 1.5, places=15)
        self.assertAlmostEqual(res_over["false_alarm_s"], 2.5, places=15)
        self.assertAlmostEqual(res_over["speaker_confusion_s"], 0.0, places=15)
        self.assertAlmostEqual(res_over["der"], 8.0 / 7.0, places=15)
        self.assertTrue(res_over["exceeds_unit_der"])

        # Large-K (K = 8 > 7) Hungarian path check
        eye8 = np.eye(8, dtype=np.int64)
        perm8 = np.roll(eye8, shift=3, axis=0)
        res8 = diarization_tsvad.compute_der_components(eye8, perm8, frame_dur_s=0.25)
        self.assertAlmostEqual(res8["der"], 0.0, places=15)

    def test_tsvad_posterior_and_gss_mvdr_distortionless(self):
        # Hand check: prior=0.5, c_sp=0.80, c_emb=0.47 -> logit = 3*(0.55) + 4*(0.02) = 1.73
        p = diarization_tsvad.tsvad_spatial_spectral_posterior(
            0.80, 0.47, prior=0.5, kappa_sp=3.0, tau_sp=0.25, kappa_emb=4.0, tau_emb=0.45
        )
        self.assertAlmostEqual(float(p), 1.0 / (1.0 + math.exp(-1.73)), places=14)

        # Rank-1 target check: Souden MVDR preserves exact reference mic response a_1[ref_mic]
        a1 = np.array([1.0 + 0.0j, np.exp(-1j * 0.9)], dtype=np.complex128)
        a2 = np.array([1.0 + 0.0j, np.exp(+1j * 1.1)], dtype=np.complex128)
        idx = np.arange(40, dtype=np.float64)
        # First 20 frames pure speaker 1, next 20 frames pure speaker 2
        s1 = (idx < 20).astype(np.float64) * np.exp(1j * 0.3 * idx)
        s2 = (idx >= 20).astype(np.float64) * np.exp(1j * 0.7 * idx)
        x_obs = np.outer(a1, s1) + np.outer(a2, s2)
        gamma = np.vstack([(idx < 20).astype(np.float64), (idx >= 20).astype(np.float64)])
        out = diarization_tsvad.gss_soft_gated_mvdr(x_obs, gamma, target_index=0, diag_loading=1e-6, ref_mic=0)
        self.assertAlmostEqual(abs(np.vdot(out["weights"], a1) - a1[0]), 0.0, places=10)
        self.assertLess(abs(np.vdot(out["weights"], a2)), 1e-5)


class AVAlignmentTests(unittest.TestCase):
    def test_cadence_lag_onset_clipping_and_confidence_gate(self):
        sync = av_alignment.map_video_to_audio_frames(15, 0.016, 25.0, av_offset_s=0.0)
        self.assertEqual(
            sync["video_frame_indices"].tolist(),
            [0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5],
        )
        self.assertEqual(sync["onset_audit"]["clipped_consonant_frames"], 0)

        lagged = av_alignment.map_video_to_audio_frames(
            30, 0.016, 25.0, av_offset_s=0.120, onset_frame=10, consonant_frames=8
        )
        self.assertEqual(lagged["onset_audit"]["clipped_consonant_frames"], 7)
        self.assertAlmostEqual(lagged["onset_audit"]["consonant_loss_fraction"], 0.875, places=15)

        m_aud = np.array([0.82, 0.82, 0.82], dtype=np.float64)
        m_av = np.array([0.96, 0.65, 0.05], dtype=np.float64)
        alpha = np.array([0.90, 0.40, 0.05], dtype=np.float64)
        fused = av_alignment.av_confidence_gated_mask(m_aud, m_av, alpha)
        np.testing.assert_allclose(fused, [0.946, 0.752, 0.7815], rtol=0.0, atol=1e-15)

    def test_av_boundary_validation_and_2d_broadcast(self):
        m2d_aud = np.full((3, 4), 0.8, dtype=np.float64)
        m2d_av = np.full((3, 4), 0.2, dtype=np.float64)
        alpha_t = np.array([0.0, 0.25, 0.5, 1.0], dtype=np.float64)
        out_2d = av_alignment.av_confidence_gated_mask(m2d_aud, m2d_av, alpha_t)
        np.testing.assert_allclose(out_2d[0], [0.8, 0.65, 0.5, 0.2], rtol=0.0, atol=1e-15)

        for bad_kwargs in (
            {"n_audio_frames": True, "audio_hop_s": 0.016, "video_fps": 25.0},
            {"n_audio_frames": np.bool_(True), "audio_hop_s": 0.016, "video_fps": 25.0},
            {"n_audio_frames": 10, "audio_hop_s": 0.016, "video_fps": 25.0, "onset_frame": -1},
            {"n_audio_frames": 10, "audio_hop_s": 0.016, "video_fps": 25.0, "consonant_frames": 0},
            {"n_audio_frames": 10, "audio_hop_s": 0.0, "video_fps": 25.0},
            {"n_audio_frames": 10, "audio_hop_s": 1.0, "video_fps": 1e20},
        ):
            with self.assertRaises(ValueError):
                av_alignment.map_video_to_audio_frames(**bad_kwargs)

        with self.assertRaises(ValueError):
            av_alignment.av_confidence_gated_mask(m2d_aud, m2d_av, 1.2)
        with self.assertRaises(ValueError):
            av_alignment.av_confidence_gated_mask(m2d_aud, m2d_av, True)


class EendAndEnrollmentRobustnessTests(unittest.TestCase):
    def test_e17_17_eend_eda_variable_stopping_vs_fixed_slot_pit_bound(self):
        w_exist = np.array([1.0, 1.0, 1.0], dtype=np.float64)
        emb_3spk = np.array(
            [
                [2.0, 2.0, 2.0, -2.0],
                [-2.0, 2.0, 2.0, -2.0],
                [-2.0, -2.0, 2.0, 2.0],
            ],
            dtype=np.float64,
        )
        attr_2spk = np.array([[2.0, 0.0, 0.1], [0.0, 2.0, 0.1], [0.0, 0.0, 0.0]], dtype=np.float64)
        stop2 = eend_enrollment_robustness.eend_eda_inference(emb_3spk, attr_2spk, w_exist, b_exist=-1.0)
        # Logits: [2-1, 2-1, 0.2-1] = [1.0, 1.0, -0.8] -> stops at s_hat = 2
        np.testing.assert_allclose(stop2["existence_logits"], [1.0, 1.0, -0.8], atol=1e-15)
        self.assertEqual(stop2["estimated_num_speakers"], 2)

        attr_3spk = np.array(
            [[2.0, 0.0, 0.0, 0.1], [0.0, 2.0, 0.0, 0.1], [0.0, 0.0, 2.0, 0.1]],
            dtype=np.float64,
        )
        stop3 = eend_enrollment_robustness.eend_eda_inference(emb_3spk, attr_3spk, w_exist, b_exist=-1.0)
        self.assertEqual(stop3["estimated_num_speakers"], 3)
        ref_3spk = np.array([[1, 1, 1, 0], [0, 1, 1, 0], [0, 0, 1, 1]], dtype=np.int64)
        np.testing.assert_array_equal(stop3["binary_activity"], ref_3spk)
        der_eda = diarization_tsvad.compute_der_components(ref_3spk, stop3["binary_activity"], 0.5)
        self.assertAlmostEqual(der_eda["der"], 0.0, places=15)

        # Fixed 2-slot PIT lower bound on 3-speaker scene [3, 2, 2] frames (total 7 frames): 2/7
        pit_bound = eend_enrollment_robustness.fixed_slot_der_lower_bound(ref_3spk, max_slots=2, frame_dur_s=0.5)
        self.assertEqual(pit_bound["min_unmatched_frames"], 2)
        self.assertEqual(pit_bound["ref_speaker_frames"], 7)
        self.assertAlmostEqual(pit_bound["der_lower_bound"], 2.0 / 7.0, places=15)
        der_oracle = diarization_tsvad.compute_der_components(ref_3spk, pit_bound["oracle_slot_hypothesis"], 0.5)
        self.assertAlmostEqual(der_oracle["der"], 2.0 / 7.0, places=15)

    def test_e17_18_enrollment_contamination_and_diarization_gated_pooling(self):
        e_tar = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        e_int_orth = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        orth = eend_enrollment_robustness.contaminated_enrollment_metrics(e_tar, e_int_orth, leakage_alpha=0.5)
        # Orthogonal: e_cont = [1, 0.5, 0] / sqrt(1.25) => cos_tar = 2/sqrt(5), cos_int = 1/sqrt(5)
        self.assertAlmostEqual(orth["target_cosine"], 2.0 / math.sqrt(5.0), places=15)
        self.assertAlmostEqual(orth["interferer_cosine"], 1.0 / math.sqrt(5.0), places=15)
        self.assertAlmostEqual(orth["clean_interferer_gate"], 1.0 / (1.0 + math.exp(3.2)), places=15)

        # Correlated rho_12 = 0.3, alpha = 0.5: ||e_1 + 0.5 e_2|| = sqrt(1 + 0.3 + 0.25) = sqrt(1.55)
        e_int_corr = np.array([0.3, math.sqrt(0.91), 0.0], dtype=np.float64)
        corr = eend_enrollment_robustness.contaminated_enrollment_metrics(e_tar, e_int_corr, leakage_alpha=0.5)
        self.assertAlmostEqual(corr["target_cosine"], 1.15 / math.sqrt(1.55), places=15)
        self.assertAlmostEqual(corr["interferer_cosine"], 0.80 / math.sqrt(1.55), places=15)

        # 6-frame enrollment buffer where frames 4, 5 are high-SNR interferer crosstalk during target pause
        enr_frames = np.column_stack([e_tar, e_tar, e_tar, e_tar, e_int_orth, e_int_orth])
        enr_snr = np.array([20.0, 20.0, 20.0, 20.0, 15.0, 15.0], dtype=np.float64)
        enr_gamma = np.array(
            [[0.95, 0.95, 0.95, 0.95, 0.05, 0.05], [0.04, 0.04, 0.04, 0.04, 0.92, 0.92]],
            dtype=np.float64,
        )
        pooled = eend_enrollment_robustness.gated_enrollment_pooling(
            enr_frames, enr_snr, enr_gamma, target_index=0, min_snr_db=10.0, min_exclusive_prob=0.6
        )
        self.assertEqual(pooled["vad_selected_indices"], [0, 1, 2, 3, 4, 5])
        self.assertEqual(pooled["joint_selected_indices"], [0, 1, 2, 3])
        self.assertEqual(pooled["rejected_crosstalk_indices"], [4, 5])
        # VAD-only mean is (4*e_tar + 2*e_int)/6 -> normalized [2/sqrt(5), 1/sqrt(5), 0]
        np.testing.assert_allclose(
            pooled["vad_only_embedding"], [2.0 / math.sqrt(5.0), 1.0 / math.sqrt(5.0), 0.0], atol=1e-15
        )
        np.testing.assert_allclose(pooled["gated_embedding"], [1.0, 0.0, 0.0], atol=1e-15)

        # Negative controls
        with self.assertRaises(ValueError):
            eend_enrollment_robustness.fixed_slot_der_lower_bound(np.ones((2, 3), dtype=int), max_slots=True)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            eend_enrollment_robustness.contaminated_enrollment_metrics([2.0, 0.0], [0.0, 1.0])
        with self.assertRaises(ValueError):
            eend_enrollment_robustness.gated_enrollment_pooling(enr_frames, enr_snr, enr_gamma, min_snr_db=50.0)


class DefensiveContractsTests(unittest.TestCase):
    def test_extreme_scale_invariance_and_nonfinite_rejection(self):
        s_ref = np.array([1.0, 1.0, -1.0, -1.0], dtype=np.float64)
        v_int = np.array([1.0, -1.0, 1.0, -1.0], dtype=np.float64)
        y_pri = np.array([1.0, 0.0, -1.0, 0.0], dtype=np.float64)
        y_fric = np.array([0.0, 0.1, 0.0, -0.1], dtype=np.float64)

        for scale in (1e-200, 1e-80, 1e80, 1e140):
            diag = conditional_tse.inspect_standard_sisdr(
                scale * (s_ref + 0.1 * v_int), scale * s_ref
            )
            self.assertFalse(diag["projection_undefined"])
            self.assertAlmostEqual(diag["standard_sisdr_db"], 20.0, places=12)

            byp = conditional_tse.zero_overlap_bypass_decision(
                scale * y_pri, scale * y_fric, threshold_db=-18.0
            )
            self.assertTrue(byp["bypass_separation"])
            self.assertAlmostEqual(byp["energy_ratio_db"], -20.0, places=12)
            np.testing.assert_allclose(
                byp["output"], scale * (y_pri + y_fric), rtol=1e-13, atol=0.0
            )

        # GSS Souden MVDR weights must remain invariant even at 1e-200 and 1e140 amplitudes
        a1 = np.array([1.0 + 0.0j, np.exp(-1j * 0.8)], dtype=np.complex128)
        a2 = np.array([1.0 + 0.0j, np.exp(+1j * 1.2)], dtype=np.complex128)
        idx = np.arange(20, dtype=np.float64)
        s1 = (idx < 10).astype(np.float64) * np.exp(1j * 0.25 * idx)
        s2 = (idx >= 10).astype(np.float64) * np.exp(1j * 0.65 * idx)
        x0 = np.outer(a1, s1) + np.outer(a2, s2)
        gamma = np.vstack([(idx < 10).astype(np.float64), (idx >= 10).astype(np.float64)])
        base = diarization_tsvad.gss_soft_gated_mvdr(x0, gamma, target_index=0, diag_loading=1e-4)
        for scale in (1e-200, 1e140):
            scaled = diarization_tsvad.gss_soft_gated_mvdr(
                scale * x0, gamma, target_index=0, diag_loading=1e-4
            )
            np.testing.assert_allclose(scaled["weights"], base["weights"], rtol=0.0, atol=1e-12)
            np.testing.assert_allclose(
                scaled["enhanced_frames"], scale * base["enhanced_frames"], rtol=1e-12, atol=1e-13 * scale
            )

        # Bool, NaN/Inf, and domain/overflow rejections
        with self.assertRaises(ValueError):
            conditional_tse.film_conditioning(np.full((2, 2), 10.0), [1e308, 1.0], [0.0, 0.0])
        with self.assertRaises(ValueError):
            conditional_tse.tse_dual_state_loss(s_ref, s_ref, s_ref, tau_floor=True)
        with self.assertRaises(ValueError):
            distance_bubble.fresnel_distance(1e200, 1e200, 343.0)
        with self.assertRaises(ValueError):
            distance_bubble.spherical_ula_steering([0.0], 0.4, 0.0, 1e308, c=1.0)
        with self.assertRaises(ValueError):
            distance_bubble.broadside_two_tap_real_fir([1.0 + 1.0j], 1e-16, 16000.0)
        with self.assertRaises(ValueError):
            diarization_tsvad.compute_der_components([[True, False]], [[1, 0]], 0.5)
        with self.assertRaises(ValueError):
            diarization_tsvad.gss_soft_gated_mvdr(x0, gamma, target_index=True)
        with self.assertRaises(ValueError):
            diarization_tsvad.tsvad_spatial_spectral_posterior(1.2, 0.5)
        with self.assertRaises(ValueError):
            tse_audio.measure_signal(np.zeros((1, 32000)), "near_target_ref", pcm="yes")  # type: ignore[arg-type]


class TSEAudioAndExercisesTests(unittest.TestCase):
    def test_exercise_runner_covers_all_18_ids_and_hand_values(self):
        payload = ex17.run_experiments()
        json_text = json.dumps(payload, allow_nan=False)
        self.assertTrue(len(json_text) > 100)
        ex = payload["exercises"]
        self.assertEqual(tuple(ex.keys()), ex17.EXERCISE_IDS)
        self.assertEqual(len(ex), 18)

        self.assertEqual(ex["E17-01"]["enrolled_modulated"], [[1.5, -1.5, 3.5], [0.5, 1.0, -0.25]])
        self.assertAlmostEqual(ex["E17-02"]["target_present"]["loss_db"], -20.0, places=12)
        self.assertAlmostEqual(ex["E17-02"]["target_absent_moderate_leak"]["loss_db"], 20.0, places=12)
        self.assertAlmostEqual(
            ex["E17-04"]["d_10cm"]["exact_extra_path_mm"],
            (math.sqrt(0.17) - 0.40) * 1000.0,
            places=12,
        )
        self.assertAlmostEqual(
            ex["E17-06"]["outer_amplitude_ratio_r0_over_r1"], 4.0 / math.sqrt(17.0), places=14
        )
        self.assertAlmostEqual(ex["E17-08"]["imperfect_2cluster"]["der"], 4.0 / 7.0, places=14)
        self.assertTrue(ex["E17-08"]["permutation_invariant_der"])
        self.assertAlmostEqual(ex["E17-09"]["theoretical_bound"]["der_lower_bound"], 2.0 / 7.0, places=14)
        self.assertAlmostEqual(ex["E17-10"]["over_triggering_3track"]["der"], 8.0 / 7.0, places=14)
        self.assertAlmostEqual(ex["E17-11"]["joint_der"], 0.0, places=14)
        self.assertEqual(ex["E17-13"]["audio_frames_per_video_exposure"], [2, 3, 2, 3, 2, 3])
        self.assertEqual(
            [row["clipped_consonant_frames"] for row in ex["E17-14"]["camera_lag_consonant_clipping_sweep"]],
            [0, 2, 5, 7],
        )
        np.testing.assert_allclose(ex["E17-15"]["mask_fused"], [0.946, 0.752, 0.7815], rtol=0.0, atol=1e-14)
        self.assertLess(ex["E17-16"]["pcm"]["bubble_d10cm"]["integer_reference_nmse"], 0.025)
        self.assertEqual(ex["E17-17"]["two_speaker_stopping_check"]["estimated_num_speakers"], 2)
        self.assertAlmostEqual(ex["E17-17"]["three_speaker_eend_eda"]["der"], 0.0, places=14)
        self.assertAlmostEqual(ex["E17-17"]["fixed_2slot_oracle_der"], 2.0 / 7.0, places=14)
        self.assertAlmostEqual(
            ex["E17-18"]["orthogonal_speakers_alpha0p5"]["target_cosine"], 2.0 / math.sqrt(5.0), places=14
        )
        self.assertAlmostEqual(ex["E17-18"]["six_frame_enrollment_remedy"]["gated_target_cosine"], 1.0, places=14)

    def test_generate_and_strict_check_tse_audio_in_tempdir(self):
        with self.assertRaises(ValueError):
            tse_gen.generate_assets(tse_gen.DEFAULT_OUTPUT, check="yes")  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            tse_gen.check_assets(Path("codes/chapters/../tse_audio"))

        manifest_mem, blobs = tse_gen.prepare_assets()
        self.assertEqual(set(blobs.keys()), set(tse_gen.MEMBERS))
        self.assertEqual(len(tse_gen.MEMBERS), 7)
        self.assertEqual(manifest_mem["common_export_gain"], 1.0)

        # Verify WAV headers and PCM headroom
        for fname, ch_count in (
            ("tse_near_target_ref.wav", 1),
            ("tse_far_interferer_ref.wav", 1),
            ("tse_array_3ch.wav", 3),
            ("tse_farfield_dsb.wav", 1),
            ("tse_bubble_d10cm.wav", 1),
            ("tse_bubble_d2cm_illcond.wav", 1),
        ):
            with wave.open(io.BytesIO(blobs[fname]), "rb") as reader:
                self.assertEqual(
                    (reader.getframerate(), reader.getnchannels(), reader.getnframes(), reader.getsampwidth()),
                    (16000, ch_count, 32000, 2),
                )

        # Check physics: 10 cm bubble massively outperforms far-field DSB and ill-conditioned 2 cm bubble
        nmse_dsb = manifest_mem["samples"]["farfield_dsb"]["float_measurements"]["reference_nmse"]
        nmse_b10 = manifest_mem["samples"]["bubble_d10cm"]["float_measurements"]["reference_nmse"]
        nmse_b02 = manifest_mem["samples"]["bubble_d2cm_illcond"]["float_measurements"]["reference_nmse"]
        self.assertGreater(nmse_dsb, 1.2)
        self.assertLess(nmse_b10, 0.03)
        self.assertGreater(nmse_b02, 2.0)

        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "tse_audio"
            with self.assertRaises(ValueError):
                tse_gen.generate_assets(out_dir, check=True)
            generated = tse_gen.generate_assets(out_dir, check=False)
            checked = tse_gen.generate_assets(out_dir, check=True)
            self.assertEqual(generated["files"], checked["files"])

            # Extra unexpected file must be rejected by strict --check
            (out_dir / "extra.txt").write_text("tamper", encoding="utf-8")
            with self.assertRaises(ValueError):
                tse_gen.generate_assets(out_dir, check=True)


if __name__ == "__main__":
    unittest.main()
