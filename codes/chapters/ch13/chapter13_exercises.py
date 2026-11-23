"""Run E13-01..21 with NumPy and in-memory PCM; never write or download.

Usage: python -m codes.chapters.ch13.chapter13_exercises [--exercise E13-01]
Metadata is outside the stable exercise-ID map. Known-model controls are not
the MATLAB DANSE chain, a theorem proof, or a measured distributed system.
"""
from __future__ import annotations

if __name__ == '__main__' and not __package__:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import argparse
import json
import math
import numpy as np
from codes.chapters.ch00.core.audio_samples import pcm16_bytes, read_pcm16
from codes.chapters.ch04.core.covariance import _load_covariance
from codes.chapters.ch13.core.distributed import (
    mwf_weights, mse_components, compressed_mwf, distributed_updates,
    run_covariance_experiments, jacobi_control, step_size_control, tree_sum_control, gevd_control, broadcast_statistics_control,
)
from codes.chapters.ch13.core.distributed_audio import (
    SAMPLE_RATE, POWER, run_experiment, measure_signal, OUTPUT_REFERENCES,
)
from codes.chapters.ch13.core.dasf_and_ptp import (
    ptp_vs_sro_phase_budget, dasf_unification_control,
)

EXERCISE_IDS = tuple('E13-'+str(i).zfill(2) for i in range(1, 22))


def _json(value):
    if isinstance(value, np.ndarray):
        return np.stack([value.real, value.imag], axis=-1).tolist() if np.iscomplexobj(value) else value.tolist()
    if isinstance(value, (complex, np.complexfloating)):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(item) for item in value]
    return value


def resource_control():
    fs = 16000; bits = 16; channels = 2; compressed = 1
    network = []
    for nodes in (2, 4, 8):
        total = nodes*channels; d = channels+(nodes-1)*compressed
        network.append({'nodes': nodes, 'channels_per_node': channels, 'compressed_streams_per_node': compressed,
                        'central_upload_PCM16_bits_per_second': total*fs*bits,
                        'compressed_broadcast_PCM16_bits_per_second': nodes*compressed*fs*bits,
                        'compressed_full_mesh_unicast_PCM16_bits_per_second': nodes*(nodes-1)*compressed*fs*bits,
                        'central_matrix_dimension': total, 'node_matrix_dimension': d,
                        'central_full_complex128_covariance_bytes': total**2*16,
                        'node_full_complex128_covariance_bytes': d**2*16,
                        'network_full_complex128_covariance_bytes': nodes*d**2*16,
                        'central_cubic_dimension_proxy': total**3,
                        'network_cubic_dimension_proxy': nodes*d**3})
    return {'sample_rate_hz': fs, 'network': network,
            'two_nodes_PCM16_raw_upload_bps': 4*fs*16,
            'two_nodes_PCM16_scalar_exchange_bps': 2*fs*16,
            'two_nodes_float32_scalar_exchange_bps': 2*fs*32,
            'two_nodes_STFT_complex_float32_exchange_bps': 2*257*64*(fs//128),
            'STFT': {'n_fft': 512, 'hop': 128, 'bins': 257, 'bits_per_complex_value': 64,
                     'new_audio_seconds_per_call': 128/fs},
            'central_second_node_output_return_PCM16_bps': fs*16,
            'scope': 'payload only; excludes framing, timestamps, headers, loss recovery and radio scheduling; cubic dimensions are not measured flops or runtime'}


def dwacd_sro_gating_control(seed: int = 13021, trials: int = 20):
    """DWACD-inspired two-window cross-spectrum slope and known activity gating control (E13-21)."""
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError('seed must be an integer')
    if isinstance(trials, bool) or not isinstance(trials, int) or trials < 1:
        raise ValueError('trials must be a positive integer')

    fs = 16000
    n_fft = 512
    hop = 256
    eps_true_ppm = 100.0
    eps_true = eps_true_ppm * 1e-6
    lag_b = 200
    seg_len = 100
    k_min, k_max = 16, 128
    bins = np.arange(k_min, k_max + 1, dtype=np.float64)
    n_bins = int(len(bins))
    tau0_samples = 1.75
    speech_power = 100.0
    active_noise_power = 1.0  # +20 dB SNR during active speech

    delta_n_true = eps_true * lag_b * hop  # 5.12 samples
    slope_true = (2.0 * np.pi / n_fft) * delta_n_true  # 0.02 * pi rad/bin
    omega = (2.0 * np.pi / n_fft) * bins
    true_unit = np.exp(1j * omega * delta_n_true)

    tau_max = 15.0
    coarse_grid = np.linspace(-tau_max, tau_max, 1201)
    coarse_basis = np.exp(-1j * np.outer(coarse_grid, omega))
    phi = (math.sqrt(5.0) - 1.0) / 2.0

    def _wacd_peak_search(cross_prod):
        obj = np.real(coarse_basis @ cross_prod)
        idx = int(np.argmax(obj))
        lo = float(coarse_grid[max(0, idx - 1)])
        hi = float(coarse_grid[min(len(coarse_grid) - 1, idx + 1)])

        def _score(tau):
            return float(np.real(np.dot(cross_prod, np.exp(-1j * omega * tau))))

        c = hi - phi * (hi - lo)
        d = lo + phi * (hi - lo)
        fc, fd = _score(c), _score(d)
        for _ in range(65):
            if fc < fd:
                lo, c, fc = c, d, fd
                d = lo + phi * (hi - lo)
                fd = _score(d)
            else:
                hi, d, fd = d, c, fc
                c = hi - phi * (hi - lo)
                fc = _score(c)
        tau_hat = 0.5 * (lo + hi)
        eps_hat_ppm = (tau_hat / (lag_b * hop)) * 1e6
        return tau_hat, eps_hat_ppm

    def _run_case(mask1, mask2, inactive_noise_power, case_seed):
        rng = np.random.default_rng(case_seed)
        joint_mask = mask1 & mask2
        n_act1 = int(np.sum(mask1))
        n_act2 = int(np.sum(mask2))
        n_joint = int(np.sum(joint_mask))
        j_joint = np.where(joint_mask)[0].astype(np.float64)
        j_act1 = np.where(mask1)[0].astype(np.float64)
        j_act2 = np.where(mask2)[0].astype(np.float64)
        centroid_lag = float((lag_b + np.mean(j_act2)) - np.mean(j_act1))
        expected_ungated_centroid_ppm = float(eps_true_ppm * (centroid_lag / lag_b))

        rho_joint_sq = np.abs(np.mean(np.exp(1j * np.outer(omega, eps_true * hop * j_joint)), axis=1)) ** 2
        rho_act1_sq = np.abs(np.mean(np.exp(1j * np.outer(omega, eps_true * hop * j_act1)), axis=1)) ** 2

        snr_gamma_gated_lin = float(np.mean(
            n_joint * (speech_power ** 2) * rho_joint_sq / ((speech_power + active_noise_power) ** 2)
        ))
        n_inact1 = seg_len - n_act1
        snr_gamma_ungated_lin = float(np.mean(
            (n_act1 ** 2) * (speech_power ** 2) * rho_act1_sq
            / (n_act1 * ((speech_power + active_noise_power) ** 2) + n_inact1 * (inactive_noise_power ** 2))
        ))

        def _synth(ell0, mask):
            x1 = np.empty((n_bins, seg_len), dtype=np.complex128)
            x2 = np.empty((n_bins, seg_len), dtype=np.complex128)
            ell = ell0 + np.arange(seg_len, dtype=np.float64)
            idx_a = np.where(mask)[0]
            if len(idx_a):
                s = math.sqrt(speech_power / 2.0) * (
                    rng.standard_normal((n_bins, len(idx_a))) + 1j * rng.standard_normal((n_bins, len(idx_a)))
                )
                n1 = math.sqrt(active_noise_power / 2.0) * (
                    rng.standard_normal((n_bins, len(idx_a))) + 1j * rng.standard_normal((n_bins, len(idx_a)))
                )
                n2 = math.sqrt(active_noise_power / 2.0) * (
                    rng.standard_normal((n_bins, len(idx_a))) + 1j * rng.standard_normal((n_bins, len(idx_a)))
                )
                drift = tau0_samples + eps_true * hop * ell[idx_a]
                phase = np.exp(-1j * np.outer(omega, drift))
                x1[:, idx_a] = s + n1
                x2[:, idx_a] = s * phase + n2
            idx_i = np.where(~mask)[0]
            if len(idx_i):
                scale = math.sqrt(inactive_noise_power / 2.0)
                x1[:, idx_i] = scale * (
                    rng.standard_normal((n_bins, len(idx_i))) + 1j * rng.standard_normal((n_bins, len(idx_i)))
                )
                x2[:, idx_i] = scale * (
                    rng.standard_normal((n_bins, len(idx_i))) + 1j * rng.standard_normal((n_bins, len(idx_i)))
                )
            return x1, x2

        gated_errs = []
        ungated_errs = []
        first_detail = None
        for t in range(trials):
            x1_1, x2_1 = _synth(0, mask1)
            x1_2, x2_2 = _synth(lag_b, mask2)

            g1_g = np.mean(x1_1[:, joint_mask] * np.conj(x2_1[:, joint_mask]), axis=1)
            g2_g = np.mean(x1_2[:, joint_mask] * np.conj(x2_2[:, joint_mask]), axis=1)
            p_g = g2_g * np.conj(g1_g)
            tau_g, eps_g = _wacd_peak_search(p_g)
            msc_g = float(np.mean(
                (np.abs(g1_g) ** 2) / (
                    np.mean(np.abs(x1_1[:, joint_mask]) ** 2, axis=1)
                    * np.mean(np.abs(x2_1[:, joint_mask]) ** 2, axis=1)
                )
            ))
            ph_rmse_g = float(np.sqrt(np.mean(np.angle(p_g * np.conj(true_unit)) ** 2)))

            g1_u = np.mean(x1_1 * np.conj(x2_1), axis=1)
            g2_u = np.mean(x1_2 * np.conj(x2_2), axis=1)
            p_u = g2_u * np.conj(g1_u)
            tau_u, eps_u = _wacd_peak_search(p_u)
            msc_u = float(np.mean(
                (np.abs(g1_u) ** 2) / (
                    np.mean(np.abs(x1_1) ** 2, axis=1)
                    * np.mean(np.abs(x2_1) ** 2, axis=1)
                )
            ))
            ph_rmse_u = float(np.sqrt(np.mean(np.angle(p_u * np.conj(true_unit)) ** 2)))

            gated_errs.append(eps_g - eps_true_ppm)
            ungated_errs.append(eps_u - eps_true_ppm)
            if t == 0:
                _, eps_rev = _wacd_peak_search(np.conj(p_g))
                first_detail = {
                    'gated': {
                        'delta_n_hat_samples': float(tau_g),
                        'eps_hat_ppm': float(eps_g),
                        'error_ppm': float(eps_g - eps_true_ppm),
                        'mean_msc': msc_g,
                        'phase_rmse_rad': ph_rmse_g,
                    },
                    'ungated': {
                        'delta_n_hat_samples': float(tau_u),
                        'eps_hat_ppm': float(eps_u),
                        'error_ppm': float(eps_u - eps_true_ppm),
                        'mean_msc': msc_u,
                        'phase_rmse_rad': ph_rmse_u,
                    },
                    'reversed_conjugation_eps_hat_ppm': float(eps_rev),
                }
        return {
            'active_frames_seg1': n_act1,
            'active_frames_seg2': n_act2,
            'jointly_active_pairs': n_joint,
            'inactive_noise_power': float(inactive_noise_power),
            'active_centroid_lag_frames': centroid_lag,
            'expected_ungated_centroid_ppm': expected_ungated_centroid_ppm,
            'analytic_cross_psd_snr_gated_db': float(10.0 * math.log10(snr_gamma_gated_lin)),
            'analytic_cross_psd_snr_ungated_db': float(10.0 * math.log10(snr_gamma_ungated_lin)),
            'single_trial': first_detail,
            'trials': trials,
            'gated_rmse_ppm': float(np.sqrt(np.mean(np.square(gated_errs)))),
            'ungated_rmse_ppm': float(np.sqrt(np.mean(np.square(ungated_errs)))),
            'ungated_mean_bias_ppm': float(np.mean(ungated_errs)),
        }

    mask_60_head = np.zeros(seg_len, dtype=bool)
    mask_60_head[:60] = True
    mask_60_tail = np.zeros(seg_len, dtype=bool)
    mask_60_tail[40:] = True
    mask_20_head = np.zeros(seg_len, dtype=bool)
    mask_20_head[:20] = True

    return {
        'parameters': {
            'fs_hz': fs,
            'n_fft': n_fft,
            'hop_samples': hop,
            'hop_ms': float(1000.0 * hop / fs),
            'eps_true_ppm': eps_true_ppm,
            'lag_frames_B': lag_b,
            'lag_seconds': float(lag_b * hop / fs),
            'segment_frames_Lw': seg_len,
            'k_min': k_min,
            'k_max': k_max,
            'n_bins': n_bins,
            'speech_snr_db': 20.0,
            'tau0_samples': tau0_samples,
            'seed': seed,
        },
        'analytic': {
            'true_drift_samples': float(delta_n_true),
            'true_phase_slope_rad_per_bin': float(slope_true),
            'phase_at_k_min_rad': float(omega[0] * delta_n_true),
            'phase_at_k_max_rad': float(omega[-1] * delta_n_true),
            'wraps_past_pi_at_k_max': bool(omega[-1] * delta_n_true > math.pi),
            'unambiguous_wacd_drift_limit_samples': float(n_fft / 2.0),
            'intra_frame_drift_samples': float(eps_true * n_fft),
        },
        'symmetric_60pct_speech_40pct_quiet_silence': _run_case(mask_60_head, mask_60_head, 1.0, seed),
        'asymmetric_60pct_speech_40pct_quiet_silence': _run_case(mask_60_head, mask_60_tail, 1.0, seed + 1),
        'symmetric_60pct_speech_40pct_local_noise_bursts': _run_case(mask_60_head, mask_60_head, 2500.0, seed + 2),
        'sparse_20pct_speech_80pct_local_noise_bursts': _run_case(mask_20_head, mask_20_head, 1000.0, seed + 3),
        'scope': 'synthetic complex-STFT two-window cross-spectrum slope control under known SRO; no original coherence normalization or recursive state; not a measured room recording or full streaming DWACD tracker',
    }


def run_experiments():
    covariance = run_covariance_experiments(); model = covariance['model']
    rs, rn, a = model['Rs'], model['correlated'], model['a']
    correlated = covariance['cases']['correlated']
    report, signals = run_experiment()
    decoded = {key: read_pcm16(pcm16_bytes(value, SAMPLE_RATE))[1] for key, value in signals.items()}
    pcm = {key: measure_signal(decoded[key], key, reference=decoded[ref], pcm=True)
           for key, ref in OUTPUT_REFERENCES.items()}
    raw_projection = np.array([[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 2., -.5]])
    complex_a = np.array([1., .5j, 2., -.5]); complex_rs = np.outer(complex_a, complex_a.conj())
    complex_w = mwf_weights(complex_rs, np.eye(4), 1)
    scale_controls = []
    for scale in (1e-300, 1., 1e300):
        projection = raw_projection.copy(); projection[2] *= scale
        result = compressed_mwf(rs, rn, projection)
        scale_controls.append({'compression_scale': scale, 'effective_weights': result['weights'],
                               'normalized_condition_number': result['condition_number']})
    raw_conditions = []
    for scale in (1e-10, 1., 1e10):
        t = raw_projection.copy(); t[2] *= scale
        raw_conditions.append({'scale': scale, 'unnormalized_covariance_condition': float(np.linalg.cond(t@(rs+rn)@t.T))})
    steps = distributed_updates(rs, rn, model['nodes'], model['references'], [a[:2], a[2:]], max_updates=6, tolerance=0.)
    multi_a = np.array([[1., 0.], [0., 1.], [1., 1.], [1., -1.]])
    multi_rs = multi_a@multi_a.T; multi_rn = np.eye(4)
    multi_full = np.column_stack([mwf_weights(multi_rs, multi_rn, r) for r in (0, 1)])
    multi_t = [[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 1.]]
    multi_compressed = [compressed_mwf(multi_rs, multi_rn, multi_t, r) for r in (0, 1)]
    zero = distributed_updates(rs, model['white'], model['nodes'], model['references'],
                               [np.zeros(2), np.zeros(2)], zero_policy='local', max_updates=1, tolerance=0.)
    rejected = {}
    for name, target, noise in [('indefinite_noise', np.eye(2), [[1., 2.], [2., 1.]]),
                                ('ill_conditioned', np.diag([1., 0.]), np.diag([1., 1e-14]))]:
        try:
            mwf_weights(target, noise)
        except (ValueError, np.linalg.LinAlgError) as error:
            rejected[name] = {'rejected': True, 'diagnostic': str(error)}
        else:
            raise AssertionError('expected explicit covariance rejection')
    ill_rs = np.diag([1., 0.]); ill_rn = np.diag([1., 1e-14]); ill_rx = ill_rs+ill_rn
    loaded = _load_covariance(ill_rx, .001, 1e12)
    load = float((loaded-ill_rx)[0, 0].real)
    loaded_w = np.linalg.solve(loaded, ill_rs[:, 0])
    exercises = {
        'E13-01': {'a': complex_a, 'reference': 1, 'cross_covariance': complex_rs[:, 1], 'weights': complex_w,
                   'response': np.vdot(complex_w, complex_a), 'desired_response': complex_a[1],
                   'components': mse_components(complex_rs, np.eye(4), complex_w, 1),
                   'node2_weights': correlated['node2_weights'], 'node2_components': correlated['node2_components']},
        'E13-02': {'normalized_scale_controls': scale_controls, 'raw_conditions': raw_conditions,
                   'scope': 'row normalization preserves a nonzero compression subspace; it cannot recover a zero or lost direction'},
        'E13-03': steps,
        'E13-04': covariance['round_robin'],
        'E13-05': {'simultaneous': covariance['simultaneous'], 'simultaneous_half': covariance['simultaneous_half'],
                   'scope': 'finite known-covariance control; relaxation applies only to broadcast local weights, with solved output proposals recorded unchanged'},
        'E13-06': {'A': multi_a, 'full_weights': multi_full, 'remote_weight_rank': int(np.linalg.matrix_rank(multi_full[2:])),
                   'compressed': multi_compressed, 'full_components': [mse_components(multi_rs, multi_rn, multi_full[:, r], r) for r in (0, 1)]},
        'E13-07': {'zero_compression_local_policy': zero, 'default_policy': 'reject zero peer compression',
                   'zero_reference_target_power': mse_components(np.diag([1., 0.]), np.eye(2), np.zeros(2), 1)},
        'E13-08': {'rejections': rejected, 'ill_Rs': ill_rs, 'ill_Rn': ill_rn,
                   'loaded_covariance': loaded, 'additive_load': load, 'loaded_weights': loaded_w,
                   'physical_components': mse_components(ill_rs, ill_rn, loaded_w),
                   'extra_objective_penalty': float(load*np.vdot(loaded_w, loaded_w).real),
                   'scope': 'loading changes the objective; invalid original noise covariance is rejected before loading'},
        'E13-09': resource_control(),
        'E13-10': {'resources': resource_control()['network'],
                   'observed_solve_counts': {key: covariance[key]['total_solve_count'] for key in ('round_robin', 'simultaneous', 'simultaneous_half')},
                   'new_audio_seconds_per_STFT_call': 128/16000,
                   'scope': 'matrix bytes and count proxies; no measured hardware runtime or power'},
        'E13-11': {'parameters': report['parameters']['clock'], 'clock': report['clock'],
                   'synchronous_float': report['float_components']['central_white'],
                   'misaligned_float': report['float_components']['clock_misaligned_white'],
                   'linear_corrected_float': report['float_components']['clock_linear_corrected_white'],
                   'PCM': {key: pcm[key] for key in ('clock_misaligned_white', 'clock_linear_corrected_white')},
                   'warning': 'interpolation attenuates frequencies and changes noise covariance; lower total MSE is not recovery of the original synchronous optimum'},
        'E13-12': {'packet': report['parameters']['packet'],
                   'analytic': {key: report['analytic_expected'][key] for key in ('packet_zerofill_white', 'packet_local_fallback_white')},
                   'float': {key: report['float_components'][key] for key in ('packet_zerofill_white', 'packet_local_fallback_white')},
                   'PCM': {key: pcm[key] for key in ('packet_zerofill_white', 'packet_local_fallback_white')}},
        'E13-13': {'transport': report['transport'], 'parameters': report['parameters']['transport'],
                   'float': report['float_components']['transport_pcm16_white'], 'PCM': pcm['transport_pcm16_white'],
                   'three_tier_audit': {
                       'target_steady_power': POWER, 'finite_actual_target_power': report['finite_steady_target_power'],
                       'float_basis_covariance': report['finite_steady_basis_covariance'],
                       'finite_source_noise_cross': report['finite_steady_source_noise_cross'],
                       'PCM': {key: pcm[key] for key in ('central_white', 'central_correlated', 'central_node2_correlated', 'stale_correlated')},
                       'scope': 'all integer scores use the genuine same-reference PCM; no fitted gain, delay or averaged dB',
                   }},
        'E13-14': {'unrelaxed': jacobi_control(), 'half': jacobi_control(.5, 40)},
        'E13-15': step_size_control(),
        'E13-16': tree_sum_control(),
        'E13-17': gevd_control(),
        'E13-18': broadcast_statistics_control(),
        'E13-19': ptp_vs_sro_phase_budget(),
        'E13-20': dasf_unification_control(),
        'E13-21': dwacd_sro_gating_control(),
    }
    if set(exercises) != set(EXERCISE_IDS):
        raise AssertionError('exercise ID set differs')
    return _json({'schema_version': 1, 'scope': 'original finite known-model teaching controls',
                  'complex_encoding': 'complex scalar -> [real,imag], complex arrays -> trailing [real,imag] axis',
                  'exercises': exercises})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exercise', choices=EXERCISE_IDS)
    args = parser.parse_args(argv)
    result = run_experiments()
    if args.exercise:
        result['exercises'] = {args.exercise: result['exercises'][args.exercise]}
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
