"""Independent-scope AEC arithmetic experiments E06-12..20, E06-26..27, and E06-29..32.

All values are dimensionless unless units are explicitly recorded. These
small cases isolate identifiability, state timing and score definitions;
they are not speech-quality, device, or convergence benchmarks.
"""
from __future__ import annotations

# Allow the documented direct-file command as well as python -m.
if __name__ == "__main__" and not __package__:
    import sys as _chapter_entry_sys
    from pathlib import Path as _ChapterEntryPath
    _chapter_entry_sys.path.insert(0, str(_ChapterEntryPath(__file__).resolve().parents[3]))

import json
import numpy as np

from codes.chapters.ch06.core.aec import NLMSState
from codes.chapters.ch06.aec_kalman_scalar_demo import scalar_kalman_step
from codes.chapters.ch06.aec_controlled_doubletalk import increment_metrics
from codes.chapters.ch06.core.reference_timing import reference_gain_order, echo_tail_activity
from codes.chapters.ch06.core.nonlinear_and_stereo import (
    hammerstein_orthogonal_basis,
    stereo_decorrelation_and_loop_deadlines,
)


def shared_reference_normalization() -> dict:
    x = np.ones(2)
    records = {}
    for label, denominator in [('joint', float(x @ x)), ('separate', 1.)]:
        w = np.zeros(2)
        errors, paths = [], []
        for _ in range(3):
            error = 1. - float(w @ x)
            w += 1.5 * error * x / denominator
            errors.append(error)
            paths.append(w.tolist())
        records[label] = {'prior_errors': errors, 'posterior_weights': paths,
                          'next_error': float(1. - w @ x)}
    return {'reference': x.tolist(), 'microphone': 1., 'step_size': 1.5,
            'initial_weights': [0., 0.], 'normalization': records,
            'scope': 'one joint two-coordinate update, not two independent echo targets'}


def reference_identifiability() -> dict:
    training = np.array([[1., 1.], [-1., -1.]])  # rows are time samples
    truth, alternative = np.array([1., 2.]), np.array([1.5, 1.5])
    holdout = np.array([1., -1.])
    covariance = training.T @ training / len(training)
    return {'training_rows': training.tolist(), 'true_weights': truth.tolist(),
            'alternative_weights': alternative.tolist(), 'correlation': covariance.tolist(),
            'correlation_eigenvalues': np.linalg.eigvalsh(covariance).tolist(),
            'true_training_echo': (training @ truth).tolist(),
            'alternative_training_echo': (training @ alternative).tolist(),
            'holdout': holdout.tolist(), 'true_holdout_echo': float(holdout @ truth),
            'alternative_holdout_echo': float(holdout @ alternative),
            'scope': 'equal programs leave a null direction; no holdout adaptation'}


def causal_delay_support() -> dict:
    physical_delay, taps = 3, 3
    cases = []
    for reference_delay in (1, 3, 4):
        remaining = physical_delay - reference_delay
        cases.append({'reference_delay_samples': reference_delay,
                      'remaining_path_delay_samples': remaining,
                      'representable_by_causal_fir': bool(0 <= remaining < taps)})
    return {'physical_delay_samples': physical_delay, 'filter_taps': taps,
            'support_samples': [0, taps - 1], 'cases': cases,
            'scope': 'pure integer delay; delaying the algorithm reference too far requires prediction'}


def clock_drift_phase() -> dict:
    fs, frequency, rho = 16000., 1000., 250e-6
    seconds = np.array([0., 1., 2.])
    delay = rho * fs * seconds
    phase = 2 * np.pi * frequency * delay / fs
    return {'sample_rate_hz': fs, 'frequency_hz': frequency, 'drift_ppm': rho * 1e6,
            'times_seconds': seconds.tolist(), 'uncompensated_delay_samples': delay.tolist(),
            'phase_error_radians': phase.tolist(),
            'residual_to_echo_power_ratio': (2 - 2 * np.cos(phase)).tolist(),
            'scope': 'declared linearly increasing required delay; instantaneous single-frequency equal-amplitude phasor score, not whole-record ERLE'}


def missing_reference_block() -> dict:
    x, h = np.arange(1., 7.), np.array([1., .5])
    available = x.copy()
    available[2:4] = 0.
    echo = np.convolve(x, h)[:x.size]
    state = NLMSState(2, initial_weights=h)
    residual, predicted = state.process(available, echo, freeze=np.ones(x.size, bool))
    return {'physical_reference': x.tolist(), 'available_reference': available.tolist(),
            'path': h.tolist(), 'echo': echo.tolist(), 'prediction': predicted.tolist(),
            'prior_residual': residual.tolist(), 'missing_interval_samples': [2, 4],
            'first_restored_reference_sample': 4, 'first_fully_recovered_output_sample': 5,
            'scope': 'known fixed path and frozen coefficients; reference history still advances'}


def prior_posterior_scoring() -> dict:
    result = scalar_kalman_step(previous_weight=0., previous_variance=1., transition=1.,
                               process_variance=0., reference=1., observation=2.,
                               observation_variance=1.)
    prior = float(result['prior_error'].real)
    posterior = float(result['posterior_error'].real)
    return {'inputs': {'prior_weight': 0, 'prior_variance': 1, 'reference': 1,
                       'observation': 2, 'observation_variance': 1},
            'gain': float(result['gain'].real), 'posterior_weight': float(result['posterior_weight'].real),
            'prior_error': prior, 'posterior_error': posterior,
            'same_observation_fit_power_difference_db': float(20 * np.log10(abs(prior / posterior))),
            'scope': 'posterior fits the scored observation; not a causal cancellation improvement'}


def block_db_aggregation() -> dict:
    before, after = np.array([1., 100.]), np.array([.01, 100.])
    scores = 10 * np.log10(before / after)
    return {'equal_length_block_input_energies': before.tolist(),
            'block_output_energies': after.tolist(), 'block_power_ratios_db': scores.tolist(),
            'arithmetic_mean_block_db': float(np.mean(scores)),
            'combined_energy_ratio_db': float(10 * np.log10(before.sum() / after.sum())),
            'scope': 'different defined summaries; no double talk and no discarded/failed blocks'}


def normalized_update_expectation() -> dict:
    x, s = np.array([1., -1., 2., -2.]), np.array([1., -1., -.5, .5])
    h = .8
    d = h * x + s
    return {'reference_states': x.tolist(), 'near_states': s.tolist(),
            'state_probabilities': [.25] * 4, 'true_path': h,
            'cross_second_moment': float(np.mean(x * s)),
            'mean_near_over_reference': float(np.mean(s / x)),
            'batch_least_squares_path': float(x @ d / (x @ x)),
            'expected_nlms_equilibrium': float(h + np.mean(s / x)),
            'scope': 'IID equal-probability state draws, epsilon=0, scalar constant-step NLMS, current draw independent of prior weight; not a four-sample online trajectory'}


def increment_orthogonal_decomposition() -> dict:
    target = np.array([1., -1., 1., -1.])
    orthogonal = .3 * np.array([1., 1., -1., -1.])
    increment = .8 * target + orthogonal
    metrics = increment_metrics(np.zeros(4), increment, target, 0, 4, sample_unit='dimensionless')
    gain = metrics['increment_projection_gain_no_delay_fit']
    remainder = increment - gain * target
    return {'target': target.tolist(), 'increment': increment.tolist(),
            'projection_gain': gain, 'orthogonal_remainder': remainder.tolist(),
            'target_remainder_inner_product': float(target @ remainder),
            'gain_error_term': float((gain - 1)**2),
            'orthogonal_error_term': float(remainder @ remainder / (target @ target)),
            'total_relative_squared_error': metrics['fixed_sample_relative_squared_error'],
            'scope': 'output increment may contain adaptive-state changes; decomposition is algebra, not isolated near-end recovery'}


def pbfdkf_double_talk_and_forgetting() -> dict:
    """E06-31: PBFDKF measurement-variance inflation under double-talk and state forgetting factor A."""
    true_weight = 1.0
    initial_weight = 0.8
    initial_error = true_weight - initial_weight
    initial_misalignment_db = float(20.0 * np.log10(abs(initial_error) / abs(true_weight)))
    prior_variance = 0.04
    regressor = 1.0

    def evaluate_regime(near_end: float, psi_e: float) -> dict:
        residual = initial_error * regressor + near_end
        kalman_gain = float(prior_variance * regressor / (prior_variance * (regressor ** 2) + psi_e))
        weight_increment = float(kalman_gain * residual)
        posterior_weight = float(initial_weight + weight_increment)
        posterior_error = float(true_weight - posterior_weight)
        posterior_variance = float((1.0 - kalman_gain * regressor) * prior_variance)
        error_change_db = float(20.0 * np.log10(abs(posterior_error) / abs(initial_error)))
        normalized_misalignment_db = float(20.0 * np.log10(abs(posterior_error) / abs(true_weight)))
        return {
            'near_end_s': float(near_end),
            'residual_e': float(residual),
            'psi_e': float(psi_e),
            'kalman_gain': kalman_gain,
            'weight_increment': weight_increment,
            'posterior_weight': posterior_weight,
            'posterior_error': posterior_error,
            'abs_posterior_error': float(abs(posterior_error)),
            'posterior_variance': posterior_variance,
            'error_change_db': error_change_db,
            'normalized_misalignment_db': normalized_misalignment_db,
        }

    dt_residual = initial_error * regressor + 2.0
    opp_residual = initial_error * regressor - 2.0
    regimes = {
        'clean_far_end_single_talk': evaluate_regime(0.0, 0.01),
        'unguarded_double_talk': evaluate_regime(2.0, 0.01),
        'adapted_pbfdkf_double_talk': evaluate_regime(2.0, float(dt_residual ** 2)),
        'adapted_pbfdkf_ml_variance': evaluate_regime(
            2.0, float(max(0.0, dt_residual ** 2 - (regressor ** 2) * prior_variance))
        ),
        'adapted_pbfdkf_opposite_polarity': evaluate_regime(-2.0, float(opp_residual ** 2)),
    }

    sigma_w_sq = 1.0
    psi_e_ss = 0.01
    steady_state = {}
    for label, a_val in [('A_0_999', 0.999), ('A_0_99', 0.99)]:
        psi_delta = float((1.0 - a_val ** 2) * sigma_w_sq)
        # Positive root of (P^-)^2 + (psi_e*(1 - A^2) - psi_delta)*P^- - psi_delta*psi_e = 0
        b_coef = psi_e_ss * (1.0 - a_val ** 2) - psi_delta
        c_coef = -psi_delta * psi_e_ss
        p_minus_ss = float((-b_coef + np.sqrt(b_coef ** 2 - 4.0 * c_coef)) / 2.0)
        k_ss = float(p_minus_ss / (p_minus_ss + psi_e_ss))
        p_post_ss = float((1.0 - k_ss) * p_minus_ss)
        steady_state[label] = {
            'state_transition_A': float(a_val),
            'process_variance_psi_delta': psi_delta,
            'time_constant_blocks': float(-1.0 / np.log(a_val)),
            'prior_variance_ss': p_minus_ss,
            'kalman_gain_ss': k_ss,
            'posterior_variance_ss': p_post_ss,
            'excess_misalignment_ss_db': float(10.0 * np.log10(p_post_ss / sigma_w_sq)),
        }

    p_999 = steady_state['A_0_999']['posterior_variance_ss']
    p_99 = steady_state['A_0_99']['posterior_variance_ss']
    k_999 = steady_state['A_0_999']['kalman_gain_ss']
    k_99 = steady_state['A_0_99']['kalman_gain_ss']

    return {
        'model': {
            'true_weight_W': float(true_weight),
            'initial_weight_W0': float(initial_weight),
            'initial_error_W0_tilde': float(initial_error),
            'initial_misalignment_db': initial_misalignment_db,
            'prior_variance_P0_minus': float(prior_variance),
            'regressor_X': float(regressor),
            'measurement_noise_v': 0.0,
        },
        'regimes': regimes,
        'steady_state_forgetting': {
            'sigma_W_sq': float(sigma_w_sq),
            'psi_e': float(psi_e_ss),
            'cases': steady_state,
            'gain_ratio_099_over_0999': float(k_99 / k_999),
            'excess_floor_increase_db': float(10.0 * np.log10(p_99 / p_999)),
        },
        'scope': (
            'Single-bin scalar frequency-domain Kalman AEC update and Riccati fixed-point '
            'comparison; illustrates variance-gated soft step control rather than full-band speech scoring.'
        ),
    }


def res_wiener_and_cng_matching() -> dict:
    """E06-32: residual echo suppression Wiener gain with floor G_min and CNG level matching."""
    g_min = 0.1
    phi_nn = 0.01

    def evaluate_frame(phi_ss: float, phi_bb: float) -> dict:
        g_raw = float(phi_ss / (phi_ss + phi_bb + phi_nn))
        g_res = float(max(g_min, g_raw))
        clamped = bool(g_raw < g_min)
        gain_db = float(20.0 * np.log10(g_res))
        echo_suppression_db = float(-gain_db)
        output_echo_psd = float((g_res ** 2) * phi_bb)
        output_noise_no_cng = float((g_res ** 2) * phi_nn)
        output_noise_no_cng_db = float(10.0 * np.log10(output_noise_no_cng))
        pumping_no_cng_db = float(10.0 * np.log10(output_noise_no_cng / phi_nn))
        phi_cng = float(max(0.0, (1.0 - g_res ** 2) * phi_nn))
        cng_amp_scale = float(np.sqrt(max(0.0, 1.0 - g_res ** 2)))
        output_noise_with_cng = float(output_noise_no_cng + phi_cng)
        pumping_with_cng_db = float(10.0 * np.log10(output_noise_with_cng / phi_nn))
        return {
            'phi_ss': float(phi_ss),
            'phi_bb': float(phi_bb),
            'phi_nn': float(phi_nn),
            'g_unclamped': g_raw,
            'g_res': g_res,
            'clamped_by_floor': clamped,
            'gain_db': gain_db,
            'speech_attenuation_db': gain_db,
            'echo_suppression_db': echo_suppression_db,
            'output_echo_psd': output_echo_psd,
            'output_noise_psd_without_cng': output_noise_no_cng,
            'output_noise_db_without_cng': output_noise_no_cng_db,
            'noise_pumping_db_without_cng': pumping_no_cng_db,
            'cng_psd': phi_cng,
            'cng_amplitude_scale': cng_amp_scale,
            'output_noise_psd_with_cng': output_noise_with_cng,
            'noise_pumping_db_with_cng': pumping_with_cng_db,
        }

    return {
        'g_min': float(g_min),
        'g_min_db': float(20.0 * np.log10(g_min)),
        'target_background_noise_psd': float(phi_nn),
        'target_background_noise_db': float(10.0 * np.log10(phi_nn)),
        'frames': {
            'far_end_single_talk': evaluate_frame(0.0, 0.25),
            'double_talk': evaluate_frame(1.0, 0.25),
        },
        'scope': (
            'Single-bin stationary uncorrelated PSD algebra for Wiener RES gain clamping '
            'and CNG power matching; assumes accurate PSD estimation and uncorrelated CNG phase.'
        ),
    }


def run_experiments() -> dict:
    functions = [
        shared_reference_normalization,
        reference_identifiability,
        causal_delay_support,
        clock_drift_phase,
        missing_reference_block,
        prior_posterior_scoring,
        block_db_aggregation,
        normalized_update_expectation,
        increment_orthogonal_decomposition,
    ]
    cases = {f'E06-{index:02d}': function() for index, function in enumerate(functions, 12)}
    cases['E06-26'] = reference_gain_order()
    cases['E06-27'] = echo_tail_activity()
    cases['E06-29'] = hammerstein_orthogonal_basis()
    cases['E06-30'] = stereo_decorrelation_and_loop_deadlines()
    cases['E06-31'] = pbfdkf_double_talk_and_forgetting()
    cases['E06-32'] = res_wiener_and_cng_matching()
    return cases


if __name__ == '__main__':
    print(json.dumps(run_experiments(), ensure_ascii=False, indent=2, allow_nan=False))
