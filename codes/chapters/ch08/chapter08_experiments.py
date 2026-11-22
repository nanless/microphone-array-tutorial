"""Original small experiments E08-10..27, without upstream imports.

Run with python -m codes.chapters.ch08.chapter08_experiments. No files are written.
Substeps and supplied source/slot models are not complete separation systems.
"""
from __future__ import annotations

# Allow the documented direct-file command as well as python -m.
if __name__ == "__main__" and not __package__:
    import sys as _chapter_entry_sys
    from pathlib import Path as _ChapterEntryPath
    _chapter_entry_sys.path.insert(0, str(_ChapterEntryPath(__file__).resolve().parents[3]))

import json
import hashlib
import math
from pathlib import Path
import struct
import wave
import numpy as np

from codes.chapters.ch00.core.audio_samples import css_overlap_case
from codes.chapters.ch08.core.css import match_two_source_overlap, overlap_application_gain
from codes.chapters.ch08.core.separation import guided_activity_posterior
from codes.chapters.ch08.core.separation_metrics_streaming import (
    ci_sdr_vs_si_sdr_experiment,
    streaming_mac_and_bypass_experiment,
)
from codes.chapters.ch00.cross_chapter.enhancement_step_exercises import ip_row

ROOT = Path(__file__).resolve().parents[3]
CSS_STEMS = ('css_overlap_reference', 'css_overlap_mixture',
             'css_overlap_naive', 'css_overlap_aligned')


def _integer_centered_si_sdr(estimate, reference):
    """Independent centered PCM score from exact integer second moments.

    Both inputs are decoded signed PCM integers. No waveform gain or delay is
    fitted apart from the scalar projection that defines SI-SDR itself.
    Centered moment numerators retain the common sample denominator exactly.
    This fixture has nonzero projection and residual; silent/perfect cases are
    rejected rather than assigned the teaching score's finite numerical caps.
    """
    n = len(reference)
    if n == 0 or len(estimate) != n:
        raise ValueError('PCM score requires equal nonempty records')
    es = n * sum(x*x for x in reference) - sum(reference)**2
    ee = n * sum(x*x for x in estimate) - sum(estimate)**2
    cross = n * sum(x*y for x, y in zip(estimate, reference)) - sum(estimate)*sum(reference)
    target = cross**2
    residual = es*ee - target
    if es <= 0 or ee <= 0 or target <= 0 or residual <= 0:
        raise ValueError('PCM fixture score is undefined, orthogonal or perfect')
    return {'si_sdr_db': 10*(math.log10(target)-math.log10(residual)),
            'reference_centered_energy_numerator': es,
            'estimate_centered_energy_numerator': ee,
            'centered_cross_numerator': cross,
            'projection_energy_ratio_numerator': target,
            'residual_energy_ratio_denominator': residual}


def css_published_audio(*, manifest_path=None, audio_directory=None):
    """Read/check the four published CSS WAVs and score their actual PCM.

    This function never generates or repairs assets. The main manifest must
    describe all current generator inputs, and each CSS file must have its
    declared SHA, PCM16 format, length, channels and common gain. Statistics
    are independently recomputed from integer samples, not copied from the
    manifest or regenerated/quantized float arrays.
    """
    from codes.chapters.ch00.examples.generate_audio_samples import INPUTS
    manifest_path = Path(manifest_path) if manifest_path is not None else ROOT/'codes/chapters/ch00/audio/MANIFEST.json'
    directory = Path(audio_directory) if audio_directory is not None else ROOT/'codes/chapters/ch08/audio'
    manifest = json.loads(manifest_path.read_text(),
                          parse_constant=lambda _: (_ for _ in ()).throw(
                              ValueError('main audio manifest contains a nonfinite JSON value')))
    sources = manifest.get('generator_inputs', {})
    if set(sources) != set(INPUTS):
        raise ValueError('main audio generator input set differs')
    for name, digest in sources.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'main audio source SHA differs: {name}')
    group = manifest['groups']['css_overlap']
    parameters = group['parameters']
    if group['common_export_gain'] != 1 or parameters['common_export_gain'] != 1:
        raise ValueError('CSS common export gain differs')
    if parameters['score_interval_samples'] != [0, 32000]:
        raise ValueError('CSS score interval differs')
    records = [item for item in manifest['files'] if item['group'] == 'css_overlap']
    if len(records) != 4 or {item['file'] for item in records} != {s+'.wav' for s in CSS_STEMS}:
        raise ValueError('CSS manifest file set differs')
    decoded, measurements = {}, {}
    for item in records:
        stem = Path(item['file']).stem
        path = directory/item['file']
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        expected_channels = 1 if stem == 'css_overlap_mixture' else 2
        if (digest != item['sha256'] or item['chapter'] != 'ch08'
                or item['common_export_gain'] != 1 or item['sample_rate_hz'] != 16000
                or item['samples'] != 32000 or item['channels'] != expected_channels):
            raise ValueError(f'CSS asset metadata or SHA differs: {stem}')
        with wave.open(str(path), 'rb') as handle:
            if (handle.getsampwidth() != 2 or handle.getcomptype() != 'NONE'
                    or handle.getframerate() != 16000 or handle.getnframes() != 32000
                    or handle.getnchannels() != expected_channels):
                raise ValueError(f'CSS PCM format differs: {stem}')
            blob = handle.readframes(32000)
        if len(blob) != 2*32000*expected_channels:
            raise ValueError(f'CSS PCM data is truncated: {stem}')
        integers = struct.unpack('<'+'h'*(32000*expected_channels), blob)
        decoded[stem] = [list(integers[channel::expected_channels]) for channel in range(expected_channels)]
        measurements[stem] = {'sha256': digest, 'channels': expected_channels,
                              'samples_per_channel': 32000}
    reference = decoded['css_overlap_reference']
    for stem in CSS_STEMS[1:]:
        channels = decoded[stem]
        scores = [_integer_centered_si_sdr(channels[min(i, len(channels)-1)], reference[i]) for i in range(2)]
        measurements[stem]['fixed_identity_scores'] = scores
        recorded = parameters['pcm_si_sdr_db'][stem]
        if (not isinstance(recorded, list) or len(recorded) != 2 or
                any(type(value) not in (int, float) or not math.isfinite(value)
                    for value in recorded)):
            raise ValueError(f'CSS recorded PCM scores require two finite numeric values: {stem}')
        if any(abs(score['si_sdr_db']-expected) > 1e-10 for score, expected in zip(scores, recorded)):
            raise ValueError(f'CSS recorded PCM score differs: {stem}')
    return {'manifest_path': str(manifest_path.relative_to(ROOT)) if manifest_path.is_relative_to(ROOT) else str(manifest_path),
            'source_sha256': sources, 'common_export_gain': 1,
            'score_interval_samples': [0, 32000], 'samples_per_channel': 32000,
            'centering': 'per record; exact integer centered second moments',
            'fixed_identity_no_pit_no_delay_fit': True,
            'actual_pcm_measurements': measurements}


def whitening_minicase():
    """Four equiprobable non-Gaussian points; whitening is not independence."""
    sources = np.array([[-1., -1.], [-1., 1.], [1., -1.], [1., 1.]]).T
    mixing = np.array([[1., .5], [.5, 1.]])
    mixture = mixing @ sources
    whitener = np.linalg.inv(mixing)  # this positive symmetric A gives a valid whitening
    whitened = whitener @ mixture
    rotation = np.array([[1., 1.], [-1., 1.]])/np.sqrt(2)
    rotated = rotation @ whitened
    # Event identities are derived from the exact integer source coordinates,
    # rather than treating a floating-point tolerance as a probability model.
    zero1 = sources[0]+sources[1] == 0
    zero2 = -sources[0]+sources[1] == 0
    return {'source_points': sources, 'mixing': mixing, 'mixture_points': mixture,
            'mixture_covariance': mixture @ mixture.T/4,
            'whitener': whitener, 'whitened_points': whitened,
            'rotation': rotation, 'rotated_points': rotated,
            'rotated_covariance': rotated @ rotated.T/4,
            'probability_first_zero': np.mean(zero1), 'probability_second_zero': np.mean(zero2),
            'probability_both_zero': np.mean(zero1 & zero2),
            'independent_product_probability': np.mean(zero1)*np.mean(zero2),
            'rank_one_covariance': np.ones((2, 2))}


def congruence_model_boundary():
    """With the first shape I, the whitened remaining shapes must commute."""
    shapes = np.array([np.eye(2), np.diag([1., 2.]), [[2., 1.], [1., 2.]]])
    commutator = shapes[1] @ shapes[2]-shapes[2] @ shapes[1]
    return {'shapes': shapes, 'shape_eigenvalues': np.linalg.eigvalsh(shapes),
            'whitened_remaining_commutator': commutator,
            'joint_congruence_exists_for_all_three': False,
            'scope': 'positive definite first shape I; pairwise existence does not imply joint existence'}


def fastmnmf_scale_objective():
    """A fixed model's physical likelihood is invariant under Q/g rescaling."""
    q = np.array([[1., -1.], [0., 1.]])
    shapes = np.array([[[3., 2.], [2., 2.]], [[4., 1.], [1., 1.]]])
    diagonal = np.array([np.diag(q @ shape @ q.T) for shape in shapes])
    observation = np.array([2., 1.])
    cases = []
    for multiplier in (1., 2.):
        transform = multiplier*q
        d = multiplier**2*diagonal.sum(0)
        z = transform @ observation
        quadratic = float(np.sum(abs(z)**2/d))
        log_variance = float(np.log(d).sum())
        jacobian = float(-2*np.linalg.slogdet(transform)[1])
        physical_covariance = np.linalg.solve(transform, np.diag(d)) @ np.linalg.inv(transform).T
        cases.append({'Q': transform, 'diagonal_variance': d, 'transformed_observation': z,
                      'quadratic_term': quadratic, 'log_variance_term': log_variance,
                      'jacobian_term': jacobian, 'full_objective': quadratic+log_variance+jacobian,
                      'objective_without_jacobian': quadratic+log_variance,
                      'physical_covariance': physical_covariance})
    return {'cases': cases, 'frames': 1, 'frequencies': 1,
            'scope': 'proper complex likelihood, fixed parameters; not an estimated FastMNMF model'}


def mixture_consistency_minicase():
    """Fixed positive-weight KKT projection; sum consistency is not separation."""
    estimate, mixture, variance = np.array([1., 1.]), 4., np.array([1., 3.])
    correction = mixture-estimate.sum()
    projected = estimate+variance/variance.sum()*correction
    equal = estimate+correction/len(estimate)
    return {'mixture': mixture, 'estimates': estimate, 'variances': variance,
            'mixture_residual': correction, 'weighted_projection': projected,
            'equal_projection': equal,
            'weighted_cost_at_weighted_projection': float(.5*np.sum((projected-estimate)**2/variance)),
            'weighted_cost_at_equal_projection': float(.5*np.sum((equal-estimate)**2/variance)),
            'kkt_correction_over_variance': (projected-estimate)/variance,
            'scope': 'one shared reference and scale; positive fixed weights, not learned confidence'}


def css_polarity_minicase():
    """Absolute-correlation matching fixes permutation, not polarity or gain."""
    previous = np.array([[1., -1., 0., 0.], [0., 0., 1., -1.]])
    cases = []
    for multiplier in (-1., -2.):
        current = multiplier*previous
        matching = match_two_source_overlap(previous, current)
        aligned = current[matching['current_indices_for_previous']]
        gain = overlap_application_gain(previous, aligned)
        signed = previous @ current.T/2/abs(multiplier)
        cases.append({'current_multiplier': multiplier, 'current': current,
                      'matching': matching, 'signed_centered_correlation': signed,
                      'naive_equal_overlap': (previous+aligned)/2,
                      'polarity_only_overlap': (previous+np.sign(multiplier)*aligned)/2,
                      'direct_application_gain': gain,
                      'corrected_equal_overlap': (previous+gain[:, None]*aligned)/2})
    return {'previous': previous, 'samples_per_overlap': 4, 'cases': cases,
            'scope': 'same sample times, nonzero clean overlap, known scalar relation; no universal identity/delay repair'}


def _auxiva_cost(w_matrix: np.ndarray, covariances: tuple[np.ndarray, ...]) -> float:
    """Evaluate fixed-V AuxIVA auxiliary cost sum_n w_n^H V_n w_n - 2 ln|det W|."""
    w = np.asarray(w_matrix, dtype=np.complex128)
    quad_sum = 0.0
    for idx, cov in enumerate(covariances):
        col = w[idx].conj()
        quad_sum += float(np.vdot(col, cov @ col).real)
    sign, logdet = np.linalg.slogdet(w)
    if sign == 0 or not math.isfinite(float(logdet)):
        raise ValueError('demixing matrix must be nonsingular')
    return float(quad_sum - 2.0 * float(logdet))


def _auxiva_ip_sweep(w_matrix: np.ndarray, covariances: tuple[np.ndarray, ...]):
    """Run one sequential Iterative Projection (IP) sweep using M linear solves."""
    w = np.asarray(w_matrix, dtype=np.complex128).copy()
    dim = w.shape[0]
    steps = []
    for src in range(dim):
        cov = covariances[src]
        e_vec = np.zeros(dim, dtype=np.complex128)
        e_vec[src] = 1.0
        direction = np.linalg.solve(w @ cov, e_vec)
        quad = float(np.vdot(direction, cov @ direction).real)
        w[src] = (direction / math.sqrt(quad)).conj()
        steps.append({
            'source_index': src + 1,
            'unnormalized_direction': direction.copy(),
            'direction_quadratic_form': quad,
            'W_after': w.copy(),
            'cost_after': _auxiva_cost(w, covariances),
        })
    return w, steps


def _auxiva_iss_sweep(w_matrix: np.ndarray, covariances: tuple[np.ndarray, ...]):
    """Run one inversion-free Iterative Source Steering (ISS) rank-1 sweep."""
    w = np.asarray(w_matrix, dtype=np.complex128).copy()
    dim = w.shape[0]
    steps = []
    for k in range(dim):
        w_k = w[k].conj()
        v_vec = np.zeros(dim, dtype=np.complex128)
        for n in range(dim):
            cov_n = covariances[n]
            denom = float(np.vdot(w_k, cov_n @ w_k).real)
            if denom <= 0.0 or not math.isfinite(denom):
                raise ValueError('ISS denominator must be strictly positive and finite')
            if n == k:
                v_vec[k] = 1.0 - 1.0 / math.sqrt(denom)
            else:
                w_n = w[n].conj()
                v_vec[n] = np.vdot(w_n, cov_n @ w_k) / denom
        w = w - np.outer(v_vec, w[k])
        steps.append({
            'steering_index': k + 1,
            'steering_vector': v_vec.copy(),
            'W_after': w.copy(),
            'cost_after': _auxiva_cost(w, covariances),
        })
    return w, steps


def auxiva_iss_vs_ip_experiment(sweeps: int = 5) -> dict:
    """Compare AuxIVA IP vs inversion-free ISS rank-1 sweeps on a 2x2 complex bin."""
    if isinstance(sweeps, bool) or not isinstance(sweeps, int) or sweeps < 1:
        raise ValueError('sweeps must be a positive integer')
    v1 = np.array([[4.0, 4.0j], [-4.0j, 20.0]], dtype=np.complex128)
    v2 = np.array([[20.0, 4.0j], [-4.0j, 4.0]], dtype=np.complex128)
    covariances = (v1, v2)
    w0 = np.eye(2, dtype=np.complex128)

    # Exact closed-form global minimum via generalized eigenvalues of (V1, V2).
    mu1 = 3.0 - 2.0 * math.sqrt(2.0)
    mu2 = 3.0 + 2.0 * math.sqrt(2.0)
    det_v2 = float(np.linalg.det(v2).real)
    exact_global_min_cost = 2.0 + math.log(mu1 * det_v2)

    q0 = _auxiva_cost(w0, covariances)
    w_ip = w0.copy()
    w_iss = w0.copy()
    ip_costs = [q0]
    iss_costs = [q0]
    sweep1_ip = None
    sweep1_iss = None
    for idx in range(sweeps):
        w_ip, ip_steps = _auxiva_ip_sweep(w_ip, covariances)
        w_iss, iss_steps = _auxiva_iss_sweep(w_iss, covariances)
        ip_costs.append(ip_steps[-1]['cost_after'])
        iss_costs.append(iss_steps[-1]['cost_after'])
        if idx == 0:
            sweep1_ip = ip_steps
            sweep1_iss = iss_steps

    return {
        'V1': v1,
        'V2': v2,
        'W0': w0,
        'sweeps': sweeps,
        'initial_cost': q0,
        'generalized_eigenvalues': [mu1, mu2],
        'exact_global_min_cost': exact_global_min_cost,
        'sweep1_ip': sweep1_ip,
        'sweep1_iss': sweep1_iss,
        'ip_cost_history': ip_costs,
        'iss_cost_history': iss_costs,
        'ip_monotone_decreasing': all(b < a for a, b in zip(ip_costs, ip_costs[1:])),
        'iss_monotone_decreasing': all(b < a for a, b in zip(iss_costs, iss_costs[1:])),
        'ip_final_W': w_ip,
        'iss_final_W': w_iss,
        'ip_final_cost': ip_costs[-1],
        'iss_final_cost': iss_costs[-1],
        'final_cost_gap_ip_vs_iss': abs(ip_costs[-1] - iss_costs[-1]),
        'ip_gap_to_optimum': abs(ip_costs[-1] - exact_global_min_cost),
        'iss_gap_to_optimum': abs(iss_costs[-1] - exact_global_min_cost),
        'ip_final_weighted_grams': [w_ip @ v1 @ w_ip.conj().T, w_ip @ v2 @ w_ip.conj().T],
        'iss_final_weighted_grams': [w_iss @ v1 @ w_iss.conj().T, w_iss @ v2 @ w_iss.conj().T],
        'ip_matrix_solves_per_sweep': 2,
        'iss_matrix_solves_per_sweep': 0,
        'scope': (
            'fixed 2x2 Hermitian positive-definite weighted covariances at one frequency bin; '
            'compares auxiliary-function IP against inversion-free ISS rank-1 sweeps, '
            'not full broadband AuxIVA'
        ),
    }


def _plain(value):
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {'real': value.real.tolist(), 'imag': value.imag.tolist()}
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def run_experiments() -> dict:
    result = {}
    covariance = np.array([[2, 1j], [-1j, 2]])
    sequential = np.eye(2, dtype=complex)
    for row in range(2):
        sequential[row] = ip_row(sequential, covariance, row).conj()
    simultaneous = np.array([ip_row(np.eye(2), covariance, row).conj() for row in range(2)])
    result['E08-10'] = {'covariance': covariance, 'sequential_rows': sequential,
                        'simultaneous_rows': simultaneous,
                        'sequential_weighted_gram': sequential @ covariance @ sequential.conj().T,
                        'simultaneous_weighted_gram': simultaneous @ covariance @ simultaneous.conj().T}

    y = np.array([[1., 0., 1.], [0., 1., 1.]])
    reference = np.array([1., 2., 3.])
    gram = y @ y.T
    joint = np.linalg.solve(gram, y @ reference)
    scalar = (y @ reference) / np.diag(gram)
    result['E08-11'] = {'outputs': y, 'reference': reference, 'gram': gram,
                        'joint_gains': joint, 'individual_gains': scalar,
                        'joint_sum': joint @ y, 'individual_sum': scalar @ y,
                        'complex_application_gain': np.array([2 + 1j]),
                        'pyroomacoustics_returned_conjugate': np.array([2 - 1j])}

    activity = np.array([[1], [0], [0]])
    densities = np.array([[9., 1.], [1., 1.], [1., 1.]])
    prior = np.array([.5, .5])
    posterior = guided_activity_posterior(prior, densities, activity)
    updated = posterior.mean(0)
    before, after = np.log(prior @ densities[0]), np.log(updated @ densities[0])
    result['E08-12'] = {'activity': activity, 'fixed_densities': densities, 'posterior': posterior,
                        'initial_priors': prior, 'mean_posterior_priors': updated,
                        'conditional_log_likelihood_before': before, 'conditional_log_likelihood_after': after,
                        'change': after - before,
                        'scope': 'shapes fixed; inactive-only background frames have conditional density one'}

    power = np.array([[1., 1.], [9., 9.]])
    basis, activation = np.ones((2, 1)), np.ones((1, 2))
    def objective(b, h):
        model = b @ h
        return float(np.sum(power / model + np.log(model)))
    costs = [objective(basis, activation)]
    model = basis @ activation
    basis *= np.sqrt(((power / model ** 2) @ activation.T) / ((1 / model) @ activation.T))
    costs.append(objective(basis, activation))
    model = basis @ activation
    activation *= np.sqrt((basis.T @ (power / model ** 2)) / (basis.T @ (1 / model)))
    costs.append(objective(basis, activation))
    result['E08-13'] = {'observed_power': power, 'updated_basis': basis, 'updated_activation': activation,
                        'updated_power': basis @ activation, 'objectives_initial_basis_activation': costs,
                        'scope': 'positive one-basis Gaussian/IS NMF substeps with demixing fixed'}

    a = np.array([.5 * np.eye(2), .25 * np.eye(2), .25 * np.eye(2)])
    b = np.array([np.diag([.6, .2]), np.diag([.2, .4]), np.diag([.2, .4])])
    observation = np.array([1., 2.])
    result['E08-14'] = {'source_scms_a': a, 'source_scms_b': b,
                        'observation_scm_a': a.sum(0), 'observation_scm_b': b.sum(0),
                        'first_wiener_image_a': a[0] @ np.linalg.solve(a.sum(0), observation),
                        'first_wiener_image_b': b[0] @ np.linalg.solve(b.sum(0), observation)}

    q = np.array([[1., -1.], [0., 1.]])
    g1, g2 = np.array([[3., 2.], [2., 2.]]), np.array([[4., 1.], [1., 1.]])
    diagonal = np.array([q @ g @ q.T for g in (g1, g2)])
    x = np.array([2., 1.])
    z = q @ x
    separated_z = np.array([np.diag(d) / np.diag(diagonal.sum(0)) * z for d in diagonal])
    images = np.array([np.linalg.solve(q, item) for item in separated_z])
    result['E08-15'] = {'Q': q, 'shapes': np.array([g1, g2]), 'congruence_diagonals': diagonal,
                        'QQH': q @ q.T, 'commutator': g1 @ g2 - g2 @ g1,
                        'observation': x, 'transformed_observation': z, 'diagonal_images': separated_z,
                        'microphone_images': images, 'sum_images': images.sum(0)}

    audio = css_overlap_case()
    same = np.tile([1., -1., 1., -1.], (2, 1))
    result['E08-16'] = {'audio_parameters': audio['parameters'], 'wav_stems': list(audio['signals']),
                        'float_matching': audio['parameters']['matching'],
                        'floating_si_sdr_db': audio['parameters']['floating_si_sdr_db'],
                        'published_audio': css_published_audio(),
                        'silence': match_two_source_overlap(np.zeros((2, 4)), np.zeros((2, 4))),
                        'identical_slots': match_two_source_overlap(same, same)}
    result['E08-17'] = whitening_minicase()
    result['E08-18'] = congruence_model_boundary()
    result['E08-19'] = fastmnmf_scale_objective()
    result['E08-20'] = mixture_consistency_minicase()

    target, other = np.array([1., 1j]), np.array([-.8, 1.])
    observed = target + other
    ideal = target / observed
    bounded = np.clip((target * observed.conj()).real / abs(observed) ** 2, 0, 1)
    unit_mask_control = {'target': target, 'interference': other, 'mixture': observed,
                         'complex_oracle_mask': ideal, 'best_bounded_real_mask': bounded,
                         'bounded_output': bounded * observed,
                         'squared_errors': abs(bounded * observed - target) ** 2}
    from codes.chapters.ch08.core.mask_representation import run_experiment
    from codes.chapters.ch08.examples.mask_representation_demo import check_assets
    mask_report, _ = run_experiment()
    result['E08-21'] = {'unit_scale_complex_control': unit_mask_control,
                        'floating_point_experiment': mask_report,
                        'published_audio': check_assets()}

    from codes.chapters.ch08.core.consistency_teaching import run_experiment as consistency_experiment
    from codes.chapters.ch08.core.overiva_teaching import run_experiment as overiva_experiment
    result['E08-22'] = consistency_experiment()
    result['E08-23'] = overiva_experiment()

    from codes.chapters.ch08.core.css_audio import (FILE_NAMES, generate_experiment,
        parameters as css_parameters, analytic_measurements, measure_signal)
    from codes.chapters.ch08.examples.generate_css_audio import check_assets as check_css_assets
    css_control = generate_experiment()
    css_signals = css_control['signals']
    result['E08-24'] = {'four_point_discrete_control': css_polarity_minicase(),
                        'parameters': css_parameters(),
                        'overlap_matching': css_control['overlap_matching'],
                        'overlap_application_gain': css_control['overlap_application_gain'],
                        'analytic': {name: analytic_measurements(name) for name in FILE_NAMES},
                        'floating_point_measurements': {name: measure_signal(value, css_signals['reference'])
                                                        for name, value in css_signals.items()},
                        'published_audio': check_css_assets()}
    result['E08-25'] = ci_sdr_vs_si_sdr_experiment()
    result['E08-26'] = streaming_mac_and_bypass_experiment()
    result['E08-27'] = auxiva_iss_vs_ip_experiment()
    return _plain(result)


if __name__ == '__main__':
    print(json.dumps(run_experiments(), ensure_ascii=False, indent=2, allow_nan=False))
