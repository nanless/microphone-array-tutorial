"""Teaching experiments E09-10..23; no writes on import/run.

These are small arithmetic/controlled-signal examples, not complete EKF, UKF,
PHD, GOSPA tracking systems or device evaluations. No upstream imports.
Run: python -m codes.chapters.ch09.chapter09_experiments
"""
from __future__ import annotations

# Allow the documented direct-file command as well as python -m.
if __name__ == "__main__" and not __package__:
    import sys as _chapter_entry_sys
    from pathlib import Path as _ChapterEntryPath
    _chapter_entry_sys.path.insert(0, str(_ChapterEntryPath(__file__).resolve().parents[3]))

import json
import math
from pathlib import Path
import numpy as np
from codes.chapters.ch02.core.conventions import finite_real_array, finite_real_scalar
from codes.chapters.ch09.core.echo_av_fusion import (
    small_set_distances,
    av_oosm_and_mot_experiment,
    echo_clutter_gating_experiment,
)
from codes.chapters.ch09.core.tracking import CircularParticleFilter, wrap_angle
from codes.chapters.ch09.core.tracking_audio import analyze_array, read_pcm16
from codes.chapters.ch09.examples.chapter09_tracking_audio import OUTPUT as TRACKING_OUTPUT, generate
from codes.chapters.ch09.core.imm_teaching import imm_missing_observation_example
from codes.chapters.ch09.examples.tracking_lifecycle_demo import run_demo as lifecycle_demo
from codes.chapters.ch09.core.spherical_tracking import (
    direction_from_angles, angles_from_direction, rotate_direction,
    normalized_direction_covariance, vmf_static_update, vmf_resultant_length,
)


def jpda_vs_gnn_crossing_case(
    predicted_bearings_deg,
    detections_deg,
    *,
    prior_variance_deg2: float = 2.0,
    measurement_variance_deg2: float = 2.0,
) -> dict:
    """Exact 2-track, 2-detection JPDA vs hard GNN update with P_D=1 and zero clutter."""
    pred = finite_real_array(predicted_bearings_deg, 'predicted_bearings_deg')
    meas = finite_real_array(detections_deg, 'detections_deg')
    p_minus = finite_real_scalar(prior_variance_deg2, 'prior_variance_deg2')
    r_meas = finite_real_scalar(measurement_variance_deg2, 'measurement_variance_deg2')
    if pred.shape != (2,) or meas.shape != (2,):
        raise ValueError('predicted_bearings_deg and detections_deg must each have length 2')
    if p_minus <= 0.0 or r_meas <= 0.0:
        raise ValueError('prior and measurement variances must be strictly positive')

    s_var = p_minus + r_meas
    k_gain = p_minus / s_var
    p_cond = (1.0 - k_gain) * p_minus

    # Row j indexes detection z_j, column t indexes track theta_t^-.
    nu = meas[:, None] - pred[None, :]
    norm_factor = 1.0 / math.sqrt(2.0 * math.pi * s_var)
    likelihood = norm_factor * np.exp(-0.5 * (nu ** 2) / s_var)

    # Joint events under P_D = 1, no clutter:
    # Theta_1 (straight): (1->1, 2->2); Theta_2 (swapped): (2->1, 1->2).
    e1 = 0.5 * float(nu[0, 0] ** 2 + nu[1, 1] ** 2) / s_var
    e2 = 0.5 * float(nu[1, 0] ** 2 + nu[0, 1] ** 2) / s_var
    w1_kernel = math.exp(-e1)
    w2_kernel = math.exp(-e2)
    log_odds = e2 - e1
    # Log-domain normalization avoids overflow for arbitrary separations.
    if log_odds >= 0.0:
        p_straight = 1.0 / (1.0 + math.exp(-log_odds))
        p_swapped = math.exp(-log_odds) * p_straight
    else:
        p_swapped = 1.0 / (1.0 + math.exp(log_odds))
        p_straight = math.exp(log_odds) * p_swapped

    beta = np.array([[p_straight, p_swapped], [p_swapped, p_straight]], dtype=float)

    # Hard GNN picks the maximum-likelihood 1-to-1 permutation.
    if p_straight >= p_swapped:
        gnn_event = 'straight'
        gnn_innovations = np.array([nu[0, 0], nu[1, 1]], dtype=float)
        gnn_assigned = np.array([meas[0], meas[1]], dtype=float)
    else:
        gnn_event = 'swapped'
        gnn_innovations = np.array([nu[1, 0], nu[0, 1]], dtype=float)
        gnn_assigned = np.array([meas[1], meas[0]], dtype=float)
    gnn_means = pred + k_gain * gnn_innovations
    gnn_vars = np.array([p_cond, p_cond], dtype=float)

    # JPDA combined innovation and spread-of-the-means covariance inflation.
    nu_bar = np.sum(beta * nu, axis=0)
    jpda_means = pred + k_gain * nu_bar
    spread_vars = (k_gain ** 2) * (np.sum(beta * (nu ** 2), axis=0) - nu_bar ** 2)
    jpda_vars = p_cond + spread_vars

    return {
        'predicted_bearings_deg': pred,
        'detections_deg': meas,
        'innovation_variance_deg2': s_var,
        'kalman_gain': k_gain,
        'conditional_posterior_variance_deg2': p_cond,
        'innovations_deg': nu,
        'likelihood_matrix_per_deg': likelihood,
        'straight_event_kernel': w1_kernel,
        'swapped_event_kernel': w2_kernel,
        'log_odds_straight_vs_swapped': log_odds,
        'joint_event_probabilities': {'straight': p_straight, 'swapped': p_swapped},
        'marginal_beta': beta,
        'gnn': {
            'selected_event': gnn_event,
            'assigned_detections_deg': gnn_assigned,
            'posterior_means_deg': gnn_means,
            'posterior_variances_deg2': gnn_vars,
            'posterior_separation_deg': float(abs(gnn_means[1] - gnn_means[0])),
        },
        'jpda': {
            'combined_innovations_deg': nu_bar,
            'posterior_means_deg': jpda_means,
            'spread_of_means_variance_deg2': spread_vars,
            'posterior_variances_deg2': jpda_vars,
            'posterior_separation_deg': float(abs(jpda_means[1] - jpda_means[0])),
        },
    }


def jpda_vs_gnn_crossing_experiment() -> dict:
    """Compare JPDA marginal weights & covariance inflation vs hard GNN on symmetric crossings."""
    return {
        'prior_variance_deg2': 2.0,
        'measurement_variance_deg2': 2.0,
        'innovation_variance_deg2': 4.0,
        'kalman_gain': 0.5,
        'detection_probability': 1.0,
        'clutter_intensity': 0.0,
        'separated_approach': jpda_vs_gnn_crossing_case([-2.0, 2.0], [-1.5, 1.5]),
        'ambiguous_crossing': jpda_vs_gnn_crossing_case([-0.5, 0.5], [-0.4, 0.4]),
        'scope': (
            '1-D scalar-azimuth 2-track, 2-detection symmetric crossing with P_D=1 and zero clutter; '
            'exact 2-event JPDA vs hard GNN single-frame update, not full multi-frame MHT or clutter gating'
        ),
    }


def _plain(value):
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def actual_tracking_audio(audio_directory=TRACKING_OUTPUT):
    """Strictly check published assets, then decode/analyze the actual WAV.

    generate(check=True) may replay in memory to validate the source contract;
    it never repairs files. The PCM result below is a new analysis of bytes
    read from the checked directory, not that in-memory replay's result.
    """
    directory = Path(audio_directory)
    manifest = generate(directory, check=True)
    array_bytes = (directory/'array_noisy.wav').read_bytes()
    analysis = analyze_array(read_pcm16(array_bytes), export_gain=manifest['common_export_gain'])
    return {'audio_directory': str(directory), 'asset_check': 'strict_read_only_passed',
        'files': manifest['files'], 'source_sha256': manifest['source_sha256'],
        'float_scores': manifest['float_analysis']['scores'], 'pcm_scores': analysis['scores'],
        'common_export_gain': manifest['common_export_gain'], 'config': manifest['analysis_config'],
        'first_state_time_s': analysis['frames']['state_time_s'][0],
        'first_available_time_s': analysis['frames']['available_time_s'][0],
        'missing_state_times_s': [t for t, valid in zip(analysis['frames']['state_time_s'],
                                                     analysis['frames']['observation_valid']) if not valid]}


def run_experiments(*, audio_directory=TRACKING_OUTPUT):
    result = {}

    # E09-10: Bearing-only EKF one-step update (book clockwise atan2(x, y)).
    position, prior = np.array([1., 2.]), np.eye(2)*.25
    observation_variance, innovation = .01, .1
    jacobian = np.array([position[1], -position[0]])/(position@position)
    innovation_variance = float(jacobian@prior@jacobian+observation_variance)
    gain = prior@jacobian/innovation_variance
    residual_map = np.eye(2)-np.outer(gain, jacobian)
    result['E09-10'] = {'position_m': position, 'prior_covariance_m2': prior,
        'bearing_rad': math.atan2(*position), 'innovation_rad': innovation,
        'jacobian_rad_per_m': jacobian, 'innovation_variance_rad2': innovation_variance,
        'gain_m_per_rad': gain, 'updated_position_m': position+gain*innovation,
        'posterior_covariance_m2': residual_map@prior@residual_map.T+observation_variance*np.outer(gain, gain)}

    # E09-11: Unscented transform vs EKF on h(x) = x^2.
    points, wm, wc = np.array([0., 1., -1.]), np.array([0., .5, .5]), np.array([2., .5, .5])
    transformed = points**2
    mean = wm@transformed
    result['E09-11'] = {'alpha': 1., 'beta': 2., 'kappa': 0., 'sigma_points': points,
        'mean_weights': wm, 'covariance_weights': wc, 'transformed': transformed,
        'mean': mean, 'variance': wc@((transformed-mean)**2),
        'state_observation_cross_covariance': wc@(points*(transformed-mean)),
        'boundary': 'squared observation cannot distinguish x from -x; cross covariance and UKF gain are zero'}

    # E09-12: Circular mean, antipodal ambiguity, and log-weight prior support.
    particles = CircularParticleFilter([179., -179.])
    undefined = CircularParticleFilter([90., -90.])
    undefined_angle = None
    status = 'not_checked'
    try:
        undefined_angle = undefined.estimate()
        status = 'unexpected_defined_antipodal_mean'
    except ValueError:
        status = 'undefined_antipodal_mean'
    supported = CircularParticleFilter([0., 10.])
    supported.update(10., .01, clutter_probability=0.)
    first = supported.weights.copy()
    first_log = supported.log_weights
    underflow = supported.weight_underflow_indices
    supported.update(0., .01, clutter_probability=0.)
    literal_zero = CircularParticleFilter([0., 10.])
    literal_zero.set_prior_weights([0., 1.])
    literal_zero.update(0., .01, clutter_probability=0.)
    result['E09-12'] = {'particles_deg': [179., -179.], 'linear_mean_deg': 0.,
        'circular_mean_deg': particles.estimate(), 'resultant_length': math.cos(math.radians(1)),
        'antipodal_status': status, 'antipodal_angle_deg': undefined_angle,
        'support_after_first': first, 'support_after_second': supported.weights,
        'log_weights_after_first': first_log, 'linear_underflow_indices_after_first': underflow,
        'log_weights_after_second': supported.log_weights,
        'explicit_zero_prior_after_update': literal_zero.weights,
        'prior_reset_method': 'set_prior_weights',
        'support_note': 'persistent finite log weights retain positive support; set_prior_weights explicitly resets support, including same-value displayed weights'}

    # E09-13 (merged old E09-14 + E09-22): PHD intensity vs cardinality variance & Bernoulli empty observation.
    distributions = {
        'A': (np.array([[1., 0.], [0., 1.]]), np.array([.5, .5])),
        'B': (np.array([[0., 0.], [1., 1.]]), np.array([.5, .5])),
    }
    phd_result = {}
    for name, (sets, probabilities) in distributions.items():
        counts = sets.sum(axis=1)
        mean_count = probabilities@counts
        phd_result[f'model_{name}_grid_mass'] = probabilities@sets
        phd_result[f'model_{name}_count_probability_0_1_2'] = [float(probabilities[counts == n].sum()) for n in range(3)]
        phd_result[f'model_{name}_count_variance'] = float(probabilities@((counts-mean_count)**2))
    mass = phd_result['model_A_grid_mass']
    detection_probability = .8
    remaining = (1-detection_probability)*mass
    prior_existence = .8
    bernoulli_cases = []
    for p_d in (.9, .1):
        empty_probability = 1-prior_existence*p_d
        bernoulli_cases.append({'detection_probability': p_d,
            'empty_observation_probability': empty_probability,
            'bernoulli_posterior_existence': prior_existence*(1-p_d)/empty_probability,
            'different_poisson_phd_posterior_mass': prior_existence*(1-p_d)})
    result['E09-13'] = {**phd_result, 'grid_mass': mass, 'expected_count': float(mass.sum()),
        'detection_probability': detection_probability, 'empty_observation_posterior_mass': remaining,
        'empty_observation_expected_count': float(remaining.sum()),
        'bernoulli_empty_observation': {'prior_existence': prior_existence, 'cases': bernoulli_cases,
            'scope': 'at most one target, no clutter/birth/survival transition in this update, known constant detection probability; not VAD or identity'}}

    # E09-14 (old E09-15): OSPA vs GOSPA normalization & p<1 boundary.
    result['E09-14'] = {'original': small_set_distances([32], [30, 70]),
        'added_common_correct_target': small_set_distances([32, 100], [30, 70, 100]),
        'both_empty': small_set_distances([], []), 'one_empty': small_set_distances([], [30, 70]),
        'identity_note': 'permutation is optimized independently; equal sets do not establish track identity'}

    # E09-15 (merged old E09-16 + E09-23): Fixed-lag RTS smoothing & OOSM late-observation replay.
    initial_mean, initial_variance, process_variance, observation_variance = 0., 1., 1., 1.
    observations = [1., 0.]
    mean, variance = initial_mean, initial_variance
    filtered_means, filtered_variances, predicted_variances = [], [], []
    for observation in observations:
        predicted_variance = variance+process_variance
        gain = predicted_variance/(predicted_variance+observation_variance)
        mean = mean+gain*(observation-mean)
        variance = (1-gain)**2*predicted_variance+gain**2*observation_variance
        predicted_variances.append(predicted_variance)
        filtered_means.append(mean)
        filtered_variances.append(variance)
    m1, m2 = filtered_means
    p1, p2 = filtered_variances
    predicted_p2 = predicted_variances[1]
    smoothing_gain = p1/predicted_p2
    observation_covariance = np.array([[3., 2.], [2., 4.]])
    final_state_observation_cross = np.array([2., 3.])
    oosm_joint_gain = np.linalg.solve(observation_covariance, final_state_observation_cross)
    replay_mean, replay_variance = 0., 1.
    replay = []
    for observation in (2., -1.):
        prior_var = replay_variance+1
        k_gain = prior_var/(prior_var+1)
        replay_mean += k_gain*(observation-replay_mean)
        replay_variance = (1-k_gain)*prior_var
        replay.append({'mean': replay_mean, 'variance': replay_variance})
    drop_mean, drop_variance = -.75, .75
    wrong_gain = drop_variance/(drop_variance+1)
    result['E09-15'] = {'initial_mean': initial_mean, 'initial_variance': initial_variance,
        'Q': process_variance, 'R': observation_variance, 'observations': observations,
        'filtered_means': [m1, m2], 'filtered_variances': [p1, p2],
        'smoothing_gain': smoothing_gain, 'smoothed_time1_mean': m1+smoothing_gain*(m2-m1),
        'smoothed_time1_variance': p1+smoothing_gain**2*(p2-predicted_p2),
        'target_time_step': 1, 'latest_used_observation_step': 2, 'extra_observation_delay_steps': 1,
        'oosm_late_observation': {'P0': 1., 'Q': 1., 'R': 1.,
            'observations_by_measurement_time': [2., -1.],
            'arrival_order_measurement_times': [2, 1], 'final_state_time': 2,
            'joint_observation_covariance': observation_covariance,
            'final_state_observation_cross_covariance': final_state_observation_cross,
            'joint_gain': oosm_joint_gain, 'joint_mean': float(oosm_joint_gain@np.array([2., -1.])),
            'joint_variance': float(3-oosm_joint_gain@final_state_observation_cross),
            'chronological_replay': replay,
            'drop_late_observation': {'mean': drop_mean, 'variance': drop_variance},
            'wrong_current_time_update': {'mean': drop_mean+wrong_gain*(2-drop_mean),
                                          'variance': (1-wrong_gain)*drop_variance},
            'scope': 'scalar Euclidean random walk; chronological replay uses original measurement times once, not a generic out-of-sequence tracker'}}

    # E09-16 (old E09-19): Published dual-mic PCM tracking_audio timing & denominators.
    result['E09-16'] = actual_tracking_audio(audio_directory)

    # E09-17 (merged old E09-20 + E09-21): Fig 56 information boundaries (correlated obs + bearing scale nullspace).
    observation_map = np.ones((2, 1))
    correlated_noise = np.array([[1., .9], [.9, 1.]])
    corr_innovation_cov = observation_map@observation_map.T+correlated_noise
    corr_joint_gain = np.linalg.solve(corr_innovation_cov, observation_map).T
    independent_covariance = observation_map@observation_map.T+np.eye(2)
    independent_gain = np.linalg.solve(independent_covariance, observation_map).T
    correlated_obs = {'prior_mean': 0., 'prior_variance': 1., 'observations': [1., 1.],
        'measurement_noise_covariance': correlated_noise, 'innovation_covariance': corr_innovation_cov,
        'joint_gain': corr_joint_gain[0], 'joint_mean': float((corr_joint_gain@np.ones(2))[0]),
        'joint_variance': float((1-corr_joint_gain@observation_map)[0, 0]),
        'wrong_independent_mean': float((independent_gain@np.ones(2))[0]),
        'wrong_independent_variance': float((1-independent_gain@observation_map)[0, 0]),
        'perfect_correlation_limit': {'mean': .5, 'variance': .5},
        'scope': 'static state and correlated zero-mean Gaussian noise; rho=1 is a singular limit'}
    initial = np.array([1., 2.]); velocity = np.array([1., 0.]); times = np.arange(3.)
    positions = initial+times[:, None]*velocity
    rows = np.column_stack((positions[:, 1], -positions[:, 0],
                            times*positions[:, 1], -times*positions[:, 0]))
    scale_null = np.r_[initial, velocity]
    scale_ambiguity = {'array_center_m': [0., 0.], 'initial_position_m': initial,
        'velocity_m_s': velocity, 'scale_factor': 2., 'times_s': times,
        'positions_m': positions, 'scaled_positions_m': 2*positions,
        'bearings_deg': np.rad2deg(np.arctan2(positions[:, 0], positions[:, 1])),
        'scaled_bearings_deg': np.rad2deg(np.arctan2(2*positions[:, 0], 2*positions[:, 1])),
        'ranges_m': np.hypot(positions[:, 0], positions[:, 1]),
        'scaled_ranges_m': 2*np.hypot(positions[:, 0], positions[:, 1]),
        'stacked_jacobian_without_positive_row_denominators': rows,
        'stacked_rank': int(np.linalg.matrix_rank(rows)), 'scale_null_direction': scale_null,
        'stacked_null_product': rows@scale_null,
        'scope': 'fixed array, unknown constant Cartesian velocity; instantaneous geometric bearing at source-state time, propagation delay ignored; not equivalence of receiver-clock retarded audio'}
    result['E09-17'] = {'correlated_observations': correlated_obs, 'bearing_scale_ambiguity': scale_ambiguity}

    # E09-18 (old E09-24): Full log-density IMM & consecutive missing-observation coasting.
    result['E09-18'] = imm_missing_observation_example()

    # E09-19 (old E09-25): Dual-mic PCM integer half-sample clock lifecycle (Fig72).
    result['E09-19'] = lifecycle_demo(audio_directory)

    # E09-20 (old E09-26): 3D attitude rotation, tangent covariance projection & spherical vMF update.
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    array_direction = direction_from_angles(30., 0.)
    world_direction = rotate_direction(array_direction, rotation)
    sphere_update = vmf_static_update([0., 1., 0.], 3., [1., 0., 0.], 4.)
    result['E09-20'] = {
        'array_angles_deg': [30., 0.], 'array_direction': array_direction,
        'array_to_world': rotation, 'world_direction': world_direction,
        'world_angles': angles_from_direction(world_direction),
        'north_pole_from_two_azimuths': [direction_from_angles(0., 90.), direction_from_angles(120., 90.)],
        'north_pole_angles': angles_from_direction([0., 0., 1.]),
        'local_normalization': normalized_direction_covariance([0., 2., 0.], np.eye(3)*.04),
        'static_vmf_posterior': sphere_update,
        'posterior_angles': angles_from_direction(sphere_update['mean_direction']),
        'posterior_resultant_length': vmf_resultant_length(5.),
        'uniform_polar_cap_10deg_probability': (1-math.cos(math.radians(10)))/2,
        'equal_theta_phi_grid_density_ratio_at_80deg': 1/math.cos(math.radians(80)),
        'scope': 'known proper rotation and one static independent vMF likelihood; not dynamic FvMFF, pose estimation or a calibrated tangent filter',
    }

    # E09-21 (old E09-27): Loudspeaker residual-echo clutter sector gating & double-talk release.
    result['E09-21'] = echo_clutter_gating_experiment()

    # E09-22 (old E09-28): Audio-visual OOSM retrodiction, FoV blind-zone & CLEAR MOT vs OSPA.
    result['E09-22'] = av_oosm_and_mot_experiment()

    # E09-23: JPDA marginal association weights & covariance inflation vs hard GNN crossing.
    result['E09-23'] = jpda_vs_gnn_crossing_experiment()

    return _plain(result)


if __name__ == '__main__':
    print(json.dumps(run_experiments(), ensure_ascii=False, indent=2, allow_nan=False))
