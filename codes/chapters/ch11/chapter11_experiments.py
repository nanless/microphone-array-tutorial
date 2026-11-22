"""E11-10..22: computed selection examples, not product benchmark scores.

Run ``python -m codes.chapters.ch11.chapter11_experiments``. No network, model or
recording is used. Tables are constructed; E11-13/18 strictly check and read
published mathematical-tone WAVs. Missing or stale assets are never regenerated.
"""
from __future__ import annotations

# Allow the documented direct-file command as well as python -m.
if __name__ == "__main__" and not __package__:
    import sys as _chapter_entry_sys
    from pathlib import Path as _ChapterEntryPath
    _chapter_entry_sys.path.insert(0, str(_ChapterEntryPath(__file__).resolve().parents[3]))
import json
import math
from fractions import Fraction
import numpy as np
from codes.chapters.ch11.core.selection import (upper_limit_verdict, pareto_minima,
                                            small_slot_word_errors)
from codes.chapters.ch00.cross_chapter.tracking_time_exercises import zero_failure_upper_bound


def pareto_anchor():
    names = ['A', 'B', 'C', 'D']
    latency, wer, power = np.array([80., 100., 90., 160.]), np.array([.12, .10, .14, .08]), np.array([.8, .8, 1., .7])
    feasible = np.flatnonzero((latency <= 150) & (power <= 1))
    indices = pareto_minima(np.column_stack([latency[feasible], wer[feasible], power[feasible]]))
    latency_difference = (latency[1] - latency[0]) / 150
    wer_difference = (wer[0] - wer[1]) / .20
    rankings = {}
    for weight in (.5, .8):
        loss = weight * wer[feasible] / .2 + (1-weight) * latency[feasible] / 150
        rankings[str(weight)] = {'losses': loss.tolist(), 'winner': names[feasible[np.argmin(loss)]]}
    return {'names': names, 'latency_ms': latency.tolist(), 'wer': wer.tolist(), 'power_w': power.tolist(),
            'feasible': [names[i] for i in feasible], 'pareto': [names[feasible[i]] for i in indices],
            'wer_scale': .2, 'latency_scale_ms': 150., 'weight_crossing': latency_difference / (latency_difference + wer_difference),
            'weighted_rankings': rankings}


def deadline_acceptance_anchor():
    """One planned response per simulated session, including absent output."""
    deadline_ms = 150
    target_rate = Fraction(1, 100)
    tail_alpha = Fraction(1, 160)  # Two candidates, two scenarios, two tails.
    candidate_responses = {
        'A': {'normal': [100] * 600, 'stress': [100] * 17 + [180, None, None]},
        'B': {'normal': [125] * 600, 'stress': [125] * 19 + [150]},
    }
    candidates = {}
    for candidate, scenarios in candidate_responses.items():
        rows = {}
        delivered = []
        for scenario, responses in scenarios.items():
            observed = [value for value in responses if value is not None]
            delivered.extend(observed)
            n = len(responses)
            late = sum(value > deadline_ms for value in observed)
            absent = n - len(observed)
            violations = late + absent
            # At the required probability, large violation counts reject an
            # acceptable risk. Zero counts use the existing exact upper limit.
            upper = tail = None
            if violations == 0:
                upper = zero_failure_upper_bound(n, 1 - float(tail_alpha))
                verdict = 'pass' if upper <= target_rate else 'undetermined'
            else:
                tail = sum(math.comb(n, k) * target_rate**k * (1 - target_rate)**(n - k)
                           for k in range(violations, n + 1))
                verdict = 'fail' if tail < tail_alpha else 'undetermined'
            rows[scenario] = {'planned': n, 'on_time': n - violations,
                              'late': late, 'absent': absent, 'violations': violations,
                              'zero_count_upper': upper,
                              'upper_tail_at_limit': float(tail) if tail is not None else None,
                              'verdict': verdict}
        planned = sum(row['planned'] for row in rows.values())
        violations = sum(row['violations'] for row in rows.values())
        verdicts = {row['verdict'] for row in rows.values()}
        candidates[candidate] = {
            'scenarios': rows,
            'pooled_violations': violations,
            'pooled_planned': planned,
            'pooled_observed_rate': violations / planned,
            'delivered_only_p99_ms': sorted(delivered)[(99 * len(delivered) + 99) // 100 - 1],
            'verdict': 'fail' if 'fail' in verdicts else
            'undetermined' if 'undetermined' in verdicts else 'pass',
        }
    return {'constructed_data': True, 'deadline_ms': deadline_ms,
            'target_violation_probability': float(target_rate),
            'per_bound_alpha': float(tail_alpha), 'candidates': candidates}


def simultaneous_risk_anchor():
    scenarios, trials, alpha, target = 3, 100, .05, .01
    per_alpha = alpha / scenarios
    minimum = math.ceil(math.log(per_alpha) / math.log1p(-target))
    return {'scenarios': scenarios, 'zero_failures_each': 0, 'trials_each': trials,
            'joint_alpha': alpha, 'allocated_alpha_each': per_alpha,
            'separate_95_upper': zero_failure_upper_bound(trials),
            'simultaneous_upper_each': zero_failure_upper_bound(trials, 1-per_alpha),
            'target_upper': target, 'minimum_zero_failure_trials_each': minimum,
            'previous_n_upper': zero_failure_upper_bound(minimum-1, 1-per_alpha),
            'minimum_n_upper': zero_failure_upper_bound(minimum, 1-per_alpha),
            'deadline_acceptance': deadline_acceptance_anchor()}


def interaction_anchor():
    errors, words = {'none': 20, 'A': 22, 'B': 23, 'AB': 10}, 100
    rates = {k: v / words for k, v in errors.items()}
    difference = errors['AB'] - errors['A'] - errors['B'] + errors['none']
    return {'errors': errors, 'reference_words_each': words, 'wer': rates,
            'interaction_percentage_points': 100 * difference / words,
            'best_configuration': min(errors, key=errors.get),
            'constructed_scores_not_asr_run': True}


def audio_anchor():
    from codes.chapters.ch11.examples.generate_selection_audio import check_main_selection_assets
    return check_main_selection_assets()


def interval_dominance_anchor():
    from codes.chapters.ch11.core.selection import interval_dominates
    candidates = {'A': [[80, 90], [.10, .12]], 'B': [[100, 110], [.14, .16]],
                  'C': [[85, 95], [.11, .13]]}
    examples = [{'A': [90, .12], 'C': [85, .11]}, {'A': [80, .10], 'C': [95, .13]}]
    for example in examples:
        # Singletons are valid interval bounds on these exact selected points.
        av = [[v, v] for v in example['A']]
        cv = [[v, v] for v in example['C']]
        example['A_dominates_C'] = interval_dominates(av, cv)
        example['C_dominates_A'] = interval_dominates(cv, av)
    return {'intervals_latency_ms_wer': candidates, 'bounds_are_deterministic_rectangles': True,
            'A_dominates_B': interval_dominates(candidates['A'], candidates['B']),
            'A_dominates_C': interval_dominates(candidates['A'], candidates['C']),
            'C_dominates_A': interval_dominates(candidates['C'], candidates['A']),
            'point_counterexamples': examples, 'latency_upper_limit_ms': 92,
            'latency_verdicts': {name: upper_limit_verdict([rows[0]], [92])
                                 for name, rows in candidates.items()}}


def poisson_exposure_anchor():
    from codes.chapters.ch11.core.selection import zero_event_poisson_upper_bound
    exposure, confidence, target = 2., .95, .1
    alpha = 1 - confidence
    minimum = -math.log(alpha) / target
    return {'exposure_hours': exposure, 'observed_events': 0, 'confidence': confidence,
            'poisson_rate_upper_per_hour': zero_event_poisson_upper_bound(exposure, confidence),
            'target_rate_per_hour': target, 'minimum_fixed_exposure_hours': minimum,
            'zero_probability_at_upper': math.exp(-zero_event_poisson_upper_bound(exposure, confidence)*exposure),
            'scope': 'fixed exposure, homogeneous Poisson process; no KWS execution'}


def utterance_boundary_anchor():
    hypotheses = [['a'], ['b']]
    cases = {}
    for label, utterances in (('whole', [['a', 'b']]), ('split', [['a'], ['b']])):
        speakers = [0] * len(utterances)
        cases[label] = {'reference_utterances': utterances, 'speaker_labels': speakers,
                        **small_slot_word_errors(utterances, speakers, hypotheses)}
    return {'hypothesis_slots': hypotheses, 'cases': cases,
            'scope': 'same words; only reference utterance boundaries change; not ASR improvement'}


def time_edit_anchor():
    from codes.chapters.ch11.core.selection import token_edit_distance, tiny_time_constrained_edit_distance
    reference, hypothesis = ['a'], ['a']
    reference_times = [[0., 1.]]
    cases = []
    for label, times, collar in (('late', [[10., 11.]], 0.),
                                ('overlapping', [[.2, .8]], 0.),
                                ('wide_collar', [[10., 11.]], 10.)):
        error = tiny_time_constrained_edit_distance(reference, hypothesis, reference_times, times, collar)
        cases.append({'name': label, 'hypothesis_intervals_s': times, 'collar_s': collar,
                      'time_constrained_errors': error, 'time_constrained_wer': error / len(reference)})
    return {'reference_tokens': reference, 'hypothesis_tokens': hypothesis,
            'reference_intervals_s': reference_times,
            'ordinary_errors': token_edit_distance(reference, hypothesis), 'reference_words': len(reference),
            'cases': cases, 'scope': 'tiny unit-cost time edit, not complete tcpWER or original upstream run'}


def scenario_audio_anchor(directory=None):
    from codes.chapters.ch11.core.selection_audio import build_fixture, analytic_results, analyze_fixture
    from codes.chapters.ch11.examples.generate_selection_audio import check_assets
    fixture = build_fixture()
    manifest = check_assets() if directory is None else check_assets(directory)
    return {'parameters': fixture['parameters'], 'analytic': analytic_results(fixture),
            'floating_point': analyze_fixture(fixture), 'published_audio': manifest}


def scoring_events_anchor():
    from codes.chapters.ch11.core.selection import score_session_events, token_edit_distance
    sessions = [{'session_id': 's0', 'reference': ['a', 'b'], 'hypothesis': ['a', 'c']},
                {'session_id': 's1', 'reference': ['d', 'e'], 'hypothesis': []},
                {'session_id': 's2', 'reference': ['f', 'g', 'h']},
                {'session_id': 's3', 'reference': ['i'], 'hypothesis': 'i'},
                {'session_id': 's4', 'reference': ['j', 'k'], 'hypothesis': ['j', 'k']}]
    def controlled_scorer(reference, hypothesis):
        if reference == ['j', 'k']:
            raise RuntimeError('controlled teaching scorer failure')
        return token_edit_distance(reference, hypothesis)
    return {'planned_inputs': sessions, **score_session_events(sessions, scorer=controlled_scorer),
            'scope': 'constructed tokens and an actual controlled callback exception; no ASR or production scorer run'}


def physical_selection_anchor():
    from codes.chapters.ch11.core.selection_physics import beamformer_selection_case
    return beamformer_selection_case()


def multizone_cabin_anchor():
    from codes.chapters.ch11.core.automotive_selection import multizone_cabin_selection_case
    return multizone_cabin_selection_case()


def icc_certification_anchor():
    from codes.chapters.ch11.core.automotive_selection import icc_and_certification_selection_case
    return icc_and_certification_selection_case()


def run_experiments(*, scenario_directory=None):
    functions = [pareto_anchor, simultaneous_risk_anchor, interaction_anchor, audio_anchor,
                 interval_dominance_anchor, poisson_exposure_anchor,
                 utterance_boundary_anchor, time_edit_anchor,
                 lambda: scenario_audio_anchor(scenario_directory), scoring_events_anchor,
                 physical_selection_anchor, multizone_cabin_anchor, icc_certification_anchor]
    return {f'E11-{i}': function() for i, function in enumerate(functions, 10)}


if __name__ == '__main__':
    print(json.dumps(run_experiments(), ensure_ascii=False, indent=2, allow_nan=False))
