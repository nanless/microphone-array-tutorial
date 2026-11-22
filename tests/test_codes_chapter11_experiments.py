"""Independent arithmetic, small counterexamples and actual PCM for chapter 11."""
import math
import json
import shutil
import struct
from fractions import Fraction
from itertools import product
from pathlib import Path
import tempfile
import unittest
import wave
import numpy as np
from codes.chapters.ch11.core.selection import (upper_limit_verdict, pareto_minima,
                                            small_slot_word_errors, token_edit_distance)
from codes.chapters.ch11.core.automotive_selection import bus_clock_drift_metrics, icc_comb_and_howling_metrics
from codes.chapters.ch00.core.audio_samples import selection_tradeoff_case
from codes.chapters.ch11.chapter11_experiments import run_experiments, audio_anchor, scenario_audio_anchor
from codes.chapters.ch11.examples.generate_selection_audio import generate, check_main_selection_assets
from codes.chapters.ch00.cross_chapter.engineering_boundary_exercises import paired_sign_test_lower_is_better
from codes.chapters.ch00.cross_chapter.tracking_time_exercises import zero_failure_upper_bound


class Chapter11ExperimentsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Unit fixtures are generated only in a temporary directory. The main
        # E11-13 branch always reads the already published four repository WAVs.
        cls.temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.scenario_directory = Path(cls.temporary.name) / 'scenarios'
        generate(cls.scenario_directory)
        cls.results = run_experiments(scenario_directory=cls.scenario_directory)
        cls.audio = selection_tradeoff_case()

    def test_identifiers_and_pareto_analytic_intersection(self):
        self.assertEqual(list(self.results), [f'E11-{i}' for i in range(10, 23)])
        x = self.results['E11-10']
        self.assertEqual(x['feasible'], ['A', 'B', 'C'])
        self.assertEqual(x['pareto'], ['A', 'B'])
        self.assertAlmostEqual(x['weight_crossing'], 4/7)
        self.assertEqual([x['weighted_rankings'][k]['winner'] for k in ['0.5', '0.8']], ['A', 'B'])
        np.testing.assert_allclose(x['weighted_rankings']['0.5']['losses'], [17/30, 7/12, 13/20])

    def test_pareto_equal_candidates_remain_and_invalid_rejected(self):
        self.assertEqual(pareto_minima([[1, 2], [1, 2], [2, 2], [3, 0]]), [0, 1, 3])
        for invalid in ([], [[np.nan]], [[1j]], [[True]]):
            with self.assertRaises(ValueError):
                pareto_minima(invalid)

    def test_selection_keeps_exact_integer_order(self):
        for lower in (2**53, int(np.iinfo(np.int64).max)-1, 10**400):
            with self.subTest(lower=lower):
                self.assertEqual(pareto_minima([[lower], [lower+1]]), [0])
                self.assertEqual(upper_limit_verdict([[lower+1, lower+1]], [lower]), 'fail')
                self.assertEqual(upper_limit_verdict([[lower, lower+1]], [lower]), 'undetermined')
                self.assertEqual(upper_limit_verdict([[lower, lower]], [lower]), 'pass')
                with self.assertRaises(ValueError):
                    upper_limit_verdict([[lower+1, lower]], [lower])
        # Both mixed columns and mixed endpoints preserve the integer operand.
        self.assertEqual(pareto_minima([[float(2**53), 1.5], [2**53+1, 1.5]]), [0])
        self.assertEqual(upper_limit_verdict([[2**53+1, 2**53+1]], [float(2**53)]), 'fail')
        self.assertEqual(upper_limit_verdict([[float(2**53), 2**53+1]], [2**53]), 'undetermined')
        for invalid in (True, 1j, float('nan'), float('inf'), '2'):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    pareto_minima([[invalid]])
                with self.assertRaises(ValueError):
                    upper_limit_verdict([[invalid, invalid]], [1])

    def test_joint_risk_coverage_uses_bonferroni_not_independence(self):
        x = self.results['E11-11']
        upper = x['simultaneous_upper_each']
        self.assertAlmostEqual((1-upper)**100, 1/60)
        self.assertEqual(x['minimum_zero_failure_trials_each'], 408)
        self.assertGreater(.99**407, 1/60)
        self.assertLessEqual(.99**408, 1/60)

    def test_deadline_counts_missing_and_late_responses_in_pressure_scenario(self):
        result = self.results['E11-11']['deadline_acceptance']
        a, b = (result['candidates'][candidate] for candidate in ('A', 'B'))
        a_stress, b_stress = a['scenarios']['stress'], b['scenarios']['stress']
        self.assertEqual((a_stress['on_time'], a_stress['late'], a_stress['absent']), (17, 1, 2))
        self.assertEqual(b_stress['on_time'], 20)  # A delivery at 150 ms is on time.
        self.assertEqual((a['pooled_violations'], a['pooled_planned']), (3, 620))
        self.assertEqual((a['delivered_only_p99_ms'], b['delivered_only_p99_ms']), (100, 125))
        self.assertEqual((a['verdict'], b['verdict']), ('fail', 'undetermined'))
        self.assertAlmostEqual((1 - a['scenarios']['normal']['zero_count_upper']) ** 600, 1/160)
        self.assertAlmostEqual((1 - b_stress['zero_count_upper']) ** 20, 1/160)
        # Independent lower-tail complement at the specification limit p=1%.
        p, q = Fraction(1, 100), Fraction(99, 100)
        complement = 1 - q**20 - 20*p*q**19 - 190*p*p*q**18
        self.assertAlmostEqual(a_stress['upper_tail_at_limit'], float(complement))
        self.assertLess(complement, Fraction(1, 160))

    def test_interval_verdict_boundaries_and_missingness(self):
        self.assertEqual(upper_limit_verdict([[150, 150]], [150]), 'pass')
        self.assertEqual(upper_limit_verdict([[150, 151]], [150]), 'undetermined')
        self.assertEqual(upper_limit_verdict([None, [2, 3]], [1, 1]), 'fail')
        for bounds, limits in (([[2, 1]], [3]), ([[0, np.inf]], [3]), ([], []), ([None], [])):
            with self.assertRaises(ValueError):
                upper_limit_verdict(bounds, limits)

    def test_interaction_is_computed_from_four_cells(self):
        x = self.results['E11-12']
        self.assertEqual(x['interaction_percentage_points'], -15)
        self.assertEqual(x['best_configuration'], 'AB')
        # At A absent, B costs +3; at A present, B changes errors by -12.
        self.assertEqual((x['errors']['B']-x['errors']['none'], x['errors']['AB']-x['errors']['A']), (3, -12))

    def test_slot_reuse_has_global_permutation_penalty(self):
        reordered = small_slot_word_errors([['a'], ['b'], ['c']], [0, 1, 2], [['b'], ['a', 'c']])
        self.assertEqual((reordered['cp_errors'], reordered['orc_errors']), (2, 0))
        identity = small_slot_word_errors([['a'], ['b'], ['c']], [0, 1, 0], [['a', 'c'], ['b']])
        self.assertEqual((identity['cp_errors'], identity['orc_errors']), (0, 0))

    def test_edit_distance_and_bounded_domain(self):
        self.assertEqual(token_edit_distance(['a', 'b', 'c'], ['a', 'x', 'c', 'd']), 2)
        self.assertEqual(token_edit_distance([], ['a', 'b']), 2)
        with self.assertRaises(ValueError):
            small_slot_word_errors([['a']]*7, [0]*7, [['a']])
        with self.assertRaises(ValueError):
            small_slot_word_errors([[]], [0], [['a']])

    def test_sign_test_preserves_big_integer_and_mixed_comparisons(self):
        for baseline, candidate in (([2**53+1], [2**53]),
                                    ([np.iinfo(np.int64).max], [np.iinfo(np.int64).max-1]),
                                    ([2**53+1, 1.5], [float(2**53), 1.5]),
                                    ([10**400], [10**400-1])):
            result = paired_sign_test_lower_is_better(baseline, candidate)
            self.assertEqual(result['wins'], 1)
            self.assertEqual(result['one_sided_pvalue'], .5)

    def test_sign_tail_reports_positive_underflow(self):
        x = paired_sign_test_lower_is_better([1]*1074, [0]*1074)
        self.assertEqual(x['one_sided_pvalue'], float.fromhex('0x0.0000000000001p-1022'))
        self.assertFalse(x['pvalue_underflow'])
        y = paired_sign_test_lower_is_better([1]*2000, [0]*2000)
        self.assertEqual(y['one_sided_pvalue'], 0)
        self.assertTrue(y['pvalue_underflow'])
        self.assertEqual(y['status'], 'positive_pvalue_underflow')
        self.assertAlmostEqual(y['log_one_sided_pvalue'], -2000*math.log(2))
        z = paired_sign_test_lower_is_better([1], [1])
        self.assertIsNone(z['log_one_sided_pvalue'])

    def test_zero_failure_extreme_integer_is_controlled(self):
        self.assertAlmostEqual(zero_failure_upper_bound(298), 0.0100024317, places=7)
        self.assertAlmostEqual(zero_failure_upper_bound(299), 0.0099691454, places=7)
        self.assertGreater(zero_failure_upper_bound(10**300), 0)
        with self.assertRaisesRegex(ValueError, 'range'):
            zero_failure_upper_bound(10**400)

    def test_filter_response_from_boxcar_trigonometric_sum(self):
        for length in (3, 9):
            name = f'selection_fir{length}'
            def response(f):
                return math.sin(length*math.pi*f/16000)/(length*math.sin(math.pi*f/16000))
            expected = abs(response(1500)/response(500))
            self.assertAlmostEqual(self.audio['analytic_response'][name]['amplitude_response']['1500'], expected)
            gain_noise = abs(response(3500)/response(500))
            # Equal 500/1500 reference energies, 500 response exactly one.
            expected_nmse = ((expected-1)**2 + gain_noise**2)/2
            self.assertAlmostEqual(self.audio['float_analysis']['candidates'][name]['aligned_total_nmse'], expected_nmse)

    def test_causal_convolution_startup_by_explicit_sample_sum(self):
        signals = self.audio['signals']
        for length in (3, 9):
            name = f'selection_fir{length}'
            coefficient = self.audio['parameters']['filters'][name]['taps'][0]
            for n in [0, 1, 7, 9, 318, 1600, 31999]:
                expected = coefficient * sum(signals['selection_mixture'][0, n-k]
                                             for k in range(min(length, n+1)))
                self.assertAlmostEqual(signals[name][0, n], expected, places=15)

    def test_pcm_readback_uses_pcm_reference_and_independent_fourier_projection(self):
        # Read repository bytes, not fresh encoder output. Integer sums and a
        # Fourier projection are independent of the entry's seven-column LS.
        directory = Path(__file__).resolve().parents[1] / 'codes/chapters/ch11/audio'
        values, integers = {}, {}
        for name in self.audio['signals']:
            with wave.open(str(directory/(name+'.wav')), 'rb') as wav:
                self.assertEqual((wav.getnchannels(), wav.getframerate(), wav.getsampwidth(), wav.getnframes()),
                                 (1, 16000, 2, 32000))
                integers[name] = struct.unpack('<32000h', wav.readframes(32000))
            values[name] = np.asarray(integers[name], dtype=float)/32768
        clean = values['selection_clean'][1600:30400]
        def amplitude(signal, f, start, stop):
            return 2*abs(np.dot(signal[start:stop], np.exp(-2j*np.pi*f*np.arange(start, stop)/16000)))/(stop-start)
        for length, expected_error in ((3, 91223542800), (9, 249713672400)):
            name, delay = f'selection_fir{length}', (length-1)//2
            output, start, stop = values[name], 1600+delay, 30400+delay
            actual = self.results['E11-13']['pcm_analysis']['candidates'][name]
            target = amplitude(output, 1500, start, stop)/amplitude(values['selection_clean'], 1500, 1600, 30400)
            noise = amplitude(output, 3500, start, stop)/amplitude(values['selection_mixture'], 3500, 1600, 30400)
            self.assertAlmostEqual(actual['target_1500_retention'], target, places=12)
            self.assertAlmostEqual(actual['noise_attenuation_db'], -20*math.log10(noise), places=10)
            denominator = sum(v*v for v in integers['selection_clean'][1600:30400])
            numerator = sum((integers[name][i+delay]-integers['selection_clean'][i])**2 for i in range(1600,30400))
            self.assertEqual((numerator,denominator), (expected_error,791595378000))
            self.assertEqual(self.results['E11-13']['integer_analysis'][name]['integer_error_squared_sum'], numerator)
            self.assertAlmostEqual(actual['aligned_total_nmse'], numerator/denominator, places=14)

    def test_interval_dominance_checks_all_vertices_and_counterexamples(self):
        x = self.results['E11-14']
        self.assertTrue(x['A_dominates_B'])
        self.assertFalse(x['A_dominates_C'])
        self.assertFalse(x['C_dominates_A'])
        for a in product((80,90),(Fraction(10,100),Fraction(12,100))):
            for b in product((100,110),(Fraction(14,100),Fraction(16,100))):
                self.assertTrue(all(av < bv for av,bv in zip(a,b)))
        self.assertEqual([(r['A_dominates_C'],r['C_dominates_A']) for r in x['point_counterexamples']],
                         [(False,True),(True,False)])
        self.assertEqual(x['latency_verdicts'],{'A':'pass','B':'fail','C':'undetermined'})

    def test_poisson_zero_events_uses_hours_not_trial_count(self):
        x = self.results['E11-15']
        rate = -math.log(.05)/2
        self.assertAlmostEqual(x['poisson_rate_upper_per_hour'],rate)
        self.assertAlmostEqual(math.exp(-2*rate),.05)
        self.assertAlmostEqual(x['minimum_fixed_exposure_hours'],-math.log(.05)/.1)
        self.assertGreater(x['poisson_rate_upper_per_hour'],.1)

    def test_whole_utterance_cannot_be_split_between_output_slots(self):
        x = self.results['E11-16']['cases']
        self.assertEqual((x['whole']['orc_errors'],x['split']['orc_errors']), (2,0))
        self.assertEqual((x['whole']['cp_errors'],x['split']['cp_errors']), (2,2))
        self.assertEqual((x['whole']['reference_words'],x['split']['reference_words']), (2,2))
        self.assertEqual(x['whole']['orc_wer'],1)

    def test_time_constraint_allows_error_rate_above_one(self):
        x = self.results['E11-17']
        self.assertEqual(x['ordinary_errors'],0)
        self.assertEqual([r['time_constrained_errors'] for r in x['cases']], [2,0,0])
        self.assertEqual(x['cases'][0]['time_constrained_wer'],2)
        from codes.chapters.ch11.core.selection import tiny_time_constrained_edit_distance
        # Enumerate paths (rather than reproduce the core's dynamic program).
        ref,hyp = ['a','b'],['a','x']
        rt,ht = [(0,1),(2,3)],[(.2,.8),(10,11)]
        def all_costs(i,j,cost):
            if i==len(ref) and j==len(hyp):
                yield cost
            if i<len(ref):yield from all_costs(i+1,j,cost+1)
            if j<len(hyp):yield from all_costs(i,j+1,cost+1)
            if i<len(ref) and j<len(hyp):
                overlap=max(rt[i][0],ht[j][0])<min(rt[i][1],ht[j][1])
                yield from all_costs(i+1,j+1,cost+((ref[i]!=hyp[j]) if overlap else 2))
        self.assertEqual(tiny_time_constrained_edit_distance(ref,hyp,rt,ht),min(all_costs(0,0,0)))

    def test_scenario_composition_has_distinct_robust_objectives(self):
        x = self.results['E11-18']
        def response(length,frequency):
            # Explicit finite cosine sum around the filter center.
            half=(length-1)//2
            return (1+2*sum(math.cos(2*math.pi*frequency*k/16000) for k in range(1,half+1)))
        expected={}
        for length in (3,9):
            a=response(length,1500)/response(length,500)
            v=response(length,3500)/response(length,500)
            single,dual=v*v,((a-1)**2+v*v)/2
            expected['fir'+str(length)]=(single,dual)
            for scene,value in [('single',single),('dual',dual)]:
                self.assertAlmostEqual(x['analytic']['candidates'][f'selection_{scene}_fir{length}.wav']['aligned_total_nmse'],value)
        for key,(single,dual) in expected.items():
            endpoints=[q*single+(1-q)*dual for q in (.25,.75)]
            self.assertAlmostEqual(x['analytic']['selection']['worst_over_q_interval'][key],max(endpoints))
            self.assertAlmostEqual(x['analytic']['selection']['worst_scene'][key],max(single,dual))
            self.assertNotEqual(max(endpoints),max(single,dual))

    def test_scenario_actual_pcm_integer_sums_and_missing_assets(self):
        manifest=self.results['E11-18']['published_audio']
        expected={'single':(395812398600,[87210066600,152494200]),
                  'dual':(791595378000,[91223542800,249713672400])}
        def read(stem):
            with wave.open(str(self.scenario_directory/(stem+'.wav')),'rb') as w:
                self.assertEqual((w.getframerate(),w.getnchannels(),w.getsampwidth(),w.getnframes()),(16000,1,2,32008))
                return struct.unpack('<32008h',w.readframes(32008))
        for scene,(den,errors) in expected.items():
            ref=read(f'selection_{scene}_target')
            self.assertEqual(sum(v*v for v in ref[1600:30400]),den)
            for length,error in zip((3,9),errors):
                stem=f'selection_{scene}_fir{length}';out=read(stem);delay=(length-1)//2
                observed=sum((out[i+delay]-ref[i])**2 for i in range(1600,30400))
                self.assertEqual(observed,error)
                row=manifest['pcm_analysis']['candidates'][stem+'.wav']
                self.assertEqual(row['integer_error_squared_sum'],observed)
                self.assertAlmostEqual(row['aligned_total_nmse'],float(Fraction(error,den)))
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):scenario_audio_anchor(d)
            self.assertEqual(list(Path(d).iterdir()),[])

    def test_main_selection_check_rejects_bad_pcm_nan_and_boolean_without_repair(self):
        root=Path(__file__).resolve().parents[1]
        metadata=json.loads((root/'codes/chapters/ch00/audio/MANIFEST.json').read_text())
        with tempfile.TemporaryDirectory() as directory:
            d=Path(directory)
            paths=[*metadata['generator_inputs'],'codes/chapters/ch00/audio/MANIFEST.json',
                   *[f'codes/chapters/ch11/audio/selection_{s}.wav' for s in ('clean','mixture','fir3','fir9')]]
            for path in paths:
                target=d/path;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(root/path,target)
            check_main_selection_assets(d)
            manifest=d/'codes/chapters/ch00/audio/MANIFEST.json';original=manifest.read_bytes()
            for invalid in (float('nan'),True):
                changed=json.loads(original);changed['groups']['selection_tradeoff']['common_export_gain']=invalid
                manifest.write_text(json.dumps(changed))
                before=(manifest.read_bytes(),manifest.stat().st_mtime_ns)
                with self.assertRaises(ValueError):check_main_selection_assets(d)
                self.assertEqual((manifest.read_bytes(),manifest.stat().st_mtime_ns),before)
            manifest.write_bytes(original)
            wav=d/'codes/chapters/ch11/audio/selection_fir3.wav'
            damaged=bytearray(wav.read_bytes());damaged[-1]^=1;wav.write_bytes(damaged)
            before=(wav.read_bytes(),wav.stat().st_mtime_ns)
            with self.assertRaises(ValueError):check_main_selection_assets(d)
            self.assertEqual((wav.read_bytes(),wav.stat().st_mtime_ns),before)

    def test_multizone_cabin_selection_independent_geometry_and_guard(self):
        res = self.results['E11-21']
        self.assertEqual(res['bus_verdicts'], {'a2b_sync': 'pass', 'bt_async_uncorrected': 'fail'})
        bt = res['bus_clocks']['bt_async_uncorrected']
        expected_drift_s = float(Fraction(60 * 50, 1_000_050))
        self.assertAlmostEqual(bt['exact_drift_s'], expected_drift_s, places=15)
        self.assertAlmostEqual(bt['max_allowed_drift_us'], 1e6 / 288_000, places=12)

        cands = res['candidates']
        # Independent rational checks on ceiling & headrest distances and weight norms.
        self.assertAlmostEqual(cands['ceiling_dsb']['driver_output_power'], float(Fraction(25, 164)), places=14)
        self.assertAlmostEqual(cands['ceiling_dsb']['weight_norm_squared'], float(Fraction(761, 1417)), places=14)
        self.assertAlmostEqual(cands['headrest_dsb']['weight_norm_squared'], float(Fraction(36, 61)), places=14)
        self.assertAlmostEqual(cands['headrest_dsb']['WNG_dB'], 10 * math.log10(61 / 36), places=12)
        self.assertLess(cands['headrest_lcmv']['WNG_dB'], -3.0)
        self.assertLess(cands['headrest_dsb_cross_cancel_unguarded']['LCR_dB'], 25.0)
        self.assertGreater(cands['headrest_dsb_cross_cancel_unguarded']['dt_target_loss_db'], 1.0)
        self.assertGreater(cands['headrest_dsb_cross_cancel_unguarded']['driver_only_self_cancel_loss_db'], 30.0)
        guarded = cands['headrest_dsb_cross_cancel_guarded']
        self.assertGreaterEqual(guarded['LCR_dB'], 25.0)
        self.assertGreaterEqual(guarded['WNG_dB'], -3.0)
        self.assertLessEqual(guarded['dt_target_loss_db'], 1.0)
        self.assertEqual(res['eligible_architectures'], ['headrest_dsb_cross_cancel_guarded'])
        self.assertEqual(res['selected_architecture'], 'headrest_dsb_cross_cancel_guarded')
        with self.assertRaises(ValueError):
            bus_clock_drift_metrics(True, 60.0, 4000.0, 5.0)

    def test_icc_comb_notches_howling_margin_and_teams_interval_verdict(self):
        res = self.results['E11-22']
        pa = res['path_A_8ms']
        pb = res['path_B_22ms']
        # Independent odd-multiple enumeration in [300, 3400] Hz.
        notches_a = [62.5 * (2 * m + 1) for m in range(100) if 300 <= 62.5 * (2 * m + 1) <= 3400]
        notches_b = [float(Fraction(250, 11) * (2 * m + 1))
                     for m in range(200) if Fraction(300) <= Fraction(250, 11) * (2 * m + 1) <= Fraction(3400)]
        self.assertEqual((pa['p1150_verdict'], pb['p1150_verdict']), ('pass', 'fail'))
        self.assertEqual((pa['notch_count'], pb['notch_count']), (25, 68))
        self.assertEqual(pa['notch_frequencies_hz'], notches_a)
        np.testing.assert_allclose(pb['notch_frequencies_hz'], notches_b, rtol=0, atol=1e-12)
        self.assertAlmostEqual(res['peak_to_notch_ripple_db'], 20 * math.log10(9), places=12)
        self.assertAlmostEqual(res['howling_max_gain_db'], 12.0, places=12)
        self.assertAlmostEqual(res['howling_max_gain_linear'], 10 ** 0.6, places=12)
        self.assertEqual(
            {k: v['verdict'] for k, v in res['howling_gain_candidates'].items()},
            {'12.0dB': 'pass', '15.0dB': 'fail', '18.0dB': 'fail'},
        )
        self.assertEqual(res['teams_certification_verdict'], 'undetermined')
        self.assertEqual(res['teams_certification']['per_metric_verdicts']['tclw_max_volume_db'], 'undetermined')
        for bad_delay in (3.0, 4.0):  # excess delay = dsp + 1.0 - 5.0 <= 0 ms
            with self.assertRaises(ValueError):
                icc_comb_and_howling_metrics(bad_delay)

    def test_computed_results_are_strict_json(self):
        json.loads(json.dumps(self.results,allow_nan=False))


if __name__ == '__main__':
    unittest.main()
