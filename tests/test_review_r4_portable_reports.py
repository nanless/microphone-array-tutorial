"""Independent arithmetic negatives for the portable publication report layer."""
from copy import deepcopy
import json
import math
from pathlib import Path
import unittest
import cmath
import hashlib
import tempfile
import wave
from unittest.mock import patch

from scripts import portable_figure_report_checks as checks
from scripts import run_portable_ci as ci

ROOT = Path(__file__).resolve().parents[1]


def report(number):
    from scripts.figure_report_contracts import REPORT_BY_NUMBER
    return json.loads((ROOT/REPORT_BY_NUMBER[number].report_relative).read_text())


class IndependentReportArithmeticTests(unittest.TestCase):
    def test_math_all_anchors_and_each_anchor_negative(self):
        for number in range(84,88):
            original=report(number)
            checks._math_report(number,original)
            for key,value in original['anchors'].items():
                altered=deepcopy(original)
                if isinstance(value,list):altered['anchors'][key][0]+=0.01
                else:altered['anchors'][key]+=1 if type(value)is int else .01
                with self.subTest(figure=number,anchor=key), self.assertRaises(ValueError):
                    checks._math_report(number,altered)

    def test_stored_trial_statistics_not_only_source_hash(self):
        original=report(13);checks._report13(original)
        for key in ('rate','low','high','median','q1','q3','successes','failures','exceptions'):
            altered=deepcopy(original);altered['statistics'][key][0][0]+=.01 if type(altered['statistics'][key][0][0])is float else 1
            with self.subTest(key=key),self.assertRaises(ValueError):checks._report13(altered)
        altered=deepcopy(original);altered['statistics']['trial_correct'][0][0][0]=1
        with self.assertRaises(ValueError):checks._report13(altered)

    def test_tracking_event_denominators(self):
        original=report(21);checks._report21(original)
        for field in ('raw_observed_rmse_deg','pf_observed_rmse_deg','pf_missing_rmse_deg','minimum_ess'):
            altered=deepcopy(original);altered['scores'][field]+=.1
            with self.subTest(field=field),self.assertRaises(ValueError):checks._report21(altered)
        altered=deepcopy(original);altered['frames']['resampled'][3]=not altered['frames']['resampled'][3]
        with self.assertRaises(ValueError):checks._report21(altered)
        altered=deepcopy(original)
        altered['phd_illustration']['intensity_targets_per_degree']=[v/2 for v in altered['phd_illustration']['intensity_targets_per_degree']]
        altered['phd_illustration']['integral_targets']/=2
        with self.assertRaises(ValueError):checks._report21(altered)

    def test_alignment_energy_oracle_and_wrong_ratio(self):
        original=report(20)
        for row in original['results'].values():
            row.update(quiet_energy=1.,total_energy=10.,quiet_energy_ratio_db=-10.,
                active_scale_aligned_nmse_db=10*math.log10(.75),active_alignment={
                'reference_energy':4.,'estimate_energy':9.,'reference_conjugate_dot_estimate':[3.,0.]})
        checks._report20(original)
        for field in ('quiet_energy_ratio_db','active_scale_aligned_nmse_db'):
            altered=deepcopy(original);altered['results']['output'][field]+=.01
            with self.subTest(field=field),self.assertRaises(ValueError):checks._report20(altered)
        altered=deepcopy(original);del altered['results']['input']['active_alignment']
        with self.assertRaises(ValueError):checks._report20(altered)

    def test_wearable_and_tse_independent_negatives(self):
        cases={76:('anchor_points','quadrature_90deg','intensity_norm'),
            80:('present_waveform_nmse_db',12),81:('radial_gain_db','bubble_d10cm',2),
            82:('der_cases','over_100','der'),83:('consonant_loss_percent',0)}
        for number,path in cases.items():
            original=report(number);function=getattr(checks,'_report'+str(number));function(original)
            altered=deepcopy(original);target=altered
            for key in path[:-1]:target=target[key]
            target[path[-1]]+=.01
            with self.subTest(figure=number),self.assertRaises(ValueError):function(altered)

    def test_online_arrival_uses_integer_equality_boundary(self):
        # A video frame arriving at20ms is unavailable at8ms, available24ms.
        indices,valid,arrival=checks._online_map(0,3)
        self.assertEqual(indices,[-1,0,0]);self.assertEqual(valid,[False,True,True])
        self.assertEqual(arrival,[0.,.02,.02])
        changed=report(83);changed['sync_video_indices'][0]=0
        with self.assertRaises(ValueError):checks._report83(changed)

    def test_tolerance_does_not_accept_bool_or_material_error(self):
        checks._close(1.+1e-12,1.)
        for value in (True,1,1.0001,float('nan'),float('inf')):
            with self.subTest(value=value),self.assertRaises(ValueError):checks._close(value,1.)


def numeric_fixture(number):
    """Construct numeric fixtures without claiming the old assets were regenerated.

    Current source hashes are supplied only to reach arithmetic in existing
    file oracles. Production identity/PNG guards are tested separately.
    """
    from scripts.figure_report_contracts import REPORT_BY_NUMBER
    data=report(number);contract=REPORT_BY_NUMBER[number]
    digest=hashlib.sha256((ROOT/contract.producer_relative).read_bytes()).hexdigest()
    for key in ('script_sha256','source_script_sha256','generator_sha256'):
        if key in data:data[key]=digest
    if isinstance(data.get('source_sha256'),dict):
        data['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in data['source_sha256']}
    elif 'source_sha256'in data:data['source_sha256']=digest
    if number in (65,67):
        data['audio_manifest_sha256']=hashlib.sha256((checks.quality.DISTRIBUTED_AUDIO_ROOT/'MANIFEST.json').read_bytes()).hexdigest()
    if number in (72,73):
        folder=checks.quality.TRACKING_AUDIO_ROOT if number==72 else checks.quality.CHANNEL_AUDIO_ROOT
        data['input_sha256']={name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in data['input_sha256']}
    if number==20:
        for row in data['results'].values():
            row.update(quiet_energy=1.,total_energy=10.,quiet_energy_ratio_db=-10.,
                active_scale_aligned_nmse_db=10*math.log10(.75),active_alignment={
                    'reference_energy':4.,'estimate_energy':9.,'reference_conjugate_dot_estimate':[3.,0.]})
    if number==77:
        case=data['circular_4khz_72dir']
        case['parameters']={'frequency_hz':4000.,'radius_m':.025,'sound_speed_m_s':343.,
            'direction_count':72,'regularization':.05,'iterations':25}
        # Fixed author-producer values independently reviewed against the4x4
        # normal equations. No function under test generates these weights.
        case['ls_weights']={'real':[-.6235175500351416,.18650516347260931,-.6235175500351418,-.049491063175824156],
            'imag':[.054637487371200825,-.27438888749044354,.054637487371200784,.08902518974400486]}
        case['magls_weights']={'real':[-.3172307462626342,.36622356750478857,-.31723074626263403,-.2519142324987458],
            'imag':[-.21698539533017808,-.9295920010638994,-.21698539533017794,-.15585879378340475]}
        kr=2*math.pi*4000*.025/343
        case['target_magnitudes']=[1+.55*math.sin(i*math.pi/36) for i in range(72)]
        for prefix in ('ls','magls'):
            weights=[complex(a,b) for a,b in zip(case[prefix+'_weights']['real'],case[prefix+'_weights']['imag'])]
            case[prefix+'_magnitudes']=[abs(sum(w.conjugate()*cmath.exp(1j*kr*math.cos(m*math.pi/2-i*math.pi/36))
                for m,w in enumerate(weights))) for i in range(72)]
    return data


class AllReportNumericControls(unittest.TestCase):
    NEGATIVES={13:('statistics','high',0,0),20:('results','output','active_scale_aligned_nmse_db'),
        21:('scores','pf_missing_rmse_deg'),45:('wola_outputs','fft8_keep4',5),
        56:('same_scalar_state','posterior_variance',0),57:('finite_tail','finite_db',0),
        59:('pcm_nmse','dual','fir9'),60:('analytic_mse','gls'),61:('reflection_squared_gain','a',0),
        62:('results','two_cell','cases','coherent','gs','history',1,0),
        63:('results','rank_one_clean_sc','relative_peak_stop_factor'),64:('results','intercept_prediction'),
        65:('results','rows','central_white','float','noise_power'),
        66:('results','trajectories','round_robin','history',0,'solutions',0,'condition_number'),
        67:('results','clock','corrected_peer_reconstruction_nmse'),
        68:('samples','flip','pcm_measurements','phasor_real_imag',0,0),69:('plotted_tail_centered_rms',0),
        70:('power_decay_control','relative_power',1),71:('overlap_fit','decoded_pcm_overlap_application_gain',0),
        72:('lifecycle','rows',99,'age_ticks'),
        73:('samples','healthy_array','float_measurements','mean_square_per_channel',0),
        74:('numerical_case','candidates',0,'DI_dB'),75:('results','plotted','scan_squared',0,0),
        76:('anchor_points','quadrature_90deg','intensity_norm'),77:('circular_4khz_72dir','objective_history',12),
        78:('analytic_reference_mse_ch0','wearable_fused_output'),79:('pcm_tone_magnitudes_clean_ref',0),
        80:('present_waveform_nmse_db',12),81:('radial_gain_db','bubble_d10cm',2),
        82:('der_cases','over_100','der'),83:('consonant_loss_percent',0),
        84:('anchors','parseval_energy'),85:('anchors','tikhonov_eps0p1_error_norm'),
        86:('anchors','whitening_sinr_gain_db'),87:('anchors','dft_dc_pseudo_variance')}

    def test_every_report_has_a_valid_numeric_fixture_and_independent_negative(self):
        from scripts.figure_report_contracts import REPORT_BY_NUMBER
        self.assertEqual(set(self.NEGATIVES),set(REPORT_BY_NUMBER))
        original_read=checks._read
        def fixture_read(path):
            value=original_read(path)
            if path.parent.name=='channel_audio' and path.name=='MANIFEST.json':
                value=deepcopy(value)
                value['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in value['source_sha256']}
            return value
        with tempfile.TemporaryDirectory() as directory,patch.object(checks,'_read',side_effect=fixture_read):
            for number,pointer in self.NEGATIVES.items():
                with self.subTest(figure=number):
                    data=numeric_fixture(number);path=Path(directory).resolve()/Path(REPORT_BY_NUMBER[number].report_relative).name
                    path.write_text(json.dumps(data))
                    checks.validate_numeric_report(number,data,path=path)
                    changed=deepcopy(data);item=changed
                    for key in pointer[:-1]:item=item[key]
                    item[pointer[-1]]+=1 if type(item[pointer[-1]])is int else .01
                    path.write_text(json.dumps(changed))
                    with self.assertRaises(ValueError):checks.validate_numeric_report(number,changed,path=path)

    def test_consistent_report_png_hash_cannot_hide_bad_arithmetic(self):
        from scripts import figure_report_contracts as registry
        from PIL import Image,PngImagePlugin
        entry=registry.REPORT_BY_NUMBER[83]
        with tempfile.TemporaryDirectory() as directory,patch.object(registry,'FIGURE_REPORTS',(entry,)):
            root=Path(directory).resolve();path=root/entry.report_relative;path.parent.mkdir(parents=True)
            data=report(83);data['consonant_loss_percent'][0]=0.
            path.write_text(json.dumps(data));figures=root/'figures';figures.mkdir()
            metadata=PngImagePlugin.PngInfo();metadata.add_text('NumericalReportDigest',hashlib.sha256(path.read_bytes()).hexdigest())
            Image.new('RGB',(4,4),'white').save(figures/entry.png_filename,pnginfo=metadata)
            self.assertEqual(registry.report_bindings_issues(root,figures),[])
            with self.assertRaises(ValueError):checks.validate_numeric_report(83,data,path=path)

    def test_binaural_measurement_schema_and_actual_four_wav_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory).resolve();data={'measurements':{}}
            for role in ('input','bmvdr','blcmv','null'):
                name='binaural_'+role+'.wav';path=root/name
                with wave.open(str(path),'wb') as stream:
                    stream.setparams((2,2,16000,32000,'NONE','not compressed'));stream.writeframes(b'\0'*128000)
                data['measurements'][role]={'file':name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
            errors=[];ci.check_binaural_measurement_bindings(root,data,errors);self.assertEqual(errors,[])
            (root/'binaural_null.wav').write_bytes(b'bad')
            errors=[];ci.check_binaural_measurement_bindings(root,data,errors);self.assertTrue(errors)
            data['measurements']['input']['file']='../escape.wav'
            errors=[];ci.check_binaural_measurement_bindings(root,data,errors);self.assertTrue(any('filename' in e for e in errors))

    def test_schema_rejects_unknown_numeric_field_and_boolean_substitution(self):
        original=numeric_fixture(66)
        changed=deepcopy(original);changed['results']['model']['a'][0]=True
        with self.assertRaisesRegex(ValueError,'numeric schema differs'):
            checks.validate_numeric_report(66,changed)

        changed=deepcopy(original);changed['results']['unchecked_future_curve']=[1.,2.]
        with self.assertRaisesRegex(ValueError,'numeric schema differs'):
            checks.validate_numeric_report(66,changed)
        changed=deepcopy(original);del changed['results']['model']['white']
        with self.assertRaisesRegex(ValueError,'numeric schema differs'):
            checks.validate_numeric_report(66,changed)

    def test_cross_platform_checkout_keeps_source_identity_bytes(self):
        source=(ROOT/'.github/workflows/ci.yml').read_text()
        self.assertIn('os: [ubuntu-latest, macos-15, windows-latest]',source)
        self.assertLess(source.index('git config --global core.autocrlf false'),source.index('actions/checkout@'))

    def test_source_roles_are_fixed_outside_the_report(self):
        from scripts.figure_report_contracts import REPORT_BY_NUMBER
        original=numeric_fixture(83)
        original['audio_manifest_sha256']=hashlib.sha256((checks.quality.TSE_AUDIO_ROOT/'MANIFEST.json').read_bytes()).hexdigest()
        checks._identity(REPORT_BY_NUMBER[83],original)
        for add in (True,False):
            changed=deepcopy(original)
            if add:changed['source_sha256']['unexpected/source.py']='0'*64
            else:changed['source_sha256'].pop('codes/chapters/ch15/core/av_alignment.py')
            with self.subTest(add=add),self.assertRaises(ValueError):checks._identity(REPORT_BY_NUMBER[83],changed)
        changed=numeric_fixture(45);changed['source_script']='scripts/unrelated.py'
        with self.assertRaises(ValueError):checks._identity(REPORT_BY_NUMBER[45],changed)
        changed=deepcopy(original);changed['exercise_ids']=['E99-99']
        with self.assertRaises(ValueError):checks._identity(REPORT_BY_NUMBER[83],changed)

    def test_strict_json_rejects_duplicates_and_nonfinite_numbers(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory).resolve()/'report.json'
            for text in ('{"figure":83,"figure":84}','{"value":NaN}','{"value":1e309}'):
                path.write_text(text)
                with self.subTest(text=text),self.assertRaises(ValueError):checks._read(path)
            path.write_text('{"label":"NaN / 1e309","value":1.0}')
            self.assertEqual(checks._read(path),{'label':'NaN / 1e309','value':1.0})

    def test_new_report_and_html_fixture_modules_are_actually_selected(self):
        flags=dict(has_downloads=False,is_exact_host=False,is_arm64_darwin=False)
        for stem in ('test_review_r4_portable_reports','test_review_r4_publication_en'):
            self.assertEqual(ci.test_selection_reasons(stem,(ROOT/'tests'/(stem+'.py')).read_text(),**flags),[])


if __name__=='__main__':unittest.main()
