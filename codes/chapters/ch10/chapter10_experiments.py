"""E10-18..28: engineering calculations and checked PCM, not hardware benchmarks.

Run ``python -m codes.chapters.ch10.chapter10_experiments``. Every answer is computed
from the stated inputs. Published audio is checked read-only, never regenerated.
"""
from __future__ import annotations

# Allow the documented direct-file command as well as python -m.
if __name__ == "__main__" and not __package__:
    import sys as _chapter_entry_sys
    from pathlib import Path as _ChapterEntryPath
    _chapter_entry_sys.path.insert(0, str(_ChapterEntryPath(__file__).resolve().parents[3]))
import json
import math
import hashlib
import struct
import wave
from fractions import Fraction
from pathlib import Path
import numpy as np
from codes.chapters.ch10.core.engineering import resample_sro_to_reference
from codes.chapters.ch10.core.noise_suppression import power_spectral_subtraction
from codes.chapters.ch00.core.audio_samples import agc_blocks_case, prepare_exports
from codes.chapters.ch00.examples.generate_audio_samples import INPUTS
from codes.chapters.ch10.core.noise_mismatch import build_fixture, analyze_fixture, analyze_pcm, STEMS
from codes.chapters.ch10.examples.generate_noise_mismatch import check_assets, OUTPUT as NOISE_OUTPUT
from codes.chapters.ch10.core.channel_selection import teaching_channel_selection
from codes.chapters.ch10.core.pdm_and_eval import pdm_cic_teaching_case, p340_p808_teaching_case


def subtraction_anchor():
    clean, noise = 1., np.array([1., -1.])
    observation = clean + noise
    raw_power = observation**2 - np.mean(noise**2)
    # The first bin-frame is a known noise-only observation; score the last two.
    enhanced, estimate = power_spectral_subtraction(np.array([[1., *observation]]),
                                                   np.array([0]), floor_ratio=0.)
    return {'clean': clean, 'noise': noise.tolist(), 'observation': observation.tolist(),
            'mean_cross_term': float(np.mean(2 * clean * noise)),
            'noise_power': float(estimate[0]), 'unclipped_power': raw_power.tolist(),
            'unclipped_mean': float(raw_power.mean()),
            'clipped_power': (np.abs(enhanced[0, 1:])**2).tolist(),
            'clipped_mean': float(np.mean(np.abs(enhanced[0, 1:])**2))}


def alias_anchor():
    input_rate, output_rate, frequency, count = 8000, 4000, 3000, 17
    signal = np.cos(2 * np.pi * frequency * np.arange(count) / input_rate)
    output = resample_sro_to_reference(signal, (input_rate / output_rate - 1) * 1e6)
    alias = np.cos(2 * np.pi * 1000 * np.arange(len(output)) / output_rate)
    return {'input_rate_hz': input_rate, 'output_rate_hz': output_rate,
            'frequency_hz': frequency, 'input_samples': count,
            'output': output.tolist(), 'alias_frequency_hz': 1000,
            'max_alias_error': float(np.max(np.abs(output - alias))),
            'anti_alias_filter_run': False}


ROOT = Path(__file__).resolve().parents[3]


def _checked_agc_pcm(case, audio_root):
    """Verify source identities, this group's manifest and four real WAVs.

    Replay only the AGC group to check its bytes, then independently score the
    stored integer samples. Other groups are not replayed or scored here.
    """
    audio_root = Path(audio_root)
    if '..' in audio_root.parts:
        raise ValueError('AGC audio root must not contain lexical parent traversal')
    manifest_path = audio_root / 'ch00/audio/MANIFEST.json'
    paths = [manifest_path] + [audio_root / 'ch10/audio' / (s + '.wav')
                               for s in case['signals']]
    for path in paths:
        aliases = {'/tmp', '/var', '/etc'}
        unsafe_link = any(p.is_symlink() and not (str(p) in aliases
                          and p.resolve() == Path('/private') / p.name) for p in [path, *path.parents])
        if not path.is_file() or unsafe_link:
            raise ValueError('AGC requires regular manifest/WAV files without symlinks')
    def invalid_constant(value):
        raise ValueError('nonfinite JSON constant: ' + value)
    def unique_object(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError('duplicate audio manifest key')
            result[key] = value
        return result
    def finite_float(value):
        result = float(value)
        if not math.isfinite(result):
            raise ValueError('JSON number exceeds finite float64 range')
        return result
    def same_json(actual, expected):
        # Canonical JSON also distinguishes bool from numeric 0/1 and rejects
        # overflowed JSON numeric literals such as 1e999 (not parse_constant).
        return json.dumps(actual, sort_keys=True, allow_nan=False) == json.dumps(
            expected, sort_keys=True, allow_nan=False)
    manifest = json.loads(manifest_path.read_text(), parse_constant=invalid_constant,
                          parse_float=finite_float, object_pairs_hook=unique_object)
    if not isinstance(manifest, dict):
        raise ValueError('main audio manifest must be an object')
    sources = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in INPUTS}
    if not same_json(manifest.get('schema_version'), 2) or manifest.get('generator_inputs') != sources:
        raise ValueError('main audio source identities are stale')
    expected, groups = prepare_exports({'agc_blocks': case})
    if not isinstance(manifest.get('groups'), dict) or not same_json(
            manifest['groups'].get('agc_blocks'), groups['agc_blocks']):
        raise ValueError('AGC manifest parameters or measurements are stale')
    records = manifest.get('files')
    if not isinstance(records, list) or any(not isinstance(r, dict) for r in records):
        raise ValueError('invalid audio record list')
    selected = [r for r in records if r.get('group') == 'agc_blocks']
    if any(not isinstance(r.get('file'), str) for r in selected):
        raise ValueError('AGC file identities must be strings')
    if len(selected) != 4 or {r.get('file') for r in selected} != set(expected):
        raise ValueError('AGC requires exactly four distinct manifest records')
    integers, hashes = {}, {}
    for record in selected:
        name = record['file']
        blob, info = expected[name]
        digest = hashlib.sha256(blob).hexdigest()
        if not same_json(record, {'file': name, 'chapter': 'ch10', 'sha256': digest, **info}):
            raise ValueError('AGC WAV manifest record differs: ' + name)
        path = audio_root / 'ch10/audio' / name
        actual = path.read_bytes()
        if actual != blob or hashlib.sha256(actual).hexdigest() != digest:
            raise ValueError('AGC WAV bytes differ: ' + name)
        with wave.open(str(path), 'rb') as reader:
            if (reader.getframerate(), reader.getnchannels(), reader.getsampwidth(),
                    reader.getnframes(), reader.getcomptype()) != (16000, 1, 2, 32000, 'NONE'):
                raise ValueError('AGC PCM format differs')
            raw = reader.readframes(32000)
        integers[name[:-4]] = struct.unpack('<32000h', raw)
        hashes[name] = digest
    reference = integers['agc_blocks_input']
    scores = {}
    for stem, samples in integers.items():
        windows = {}
        for name, (start, stop) in case['parameters']['score_windows_samples_half_open'].items():
            denominator = sum(v * v for v in reference[start:stop])
            numerator = sum(v * v for v in samples[start:stop])
            if denominator <= 0:
                raise ValueError('AGC score reference has zero energy')
            windows[name] = {'sample_count': stop-start,
                             'input_integer_squared_sum': denominator,
                             'output_integer_squared_sum': numerator,
                             'input_rms': math.sqrt(denominator/(stop-start))/32768,
                             'output_rms': math.sqrt(numerator/(stop-start))/32768,
                             'rms_ratio': math.sqrt(numerator/denominator)}
        scores[stem] = {'peak': max(abs(v) for v in samples)/32768, 'windows': windows}
    return scores, {'checked_read_only': True, 'generator_inputs_checked': len(sources),
                    'manifest_sha256': hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                    'wav_sha256': hashes, 'scored_from_stored_integer_pcm': True}


def agc_anchor(audio_root=None):
    case = agc_blocks_case()
    scores, validation = _checked_agc_pcm(case, ROOT / 'codes/chapters' if audio_root is None else audio_root)
    return {key: case[key] for key in ('parameters', 'float_analysis', 'limits')} | {
        'pcm_analysis': scores, 'pcm_validation': validation,
        'common_export_gain': case['export_gain_override'],
        'files': [name + '.wav' for name in case['signals']],
        'availability': {name: {'first_block_end_s': record['blocks'][0]['available_time_s'],
                               'attack_alpha': record['attack_alpha'],
                               'release_alpha': record['release_alpha']}
                         for name, record in case['block_records'].items()}}


def noise_update_anchor():
    probability, noise, ungated = Fraction(0), Fraction(1), Fraction(1)
    rows = []
    for ratio in (9, 9, 1):
        indicator = int(ratio > 2)
        probability = probability/2 + Fraction(indicator, 2)
        retention = Fraction(4, 5) + probability/5
        noise = retention*noise + (1-retention)*9
        ungated = Fraction(4, 5)*ungated + Fraction(9, 5)
        rows.append({'ratio': ratio, 'indicator': indicator,
                     'smoothed_control_probability': float(probability),
                     'retention': float(retention), 'noise_power': float(noise),
                     'noise_power_fraction': str(noise), 'ungated_noise_power': float(ungated)})
    return {'observed_power': 9, 'initial_noise_power': 1, 'rows': rows,
            'local_statistics_given': True, 'complete_mcra_run': False}


def finite_rir_anchor():
    powers = [Fraction(1, 2**n) for n in range(4)]
    energy = [sum(powers[n:]) for n in range(4)]
    levels = [10*math.log10(float(e/energy[0])) for e in energy]
    return {'sample_rate_hz': 10, 'retained_squared_samples': [float(p) for p in powers],
            'retained_tail_energy_fractions': [str(e) for e in energy],
            'retained_decay_db': levels, 'infinite_tail_slope_db_per_s': -100*math.log10(2),
            'infinite_tail_t60_s': .6/math.log10(2),
            'finite_endpoint_slope_db_per_s': levels[-1]/.3,
            'finite_endpoint_extrapolated_t60_s': 1.8/math.log10(15),
            'endpoint_interval_s': [0, .3], 'rt20_measurement': False}


def blocking_anchor():
    tasks = {'A': {'period_ms': 10, 'service_ms': 2, 'deadline_ms': 10},
             'B': {'period_ms': 100, 'service_ms': 15}}
    start_b, release_a = -1, 0
    finish_b = start_b + tasks['B']['service_ms']
    start_a = max(release_a, finish_b)
    finish_a = start_a + tasks['A']['service_ms']
    response = finish_a - release_a
    utilization = sum(Fraction(job['service_ms'], job['period_ms']) for job in tasks.values())
    # No higher-priority tasks: an already running B may contribute up to its
    # full service time, followed by A's own service, as B's start approaches
    # A's release from below. This bound is not this trajectory's response.
    upper_bound = tasks['B']['service_ms'] + tasks['A']['service_ms']
    return {'tasks': tasks, 'utilization': float(utilization),
            'nonpreemptive': {'B_start_ms': start_b, 'B_finish_ms': finish_b,
                             'A_release_ms': release_a, 'A_start_ms': start_a,
                             'A_finish_ms': finish_a, 'response_ms': response,
                             'deadline_missed': response > tasks['A']['deadline_ms']},
            'response_upper_bound_ms': upper_bound,
            'preemptive_A_response_ms': tasks['A']['service_ms'],
            'scope': 'one worker, fixed priorities A above B, no overhead'}


def clock_identifiability_anchor():
    scenarios = [(Fraction(1, 10000), Fraction(0)),
                 (Fraction(20001, 100000000), Fraction(1, 10000))]
    rows = [{'device_error_ppm': float(e*1000000), 'timer_error_ppm': float(h*1000000),
             'observed_rate_ratio': float((1+e)/(1+h)),
             'observed_rate_ratio_fraction': str((1+e)/(1+h))} for e, h in scenarios]
    ratio = (1+scenarios[0][0])/(1+scenarios[0][1])
    return {'scope': 'one actual device count and external time; reference rate is nominal only',
            'scenarios': rows, 'timestamp_slope_fraction': str(1-1/ratio),
            'first_order_ppm': float((1-1/ratio)*1000000),
            'two_actual_devices_shared_timer_ratio': str((1+Fraction(1, 10000))/(1+Fraction(0))),
            'common_timer_scale_cancels_for_two_actual_devices': True}


def noise_mismatch_anchor(directory=None):
    """Check the published six-WAV experiment, then score stored PCM again.

    ``directory`` is an explicit ordinary asset directory for an independent
    temporary check. The chapter command defaults to the published directory.
    Missing/stale assets fail; this function never generates or repairs them.
    """
    directory = NOISE_OUTPUT if directory is None else Path(directory)
    manifest = check_assets(directory)
    fixture = build_fixture()
    buffers = {stem+'.wav': (directory/(stem+'.wav')).read_bytes() for stem in STEMS}
    return {'parameters': fixture['parameters'],
            'floating_point': analyze_fixture(fixture), 'pcm_analysis': analyze_pcm(buffers),
            'pcm_validation': {'checked_read_only': True, 'files': manifest['files'],
                               'source_sha256': manifest['source_sha256'],
                               'manifest_sha256': hashlib.sha256((directory/'MANIFEST.json').read_bytes()).hexdigest(),
                               'scored_from_stored_integer_pcm': True},
            'scope': 'fixed/polluted/known expected noise power; offline centered WOLA; no adaptive MCRA or speech benchmark'}


def channel_selection_anchor(directory=None):
    """E10-26 hand control plus strict read-only check of actual six PCM files.

    No source estimation or repair: missing, changed, or stale assets fail.
    Recompute PCM integer error from actual stored samples, separately from
    the given-statistics solution and unquantized component measurements.
    """
    from codes.chapters.ch10.examples.generate_channel_audio import check_assets as check_channel_assets, DEFAULT_OUTPUT
    from codes.chapters.ch10.core.channel_audio import FILE_NAMES, SCORING_INTERVAL
    directory = DEFAULT_OUTPUT if directory is None else Path(directory)
    manifest = check_channel_assets(directory)
    integers, hashes = {}, {}
    for role, filename in FILE_NAMES.items():
        path = directory/filename
        hashes[filename] = hashlib.sha256(path.read_bytes()).hexdigest()
        channels = 3 if role in ('healthy_array', 'faulty_array') else 1
        with wave.open(str(path), 'rb') as reader:
            if (reader.getframerate(), reader.getnchannels(), reader.getnframes(),
                    reader.getsampwidth(), reader.getcomptype()) != (16000, channels, 32000, 2, 'NONE'):
                raise ValueError('fixed channel selection PCM format required')
            raw = reader.readframes(32000)
        integers[role] = struct.unpack('<'+'h'*(32000*channels), raw)
    lo, hi = SCORING_INTERVAL
    reference = integers['reference'][lo:hi]
    reference_energy = sum(value*value for value in reference)
    if reference_energy <= 0:
        raise ValueError('channel selection PCM reference must have positive energy')
    denominator = (hi-lo)*32768**2
    scores = {}
    for role in ('healthy_output', 'stale_output', 'recomputed_output'):
        values = integers[role][lo:hi]
        error_energy = sum((value-truth)**2 for value, truth in zip(values, reference))
        power_energy = sum(value*value for value in values)
        scores[role] = {'error_squared_sum_pcm_integer': error_energy,
                       'power_squared_sum_pcm_integer': power_energy,
                       'mse': error_energy/denominator,
                       'nmse': error_energy/reference_energy}
    return {'mathematical_control': teaching_channel_selection(),
            'parameters': manifest['parameters'], 'analytic': manifest['analytic'],
            'pcm_analysis': {'score_interval_samples': [lo, hi], 'samples': hi-lo,
                             'reference_squared_sum_pcm_integer': reference_energy,
                             'mse_integer_denominator': denominator, 'scores': scores},
            'pcm_validation': {'checked_read_only': True, 'wav_sha256': hashes,
                               'manifest_sha256': hashlib.sha256((directory/'MANIFEST.json').read_bytes()).hexdigest(),
                               'source_sha256': manifest['source_sha256'],
                               'scored_from_stored_integer_pcm': True},
            'unquantized_measurements': {role: sample['float_measurements']
                                         for role, sample in manifest['samples'].items()},
            'scope': manifest['limits']}


def pdm_cic_anchor():
    return pdm_cic_teaching_case()


def p340_p808_anchor():
    return p340_p808_teaching_case()


def run_experiments():
    functions = [subtraction_anchor, alias_anchor, agc_anchor, noise_update_anchor,
                 finite_rir_anchor, blocking_anchor, clock_identifiability_anchor,
                 noise_mismatch_anchor, channel_selection_anchor, pdm_cic_anchor,
                 p340_p808_anchor]
    return {f'E10-{number:02d}': function() for number, function in enumerate(functions, 18)}


if __name__ == '__main__':
    print(json.dumps(run_experiments(), ensure_ascii=False, indent=2, allow_nan=False))
