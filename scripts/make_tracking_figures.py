"""Generate figure77 from independently re-read existing Chapter09 PCM.

Only this new figure/report are written. No WAV or historical report is
regenerated. Paths are preflighted before numerical work; no race-proof or
crash-persistence guarantee is implied.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from codes.chapters.ch00.io_contracts import validate_report_destination, write_json_report
from codes.chapters.ch09.examples.tracking_lifecycle_demo import run_demo, SOURCE_PATHS as DEMO_SOURCES
from codes.chapters.ch09.examples.chapter09_tracking_audio import OUTPUT

FIGURE = ROOT/'figures/fig77_tracking_lifecycle.png'
REPORT = ROOT/'codes/chapters/ch09/reports/figure77_tracking_lifecycle.json'
SOURCE_PATHS = ('scripts/make_tracking_figures.py', *DEMO_SOURCES)


def _destination(path, official, suffix):
    path = validate_report_destination(path)
    if path.suffix != suffix or (path.resolve().is_relative_to(ROOT) and path.resolve() != official):
        raise ValueError('only the current figure/report or a normal external path is allowed')
    return path


def generate_figure(figure_path=FIGURE, report_path=REPORT, *, audio_directory=OUTPUT):
    """Read actual PCM first, then draw event and time-contract diagnostics."""
    figure_path = _destination(figure_path, FIGURE, '.png')
    report_path = _destination(report_path, REPORT, '.json')
    data = run_demo(audio_directory)
    frames, lifecycle = data['analysis_frames'], data['lifecycle']
    denominator = lifecycle['parameters']['tick_denominator_hz']
    times = np.array(frames['state_time_s'])
    publication = np.array(frames['available_time_s'])
    valid = np.array(frames['observation_valid'])
    # Ages come from exact support ticks, not subtracting rounded seconds.
    last_ticks, age_state, age_available = [], [], []
    last_tick = None
    for start, is_valid in zip(frames['start_sample'], valid):
        state_tick, available_tick = 2*start+511, 2*(start+512)
        if is_valid:
            last_tick = state_tick
        last_ticks.append(last_tick)
        age_state.append(None if last_tick is None else (state_tick-last_tick)/denominator)
        age_available.append(None if last_tick is None else (available_tick-last_tick)/denominator)
    phase_codes = [0 if row['phase']=='absent' else
                   (2*row['track_id']-1 if row['phase']=='tentative' else 2*row['track_id'])
                   for row in lifecycle['rows']]
    report = {'schema_version': 1, 'exercise_id': 'E09-25',
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                        'matplotlib': matplotlib.__version__},
        'input_sha256': data['input_sha256'], 'analysis_scores': data['analysis_scores'],
        'parameters': lifecycle['parameters'], 'lifecycle': lifecycle,
        'state_age_comparison': data['state_age_comparison'],
        'plot_data': {'state_time_s': times.tolist(), 'available_time_s': publication.tolist(),
            'observation_angle_deg': frames['observation_angle_deg'],
            'filtered_angle_deg': frames['filtered_angle_deg'],
            'valid': valid.tolist(), 'last_valid_measurement_ticks': last_ticks,
            'state_age_s': age_state, 'available_age_s': age_available, 'lifecycle_phase_code': phase_codes},
        'limits': data['limits']}
    json.dumps(report, allow_nan=False)
    plt.rcParams.update({'font.family': 'sans-serif',
        'font.sans-serif': ['PingFang SC', 'Hiragino Sans GB', 'Noto Sans CJK SC',
                           'Microsoft YaHei', 'Arial Unicode MS', 'DejaVu Sans'],
        'axes.unicode_minus': False, 'font.size': 11, 'axes.titlesize': 12,
        'axes.labelsize': 11, 'legend.fontsize': 11, 'xtick.labelsize': 10, 'ytick.labelsize': 10,
        'text.color': '#1E242B', 'axes.labelcolor': '#1E242B', 'axes.edgecolor': '#4A5058',
        'axes.linewidth': 0.9, 'axes.spines.top': False, 'axes.spines.right': False,
        'grid.color': '#DFDCD4', 'grid.linestyle': '--', 'grid.linewidth': 0.65,
        'legend.facecolor': '#FAF9F6', 'legend.edgecolor': '#DDD8CE', 'legend.framealpha': 0.92})
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 9.8), gridspec_kw={'height_ratios':[1,1,1.1]})
    axes[0].plot(times[valid], np.array(frames['observation_angle_deg'],float)[valid], '.',
                 color='#1F5A92', markersize=4.5, label='实际 PCM → GCC 方位观测')
    axes[0].plot(times, frames['filtered_angle_deg'], '--', color='#B84A2E',
                 label='原连续 KF 的独立诊断（不因 ID 重置）')
    axes[0].axvspan(times[~valid][0], times[~valid][-1], color='#C9C2B6', alpha=.28,
                    label='24 个 RMS 门控缺测帧')
    axes[0].set(xlabel='状态所属时刻 / 接收窗中心 (s)',ylabel='方位角 (°)',
                title='(a) 同一实际 PCM 的方位观测与连续卡尔曼状态轨迹')
    axes[0].legend(loc='upper left')
    axes[1].step(publication, phase_codes, where='post', color='#1F5A92', linewidth=2)
    axes[1].set(yticks=range(5),yticklabels=['无轨迹','候选 ID 1','确认 ID 1','候选 ID 2','确认 ID 2'],
                ylim=(-.4,4.6),xlabel='事件发布时刻 / 接收块到齐 (s)',
                title='(b) 连续 3 帧确认；超时退役后分配新 ID 并保持因果发布')
    first_confirm=lifecycle['events'][1]
    axes[1].annotate('第 3 帧到齐才确认：0.052 s',
        xy=(first_confirm['publication_time_s'],2),xytext=(.18,3.7),
        arrowprops={'arrowstyle':'->','color':'#1F5A92'}, fontsize=11)
    axes[1].annotate('重新确认 ID 2：1.112 s',xy=(1.112,4),xytext=(1.3,2.75),
        arrowprops={'arrowstyle':'->','color':'#1F5A92'}, fontsize=11)
    axes[2].plot(publication, age_available, color='#B84A2E',label='到齐时刻 − 最后有效观测中心')
    axes[2].plot(publication, age_state, '--',color='#1F5A92',label='状态中心 − 最后有效观测中心（对照）')
    axes[2].axhline(.20,linestyle=':',color='black',label='指定上限 0.20 s；等号保留')
    for event,color,y,text in [(lifecycle['events'][2],'#B84A2E',.20603125,'到齐龄期：1.032 s 过期'),
            (data['state_age_comparison']['events'][2],'#1F5A92',.21,'状态龄期对照：1.052 s 才过期')]:
        axes[2].plot(event['publication_time_s'],y,'o',color=color)
        offset=(.815,.268) if color=='#B84A2E' else (.915,.042)
        axes[2].annotate(text,xy=(event['publication_time_s'],y),xytext=offset,
                         bbox=dict(fc='white',ec='none',alpha=0.90,pad=1.2),
                         arrowprops={'arrowstyle':'->','color':color},fontsize=11,color=color)
    axes[2].set(xlim=(.80,1.14),ylim=(-.006,.29),xlabel='本帧接收块到齐时刻 (s)',ylabel='观测龄期 (s)',
                title='(c) 两种龄期的差等于固定窗可用延迟 16.03125 ms')
    handles, labels = axes[2].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5,.035), fontsize=11)
    for ax in axes:
        ax.grid(linestyle=':',alpha=.45)
    fig.suptitle('图77  PCM 观测 → 轨迹确认、静默保持与超时退役',fontsize=14)
    fig.tight_layout(rect=(0,.15,1,.965),h_pad=2.5)
    fig.text(.5,.008,'注：本例按固定确认与过期门限演示单轨迹槽位的建立、缺测保持与超时退役过程。',
             ha='center',fontsize=11)
    figure_path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(figure_path,dpi=150,bbox_inches='tight',facecolor='white',metadata={
        'SourceScript':'scripts/make_tracking_figures.py',
        'SourceScriptDigest':report['script_sha256'],
        'AudioManifestDigest':data['input_sha256']['MANIFEST.json'],
        'SourceClosureDigest':hashlib.sha256(json.dumps(report['source_sha256'],sort_keys=True).encode()).hexdigest()})
    plt.close(fig)
    write_json_report(report_path,report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure',type=Path,default=FIGURE)
    parser.add_argument('--report',type=Path,default=REPORT)
    parser.add_argument('--audio-directory',type=Path,default=OUTPUT)
    args=parser.parse_args()
    report=generate_figure(args.figure,args.report,audio_directory=args.audio_directory)
    print(json.dumps({'figure':str(args.figure),'report':str(args.report),'events':report['lifecycle']['events']},
                     ensure_ascii=False,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
