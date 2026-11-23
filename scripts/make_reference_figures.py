"""Figure 69: known playback-gain order and actual PCM echo-tail activity.

Run after generating the six independent chapter 6 reference-location WAVs.
This figure uses frozen known paths; it runs no gain estimator or device AEC.
"""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.make_figures import plt, finalize_figure, C_MAIN, C_BLUE, C_RED, C_GREEN
from codes.chapters.ch00.core.audio_samples import read_pcm16
from codes.chapters.ch06.examples.generate_reference_audio import DEFAULT_OUTPUT, check_assets
import numpy as np


def measurements():
    manifest = check_assets(DEFAULT_OUTPUT)
    sources = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
               for name in (*manifest['source_sha256'], 'scripts/make_figures.py',
                            'scripts/make_reference_figures.py')}
    signals = {name: read_pcm16((DEFAULT_OUTPUT/row['file']).read_bytes())[1][0]
               for name, row in manifest['samples'].items()}
    frames = signals['echo'][32000:33600].reshape(10, 160)
    centered_rms = np.sqrt(np.mean((frames-frames.mean(axis=1, keepdims=True))**2, axis=1))
    report = {'figure': 69, 'exercise_id': 'E06-28', 'source_sha256': sources,
              'audio_manifest_sha256': hashlib.sha256((DEFAULT_OUTPUT/'MANIFEST.json').read_bytes()).hexdigest(),
              'wav_sha256': {n: r['sha256'] for n, r in manifest['files'].items()},
              'parameters': manifest['parameters'], 'samples': manifest['samples'],
              'pcm_tail_activity': manifest['pcm_tail_activity'],
              'plotted_tail_centered_rms': centered_rms.tolist(),
              'limits': manifest['limits'],
              'plot_scope': 'declared gain schedule, actual PCM excerpt/power and actual PCM tail activity; no gain fitting or training'}
    return manifest, report, signals


def main():
    manifest, report, signals = measurements()
    report_path = ROOT/'codes/chapters/ch06/reports/figure69_reference_timing.json'
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    fig, axes = plt.subplots(2, 2, figsize=(10.3, 7.8), constrained_layout=True)
    time = np.arange(14400, 19201)/16000
    axes[0, 0].step(time, np.where(time >= 1., 2., 1.), where='post', color=C_BLUE,
                    label='当前增益 g[n]')
    axes[0, 0].step(time, np.where(time >= 1.1, 2., 1.), where='post', color=C_RED,
                    linestyle='--', label='延迟项增益 g[n−1600]')
    axes[0, 0].axvspan(1., 1.1, color=C_RED, alpha=.10)
    axes[0, 0].set(title='(a) 已知播放增益：历史项晚 100 ms 改变',
                   xlabel='时间 / s', ylabel='无量纲增益', ylim=(.8, 2.5))
    axes[0, 0].legend(loc='upper left', fontsize=11)
    methods = ('early_residual', 'wrong_gain_residual', 'late_residual')
    labels = ('早参考固定路径', '预测后乘当前增益', '晚参考固定路径')
    colors = (C_BLUE, C_RED, C_GREEN)
    linestyles = ('-', '--', '-.')
    markers = ('o', 's', '^')
    hatches = ('', '///', 'xxx')
    count = 80
    for idx, (method, label, color, ls, marker) in enumerate(zip(methods, labels, colors, linestyles, markers)):
        axes[0, 1].plot(np.arange(count)/16000*1000, signals[method][16000:16000+count],
                        color=color, linestyle=ls, marker=marker, markevery=(idx*4, 12),
                        markersize=5.0, label=label, linewidth=1.6)
    axes[0, 1].set(title='(b) 增益切换后 5 ms：实际 PCM 残差',
                   xlabel='从第 16000 点起 / ms', ylabel='解码幅度', ylim=(-.12, .195))
    axes[0, 1].legend(loc='upper right', fontsize=9.6, framealpha=0.96)
    x = np.arange(2)
    for i, (method, label, color, hatch) in enumerate(zip(methods, labels, colors, hatches)):
        row = manifest['samples'][method]
        values = [row['pcm_integer_measurements'][w]['mean_square_all_channels'] for w in ('step', 'stable')]
        position = x+(i-1)*.24
        axes[1, 0].bar(position, np.asarray(values)*1000, width=.22, color=color,
                       edgecolor=C_MAIN, linewidth=0.85, hatch=hatch, alpha=.82, label=label)
        analytic = np.asarray([row['analytic'][w]['mean_square'] for w in ('step', 'stable')]) * 1000
        for px, ay in zip(position, analytic):
            axes[1, 0].plot([px - .135, px + .135], [ay, ay], color='white', linewidth=4.2,
                            solid_capstyle='butt', zorder=4)
            axes[1, 0].plot([px - .135, px + .135], [ay, ay], color=C_MAIN, linewidth=2.3,
                            solid_capstyle='butt', zorder=5)
            axes[1, 0].plot(px, ay, marker='D', markersize=4.8, markerfacecolor='white',
                            markeredgecolor=C_MAIN, markeredgewidth=1.4, zorder=6)
    axes[1, 0].set(title='(c) 实际 PCM 功率与解析基准', ylabel=r'均方幅度 / $10^{-3}$',
                   xticks=x, xticklabels=['切换窗：1600 点', '稳定窗：9600 点'], ylim=(-.25, 8.8))
    axes[1, 0].plot([], [], color=C_MAIN, linewidth=2.3, marker='D', markersize=4.8,
                    markerfacecolor='white', markeredgecolor=C_MAIN, markeredgewidth=1.4, label='解析均方')
    axes[1, 0].legend(loc='upper center', ncol=2, fontsize=9.4, framealpha=0.96, columnspacing=0.9)
    centers = np.arange(10)*.01+.005
    axes[1, 1].plot(centers*1000, report['plotted_tail_centered_rms'], 'o-', color=C_RED,
                    label='纯回声：帧去均值 RMS')
    axes[1, 1].plot(centers*1000, np.zeros(10), 's--', color=C_BLUE, label='当前晚参考：全零')
    axes[1, 1].axhline(.001, color='.4', linestyle=':', label='活动门限 0.001')
    axes[1, 1].set(title='(d) 播放已停止：十帧标签均 near_only',
                   xlabel='从第 32000 点起的帧中心 / ms', ylabel='去均值 RMS', ylim=(-.004, .066))
    axes[1, 1].legend(loc='lower left', bbox_to_anchor=(0, .09), fontsize=10)
    for ax in axes.flat:
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
    fig.suptitle('图69  参考处理的位置与历史：已知路径冻结的数学基准')
    finalize_figure(fig)
    metadata = {'SourceScript': 'scripts/make_reference_figures.py',
                'SourceScriptDigest': report['source_sha256']['scripts/make_reference_figures.py'],
                'GeneratorInputs': json.dumps(report['source_sha256'], sort_keys=True),
                'AudioManifestDigest': report['audio_manifest_sha256'],
                'NumericalReportDigest': hashlib.sha256(report_path.read_bytes()).hexdigest()}
    destination = ROOT/'figures/fig69_reference_timing.png'
    metadata['FigureFilename'] = destination.name
    fig.savefig(destination, dpi=150, bbox_inches='tight', facecolor='white', metadata=metadata)
    plt.close(fig)
    print(destination.relative_to(ROOT))


if __name__ == '__main__':
    main()
