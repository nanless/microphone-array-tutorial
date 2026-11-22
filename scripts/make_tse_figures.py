"""Generate Chapter 15 figures (fig80..fig83) and numerical JSON reports.

Run from the repository root after generating tse_audio:
    .venv/bin/python scripts/make_tse_figures.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from codes.chapters.ch00.core.audio_samples import read_pcm16
from codes.chapters.ch00.io_contracts import (
    validate_parent_chain,
    validate_report_destination,
    write_json_report,
)
from codes.chapters.ch15.core.av_alignment import (
    av_confidence_gated_mask,
    map_video_to_audio_frames,
)
from codes.chapters.ch15.core.conditional_tse import (
    film_conditioning,
    tse_dual_state_loss,
    zero_overlap_bypass_decision,
)
from codes.chapters.ch15.core.diarization_tsvad import (
    compute_der_components,
    gss_soft_gated_mvdr,
    single_label_overlap_der_bound,
    tsvad_spatial_spectral_posterior,
)
from codes.chapters.ch15.core.distance_bubble import (
    broadside_distance_bubble_weights,
    fresnel_distance,
    spherical_ula_steering,
)
from codes.chapters.ch15.examples.generate_tse_audio import DEFAULT_OUTPUT, check_assets
from scripts.make_figures import (
    C_BLUE,
    C_GREEN,
    C_MAIN,
    C_ORANGE,
    C_PURPLE,
    C_RED,
    finalize_figure,
    plt,
)

REPORT_DIR = ROOT / "codes/chapters/ch15/reports"
FIG_DIR = ROOT / "figures"

SCRIPT_INPUT_PATHS = (
    "codes/chapters/ch15/core/conditional_tse.py",
    "codes/chapters/ch15/core/distance_bubble.py",
    "codes/chapters/ch15/core/diarization_tsvad.py",
    "codes/chapters/ch15/core/av_alignment.py",
    "codes/chapters/ch15/core/tse_audio.py",
    "codes/chapters/ch15/examples/generate_tse_audio.py",
    "codes/chapters/ch02/core/conventions.py",
    "codes/chapters/ch00/core/audio_samples.py",
    "codes/chapters/ch00/io_contracts.py",
    "scripts/make_figures.py",
    "scripts/make_tse_figures.py",
)


def _sha_file(rel_path: str) -> str:
    source = validate_parent_chain(ROOT / rel_path)
    if not source.is_file():
        raise ValueError(f"missing ordinary figure source: {rel_path}")
    return hashlib.sha256(source.read_bytes()).hexdigest()


def _save_report_and_png(fig, fig_name: str, report_name: str, report_payload: dict, manifest_sha: str) -> None:
    report_path = validate_report_destination(
        REPORT_DIR / report_name, forbidden_roots=(DEFAULT_OUTPUT,)
    )
    destination = validate_report_destination(
        FIG_DIR / fig_name, forbidden_roots=(DEFAULT_OUTPUT,)
    )
    write_json_report(report_path, report_payload, forbidden_roots=(DEFAULT_OUTPUT,))
    report_sha = hashlib.sha256(report_path.read_bytes()).hexdigest()

    if fig.get_constrained_layout():
        fig.get_layout_engine().set(w_pad=0.06, h_pad=0.09, wspace=0.04, hspace=0.06)
    finalize_figure(fig)
    metadata = {
        "SourceScript": "scripts/make_tse_figures.py",
        "SourceScriptDigest": report_payload["source_sha256"]["scripts/make_tse_figures.py"],
        "GeneratorInputs": json.dumps(report_payload["source_sha256"], sort_keys=True),
        "AudioManifestDigest": manifest_sha,
        "NumericalReportDigest": report_sha,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    validate_report_destination(destination, forbidden_roots=(DEFAULT_OUTPUT,))
    fig.savefig(destination, dpi=150, bbox_inches="tight", facecolor="white", metadata=metadata)
    plt.close(fig)
    print(destination.relative_to(ROOT))


def _build_fig85(sources: dict[str, str], manifest_sha: str) -> None:
    h_feat = np.array([[1.0, -0.5, 2.0], [0.5, 1.5, -1.0]], dtype=np.float64)
    z_mod = film_conditioning(h_feat, [2.0, 0.5], [-0.5, 0.25])

    v_int = np.array([1.0, -1.0, 1.0, -1.0], dtype=np.float64)
    s_zero = np.zeros(4, dtype=np.float64)
    leak_db_grid = np.linspace(-45.0, 0.0, 91, dtype=np.float64)
    absent_curves = {}
    for tau in (1e-4, 1e-3, 1e-2):
        losses = []
        for l_db in leak_db_grid:
            amp = float(10.0 ** (l_db / 20.0))
            res = tse_dual_state_loss(amp * v_int, s_zero, v_int, tau_floor=tau)
            losses.append(res["loss_db"])
        absent_curves[str(tau)] = losses

    s_ref = np.array([1.0, 1.0, -1.0, -1.0], dtype=np.float64)
    x_mix = s_ref + v_int
    s_base = s_ref + 0.1 * v_int
    alphas = np.linspace(0.2, 2.2, 41, dtype=np.float64)
    neg_sisdr_vals = [
        tse_dual_state_loss(float(a) * s_base, s_ref, x_mix)["loss_db"] for a in alphas
    ]
    waveform_err_db = [
        float(10.0 * np.log10(np.sum((a * s_base - s_ref) ** 2) / np.sum(s_ref ** 2)))
        for a in alphas
    ]

    y_pri = np.array([1.0, 0.0, -1.0, 0.0], dtype=np.float64)
    bypass_single = zero_overlap_bypass_decision(
        y_pri, np.array([0.0, 0.1, 0.0, -0.1]), threshold_db=-18.0
    )
    keep_overlap = zero_overlap_bypass_decision(
        y_pri, np.array([0.0, 0.5, 0.0, -0.5]), threshold_db=-18.0
    )

    report = {
        "figure": 80,
        "exercise_ids": ["E15-01", "E15-02"],
        "source_sha256": sources,
        "audio_manifest_sha256": manifest_sha,
        "film_input": h_feat.tolist(),
        "film_output": z_mod.tolist(),
        "leak_ratio_db_grid": leak_db_grid.tolist(),
        "absent_loss_curves_db": absent_curves,
        "scale_alphas": alphas.tolist(),
        "present_neg_sisdr_db": neg_sisdr_vals,
        "present_waveform_nmse_db": waveform_err_db,
        "bypass_examples": {
            "single_fricative_ratio_db": bypass_single["energy_ratio_db"],
            "single_fricative_bypass": bypass_single["bypass_separation"],
            "dual_overlap_ratio_db": keep_overlap["energy_ratio_db"],
            "dual_overlap_bypass": keep_overlap["bypass_separation"],
        },
        "plot_scope": "deterministic FiLM modulation, dual-state TSE loss curves, scale-invariance check, and single-speaker fricative bypass boundary",
    }

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.8), constrained_layout=True)
    frames = np.arange(3)
    axes[0, 0].plot(frames, h_feat[0], "o--", color=C_BLUE, markersize=6.5, linewidth=1.7, label="原始通道 1")
    axes[0, 0].plot(frames, z_mod[0], "s-", color=C_RED, markersize=6.5, linewidth=1.9, label="调制后通道 1")
    axes[0, 0].plot(frames, h_feat[1], "^:", color=C_GREEN, markersize=7.0, linewidth=1.7, label="原始通道 2")
    axes[0, 0].plot(frames, z_mod[1], "d-.", color=C_PURPLE, markersize=5.8, markerfacecolor="white", markeredgewidth=1.6, linewidth=1.8, label="调制后通道 2")
    axes[0, 0].set(
        title="(a) FiLM 通道仿射调制前后特征",
        xlabel="特征帧索引 t",
        ylabel="特征激活值（无量纲）",
        xticks=frames,
        xlim=(-0.22, 2.22),
        ylim=(-2.2, 4.9),
    )
    axes[0, 0].legend(ncol=2, loc="upper left", framealpha=0.96)

    for idx_t, (tau_str, label, color, style, mk) in enumerate((
        ("0.0001", r"地板系数 $\tau_{\mathrm{floor}}=10^{-4}$", C_RED, "-", "o"),
        ("0.001", r"地板系数 $\tau_{\mathrm{floor}}=10^{-3}$", C_BLUE, "--", "s"),
        ("0.01", r"地板系数 $\tau_{\mathrm{floor}}=10^{-2}$", C_GREEN, "-.", "^"),
    )):
        axes[0, 1].plot(
            leak_db_grid, absent_curves[tau_str], style, color=color,
            marker=mk, markevery=(idx_t * 5, 16), markersize=4.5, linewidth=1.8, label=label,
        )
    axes[0, 1].axhline(0.0, color="black", linestyle=":", linewidth=1.0)
    axes[0, 1].set(
        title="(b) 目标沉默段双态能量地板惩罚",
        xlabel="残余泄漏相对混合能量比 / dB",
        ylabel="沉默态损失 / dB",
        ylim=(-2.0, 44.0),
    )
    axes[0, 1].legend(loc="upper left", framealpha=0.96)

    axes[1, 0].plot(alphas, neg_sisdr_vals, "-", color=C_BLUE, marker="o", markevery=8, markersize=4.5, linewidth=2.0, label="尺度不变 -SI-SDR")
    axes[1, 0].plot(alphas, waveform_err_db, "--", color=C_ORANGE, marker="s", markevery=(4, 8), markersize=4.5, linewidth=1.8, label="逐点波形相对误差")
    axes[1, 0].axvline(1.0, color="gray", linestyle=":", linewidth=1.0, zorder=1)
    axes[1, 0].set(
        title="(c) 目标存在段：幅度缩放不变性对比",
        xlabel="输出幅度缩放系数 α",
        ylabel="损失或相对误差 / dB",
        ylim=(-24.0, 13.0),
    )
    axes[1, 0].legend(loc="upper left", framealpha=0.96)

    ratio_grid = np.linspace(-30.0, 0.0, 121, dtype=np.float64)
    sec_gate = (ratio_grid >= -18.0).astype(np.float64)
    axes[1, 1].step(ratio_grid, sec_gate, where="post", color=C_MAIN, linewidth=2.0, label="次槽独立保留标志")
    axes[1, 1].plot([-18.0, -18.0], [-0.10, 1.06], color=C_RED, linestyle="--", linewidth=1.6, zorder=2, label="旁路门限 -18 dB")
    axes[1, 1].plot([-20.0], [0.0], "o", color=C_GREEN, markersize=8, zorder=5, label="单人擦音合并 (-20 dB)")
    axes[1, 1].plot(
        [keep_overlap["energy_ratio_db"]],
        [1.0],
        "s",
        color=C_BLUE,
        markersize=8,
        zorder=5,
        label=f"双人重叠保留 ({keep_overlap['energy_ratio_db']:.2f} dB)",
    )
    axes[1, 1].set(
        title="(d) 零重叠旁路门限与擦音保护",
        xlabel="次槽对主槽能量比 / dB",
        ylabel="次槽激活状态 (0=合并, 1=分离)",
        yticks=[0, 1],
        ylim=(-0.18, 2.05),
    )
    axes[1, 1].legend(loc="upper left", handlelength=1.5, framealpha=0.96)

    for ax in axes.flat:
        ax.grid(True, alpha=0.22)
        ax.set_axisbelow(True)
    fig.suptitle("图80  条件目标说话人提取：特征调制、沉默态损失与单人旁路保护")
    _save_report_and_png(
        fig,
        "fig80_tse_conditioning_modes.png",
        "figure80_tse_conditioning_modes.json",
        report,
        manifest_sha,
    )


def _build_fig86(manifest: dict, sources: dict[str, str], manifest_sha: str) -> None:
    audio_dir = validate_parent_chain(DEFAULT_OUTPUT)
    radii = np.linspace(0.15, 5.00, 195, dtype=np.float64)
    bub_10 = broadside_distance_bubble_weights(0.10, 0.40, 3430.0, c=343.0, radii_m=radii)
    bub_02_reg = broadside_distance_bubble_weights(0.02, 0.40, 3430.0, c=343.0, reg=5e-4, radii_m=radii)

    w_dsb = np.full(3, 1.0 / 3.0, dtype=np.complex128)
    dsb_rad = [
        float(abs(np.vdot(w_dsb, spherical_ula_steering([-0.10, 0.0, 0.10], float(r), 0.0, 3430.0, c=343.0))) ** 2)
        for r in radii
    ]

    spacings = np.linspace(0.01, 0.15, 29, dtype=np.float64)
    wng_curve = [
        broadside_distance_bubble_weights(float(d), 0.40, 3430.0, c=343.0)["wng_db"]
        for d in spacings
    ]
    fresnel_curve = [fresnel_distance(float(2.0 * d), 3430.0, c=343.0) for d in spacings]

    signals = {
        key: read_pcm16((audio_dir / row["file"]).read_bytes())[1][0]
        for key, row in manifest["samples"].items()
        if row["file"] in ("tse_near_target_ref.wav", "tse_farfield_dsb.wav", "tse_bubble_d10cm.wav", "tse_bubble_d2cm_illcond.wav")
    }
    pcm_nmse = {
        key: manifest["samples"][key]["pcm_measurements"]["pcm_integer_measurements"]["integer_reference_nmse"]
        for key in ("farfield_dsb", "bubble_d10cm", "bubble_d2cm_illcond")
    }

    report = {
        "figure": 81,
        "exercise_ids": ["E15-03", "E15-04", "E15-05", "E15-12"],
        "source_sha256": sources,
        "audio_manifest_sha256": manifest_sha,
        "wav_sha256": {name: rec["sha256"] for name, rec in manifest["files"].items()},
        "radii_m": radii.tolist(),
        "radial_gain_db": {
            "dsb_d10cm": [float(10.0 * np.log10(max(v, 1e-6))) for v in dsb_rad],
            "bubble_d10cm": [float(10.0 * np.log10(max(v, 1e-6))) for v in bub_10["radial_power_response"]],
            "bubble_d2cm_reg5e4": [float(10.0 * np.log10(max(v, 1e-6))) for v in bub_02_reg["radial_power_response"]],
        },
        "spacings_m": spacings.tolist(),
        "wng_db_vs_spacing": wng_curve,
        "fresnel_m_vs_spacing": fresnel_curve,
        "pcm_integer_nmse": pcm_nmse,
        "limits": manifest["limits"],
    }

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.6), constrained_layout=True)
    dsb_db = report["radial_gain_db"]["dsb_d10cm"]
    axes[0, 0].plot(radii, dsb_db, "--", color=C_RED, marker="s", markevery=(10, 30), markersize=4.5, linewidth=1.8, label="远场 DSB (10 cm)")
    axes[0, 0].plot(radii, report["radial_gain_db"]["bubble_d10cm"], "-", color=C_BLUE, marker="o", markevery=(0, 30), markersize=4.5, linewidth=2.0, label="距离气泡 (10 cm)")
    axes[0, 0].plot(radii, report["radial_gain_db"]["bubble_d2cm_reg5e4"], "-.", color=C_ORANGE, marker="^", markevery=(20, 30), markersize=4.8, linewidth=1.8, label=r"正则气泡 (2 cm, $\delta=5\times 10^{-4}$)")
    axes[0, 0].axvline(0.40, color=C_GREEN, linestyle=":", linewidth=1.4, zorder=1, label="目标 0.4 m")
    axes[0, 0].axvline(4.00, color=C_PURPLE, linestyle=":", linewidth=1.4, zorder=1, label="干扰 4.0 m")
    axes[0, 0].plot(radii[0], dsb_db[0], "s", color=C_RED, markersize=5.5, zorder=5)
    axes[0, 0].annotate(
        f"r=0.15 m: {dsb_db[0]:.1f} dB",
        xy=(radii[0], dsb_db[0]),
        xytext=(0.75, -18.5),
        color=C_RED,
        arrowprops=dict(arrowstyle="->", color=C_RED, lw=1.1, connectionstyle="arc3,rad=-0.15"),
    )
    axes[0, 0].set_title("(a) 正横径向功率响应 (3430 Hz)")
    axes[0, 0].set(
        xlabel="声源径向距离 r / m",
        ylabel="功率增益 / dB",
        xlim=(0.05, 5.10),
        ylim=(-25.0, 23.0),
        yticks=[-20, -10, 0, 10, 20],
    )
    axes[0, 0].legend(
        ncol=2, loc="upper right", fontsize=8.2,
        framealpha=0.96, handlelength=1.25, handletextpad=0.35, columnspacing=0.6,
    )

    wng_10 = broadside_distance_bubble_weights(0.10, 0.40, 3430.0, c=343.0)["wng_db"]
    wng_02 = broadside_distance_bubble_weights(0.02, 0.40, 3430.0, c=343.0)["wng_db"]
    axes[0, 1].plot(spacings * 100.0, wng_curve, "o-", color=C_BLUE, markersize=4, linewidth=1.8, label="无正则 WNG")
    axes[0, 1].plot([10.0, 10.0], [-46.5, wng_10], color=C_GREEN, linestyle="--", linewidth=1.6, label=f"d=10 cm ({wng_10:.2f} dB)")
    axes[0, 1].plot([2.0, 2.0], [-46.5, wng_02], color=C_RED, linestyle="-.", linewidth=1.6, label=f"d=2 cm ({wng_02:.2f} dB)")
    axes[0, 1].plot([10.0, 2.0], [wng_10, wng_02], "D", color=C_MAIN, markerfacecolor="white", markeredgewidth=1.5, markersize=6, zorder=5)
    axes[0, 1].axhline(-10.0, color="0.55", linestyle=":", linewidth=1.1, zorder=1)
    axes[0, 1].text(15.2, -11.6, "−10 dB 参考线", ha="right", va="top", color="0.35")
    axes[0, 1].set(
        title="(b) 阵元间距与白噪声增益塌缩",
        xlabel="相邻阵元间距 d / cm",
        ylabel="白噪声增益 WNG / dB",
        xlim=(0.4, 15.6),
        ylim=(-46.5, 10.5),
    )
    axes[0, 1].legend(loc="lower right", framealpha=0.96)

    s0, count = 1600, 48
    t_ms = np.arange(count) / 16000.0 * 1000.0
    axes[1, 0].plot(t_ms, signals["near_target_ref"][s0:s0 + count], "-", color=C_GREEN, linewidth=2.0, label="近场目标参考")
    axes[1, 0].plot(t_ms, signals["farfield_dsb"][s0:s0 + count], "--", color=C_RED, linewidth=1.4, label="远场 DSB")
    axes[1, 0].plot(t_ms, signals["bubble_d10cm"][s0:s0 + count], "-.", color=C_BLUE, linewidth=1.6, label="10 cm 距离气泡")
    axes[1, 0].set(
        title="(c) 实际 PCM 波形对照（第 1600 点起）",
        xlabel="相对时间 / ms",
        ylabel="PCM 解码幅度",
        xlim=(-0.06, t_ms[-1] + 0.08),
        ylim=(-0.35, 0.66),
    )
    axes[1, 0].legend(loc="upper right", framealpha=0.96)

    bar_keys = ("farfield_dsb", "bubble_d10cm", "bubble_d2cm_illcond")
    bar_labels = ["远场 DSB\n(d=10 cm)", "距离气泡\n(d=10 cm)", "小孔径气泡\n(d=2 cm)"]
    bar_vals = [pcm_nmse[k] for k in bar_keys]
    bar_colors = [C_RED, C_BLUE, C_ORANGE]
    bar_hatches = ["///", "", "xxx"]
    for idx, (val, color, hatch) in enumerate(zip(bar_vals, bar_colors, bar_hatches)):
        axes[1, 1].bar(idx, val, color=color, edgecolor=C_MAIN, linewidth=0.85, hatch=hatch, width=0.55, alpha=0.85)
        axes[1, 1].text(idx, val + 0.08, f"{val:.4f}", ha="center")
    axes[1, 1].text(
        0.03,
        0.94,
        r"合成音频：$3200\ \mathrm{Hz}$，$\delta=2\times 10^{-4}$",
        transform=axes[1, 1].transAxes,
        ha="left",
        va="top",
        bbox=dict(fc="white", ec="0.75", alpha=0.92, boxstyle="round,pad=0.25"),
    )
    axes[1, 1].set(
        title="(d) 稳窗实际 PCM16 误差比 (3200 Hz)",
        xlabel="波束形成方案（共同导出增益 1.0）",
        ylabel=r"$E_{\mathrm{PCM}} / P_{\mathrm{ref}}$（无量纲）",
        xticks=np.arange(3),
        xticklabels=bar_labels,
        ylim=(0.0, 3.05),
    )

    for ax in axes.flat:
        ax.grid(True, alpha=0.22)
        ax.set_axisbelow(True)
    fig.suptitle("图81  近场球面波距离气泡：同向远场抑制与小孔径底噪边界")
    _save_report_and_png(
        fig,
        "fig81_nearfield_distance_bubble.png",
        "figure81_nearfield_distance_bubble.json",
        report,
        manifest_sha,
    )


def _build_fig87(sources: dict[str, str], manifest_sha: str) -> None:
    ref_act = np.array([[1, 1, 1, 1, 0, 0], [0, 0, 1, 1, 1, 0]], dtype=np.int64)
    hyp_oracle_single = np.array([[1, 1, 1, 1, 0, 0], [0, 0, 0, 0, 1, 0]], dtype=np.int64)
    hyp_imperfect = np.array([[1, 1, 1, 0, 1, 0], [0, 0, 0, 1, 0, 1]], dtype=np.int64)
    hyp_over100 = np.array(
        [[0, 0, 0, 0, 1, 1], [0, 0, 1, 0, 1, 1], [1, 1, 0, 0, 1, 1]], dtype=np.int64
    )

    obs_ipd_cos = np.array(
        [[0.95, 0.92, 0.82, 0.80, -0.40, 0.00], [-0.45, -0.40, 0.78, 0.81, 0.94, -0.05]],
        dtype=np.float64,
    )
    emb_sim = np.array(
        [[0.88, 0.86, 0.48, 0.47, 0.12, 0.08], [0.10, 0.14, 0.46, 0.49, 0.89, 0.06]],
        dtype=np.float64,
    )
    post_spec = tsvad_spatial_spectral_posterior(
        obs_ipd_cos, emb_sim, prior=0.5, kappa_sp=0.0, kappa_emb=4.0, tau_emb=0.45
    )
    post_joint = tsvad_spatial_spectral_posterior(
        obs_ipd_cos, emb_sim, prior=0.5, kappa_sp=3.0, tau_sp=0.25, kappa_emb=4.0, tau_emb=0.45
    )

    der_cases = {
        "single_bound": compute_der_components(ref_act, hyp_oracle_single, 0.5),
        "imperfect_2c": compute_der_components(ref_act, hyp_imperfect, 0.5),
        "over_100": compute_der_components(ref_act, hyp_over100, 0.5),
        "tsvad_joint": compute_der_components(ref_act, (post_joint >= 0.5).astype(np.int64), 0.5),
    }
    bound_check = single_label_overlap_der_bound(ref_act, 0.5)

    m_tar_excl = post_joint[0] * (1.0 - post_joint[1])
    m_interf = 1.0 - post_joint[0]

    report = {
        "figure": 82,
        "exercise_ids": ["E15-06", "E15-07", "E15-08", "E15-09", "E15-10"],
        "source_sha256": sources,
        "audio_manifest_sha256": manifest_sha,
        "ref_activity": ref_act.tolist(),
        "single_label_bound": bound_check,
        "der_cases": der_cases,
        "posteriors_spectral_only": post_spec.tolist(),
        "posteriors_joint": post_joint.tolist(),
        "gss_target_exclusive_mask": m_tar_excl.tolist(),
        "gss_interference_mask": m_interf.tolist(),
    }

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.6), constrained_layout=True)
    t_centers = np.arange(6) * 0.5 + 0.25
    axes[0, 0].bar(t_centers - 0.09, ref_act[0], width=0.16, color=C_BLUE, edgecolor=C_MAIN, linewidth=0.75, label="说话人 1 真值")
    axes[0, 0].bar(t_centers + 0.09, ref_act[1], width=0.16, color=C_ORANGE, edgecolor=C_MAIN, linewidth=0.75, hatch="///", label="说话人 2 真值")
    axes[0, 0].fill_between([1.0, 2.0], 0.0, 1.08, color=C_RED, alpha=0.14, zorder=0, label="双讲重叠区 (必漏 1.0 s)")
    axes[0, 0].set(
        title="(a) 六帧重叠会议真值活动（总参考 3.5 s）",
        xlabel="会议时间 / s（每帧 0.5 s）",
        ylabel="说话活动标志 (0/1)",
        xticks=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        yticks=[0, 1],
        ylim=(0.0, 1.72),
    )
    axes[0, 0].legend(loc="upper right", framealpha=0.96)

    case_keys = ("single_bound", "imperfect_2c", "over_100", "tsvad_joint")
    x_pos = np.arange(4)
    miss_pct = [100.0 * der_cases[k]["miss_s"] / 3.5 for k in case_keys]
    fa_pct = [100.0 * der_cases[k]["false_alarm_s"] / 3.5 for k in case_keys]
    conf_pct = [100.0 * der_cases[k]["speaker_confusion_s"] / 3.5 for k in case_keys]
    axes[0, 1].bar(x_pos, miss_pct, width=0.52, color=C_BLUE, edgecolor=C_MAIN, linewidth=0.75, label="漏检 Miss")
    axes[0, 1].bar(x_pos, fa_pct, width=0.52, bottom=miss_pct, color=C_RED, edgecolor=C_MAIN, linewidth=0.75, hatch="///", label="虚警 FA")
    bottom_conf = [m + f for m, f in zip(miss_pct, fa_pct)]
    axes[0, 1].bar(x_pos, conf_pct, width=0.52, bottom=bottom_conf, color=C_ORANGE, edgecolor=C_MAIN, linewidth=0.75, hatch="xxx", label="混淆 Conf")
    axes[0, 1].axhline(100.0, color="black", linestyle="--", linewidth=1.1, label="100% 基准线")
    for idx, k in enumerate(case_keys):
        tot = 100.0 * der_cases[k]["der"]
        axes[0, 1].text(idx, tot + 3.5, f"{tot:.1f}%", ha="center")
    axes[0, 1].set(
        title="(b) 四种输出假设的 DER 三分量分解",
        xlabel="日志假设类型",
        ylabel="DER 分量占比 / %",
        xticks=x_pos,
        xticklabels=["单标签\n下界", "含混淆\n双簇", "过触发\n(>100%)", "空谱联合\nTS-VAD"],
        ylim=(0.0, 178.0),
    )
    axes[0, 1].legend(ncol=2, loc="upper left", framealpha=0.96)

    frames = np.arange(6)
    axes[1, 0].plot(frames, post_spec[0], "o--", color=C_RED, linewidth=1.6, label="仅声纹 Spk1")
    axes[1, 0].plot(frames, post_joint[0], "s-", color=C_BLUE, linewidth=2.0, label="空谱联合 Spk1")
    axes[1, 0].plot(frames, post_joint[1], "^-.", color=C_GREEN, linewidth=1.8, label="空谱联合 Spk2")
    axes[1, 0].axhline(0.5, color="gray", linestyle=":", linewidth=1.0)
    axes[1, 0].set_title("(c) 重叠帧上 TS-VAD 后验概率恢复", pad=52)
    axes[1, 0].set(
        xlabel="帧索引 $t$",
        ylabel=r"活动后验概率 $\gamma$",
        xticks=frames,
        yticks=[0.0, 0.25, 0.5, 0.75, 1.0],
        xlim=(-0.25, 5.25),
        ylim=(-0.05, 1.06),
    )
    axes[1, 0].legend(
        ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.02),
        framealpha=0.96, handlelength=1.5, columnspacing=0.9,
    )

    axes[1, 1].plot(frames, m_tar_excl, "s-", color=C_BLUE, linewidth=2.0, label=r"目标独占 $m_{\mathrm{tar}}$")
    axes[1, 1].plot(frames, m_interf, "o--", color=C_RED, linewidth=1.8, label=r"干扰加噪 $m_{\mathrm{int}}$")
    axes[1, 1].plot(frames, post_joint[0], "^:", color=C_PURPLE, linewidth=1.5, label=r"未排他 $\gamma_1$")
    axes[1, 1].set_title("(d) GSS 协方差门控权重（提取说话人 1）", pad=52)
    axes[1, 1].set(
        xlabel="帧索引 $t$",
        ylabel="协方差累加帧权重",
        xticks=frames,
        yticks=[0.0, 0.25, 0.5, 0.75, 1.0],
        xlim=(-0.25, 5.25),
        ylim=(-0.05, 1.06),
    )
    axes[1, 1].legend(
        ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.02),
        framealpha=0.96, handlelength=1.5, columnspacing=0.9,
    )

    for ax in axes.flat:
        ax.grid(True, alpha=0.22)
        ax.set_axisbelow(True)
    fig.suptitle("图82  多人重叠会议日志：DER 下界、空谱 TS-VAD 与 GSS 掩码")
    _save_report_and_png(
        fig,
        "fig82_diarization_tsvad_gss.png",
        "figure82_diarization_tsvad_gss.json",
        report,
        manifest_sha,
    )


def _build_fig88(sources: dict[str, str], manifest_sha: str) -> None:
    sync_map = map_video_to_audio_frames(15, audio_hop_s=0.016, video_fps=25.0, av_offset_s=0.0)
    lag_map = map_video_to_audio_frames(15, audio_hop_s=0.016, video_fps=25.0, av_offset_s=0.120)

    lags_ms = np.linspace(0.0, 160.0, 33, dtype=np.float64)
    loss_pct = []
    for lag_ms in lags_ms:
        res = map_video_to_audio_frames(
            30, audio_hop_s=0.016, video_fps=25.0, av_offset_s=float(lag_ms / 1000.0), onset_frame=10, consonant_frames=8
        )
        loss_pct.append(100.0 * res["onset_audit"]["consonant_loss_fraction"])

    t_frames = np.arange(12)
    alpha_traj = np.array([0.92, 0.90, 0.88, 0.85, 0.55, 0.35, 0.12, 0.05, 0.05, 0.08, 0.65, 0.88])
    m_aud_traj = np.full(12, 0.82, dtype=np.float64)
    m_av_traj = np.array([0.96, 0.95, 0.94, 0.93, 0.72, 0.58, 0.15, 0.04, 0.05, 0.10, 0.78, 0.94])
    m_fused_traj = av_confidence_gated_mask(m_aud_traj, m_av_traj, alpha_traj)

    arch_names = ["流式条件 TSE\n(声纹前向)", "联合 SOT/DiCoW\n(分块推理)", "离线级联\nTS-VAD+GSS+ASR"]
    arch_latency_s = [0.032, 1.20, 4.50]

    report = {
        "figure": 83,
        "exercise_ids": ["E15-11"],
        "source_sha256": sources,
        "audio_manifest_sha256": manifest_sha,
        "sync_video_indices": sync_map["video_frame_indices"].tolist(),
        "lag120ms_video_indices": lag_map["video_frame_indices"].tolist(),
        "lags_ms": lags_ms.tolist(),
        "consonant_loss_percent": loss_pct,
        "occlusion_alpha_v": alpha_traj.tolist(),
        "occlusion_mask_av": m_av_traj.tolist(),
        "occlusion_mask_fused": m_fused_traj.tolist(),
        "architecture_chunk_latency_s": dict(zip(["streaming_tse", "chunked_joint", "offline_gss"], arch_latency_s)),
    }

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.5), constrained_layout=True)
    aud_idx = np.arange(15)
    axes[0, 0].step(
        aud_idx, sync_map["video_frame_indices"], where="mid",
        marker="o", markersize=4.5, color=C_BLUE, linewidth=2.0, label="零滞后 (0 ms)",
    )
    axes[0, 0].step(
        aud_idx, lag_map["video_frame_indices"], where="mid",
        marker="s", markersize=4.5, linestyle="--", color=C_RED, linewidth=1.8, label="相机滞后 (+120 ms)",
    )
    axes[0, 0].set(
        title="(a) 16 ms 音频帧与 25 fps 视频曝光映射",
        xlabel=r"音频 STFT 帧索引 $\ell$",
        ylabel="关联视频曝光帧索引 $m$",
        yticks=range(6),
        ylim=(-0.4, 5.8),
    )
    axes[0, 0].legend(loc="upper left", framealpha=0.96)

    axes[0, 1].plot(lags_ms, loss_pct, "o-", color=C_RED, markersize=4, linewidth=1.8, label="词首辅音截断率")
    axes[0, 1].axvline(120.0, color=C_BLUE, linestyle="--", linewidth=1.5, label="+120 ms 滞后 (87.5%)")
    axes[0, 1].set(
        title="(b) 视觉正滞后引起的词首辅音丢失",
        xlabel=r"视频相对音频滞后 $\Delta t_{\mathrm{av}}$ / ms",
        ylabel="128 ms 词首辅音被抑制比例 / %",
        ylim=(-5.0, 108.0),
    )
    axes[0, 1].legend(loc="upper left", framealpha=0.96)

    axes[1, 0].plot(t_frames, alpha_traj, "^:", color=C_PURPLE, linewidth=1.6, label=r"视觉置信度 $\alpha_\ell$")
    axes[1, 0].plot(t_frames, m_av_traj, "o--", color=C_RED, linewidth=1.7, label=r"未门控 $M_{\mathrm{av}}$")
    axes[1, 0].plot(t_frames, m_fused_traj, "s-", color=C_BLUE, linewidth=2.0, label=r"门控融合 $M_{\mathrm{fused}}$")
    axes[1, 0].axvspan(6.0, 9.0, color=C_ORANGE, alpha=0.15, zorder=0, label="手部遮挡嘴唇区")
    axes[1, 0].set_title("(c) 唇部遮挡下的置信度门控平滑回退", pad=12)
    axes[1, 0].set(
        xlabel=r"时间块索引 $\ell$",
        ylabel="掩码增益或置信度",
        yticks=[0.0, 0.25, 0.5, 0.75, 1.0],
        ylim=(-0.04, 1.44),
    )
    axes[1, 0].legend(
        ncol=2,
        loc="upper center",
        framealpha=0.96,
        handlelength=1.4,
        handletextpad=0.35,
        columnspacing=0.75,
        borderpad=0.3,
    )

    y_pos = np.arange(3)
    arch_colors = [C_GREEN, C_BLUE, C_ORANGE]
    arch_hatches = ["", "///", "xxx"]
    for idx, (val, color, hatch) in enumerate(zip(arch_latency_s, arch_colors, arch_hatches)):
        axes[1, 1].barh(
            idx, val, color=color, edgecolor=C_MAIN, linewidth=0.85,
            hatch=hatch, height=0.52, alpha=0.85,
        )
        axes[1, 1].text(val * 1.18, idx, f"{val:.3g} s", va="center")
    axes[1, 1].set_xscale("log")
    axes[1, 1].set_title("(d) 三类前端算法等待窗口（典型量级示意）", pad=12)
    axes[1, 1].set(
        xlabel="典型最小上下文等待窗口 / s（对数轴）",
        yticks=y_pos,
        yticklabels=arch_names,
        xlim=(0.01, 20.0),
    )

    for ax in axes.flat:
        ax.grid(True, alpha=0.22)
        ax.set_axisbelow(True)
    fig.suptitle("图83  远场视听多模态与多说话人架构：异构帧率、滞后截断与遮挡回退")
    _save_report_and_png(
        fig,
        "fig83_av_asynchrony_occlusion.png",
        "figure83_av_asynchrony_occlusion.json",
        report,
        manifest_sha,
    )


def main() -> int:
    audio_dir = validate_parent_chain(DEFAULT_OUTPUT)
    manifest = check_assets(audio_dir)
    sources = {path: _sha_file(path) for path in SCRIPT_INPUT_PATHS}
    manifest_sha = hashlib.sha256((audio_dir / "MANIFEST.json").read_bytes()).hexdigest()

    _build_fig85(sources, manifest_sha)
    _build_fig86(manifest, sources, manifest_sha)
    _build_fig87(sources, manifest_sha)
    _build_fig88(sources, manifest_sha)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
