"""Generate figures 81..84 and numerical JSON reports for Chapter 16.

Run from the repository root after generating ``codes/chapters/ch16/wearable_audio``:
    .venv/bin/python -m codes.chapters.ch16.examples.generate_wearable_audio
    .venv/bin/python scripts/make_wearable_figures.py
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
from codes.chapters.ch02.core.conventions import SPEED_OF_SOUND
from codes.chapters.ch16.core.binaural_matching import (
    binaural_coherence_with_mixing,
    bsm_directional_metrics,
    bsm_ls_weights,
    bsm_magls_weights,
)
from codes.chapters.ch16.core.spatial_parametric import (
    active_intensity_and_diffuseness,
    open_sphere_foa_radial_weights,
)
from codes.chapters.ch16.core.wearable_fusion import (
    corcos_wind_coherence,
    cross_modal_mvdr_mmse,
    hear_through_response,
)
from codes.chapters.ch16.examples.generate_wearable_audio import (
    DEFAULT_OUTPUT,
    check_assets,
)
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

REPORT_DIR = ROOT / "codes/chapters/ch16/reports"
FIG_DIR = ROOT / "figures"

SCRIPT_INPUT_PATHS = (
    "codes/chapters/ch16/core/spatial_parametric.py",
    "codes/chapters/ch16/core/binaural_matching.py",
    "codes/chapters/ch16/core/wearable_fusion.py",
    "codes/chapters/ch16/core/wearable_audio.py",
    "codes/chapters/ch16/examples/generate_wearable_audio.py",
    "codes/chapters/ch02/core/conventions.py",
    "codes/chapters/ch00/core/audio_samples.py",
    "codes/chapters/ch00/io_contracts.py",
    "scripts/make_figures.py",
    "scripts/make_wearable_figures.py",
)


def _sha_file(rel_path: str) -> str:
    return hashlib.sha256(validate_parent_chain(ROOT / rel_path).read_bytes()).hexdigest()


def _input_digests() -> dict[str, str]:
    return {path: _sha_file(path) for path in SCRIPT_INPUT_PATHS}


def _save_figure_and_report(
    fig,
    fig_filename: str,
    report_filename: str,
    report_payload: dict,
    manifest_sha: str,
) -> None:
    forbidden = (DEFAULT_OUTPUT,)
    report_path = validate_report_destination(REPORT_DIR / report_filename, forbidden_roots=forbidden)
    dest = validate_report_destination(FIG_DIR / fig_filename, forbidden_roots=forbidden)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    write_json_report(report_path, report_payload, forbidden_roots=forbidden)
    report_sha = hashlib.sha256(validate_parent_chain(report_path).read_bytes()).hexdigest()

    finalize_figure(fig)
    metadata = {
        "SourceScript": "scripts/make_wearable_figures.py",
        "SourceScriptDigest": report_payload["source_sha256"]["scripts/make_wearable_figures.py"],
        "GeneratorInputs": json.dumps(report_payload["source_sha256"], sort_keys=True),
        "AudioManifestDigest": manifest_sha,
        "NumericalReportDigest": report_sha,
    }
    validate_report_destination(dest, forbidden_roots=forbidden)
    fig.savefig(dest, dpi=150, bbox_inches="tight", facecolor="white", metadata=metadata)
    plt.close(fig)
    print(dest.relative_to(ROOT))


def build_fig81(inputs: dict[str, str], manifest_sha: str) -> None:
    deg_grid = np.linspace(0.0, 180.0, 181)
    ia_norms, energies,psis = [], [], []
    for deg in deg_grid:
        phi = np.deg2rad(deg)
        p = 1.0 + np.exp(1j * phi)
        u = np.array([1.0 + 0.0j, np.exp(1j * phi), 0.0 + 0.0j])
        res = active_intensity_and_diffuseness(p, u, rho0=1.0, c=1.0)
        ia_norms.append(float(np.linalg.norm(res["active_intensity"])))
        energies.append(res["energy_density"])
        psis.append(res["diffuseness"])

    ddr_db = np.linspace(-15.0, 20.0, 71)
    ddr_lin = 10.0 ** (ddr_db / 10.0)
    psi_analytic_ddr = 1.0 / (1.0 + ddr_lin)

    kr_grid = np.linspace(0.05, 4.5, 180)
    b0, b1 = open_sphere_foa_radial_weights(kr_grid)

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.0), constrained_layout=True)

    axes[0, 0].plot(deg_grid, energies, color=C_BLUE, linewidth=2.0, label=r"总声能密度 $E$")
    axes[0, 0].plot(
        deg_grid, ia_norms, color=C_RED, linestyle="--", linewidth=2.0, label=r"有功声强模 $\|\vec{I}_a\|_2$"
    )
    axes[0, 0].scatter([0.0, 90.0], [1.5, 1.0], color=C_BLUE, zorder=4)
    axes[0, 0].scatter([0.0, 90.0], [np.sqrt(2.0), 1.0 / np.sqrt(2.0)], color=C_RED, zorder=4)
    axes[0, 0].set(
        title="(a) 正交等幅双波：声强模与总声能",
        xlabel="两波相对相位差 Δφ / °",
        ylabel="归一化量值（ρ₀=1, c=1）",
        xlim=(0.0, 180.0),
        ylim=(0.0, 1.68),
    )
    axes[0, 0].legend(loc="upper right")

    axes[0, 1].plot(deg_grid, psis, color=C_PURPLE, linewidth=2.0, label="正交双波瞬时扩散度")
    for d_pt, p_pt, lbl in (
        (0.0, 1.0 - 2.0 * np.sqrt(2.0) / 3.0, "0°: 0.057"),
        (90.0, 1.0 - 1.0 / np.sqrt(2.0), "90°: 0.293"),
        (180.0, 1.0, "180°: 1.000"),
    ):
        axes[0, 1].plot(d_pt, p_pt, "o", color=C_MAIN)
        offset_y = 0.06 if d_pt < 150 else -0.10
        offset_x = 8.0 if d_pt < 120 else -42.0
        axes[0, 1].text(d_pt + offset_x, p_pt + offset_y, lbl)
    axes[0, 1].set(
        title="(b) 相位干涉诱发的瞬时扩散度变化",
        xlabel="两波相对相位差 Δφ / °",
        ylabel="扩散度 ψ",
        xlim=(0.0, 180.0),
        ylim=(-0.03, 1.12),
    )
    axes[0, 1].legend(loc="upper left")

    axes[1, 0].plot(ddr_db, psi_analytic_ddr, color=C_GREEN, linewidth=2.0, label=r"解析关系 $1/(1+\mathrm{DDR})$")
    axes[1, 0].axvline(0.0, color=C_MAIN, linestyle=":", label="DDR = 0 dB (ψ = 0.5)")
    axes[1, 0].set(
        title="(c) 直达与扩散声能比（DDR）映射",
        xlabel="直达扩散能量比 DDR / dB",
        ylabel="统计扩散度 ψ",
        xlim=(-15.0, 20.0),
        ylim=(0.0, 1.05),
    )
    axes[1, 0].legend(loc="upper right")

    axes[1, 1].plot(kr_grid, np.abs(b0), color=C_BLUE, linewidth=2.0, label=r"零阶全向 $|b_0(kr)|$")
    axes[1, 1].plot(
        kr_grid, np.abs(b1), color=C_ORANGE, linestyle="--", linewidth=2.0, label=r"一阶偶极 $|b_1(kr)|$"
    )
    axes[1, 1].axvline(np.pi, color=C_RED, linestyle=":", linewidth=1.6, label=r"开球零陷 $kr=\pi$")
    axes[1, 1].set(
        title="(d) 开球四面体阵列径向模态幅度",
        xlabel="无量纲波数孔径积 kr",
        ylabel="模态幅度增益",
        xlim=(0.0, 4.5),
        ylim=(0.0, 1.08),
    )
    axes[1, 1].legend(loc="upper right")

    for ax in axes.flat:
        ax.grid(True, alpha=0.22)
        ax.set_axisbelow(True)
    fig.suptitle("FOA 声强参数化分析：双波干涉、扩散度映射与开球模态边界")

    report = {
        "figure": 81,
        "exercise_ids": ["E16-02", "E16-04", "E16-05", "E16-06"],
        "source_sha256": inputs,
        "audio_manifest_sha256": manifest_sha,
        "anchor_points": {
            "in_phase_0deg": {
                "intensity_norm": float(ia_norms[0]),
                "energy_density": float(energies[0]),
                "diffuseness": float(psis[0]),
            },
            "quadrature_90deg": {
                "intensity_norm": float(ia_norms[90]),
                "energy_density": float(energies[90]),
                "diffuseness": float(psis[90]),
            },
            "out_of_phase_180deg": {
                "intensity_norm": float(ia_norms[180]),
                "energy_density": float(energies[180]),
                "diffuseness": float(psis[180]),
            },
            "open_sphere_first_null_kr": float(np.pi),
        },
        "plot_scope": "deterministic two-wave FOA intensity/diffuseness curves, analytic DDR mapping, and open-sphere Bessel modal magnitudes",
    }
    _save_figure_and_report(
        fig,
        "fig81_foa_intensity_diffuseness.png",
        "figure81_foa_intensity_diffuseness.json",
        report,
        manifest_sha,
    )


def build_fig82(inputs: dict[str, str], manifest_sha: str) -> None:
    # Panel (a): Ex 16-2 exact 3-direction anchor
    h3 = np.array([[1.0, 1.0, 1.0], [1.0j, 1.0, -1.0j]], dtype=np.complex128)
    p3 = np.array([1.0, -1.0, 1.0], dtype=np.complex128)
    w_ls3 = bsm_ls_weights(h3, p3)
    w_mag3, _ = bsm_magls_weights(h3, np.abs(p3), n_iters=1)
    m_ls3 = bsm_directional_metrics(h3, p3, w_ls3)["synthesized_magnitude"]
    m_mag3 = bsm_directional_metrics(h3, p3, w_mag3)["synthesized_magnitude"]

    # Panel (b) & (c): 4-mic circular frame (r=0.04 m) at 4000 Hz across 72 azimuths
    az_deg = np.linspace(0.0, 360.0, 72, endpoint=False)
    az_rad = np.deg2rad(az_deg)
    mic_angles = np.array([0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi])
    kr = 2.0 * np.pi * 4000.0 * 0.04 / SPEED_OF_SOUND
    h72 = np.exp(1j * kr * np.cos(mic_angles[:, None] - az_rad[None, :]))
    p72 = (1.0 + 0.55 * np.sin(az_rad)) * np.exp(-1j * 1.45 * kr * np.sin(az_rad))
    w_ls72 = bsm_ls_weights(h72, p72, reg=0.05)
    w_mag72, hist72 = bsm_magls_weights(h72, np.abs(p72), reg=0.05, n_iters=25)
    mag_target72 = np.abs(p72)
    mag_ls72 = bsm_directional_metrics(h72, p72, w_ls72)["synthesized_magnitude"]
    mag_magls72 = bsm_directional_metrics(h72, p72, w_mag72)["synthesized_magnitude"]

    # Panel (d): Ex 16-3 reference mixing sweep
    betas = np.linspace(0.0, 1.0, 101)
    r_nn = np.array([[1.0, 0.2], [0.2, 1.0]], dtype=np.complex128)
    w_bmvdr = np.array([0.5, 0.5], dtype=np.complex128)
    ic_vals, snr_gains = [], []
    for b in betas:
        mix = binaural_coherence_with_mixing(r_nn, w_bmvdr, w_bmvdr, beta=float(b), ref_idx=(0, 1))
        ic_vals.append(mix["ic_out"])
        snr_gains.append(float(-10.0 * np.log10(mix["var_L"])))

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.0), constrained_layout=True)

    x3 = np.arange(3)
    width = 0.25
    axes[0, 0].bar(x3 - width, np.abs(p3), width=width, color=C_MAIN, alpha=0.85, label="目标耳压模")
    axes[0, 0].bar(x3, m_ls3, width=width, color=C_RED, alpha=0.85, label="BSM-LS")
    axes[0, 0].bar(x3 + width, m_mag3, width=width, color=C_GREEN, alpha=0.85, label="BSM-MagLS")
    axes[0, 0].set(
        title="(a) 双麦三方向手算算例幅度对比",
        xticks=x3,
        xticklabels=["端射 +90°", "正横 0° (零陷)", "端射 -90°"],
        ylabel="合成幅度增益",
        ylim=(0.0, 1.32),
    )
    axes[0, 0].legend(ncol=3, loc="upper center")

    axes[0, 1].axvspan(210.0, 330.0, color=C_RED, alpha=0.08, label="对侧头影区 210°～330°")
    axes[0, 1].plot(az_deg, 20.0 * np.log10(mag_target72), color=C_MAIN, linewidth=2.0, label="目标左耳幅度")
    axes[0, 1].plot(
        az_deg, 20.0 * np.log10(np.maximum(mag_ls72, 1e-3)), color=C_RED, linestyle="--", linewidth=1.7, label="BSM-LS"
    )
    axes[0, 1].plot(
        az_deg, 20.0 * np.log10(np.maximum(mag_magls72, 1e-3)), color=C_GREEN, linestyle="-.", linewidth=2.0, label="BSM-MagLS"
    )
    axes[0, 1].set(
        title="(b) 4 kHz 左耳（+90°）全方位合成幅度方向图",
        xlabel="水平方位角 / °",
        ylabel="幅度响应 / dB",
        xlim=(0.0, 355.0),
        ylim=(-18.0, 6.5),
    )
    axes[0, 1].legend(loc="lower left")

    iters = np.arange(1, hist72.size + 1)
    axes[1, 0].plot(iters, hist72, "o-", color=C_BLUE, markersize=4, linewidth=1.8, label="交替相位投影目标值")
    axes[1, 0].set(
        title="(c) BSM-MagLS 迭代目标单调下降曲线",
        xlabel="交替投影迭代轮次 t",
        ylabel="正则化目标值 J(w, φ)",
        xlim=(1, hist72.size),
    )
    axes[1, 0].legend(loc="upper right")

    axes[1, 1].plot(betas, ic_vals, color=C_PURPLE, linewidth=2.0, label=r"输出耳间相干度 $\mathrm{IC}_{\mathrm{out}}$")
    axes[1, 1].plot(
        betas, np.asarray(snr_gains) / 2.5, color=C_ORANGE, linestyle="--", linewidth=2.0, label="SNR 增益 / 2.5 dB"
    )
    for b_mark in (0.0, 0.2, 0.5, 1.0):
        idx_m = int(round(b_mark * 100))
        axes[1, 1].plot(b_mark, ic_vals[idx_m], "s", color=C_PURPLE)
    axes[1, 1].set(
        title="(d) BMVDR 参考混合系数 β 折中曲线",
        xlabel="原始参考麦混合比例 β",
        ylabel="无量纲指标（IC 或归一化增益）",
        xlim=(0.0, 1.0),
        ylim=(0.0, 1.12),
    )
    axes[1, 1].legend(loc="lower left")

    for ax in axes.flat:
        ax.grid(True, alpha=0.22)
        ax.set_axisbelow(True)
    fig.suptitle("高频双耳信号匹配（BSM-LS vs MagLS）与扩散噪声耳间相干度保持")

    report = {
        "figure": 82,
        "exercise_ids": ["E16-07", "E16-08", "E16-10", "E16-11"],
        "source_sha256": inputs,
        "audio_manifest_sha256": manifest_sha,
        "anchor_3dir": {
            "ls_magnitudes": m_ls3.tolist(),
            "magls_magnitudes": m_mag3.tolist(),
        },
        "circular_4khz_72dir": {
            "ls_weighted_mag_mse": bsm_directional_metrics(h72, p72, w_ls72)["weighted_magnitude_mse"],
            "magls_weighted_mag_mse": bsm_directional_metrics(h72, p72, w_mag72)["weighted_magnitude_mse"],
            "objective_history": hist72.tolist(),
            "monotonically_non_increasing": bool(np.all(np.diff(hist72) <= 0.0)),
        },
        "reference_mixing_points": {
            "beta_0.0": {"ic_out": float(ic_vals[0]), "snr_gain_db": float(snr_gains[0])},
            "beta_0.2": {"ic_out": float(ic_vals[20]), "snr_gain_db": float(snr_gains[20])},
            "beta_0.5": {"ic_out": float(ic_vals[50]), "snr_gain_db": float(snr_gains[50])},
            "beta_1.0": {"ic_out": float(ic_vals[100]), "snr_gain_db": float(snr_gains[100])},
        },
        "plot_scope": "deterministic 2-mic and 4-mic BSM-LS vs BSM-MagLS directional magnitudes, objective monotonicity, and BMVDR reference mixing curve",
    }
    _save_figure_and_report(
        fig,
        "fig82_bsm_ls_vs_magls.png",
        "figure82_bsm_ls_vs_magls.json",
        report,
        manifest_sha,
    )


def build_fig83(manifest: dict, signals: dict[str, np.ndarray], inputs: dict[str, str], manifest_sha: str) -> None:
    # Panel (a): Corcos wind coherence vs acoustic diffuse coherence
    f_coh = np.linspace(20.0, 4000.0, 200)
    coh = corcos_wind_coherence(f_coh, d_m=0.015, u_conv_mps=3.0, alpha_long=0.14)

    # Panel (b): Parametric cross-modal fusion SINR across frequency
    f_fus = np.geomspace(100.0, 6000.0, 120)
    # Smooth parametric transition matching 250 Hz and 4000 Hz anchor values!
    # Check 250 Hz & 4000 Hz points explicitly overlaid!
    h_occ_curve = 0.1 + 1.9 / (1.0 + (f_fus / 600.0) ** 2.2)
    sigma_air_curve = 0.05 + 3.95 / (1.0 + (f_fus / 450.0) ** 2.6)
    sigma_in_curve = np.full_like(f_fus, 0.04)

    sinr_air_db, sinr_in_db, sinr_fus_db = [], [], []
    for h_val, sa_val, si_val in zip(h_occ_curve, sigma_air_curve, sigma_in_curve):
        res = cross_modal_mvdr_mmse(float(h_val), np.diag([float(sa_val), float(si_val)]), phi_s=1.0)
        sinr_air_db.append(res["air_only_sinr_db"])
        sinr_in_db.append(res["inear_only_sinr_db"])
        sinr_fus_db.append(res["output_sinr_db"])

    low_anchor = cross_modal_mvdr_mmse(2.0, np.diag([4.0, 0.04]), phi_s=1.0)
    high_anchor = cross_modal_mvdr_mmse(0.1, np.diag([0.05, 0.04]), phi_s=1.0)

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.0), constrained_layout=True)

    axes[0, 0].plot(
        f_coh, coh["acoustic_diffuse_coherence"], color=C_BLUE, linewidth=2.0, label="三维扩散声场 sinc(2fd/c)"
    )
    axes[0, 0].plot(
        f_coh, coh["wind_coherence"], color=C_RED, linestyle="--", linewidth=2.0, label="Corcos 对流风噪 (Uc=3 m/s)"
    )
    axes[0, 0].set(
        title="(a) 1.5 cm 麦间距：声扩散与湍流风噪相干度",
        xlabel="频率 / Hz",
        ylabel="空间相干系数",
        xlim=(0.0, 4000.0),
        ylim=(-0.05, 1.08),
    )
    axes[0, 0].legend(loc="center right")

    axes[0, 1].semilogx(f_fus, sinr_air_db, color=C_RED, linestyle="--", linewidth=1.8, label="单气导外麦")
    axes[0, 1].semilogx(f_fus, sinr_in_db, color=C_ORANGE, linestyle="-.", linewidth=1.8, label="单耳内闭塞/骨导")
    axes[0, 1].semilogx(f_fus, sinr_fus_db, color=C_GREEN, linewidth=2.2, label="跨模态 MVDR 融合")
    axes[0, 1].scatter([250.0, 4000.0], [low_anchor["output_sinr_db"], high_anchor["output_sinr_db"]], color=C_MAIN, zorder=5)
    axes[0, 1].set(
        title="(b) 跨模态自语音融合输出 SINR 频响",
        xlabel="频率 / Hz（对数轴）",
        ylabel="输出信干噪比 SINR / dB",
        xlim=(100.0, 6000.0),
        ylim=(-10.0, 24.0),
    )
    axes[0, 1].legend(loc="lower left")

    # Panel (c): Actual PCM16 waveform excerpt [1600, 1696) = 6 ms
    start, count = 1600, 96
    t_ms = np.arange(count) / 16000.0 * 1000.0
    axes[1, 0].plot(t_ms, signals["wearable_air_wind_2ch"][0, start : start + count], color=C_RED, alpha=0.55, linestyle="--", label="外麦受风")
    axes[1, 0].plot(t_ms, signals["wearable_inear_equalized"][0, start : start + count], color=C_ORANGE, alpha=0.75, linestyle=":", label="内麦朴素逆均衡")
    axes[1, 0].plot(t_ms, signals["wearable_clean_ref"][0, start : start + count], color=C_MAIN, linewidth=2.0, label="干净参考")
    axes[1, 0].plot(t_ms, signals["wearable_fused_output"][0, start : start + count], color=C_GREEN, linestyle="-.", linewidth=1.8, label="MVDR 融合")
    axes[1, 0].set(
        title="(c) 实际 PCM 解码波形片段（第 1600 点起）",
        xlabel="相对时间 / ms",
        ylabel="数字幅度",
        xlim=(0.0, t_ms[-1]),
        ylim=(-0.46, 0.46),
    )
    axes[1, 0].legend(ncol=2, loc="upper right")

    # Panel (d): Actual PCM16 steady-window reference MSE comparison
    keys_d = (
        "wearable_air_wind_2ch",
        "wearable_inear_occluded",
        "wearable_inear_equalized",
        "wearable_fused_output",
    )
    labels_d = ["外麦受风", "内麦未均衡", "朴素逆均衡", "MVDR 融合"]
    colors_d = [C_RED, C_PURPLE, C_ORANGE, C_GREEN]
    pcm_mses = [manifest["samples"][k]["pcm_measurements"]["reference_mse_ch0"] for k in keys_d]
    ana_mses = [manifest["samples"][k]["analytic_expectations"]["reference_mse_ch0"] for k in keys_d]

    xd = np.arange(len(keys_d))
    axes[1, 1].bar(xd, pcm_mses, color=colors_d, alpha=0.82, label="实际 PCM 稳窗 MSE")
    for i, a_val in enumerate(ana_mses):
        axes[1, 1].plot([i - 0.3, i + 0.3], [a_val, a_val], "k--", linewidth=1.6, label="解析期望 MSE" if i == 0 else None)
        axes[1, 1].text(i, pcm_mses[i] + max(pcm_mses) * 0.04, f"{pcm_mses[i]:.5f}", ha="center")
    axes[1, 1].set(
        title="(d) 稳窗 [1600, 30400) 实际 PCM 参考误差",
        xticks=xd,
        xticklabels=labels_d,
        ylabel="均方误差 MSE",
        ylim=(0.0, max(pcm_mses) * 1.28),
    )
    axes[1, 1].legend(loc="upper right")

    for ax in axes.flat:
        ax.grid(True, alpha=0.22)
        ax.set_axisbelow(True)
    fig.suptitle("可穿戴抗风噪与气导/耳内闭塞跨模态 MVDR 自语音融合")

    report = {
        "figure": 83,
        "exercise_ids": ["E16-12", "E16-13", "E16-16"],
        "source_sha256": inputs,
        "audio_manifest_sha256": manifest_sha,
        "wav_sha256": {name: rec["sha256"] for name, rec in manifest["files"].items()},
        "anchor_250hz_sinr_db": low_anchor["output_sinr_db"],
        "anchor_4000hz_sinr_db": high_anchor["output_sinr_db"],
        "pcm_reference_mse_ch0": dict(zip(keys_d, pcm_mses)),
        "analytic_reference_mse_ch0": dict(zip(keys_d, ana_mses)),
        "plot_scope": "Corcos wind vs diffuse coherence, parametric cross-modal SINR, and actual PCM16 excerpt/MSE from checked wearable_audio assets",
    }
    _save_figure_and_report(
        fig,
        "fig83_own_voice_fusion.png",
        "figure83_own_voice_fusion.json",
        report,
        manifest_sha,
    )


def build_fig84(manifest: dict, signals: dict[str, np.ndarray], inputs: dict[str, str], manifest_sha: str) -> None:
    freqs = np.linspace(0.0, 4500.0, 901)
    ht_2ms = hear_through_response(freqs, tau_proc_s=0.002, g_leak=1.0, g_ht=1.0)
    ht_100us = hear_through_response(freqs, tau_proc_s=0.0001, g_leak=1.0, g_ht=1.0)

    delays_ms = np.geomspace(0.05, 8.0, 160)
    first_notches_hz = 1000.0 / (2.0 * delays_ms)

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.0), constrained_layout=True)

    gain_2ms_db = 10.0 * np.log10(np.maximum(ht_2ms["power_gain"], 1e-3))
    gain_100us_db = 10.0 * np.log10(np.maximum(ht_100us["power_gain"], 1e-3))
    axes[0, 0].plot(freqs, gain_2ms_db, color=C_RED, linewidth=1.8, label="处理延迟 τ = 2.0 ms")
    axes[0, 0].plot(freqs, gain_100us_db, color=C_BLUE, linestyle="--", linewidth=2.0, label="超低延迟 τ = 0.1 ms")
    axes[0, 0].axvline(250.0, color=C_MAIN, linestyle=":", label="首陷波 250 Hz")
    axes[0, 0].set(
        title="(a) 通透模式漏声叠加功率增益频响",
        xlabel="频率 / Hz",
        ylabel="叠加功率增益 / dB",
        xlim=(0.0, 4500.0),
        ylim=(-30.0, 8.5),
    )
    axes[0, 0].legend(loc="lower right")

    axes[0, 1].loglog(delays_ms, first_notches_hz, color=C_PURPLE, linewidth=2.2, label=r"首梳状陷波 $f_0 = 1/(2\tau)$")
    axes[0, 1].scatter([0.1, 2.0], [5000.0, 250.0], color=C_MAIN, zorder=5)
    axes[0, 1].axhline(4000.0, color=C_GREEN, linestyle="--", label="语音主频带上限 4 kHz")
    axes[0, 1].set(
        title="(b) 处理延迟 τ 与首个相消陷波位置",
        xlabel="端到端处理延迟 τ / ms（对数轴）",
        ylabel="首个相消陷波频率 f₀ / Hz",
        xlim=(0.05, 8.0),
        ylim=(50.0, 15000.0),
    )
    axes[0, 1].legend(loc="upper right")

    # Panel (c): Actual PCM16 waveform excerpt [1632, 1760) = 8 ms (2 full cycles of 250 Hz!)
    start, count = 1632, 128
    t_ms = np.arange(count) / 16000.0 * 1000.0
    axes[1, 0].plot(t_ms, signals["wearable_clean_ref"][0, start : start + count], color=C_BLUE, linewidth=1.9, label="干净直达参考")
    axes[1, 0].plot(
        t_ms, signals["wearable_hearthrough_comb"][0, start : start + count], color=C_RED, linestyle="--", linewidth=1.8, label="直达 + 2 ms 延迟叠加"
    )
    axes[1, 0].set(
        title="(c) 实际 PCM 波形：250 Hz 相消与 4 kHz 倍增",
        xlabel="自第 1632 点起的时间 / ms",
        ylabel="数字幅度",
        xlim=(0.0, t_ms[-1]),
        ylim=(-0.28, 0.28),
    )
    axes[1, 0].legend(loc="upper right")

    # Panel (d): Actual PCM16 tone projection magnitudes at 250 Hz and 4000 Hz
    ref_proj = np.asarray(manifest["samples"]["wearable_clean_ref"]["pcm_measurements"]["tone_projections_ch0_real_imag"])
    ht_proj = np.asarray(manifest["samples"]["wearable_hearthrough_comb"]["pcm_measurements"]["tone_projections_ch0_real_imag"])
    ref_mag = np.hypot(ref_proj[:, 0], ref_proj[:, 1])
    ht_mag = np.hypot(ht_proj[:, 0], ht_proj[:, 1])

    x2 = np.arange(2)
    w_bar = 0.28
    axes[1, 1].bar(x2 - w_bar / 2, ref_mag, width=w_bar, color=C_BLUE, alpha=0.85, label="直达参考 PCM 幅度")
    axes[1, 1].bar(x2 + w_bar / 2, ht_mag, width=w_bar, color=C_RED, alpha=0.85, label="2 ms 叠加 PCM 幅度")
    for idx_b, val in enumerate(ht_mag):
        axes[1, 1].text(idx_b + w_bar / 2, val + 0.008, f"{val:.4f}", ha="center")
    axes[1, 1].set(
        title="(d) 稳窗实际 PCM 双频投影峰值幅度",
        xticks=x2,
        xticklabels=["250 Hz（半波长相消）", "4000 Hz（8 周期同相）"],
        ylabel="正弦峰值幅度",
        ylim=(0.0, 0.25),
    )
    axes[1, 1].legend(loc="upper left")

    for ax in axes.flat:
        ax.grid(True, alpha=0.22)
        ax.set_axisbelow(True)
    fig.suptitle("TWS 通透模式（Hear-Through）处理延迟与耳道漏声梳状滤波效应")

    report = {
        "figure": 84,
        "exercise_ids": ["E16-14", "E16-15", "E16-16"],
        "source_sha256": inputs,
        "audio_manifest_sha256": manifest_sha,
        "wav_sha256": {
            "wearable_clean_ref.wav": manifest["files"]["wearable_clean_ref.wav"]["sha256"],
            "wearable_hearthrough_comb.wav": manifest["files"]["wearable_hearthrough_comb.wav"]["sha256"],
        },
        "comb_notches_2ms_up_to_4500hz": ht_2ms["comb_notch_frequencies_hz"].tolist(),
        "pcm_tone_magnitudes_clean_ref": ref_mag.tolist(),
        "pcm_tone_magnitudes_hearthrough_comb": ht_mag.tolist(),
        "plot_scope": "analytic hear-through comb transfer gain, first-notch latency law, and actual PCM16 waveform/tone projections",
    }
    _save_figure_and_report(
        fig,
        "fig84_hear_through_comb.png",
        "figure84_hear_through_comb.json",
        report,
        manifest_sha,
    )


def main() -> int:
    manifest = check_assets(DEFAULT_OUTPUT)
    inputs = _input_digests()
    manifest_path = validate_parent_chain(DEFAULT_OUTPUT / "MANIFEST.json")
    manifest_sha = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    signals = {
        name: read_pcm16(validate_parent_chain(DEFAULT_OUTPUT / rec["file"]).read_bytes())[1]
        for name, rec in manifest["samples"].items()
    }

    build_fig81(inputs, manifest_sha)
    build_fig82(inputs, manifest_sha)
    build_fig83(manifest, signals, inputs, manifest_sha)
    build_fig84(manifest, signals, inputs, manifest_sha)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
