"""Generate Chapter 0 figures (fig89..fig92) and numerical JSON reports.

Run from the repository root:
    .venv/bin/python scripts/make_math_figures.py
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from codes.chapters.ch00.core.chapter00_math import (
    dft_parseval_and_zero_padded_conv_example,
    dft_real_white_noise_pseudo_covariance,
    kalman_joseph_step,
    orthogonal_projection_and_four_subspaces_example,
    periodic_hann_window,
    phasor_delay_and_decibel_budget_example,
    proper_complex_gaussian_and_crlb_example,
    rayleigh_quotient_and_gev_diagonalization_example,
    solve_hermitian_gev,
    stft_cola_window_reconstruction_example,
    stft_wola_roundtrip,
    svd_condition_and_tikhonov_loading_example,
    tikhonov_svd_solve,
    vec_kronecker_and_khatri_rao_example,
    wiener_orthogonality_and_kalman_joseph_example,
    wirtinger_lagrange_mvdr_and_lcmv_example,
    woodbury_and_schur_complement_example,
)
from codes.chapters.ch00.io_contracts import (
    validate_parent_chain,
    validate_report_destination,
    write_json_report,
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

REPORT_DIR = ROOT / "codes/chapters/ch00/reports"
FIG_DIR = ROOT / "figures"

SCRIPT_INPUT_PATHS = (
    "codes/chapters/ch00/core/chapter00_math.py",
    "codes/chapters/ch00/io_contracts.py",
    "scripts/make_figures.py",
    "scripts/make_math_figures.py",
)

FIGURE_SEED_GAUSSIAN = 20260089


def _sha_file(rel_path: str) -> str:
    source = validate_parent_chain(ROOT / rel_path)
    if not source.is_file():
        raise ValueError(f"missing ordinary figure source: {rel_path}")
    return hashlib.sha256(source.read_bytes()).hexdigest()


def _input_digests() -> dict[str, str]:
    return {path: _sha_file(path) for path in SCRIPT_INPUT_PATHS}


def _save_report_and_png(
    fig,
    fig_filename: str,
    report_filename: str,
    report_payload: dict,
) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    report_path = validate_report_destination(REPORT_DIR / report_filename)
    dest = validate_report_destination(FIG_DIR / fig_filename)

    write_json_report(report_path, report_payload)
    os.chmod(report_path, 0o644)
    report_sha = hashlib.sha256(validate_parent_chain(report_path).read_bytes()).hexdigest()

    if fig.get_constrained_layout():
        fig.get_layout_engine().set(w_pad=0.06, h_pad=0.09, wspace=0.04, hspace=0.06)
    finalize_figure(fig)
    metadata = {
        "SourceScript": "scripts/make_math_figures.py",
        "SourceScriptDigest": report_payload["source_sha256"]["scripts/make_math_figures.py"],
        "GeneratorInputs": json.dumps(report_payload["source_sha256"], sort_keys=True),
        "NumericalReportDigest": report_sha,
    }
    validate_report_destination(dest)
    fig.savefig(dest, dpi=150, bbox_inches="tight", facecolor="white", metadata=metadata)
    os.chmod(dest, 0o644)
    plt.close(fig)
    print(dest.relative_to(ROOT))


def _build_fig89(sources: dict[str, str]) -> None:
    ex03 = phasor_delay_and_decibel_budget_example()
    ex04 = dft_parseval_and_zero_padded_conv_example()
    ex05 = stft_cola_window_reconstruction_example()

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.8), constrained_layout=True)

    # (a) Complex phasor sum on unit circle
    ax = axes[0, 0]
    theta_grid = np.linspace(0.0, 2.0 * np.pi, 240)
    ax.plot(np.cos(theta_grid), np.sin(theta_grid), color="#9aa5b1", linestyle=":", linewidth=1.3, label="单位圆 $|z|=1$")
    ax.axhline(0.0, color="#cbd2d9", linewidth=0.9)
    ax.axvline(0.0, color="#cbd2d9", linewidth=0.9)
    ax.annotate(
        "",
        xy=(1.0, 0.0),
        xytext=(0.0, 0.0),
        arrowprops=dict(arrowstyle="->", color=C_BLUE, lw=2.2),
    )
    ax.annotate(
        "",
        xy=(1.0, -1.0),
        xytext=(1.0, 0.0),
        arrowprops=dict(arrowstyle="->", color=C_ORANGE, lw=2.2, linestyle="--"),
    )
    ax.annotate(
        "",
        xy=(1.0, -1.0),
        xytext=(0.0, 0.0),
        arrowprops=dict(arrowstyle="->", color=C_RED, lw=2.4),
    )
    ax.plot([0.0, 1.0, 1.0], [0.0, 0.0, -1.0], "o", color=C_MAIN, markersize=5)
    ax.text(1.07, 0.12, r"$Z_1=1$", color=C_BLUE)
    ax.text(1.08, -0.55, r"$Z_2=-\mathrm{j}$", color=C_ORANGE)
    ax.annotate(
        r"$Z=1-\mathrm{j}$" + "\n" + r"$(|Z|=\sqrt{2},\ -45^\circ)$",
        xy=(0.56, -0.56),
        xytext=(-1.52, -1.44),
        bbox=dict(fc="white", ec="none", alpha=0.90, pad=1.0),
        arrowprops=dict(arrowstyle="->", color=C_RED, lw=1.3),
        color=C_RED,
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set(
        title="(a) 复平面相量合成（1 kHz，Δτ=250 μs）",
        xlabel=r"实部 $\mathrm{Re}\{Z\}$",
        ylabel=r"虚部 $\mathrm{Im}\{Z\}$",
        xlim=(-1.60, 1.92),
        ylim=(-1.62, 1.48),
    )
    ax.legend(loc="upper left")

    # (b) Linear phase delay slope across frequencies
    ax = axes[0, 1]
    f_khz = np.linspace(0.0, 4.0, 161)
    f_hz = f_khz * 1000.0
    for tau_us, color, ls, lbl in (
        (62.5, C_GREEN, "-.", r"$\tau_0 = 62.5\ \mu\mathrm{s}$（1 点）"),
        (125.0, C_BLUE, "-", r"$\tau_0 = 125\ \mu\mathrm{s}$（2 点）"),
        (250.0, C_RED, "--", r"$\tau_0 = 250\ \mu\mathrm{s}$（4 点）"),
    ):
        phase_rad = -2.0 * np.pi * f_hz * (tau_us * 1e-6)
        ax.plot(f_khz, phase_rad, color=color, linestyle=ls, linewidth=2.0, label=lbl)
    ax.scatter([0.5, 1.5, 2.5], ex03["unwrapped_phases_rad"], color=C_BLUE, zorder=4, label="E00-03 三频点采样")
    ax.set(
        title="(b) 纯时延传递函数的线性相位斜率",
        xlabel=r"频率 $f$ / kHz",
        ylabel=r"解缠绕相位 $\phi(f)$ / rad",
        xlim=(0.0, 4.0),
        ylim=(-7.35, 0.45),
    )
    ax.legend(loc="lower left")

    # (c) Circular aliasing vs zero-padded linear convolution
    ax = axes[1, 0]
    y_lin = ex04["linear_convolution"]
    y_zp = ex04["zero_padded_fft_convolution"]
    y_circ_pad = np.pad(ex04["circular_convolution_no_pad"], (0, 2), constant_values=np.nan)
    n_idx = np.arange(y_lin.size)
    ax.plot(n_idx, y_lin, "o-", color=C_BLUE, linewidth=2.1, markersize=6, label="时域线性卷积真值（6 点）")
    ax.plot(n_idx, y_zp, "s--", color=C_GREEN, linewidth=1.8, markersize=5, label=r"补零 FFT 乘积（$N_{\mathrm{pad}}=8$）")
    ax.plot(np.arange(4), y_circ_pad[:4], "^:", color=C_RED, linewidth=2.0, markersize=7, label=r"未补零循环卷积（$N=4$）")
    ax.annotate(
        "尾部 y[4]=0.5\n折回混叠至 n=0",
        xy=(0.05, 1.50),
        xytext=(1.55, 1.85),
        arrowprops=dict(arrowstyle="->", color=C_RED, lw=1.4),
        color=C_RED,
    )
    ax.set(
        title="(c) 循环卷积尾部缠绕与补零线性卷积",
        xlabel=r"离散采样索引 $n$",
        ylabel=r"序列幅值 $y[n]$",
        xlim=(-0.35, 5.35),
        ylim=(-3.90, 2.85),
    )
    ax.legend(loc="lower left")

    # (d) COLA / WOLA steady-state denominator comparison
    ax = axes[1, 1]
    n_win = 32
    sig = np.ones(128, dtype=np.float64)
    w_hann = periodic_hann_window(n_win)
    w_sqrt = np.sqrt(w_hann)
    w_rect = np.ones(n_win, dtype=np.float64)
    ola_50 = stft_wola_roundtrip(sig, w_hann, w_rect, n_win // 2)["cola_denominator"][32:96]
    wola_50 = stft_wola_roundtrip(sig, w_sqrt, w_sqrt, n_win // 2)["cola_denominator"][32:96]
    sq_50 = stft_wola_roundtrip(sig, w_hann, w_hann, n_win // 2)["cola_denominator"][32:96]
    sq_75 = stft_wola_roundtrip(sig, w_hann, w_hann, n_win // 4)["cola_denominator"][32:96]
    sample_axis = np.arange(64)
    ax.plot(
        sample_axis, ola_50, color=C_BLUE, linewidth=2.2,
        marker="o", markevery=(0, 8), markerfacecolor="none", markersize=6.0,
        label="50% 重叠周期汉宁 OLA（恒为 1.0）",
    )
    ax.plot(
        sample_axis, wola_50, color=C_GREEN, linestyle="--", linewidth=2.0,
        marker="x", markevery=(4, 8), markersize=6.0,
        label="50% 重叠开方汉宁 WOLA（恒为 1.0）",
    )
    ax.plot(sample_axis, sq_50, color=C_RED, linestyle="-.", linewidth=1.9, label="50% 重叠平方汉宁（0.50 至 1.00 波纹）")
    ax.plot(sample_axis, sq_75, color=C_PURPLE, linestyle=":", linewidth=2.0, label="75% 重叠平方汉宁（恒为 1.5）")
    ax.set(
        title="(d) STFT 重叠相加合成窗分母包络",
        xlabel=r"稳定窗内局部采样索引 $n$",
        ylabel=r"叠加分母增益 $g[n]$",
        xlim=(0, 63),
        ylim=(0.35, 2.45),
    )
    ax.legend(loc="upper right")

    for ax_i in axes.flat:
        ax_i.grid(True, alpha=0.22)
        ax_i.set_axisbelow(True)
    fig.suptitle("图89  复数相量旋转、时延线性相位斜率、补零线性卷积与 COLA 重叠相加重建")

    report = {
        "figure": 89,
        "exercise_ids": ["E00-03", "E00-04", "E00-05"],
        "source_sha256": sources,
        "anchors": {
            "quadrature_phasor_magnitude": ex03["quadrature_magnitude"],
            "quadrature_phasor_phase_deg": ex03["quadrature_phase_deg"],
            "recovered_group_delay_us": ex03["recovered_group_delay_us"],
            "parseval_energy": ex04["time_energy"],
            "circular_aliasing_l2_error": ex04["circular_aliasing_l2_error"],
            "zero_padded_max_abs_error": ex04["zero_padded_max_abs_error"],
            "hann_50pct_max_abs_error": ex05["hann_ola_50pct"]["steady_max_abs_error"],
            "squared_hann_50pct_ripple": ex05["unnormalized_squared_hann_50pct_counterexample"]["steady_cola_ripple"],
        },
        "scope": "Chapter 0 complex phasor, linear-phase delay, zero-padded FFT convolution, and COLA/WOLA curves",
    }
    _save_report_and_png(
        fig,
        "fig89_complex_phasor_fourier.png",
        "figure89_complex_phasor_fourier.json",
        report,
    )


def _build_fig90(sources: dict[str, str]) -> None:
    ex06 = orthogonal_projection_and_four_subspaces_example()
    ex07 = svd_condition_and_tikhonov_loading_example()
    ex08 = woodbury_and_schur_complement_example()

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.8), constrained_layout=True)

    # (a) Column-space orthogonal decomposition & Pythagoras check
    ax = axes[0, 0]
    u_dir = np.array([2.0, 1.0]) / math.sqrt(5.0)
    b_demo = np.array([1.2, 1.6])
    b_par_demo = float(np.dot(b_demo, u_dir)) * u_dir
    t_line = np.linspace(-0.4, 2.3, 50)
    ax.plot(t_line * u_dir[0], t_line * u_dir[1], color=C_BLUE, linestyle="--", linewidth=1.6, label=r"列空间 $\mathcal{R}(\mathbf{A})$")
    ax.annotate("", xy=b_demo, xytext=(0.0, 0.0), arrowprops=dict(arrowstyle="->", color=C_MAIN, lw=2.3))
    ax.annotate("", xy=b_par_demo, xytext=(0.0, 0.0), arrowprops=dict(arrowstyle="->", color=C_GREEN, lw=2.3))
    ax.annotate("", xy=b_demo, xytext=b_par_demo, arrowprops=dict(arrowstyle="->", color=C_RED, lw=2.2))
    # Right-angle square marker at foot of perpendicular
    u_perp = np.array([-u_dir[1], u_dir[0]])
    sq_sz = 0.11
    ra_pts = np.array(
        [
            b_par_demo - sq_sz * u_dir,
            b_par_demo - sq_sz * u_dir + sq_sz * u_perp,
            b_par_demo + sq_sz * u_perp,
        ]
    )
    ax.plot(ra_pts[:, 0], ra_pts[:, 1], color=C_MAIN, linewidth=1.2)
    ax.text(0.36, 1.15, r"$\vec{b}$", color=C_MAIN)
    ax.text(0.84, 0.20, r"$\vec{b}_\parallel = \mathbf{P}_{\mathcal{R}(\mathbf{A})}\vec{b}$", color=C_GREEN)
    ax.text(1.50, 1.20, r"$\vec{b}_\perp \in \mathcal{N}(\mathbf{A}^H)$", color=C_RED)
    ax.text(
        -0.30,
        1.76,
        f"E00-06 复投影勾股核验：\n"
        f"$\\|\\vec{{b}}_\\parallel\\|_2^2={ex06['norm_sq_parallel']:.3f},\\ "
        f"\\|\\vec{{b}}_\\perp\\|_2^2={ex06['norm_sq_orthogonal']:.3f}$",
        color=C_MAIN,
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set(
        title="(a) 列空间正交投影与左零空间正交分解",
        xlabel="第一正交基坐标",
        ylabel="第二正交基坐标",
        xlim=(-0.40, 2.70),
        ylim=(-0.28, 2.15),
    )
    ax.legend(loc="lower right")

    # (b) Tikhonov spectral shrinkage factor sigma^2 / (sigma^2 + eps)
    ax = axes[0, 1]
    sig_grid = np.logspace(np.log10(5e-4), np.log10(4e2), 260)
    for eps_val, color, ls in (
        (1e-4, C_BLUE, "-"),
        (1e-2, C_GREEN, "--"),
        (1e-1, C_ORANGE, "-."),
        (1.0, C_RED, ":"),
    ):
        shrink = (sig_grid ** 2) / (sig_grid ** 2 + eps_val)
        ax.semilogx(sig_grid, shrink, color=color, linestyle=ls, linewidth=2.0, label=rf"$\varepsilon = {eps_val:g}$")
    ax.axvline(0.01, color=C_PURPLE, linestyle=":", linewidth=1.4)
    ax.axvline(10.0, color=C_MAIN, linestyle=":", linewidth=1.4)
    ax.text(7.8e-3, 0.85, r"$\sigma_2=0.01$", color=C_PURPLE, ha="right",
            bbox=dict(fc="white", ec="none", alpha=0.88, pad=0.8))
    ax.text(8.0, 0.16, r"$\sigma_1=10$", color=C_MAIN, ha="right",
            bbox=dict(fc="white", ec="none", alpha=0.88, pad=0.8))
    ax.set(
        title=r"(b) Tikhonov 谱收缩因子 $\sigma^2/(\sigma^2+\varepsilon)$",
        xlabel=r"奇异值 $\sigma$（对数坐标）",
        ylabel="模态保留权重",
        xlim=(5e-4, 4e2),
        ylim=(-0.03, 1.06),
    )
    ax.legend(loc="center right")

    # (c) Bias-variance U-curve under ill-conditioned inverse (kappa_2 = 1000)
    ax = axes[1, 0]
    g_mat = ex07["g_matrix"]
    x_true = np.array([1.0, 0.0], dtype=np.float64)
    y_clean = g_mat @ x_true
    u_mat, _, _ = np.linalg.svd(g_mat)
    delta_y = 0.10 * u_mat[:, 1]
    eps_grid = np.logspace(-6, 2, 161)
    bias_norms, noise_norms, tot_norms = [], [], []
    for eps_v in eps_grid:
        x_clean_est = tikhonov_svd_solve(g_mat, y_clean, float(eps_v))["solution"].real
        x_noisy_est = tikhonov_svd_solve(g_mat, y_clean + delta_y, float(eps_v))["solution"].real
        bias_norms.append(float(np.linalg.norm(x_clean_est - x_true)))
        noise_norms.append(float(np.linalg.norm(x_noisy_est - x_clean_est)))
        tot_norms.append(float(np.linalg.norm(x_noisy_est - x_true)))

    ax.loglog(eps_grid, tot_norms, color=C_MAIN, linewidth=2.3, label="总估计误差范数")
    ax.loglog(eps_grid, noise_norms, color=C_RED, linestyle="--", linewidth=1.9, label="弱模态噪声放大项")
    ax.loglog(eps_grid, bias_norms, color=C_BLUE, linestyle="-.", linewidth=1.9, label="正则化收缩偏差项")
    ax.scatter([0.1], [ex07["tikhonov_error_norm"]], color=C_GREEN, s=48, zorder=5, label=r"$\varepsilon=0.1$ 选点")
    ax.set(
        title=r"(c) 病态逆（$\kappa_2=1000$）偏差—方差折中曲线",
        xlabel=r"对角加载参数 $\varepsilon$（对数坐标）",
        ylabel="二范数误差",
        xlim=(1e-6, 1e2),
        ylim=(5e-4, 2e1),
    )
    ax.legend(loc="lower left")

    # (d) Woodbury rank-1 interference nulling vs INR
    ax = axes[1, 1]
    inr_db = np.linspace(-10.0, 30.0, 81)
    inr_lin = 10.0 ** (inr_db / 10.0)
    analytic_atten_db = -20.0 * np.log10(1.0 + inr_lin)
    sample_idx = np.linspace(0, 80, 11, dtype=int)
    u_int = np.array([1.0, 1.0j, -1.0], dtype=np.complex128)
    direct_atten_db = []
    for idx in sample_idx:
        sigma_i2 = inr_lin[idx] / 3.0
        r_m = np.eye(3, dtype=np.complex128) + sigma_i2 * np.outer(u_int, u_int.conj())
        ratio = float(np.linalg.norm(np.linalg.solve(r_m, u_int)) / np.linalg.norm(u_int))
        direct_atten_db.append(20.0 * math.log10(ratio))

    ax.plot(inr_db, analytic_atten_db, color=C_BLUE, linewidth=2.1, label=r"Woodbury 解析式 $-20\log_{10}(1+\mathrm{INR})$")
    ax.plot(inr_db[sample_idx], direct_atten_db, "o", color=C_RED, markersize=5, label="直接矩阵求逆采样核验")
    anchor_inr_db = 10.0 * math.log10(6.0)
    anchor_att_db = 20.0 * math.log10(ex08["interferer_attenuation_ratio"])
    ax.scatter([anchor_inr_db], [anchor_att_db], color=C_GREEN, marker="D", s=55, zorder=5)
    ax.annotate(
        "E00-08 锚点\n(INR=7.78 dB, -16.90 dB)",
        xy=(anchor_inr_db, anchor_att_db),
        xytext=(10.5, -8.5),
        arrowprops=dict(arrowstyle="->", color=C_GREEN, lw=1.4),
        color=C_GREEN,
    )
    ax.set(
        title="(d) Woodbury 秩一更新对干扰方向的抑制深度",
        xlabel="干扰噪声功率比 INR / dB",
        ylabel=r"$\|\mathbf{R}^{-1}\vec{u}\|_2 / \|\vec{u}\|_2$ / dB",
        xlim=(-10.0, 30.0),
        ylim=(-62.0, 6.0),
    )
    ax.legend(loc="lower left")

    for ax_i in axes.flat:
        ax_i.grid(True, alpha=0.22)
        ax_i.set_axisbelow(True)
    fig.suptitle("图90  列空间正交投影分解、奇异值谱、Tikhonov 对角加载与 Woodbury 秩一更新")

    report = {
        "figure": 90,
        "exercise_ids": ["E00-06", "E00-07", "E00-08"],
        "source_sha256": sources,
        "anchors": {
            "norm_sq_b": ex06["norm_sq_b"],
            "norm_sq_parallel": ex06["norm_sq_parallel"],
            "norm_sq_orthogonal": ex06["norm_sq_orthogonal"],
            "condition_number_2": ex07["condition_number_2"],
            "unregularized_error_norm": ex07["unregularized_error_norm"],
            "tikhonov_eps0p1_error_norm": ex07["tikhonov_error_norm"],
            "woodbury_interferer_attenuation_ratio": ex08["interferer_attenuation_ratio"],
            "schur_determinant": ex08["det_via_schur"],
        },
        "scope": "Chapter 0 column-space projection, Tikhonov SVD spectral filtering, and Woodbury rank-1 update",
    }
    _save_report_and_png(
        fig,
        "fig90_subspace_projection_svd.png",
        "figure90_subspace_projection_svd.json",
        report,
    )


def _build_fig91(sources: dict[str, str]) -> None:
    ex09 = rayleigh_quotient_and_gev_diagonalization_example()
    ex10 = vec_kronecker_and_khatri_rao_example()

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.8), constrained_layout=True)

    # (a) Standard Hermitian Rayleigh quotient over 2D unit vectors
    ax = axes[0, 0]
    psi_deg = np.linspace(-90.0, 90.0, 241)
    psi_rad = np.deg2rad(psi_deg)
    # R_nn = [[2, 1], [1, 2]] -> w(psi) = [cos psi, sin psi]^T -> w^T R_nn w = 2 + sin(2 psi)
    rq_std = 2.0 + np.sin(2.0 * psi_rad)
    ax.plot(psi_deg, rq_std, color=C_BLUE, linewidth=2.2, label=r"瑞利商 $\vec{w}^H\mathbf{R}_{nn}\vec{w}/\vec{w}^H\vec{w}$")
    ax.axhline(3.0, color=C_RED, linestyle="--", linewidth=1.6, label=r"$\lambda_{\max}=3.0$（$\psi=+45^\circ$）")
    ax.axhline(1.0, color=C_GREEN, linestyle="-.", linewidth=1.6, label=r"$\lambda_{\min}=1.0$（$\psi=-45^\circ$）")
    ax.scatter([-45.0, 45.0], [1.0, 3.0], color=[C_GREEN, C_RED], zorder=5, s=45)
    ax.set(
        title="(a) 厄米矩阵标准瑞利商的极值夹逼轨迹",
        xlabel="实单位权重角 ψ / °",
        ylabel="瑞利商数值",
        xlim=(-90.0, 90.0),
        ylim=(0.65, 4.55),
    )
    ax.legend(loc="upper left")

    # (b) Generalized Rayleigh quotient (GEV SINR) vs complex weight phase/amplitude ratio
    ax = axes[0, 1]
    r_nn = ex09["r_nn"]
    r_ss = ex09["r_ss"]
    # Parameterize w(phi) = [1, exp(j phi)]^T vs optimal GEV w_gev = R_nn^{-1} a_tar
    phi_deg = np.linspace(-180.0, 180.0, 361)
    gev_eq_amp = []
    for deg in phi_deg:
        w_ph = np.array([1.0, np.exp(1j * math.radians(deg))], dtype=np.complex128)
        num = float(np.vdot(w_ph, r_ss @ w_ph).real)
        den = float(np.vdot(w_ph, r_nn @ w_ph).real)
        gev_eq_amp.append(num / den)
    ax.plot(phi_deg, gev_eq_amp, color=C_PURPLE, linewidth=2.1, label=r"等幅相位扫描 $[1,\ e^{\mathrm{j}\phi}]^\top$")
    ax.axhline(ex09["max_gev_sinr"], color=C_RED, linestyle="--", linewidth=2.0, label=r"GEV 全局最优 $\lambda_{\max}=4/3$")
    ax.scatter([90.0], [ex09["unwhitened_matched_sinr"]], color=C_BLUE, s=52, zorder=5, label=r"未白化匹配滤波（$\phi=90^\circ,\ 1.00$）")
    ax.set(
        title="(b) 广义瑞利商极值与未白化匹配滤波差距",
        xlabel="第二阵元相对相位 φ / °",
        ylabel="广义瑞利商 SINR",
        xlim=(-180.0, 180.0),
        ylim=(-0.05, 2.10),
    )
    ax.legend(loc="upper left")

    # (c) Cholesky whitening geometric ellipse -> unit circle
    ax = axes[1, 0]
    ang = np.linspace(0.0, 2.0 * np.pi, 240)
    unit_circle = np.vstack([np.cos(ang), np.sin(ang)])
    # Ellipse z^T R_nn z = 1 <=> R_nn^{1/2} z on unit circle => z = R_nn^{-1/2} u
    eig_v, eig_q = np.linalg.eigh(r_nn.real)
    r_inv_sqrt = eig_q @ np.diag(1.0 / np.sqrt(eig_v)) @ eig_q.T
    ellipse_pts = r_inv_sqrt @ unit_circle
    ax.plot(ellipse_pts[0], ellipse_pts[1], color=C_ORANGE, linewidth=2.1, label=r"原空间噪声等值椭圆 $\vec{z}^\top\mathbf{R}_{nn}\vec{z}=1$")
    ax.plot(
        unit_circle[0], unit_circle[1], color=C_BLUE, linestyle="--", linewidth=2.0,
        label=r"白化后各向同性圆 $\vec{u}^\top\vec{u}=1$（$\vec{u}=\mathbf{L}_B^\top\vec{z}$）",
    )
    gev_real = solve_hermitian_gev(np.array([[2.0, -1.5], [-1.5, 2.0]]), r_nn.real)
    v_opt = gev_real["joint_diagonalizer_phi"][:, -1].real
    v_opt = v_opt / np.linalg.norm(v_opt)
    ax.annotate("", xy=0.95 * v_opt, xytext=(0.0, 0.0), arrowprops=dict(arrowstyle="->", color=C_GREEN, lw=2.2))
    ax.annotate(
        "实对称截面：\n白化消除方向偏心",
        xy=(0.55 * v_opt[0], 0.55 * v_opt[1]),
        xytext=(-2.22, 0.82),
        bbox=dict(fc="white", ec="none", alpha=0.90, pad=1.0),
        arrowprops=dict(arrowstyle="->", color=C_GREEN, lw=1.3),
        color=C_GREEN,
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.set(
        title="(c) 白化合同变换：相关噪声椭圆化为单位圆",
        xlabel="第一分量坐标",
        ylabel="第二分量坐标",
        xlim=(-2.45, 2.35),
        ylim=(-2.38, 1.52),
    )
    ax.legend(loc="lower center")

    # (d) Khatri-Rao difference coarray virtual lag multiplicity: sparse {0,1,3} vs ULA {0,1,2}
    ax = axes[1, 1]
    lags_grid = np.arange(-3, 4)
    sparse_pos = np.array([0, 1, 3])
    ula_pos = np.array([0, 1, 2])
    sparse_counts = [sum(1 for m in sparse_pos for n in sparse_pos if m - n == lg) for lg in lags_grid]
    ula_counts = [sum(1 for m in ula_pos for n in ula_pos if m - n == lg) for lg in lags_grid]
    width = 0.36
    ax.bar(lags_grid - 0.5 * width, sparse_counts, width=width, color=C_BLUE, label="稀疏阵 {0,1,3}d（7 个连续滞后）")
    ax.bar(lags_grid + 0.5 * width, ula_counts, width=width, color=C_ORANGE, alpha=0.88, label="等距阵 {0,1,2}d（仅 5 个滞后）")
    ax.set(
        title=r"(d) $\mathbf{A}^*\odot\mathbf{A}$ 差分协同阵虚拟孔径对比",
        xlabel=r"归一化整数空间滞后 $(x_m - x_n) / d$",
        ylabel="协方差元素重数",
        xticks=lags_grid,
        ylim=(0, 4.35),
    )
    ax.legend(loc="upper right")

    for ax_i in axes.flat:
        ax_i.grid(True, alpha=0.22)
        ax_i.set_axisbelow(True)
    fig.suptitle("图91  厄米瑞利商极值轨迹、广义特征值联合对角化与 Khatri–Rao 虚拟孔径展开")

    report = {
        "figure": 91,
        "exercise_ids": ["E00-09", "E00-10"],
        "source_sha256": sources,
        "anchors": {
            "r_nn_eigenvalues": ex09["r_nn_eigenvalues"].tolist(),
            "max_gev_sinr": ex09["max_gev_sinr"],
            "unwhitened_matched_sinr": ex09["unwhitened_matched_sinr"],
            "whitening_sinr_gain_db": ex09["whitening_sinr_gain_db"],
            "sparse_unique_virtual_lags": ex10["num_unique_virtual_lags"],
            "khatri_rao_vec_max_abs_error": ex10["khatri_rao_vec_max_abs_error"],
        },
        "scope": "Chapter 0 Hermitian Rayleigh bounds, GEV Cholesky whitening, and Khatri-Rao difference coarray",
    }
    _save_report_and_png(
        fig,
        "fig91_hermitian_rayleigh_gev.png",
        "figure91_hermitian_rayleigh_gev.json",
        report,
    )


def _build_fig92(sources: dict[str, str]) -> None:
    ex11 = wirtinger_lagrange_mvdr_and_lcmv_example()
    ex13 = proper_complex_gaussian_and_crlb_example()
    ex14 = wiener_orthogonality_and_kalman_joseph_example()

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.8), constrained_layout=True)

    # (a) Wirtinger projected conjugate-gradient descent to MVDR optimum (Frost / PGD)
    ax = axes[0, 0]
    a_tar = np.array([1.0, 1.0, 1.0], dtype=np.complex128)
    a_int = np.array([1.0, 1.0j, -1.0], dtype=np.complex128)
    r_cov = np.eye(3, dtype=np.complex128) + 2.0 * np.outer(a_int, a_int.conj())
    p_perp = np.eye(3, dtype=np.complex128) - np.outer(a_tar, a_tar.conj()) / 3.0
    w_dsb = a_tar / 3.0
    for mu_step, color, ls, mk, lbl in (
        (0.08, C_BLUE, "-", "o", r"步长 $\mu = 0.08$"),
        (0.16, C_GREEN, "--", "s", r"步长 $\mu = 0.16$"),
    ):
        w_cur = w_dsb.copy()
        powers = [float(np.vdot(w_cur, r_cov @ w_cur).real)]
        for _ in range(10):
            # Frost projected Wirtinger gradient step: w <- P_perp(w - mu R w) + w_dsb
            w_cur = p_perp @ (w_cur - mu_step * (r_cov @ w_cur)) + w_dsb
            powers.append(float(np.vdot(w_cur, r_cov @ w_cur).real))
        ax.plot(np.arange(11), powers, color=color, linestyle=ls, marker=mk, markersize=4.5, linewidth=2.0, label=lbl)

    ax.axhline(ex11["mvdr_output_power"], color=C_RED, linestyle=":", linewidth=2.0, label="MVDR 闭式下界 7/19 ≈ 0.368")
    ax.axhline(ex11["lcmv_output_power"], color=C_PURPLE, linestyle="-.", linewidth=1.6, label="LCMV 硬零陷功率 3/8 = 0.375")
    ax.set(
        title=r"(a) Wirtinger 共轭梯度投影收敛轨迹",
        xlabel=r"迭代步数 $t$",
        ylabel=r"输出功率 $\vec{w}^H\mathbf{R}\vec{w}$",
        xticks=range(0, 11, 2),
        xlim=(0, 10),
        ylim=(0.33, 0.60),
    )
    ax.legend(loc="upper right")

    # (b) Proper circular vs Improper elliptical complex Gaussian scatter & 2-sigma contour
    ax = axes[0, 1]
    rng = np.random.default_rng(FIGURE_SEED_GAUSSIAN)
    n_pts = 260
    z_prop = rng.normal(0.0, math.sqrt(0.5), size=(n_pts, 2))
    cov_imp = np.array([[0.8, 0.2], [0.2, 0.2]], dtype=np.float64)
    chol_imp = np.linalg.cholesky(cov_imp)
    z_imp = (chol_imp @ rng.normal(0.0, 1.0, size=(2, n_pts))).T

    ax.scatter(z_prop[:, 0], z_prop[:, 1], s=12, alpha=0.38, color=C_BLUE, label=r"圆对称 Proper（$|C_{zz}|=0$）")
    ax.scatter(z_imp[:, 0], z_imp[:, 1], s=12, alpha=0.42, color=C_RED, marker="^", label=r"非圆对称（$|C_{zz}|=0.721$）")
    ang = np.linspace(0.0, 2.0 * np.pi, 200)
    unit_c = np.vstack([np.cos(ang), np.sin(ang)])
    c_prop = 2.0 * math.sqrt(0.5) * unit_c
    c_imp = 2.0 * (chol_imp @ unit_c)
    ax.plot(c_prop[0], c_prop[1], color=C_BLUE, linewidth=2.0)
    ax.plot(c_imp[0], c_imp[1], color=C_RED, linestyle="--", linewidth=2.0)
    ax.set_aspect("equal", adjustable="box")
    ax.set(
        title="(b) 圆对称与非圆对称复高斯分布对比",
        xlabel=r"实部 $\mathrm{Re}\{z\}$",
        ylabel=r"虚部 $\mathrm{Im}\{z\}$",
        xlim=(-2.75, 2.75),
        ylim=(-2.15, 3.75),
    )
    ax.legend(loc="upper left")

    # (c) 8-point real white-noise DFT variance vs pseudo-variance (DC & Nyquist exception)
    ax = axes[1, 0]
    dft_prop = dft_real_white_noise_pseudo_covariance(n_fft=8, sigma2=1.0)
    bins = np.arange(8)
    w_bar = 0.36
    ax.bar(bins - 0.5 * w_bar, dft_prop["variance_per_bin"], width=w_bar, color=C_BLUE, label=r"方差 $\mathbb{E}[|X[k]|^2]$")
    ax.bar(
        bins + 0.5 * w_bar,
        np.abs(dft_prop["pseudo_variance_per_bin"]),
        width=w_bar,
        color=C_ORANGE,
        label=r"伪协方差模 $|\mathbb{E}[X[k]^2]|$",
    )
    ax.annotate("DC 实频点", xy=(0.18, 1.01), xytext=(0.45, 1.25), arrowprops=dict(arrowstyle="->", color=C_ORANGE, lw=1.3), color=C_ORANGE)
    ax.annotate("Nyquist 实频点", xy=(4.18, 1.01), xytext=(4.45, 1.16), arrowprops=dict(arrowstyle="->", color=C_ORANGE, lw=1.3), color=C_ORANGE)
    ax.set(
        title="(c) 实白噪声 8 点酉 DFT 的实频点非圆对称性",
        xlabel="DFT 频点索引 k",
        ylabel="二阶统计矩幅值",
        xticks=bins,
        ylim=(0.0, 1.95),
    )
    ax.legend(loc="upper right")

    # (d) Multi-step Kalman predict vs Joseph posterior variance contraction
    ax = axes[1, 1]
    f_mat = np.array([[1.0, 1.0], [0.0, 1.0]], dtype=np.float64)
    h_mat = np.array([[1.0, 0.0]], dtype=np.float64)
    q_cov = np.diag([0.0, 1.0])
    r_cov_kf = np.array([[3.0]], dtype=np.float64)
    s_cur = np.array([10.0, 2.0], dtype=np.float64)
    p_cur = np.diag([4.0, 1.0])
    obs_seq = [13.5, 14.2, 16.4, 18.1, 20.3, 22.0, 24.1, 26.2]
    pred_vars, post_vars = [], []
    for z_val in obs_seq:
        step_res = kalman_joseph_step(f_mat, h_mat, q_cov, r_cov_kf, s_cur, p_cur, np.array([z_val]))
        pred_vars.append(float(step_res["predicted_cov"][0, 0]))
        post_vars.append(float(step_res["posterior_cov_joseph"][0, 0]))
        s_cur = step_res["posterior_state"]
        p_cur = step_res["posterior_cov_joseph"]

    steps_x = np.arange(1, len(obs_seq) + 1)
    ax.plot(steps_x, pred_vars, "s--", color=C_ORANGE, linewidth=2.0, markersize=5, label=r"先验预测方差 $P_{t|t-1}[0,0]$")
    ax.plot(steps_x, post_vars, "o-", color=C_GREEN, linewidth=2.2, markersize=6, label=r"Joseph 后验方差 $P_{t|t}[0,0]$")
    ax.axhline(3.0, color=C_BLUE, linestyle=":", linewidth=1.7, label="单帧观测噪声方差 R = 3.0")
    ax.annotate(
        "E00-14 第 1 步\n5.000 → 1.875",
        xy=(1.05, 1.90),
        xytext=(3.15, 3.45),
        arrowprops=dict(arrowstyle="->", color=C_GREEN, lw=1.4),
        color=C_GREEN,
    )
    ax.set(
        title="(d) 卡尔曼预测—新息—Joseph 协方差递推",
        xlabel="递推帧步 t",
        ylabel="位置状态误差方差",
        xticks=steps_x,
        ylim=(0.45, 9.20),
    )
    ax.legend(loc="upper left")

    for ax_i in axes.flat:
        ax_i.grid(True, alpha=0.22)
        ax_i.set_axisbelow(True)
    fig.suptitle("图92  Wirtinger 共轭梯度下降、圆对称复高斯伪协方差判据、维纳正交与卡尔曼更新")

    report = {
        "figure": 92,
        "exercise_ids": ["E00-11", "E00-13", "E00-14"],
        "source_sha256": sources,
        "seed": FIGURE_SEED_GAUSSIAN,
        "anchors": {
            "mvdr_output_power": ex11["mvdr_output_power"],
            "lcmv_output_power": ex11["lcmv_output_power"],
            "improper_pseudo_variance_abs": ex13["improper_elliptical_case"]["pseudo_variance_abs"],
            "dft_dc_pseudo_variance": ex13["real_dft_8pt_propriety"]["dc_pseudo_variance"],
            "dft_nyquist_pseudo_variance": ex13["real_dft_8pt_propriety"]["nyquist_pseudo_variance"],
            "wiener_mmse": ex14["wiener"]["mmse"],
            "kalman_step1_pred_var0": float(ex14["kalman"]["predicted_cov"][0, 0]),
            "kalman_step1_post_var0": float(ex14["kalman"]["posterior_cov_joseph"][0, 0]),
        },
        "scope": "Chapter 0 Wirtinger projected gradient, proper/improper complex Gaussian & DFT propriety, and Kalman Joseph covariance",
    }
    _save_report_and_png(
        fig,
        "fig92_wirtinger_probability_bayesian.png",
        "figure92_wirtinger_probability_bayesian.json",
        report,
    )


def main() -> int:
    sources = _input_digests()
    _build_fig89(sources)
    _build_fig90(sources)
    _build_fig91(sources)
    _build_fig92(sources)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
