"""AES67/PTPv2 clock phase budget and KU Leuven DASF teaching kernels for Ch 13.

Conventions match ``codes.chapters.ch13.core.distributed``:
- Inner-product beamformer output ``y = w^H x``;
- Linear constraints ``C^H w = f`` (with ``C`` of shape ``(M, P)``);
- Node ``q`` observes local ``M_q`` channels plus ``R`` compressed channels
  ``z_v = F_v^H x_v`` from each peer ``v != q`` (total compressed dimension
  ``L_q = M_q + (U - 1) * R``).

Pure NumPy implementation (zero SciPy imports, no file I/O).
"""
from __future__ import annotations

import math
import numpy as np


def _validate_real_scalar(value: object, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise ValueError(f"{name} must be a real numeric scalar, not bool")
    scalar = float(value)
    if not math.isfinite(scalar):
        raise ValueError(f"{name} must be finite")
    return scalar


def _validate_positive_int(value: object, name: str, *, minimum: int = 1) -> int:
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        or int(value) < minimum
    ):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _validate_complex_array(
    value: object, name: str, *, ndim: int | None = None
) -> np.ndarray:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be numeric, not bool")
    raw = np.asarray(value)
    if raw.dtype.kind not in "iufc":
        raise ValueError(f"{name} must be numeric, not bool or text")
    if isinstance(value, (list, tuple)):
        def _check_no_bool(items: object) -> None:
            if isinstance(items, (list, tuple)):
                for item in items:
                    _check_no_bool(item)
            elif isinstance(items, (bool, np.bool_)):
                raise ValueError(f"{name} must not contain bool")

        _check_no_bool(value)
    with np.errstate(over="ignore", invalid="ignore"):
        out = np.asarray(raw, dtype=np.complex128)
    if ndim is not None and out.ndim != ndim:
        raise ValueError(f"{name} must have ndim == {ndim}")
    if out.size == 0 or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be non-empty and finite")
    return out


def _validate_spd_covariance(
    value: object,
    name: str = "r_nn",
    *,
    tol: float = 1e-12,
    condition_limit: float = 1e12,
) -> np.ndarray:
    mat = _validate_complex_array(value, name, ndim=2)
    if mat.shape[0] != mat.shape[1]:
        raise ValueError(f"{name} must be square")
    scale = float(np.max(np.abs(mat)))
    if scale <= 0.0:
        raise ValueError(f"{name} must be strictly positive definite")
    norm = mat.real / scale + 1j * (mat.imag / scale)
    if float(np.max(np.abs(norm - norm.conj().T))) > tol:
        raise ValueError(f"{name} must be Hermitian")
    sym = 0.5 * (mat + mat.conj().T)
    eigs = np.linalg.eigvalsh(sym)
    if float(eigs[0]) <= 0.0:
        raise ValueError(f"{name} must be strictly positive definite")
    cond = float(eigs[-1] / eigs[0])
    if not math.isfinite(cond) or cond > condition_limit:
        raise ValueError(f"{name} condition number exceeds limit")
    return sym


def _orthonormal_columns(mat: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    q_mat, r_mat = np.linalg.qr(mat, mode="reduced")
    diag = np.abs(np.diag(r_mat))
    if diag.size == 0 or float(np.max(diag)) <= 0.0 or float(np.min(diag)) <= tol * float(np.max(diag)):
        raise ValueError("compressor columns are rank-deficient")
    return q_mat


def _complete_compressor_columns(
    primary_vec: np.ndarray,
    fallback: np.ndarray,
    num_filters: int,
    *,
    tol: float = 1e-12,
    ref_norm: float = 1.0,
) -> np.ndarray:
    """Build an M_q x R orthonormal compressor whose span contains primary_vec."""
    m_q = primary_vec.shape[0]
    norm_p = float(np.linalg.norm(primary_vec))
    candidates: list[np.ndarray] = []
    if norm_p > tol * max(ref_norm, 1e-300):
        candidates.append((primary_vec / norm_p)[:, None])
    candidates.append(fallback)
    candidates.append(np.eye(m_q, dtype=np.complex128))
    stacked = np.hstack(candidates)
    cols: list[np.ndarray] = []
    for col_idx in range(stacked.shape[1]):
        vec = stacked[:, col_idx].copy()
        for basis in cols:
            vec = vec - np.vdot(basis, vec) * basis
        vec_norm = float(np.linalg.norm(vec))
        if vec_norm > tol:
            cols.append(vec / vec_norm)
            if len(cols) == num_filters:
                break
    if len(cols) < num_filters:
        raise ValueError("unable to construct full-rank local compressor")
    return np.column_stack(cols)


def ptp_vs_sro_phase_budget(
    fs_hz: float = 48000.0,
    freq_hz: float = 2400.0,
    ptp_bound_s: float = 1e-6,
    sro_ppm: float = 100.0,
    n_fft: int = 1024,
    hop: int = 512,
) -> dict[str, object]:
    """Compare AES67/IEEE-1588 PTPv2 hardware sync vs free-running SRO drift."""
    fs = _validate_real_scalar(fs_hz, "fs_hz")
    freq = _validate_real_scalar(freq_hz, "freq_hz")
    ptp_s = _validate_real_scalar(ptp_bound_s, "ptp_bound_s")
    ppm = _validate_real_scalar(sro_ppm, "sro_ppm")
    n_fft_i = _validate_positive_int(n_fft, "n_fft", minimum=1)
    hop_i = _validate_positive_int(hop, "hop", minimum=1)
    if fs <= 0.0 or freq <= 0.0 or ptp_s < 0.0 or ppm <= 0.0:
        raise ValueError("fs_hz, freq_hz and sro_ppm must be > 0 and ptp_bound_s >= 0")

    sample_period_s = 1.0 / fs
    single_samples = ptp_s * fs
    single_phase_rad = 2.0 * math.pi * freq * ptp_s
    single_phase_deg = 360.0 * freq * ptp_s

    diff_s = 2.0 * ptp_s
    diff_samples = diff_s * fs
    diff_phase_rad = 2.0 * math.pi * freq * diff_s
    diff_phase_deg = 360.0 * freq * diff_s

    # Equal-weight 2-node coherent sum H = (1 + exp(-j * diff_phase)) / 2:
    # |H|^2 = cos^2(diff_phase / 2), |H - 1|^2 = sin^2(diff_phase / 2).
    half_phase = 0.5 * diff_phase_rad
    coherent_power_gain = math.cos(half_phase) ** 2
    coherent_attenuation_db = -10.0 * math.log10(coherent_power_gain)
    waveform_phase_mse = math.sin(half_phase) ** 2
    waveform_phase_mse_db = (
        10.0 * math.log10(waveform_phase_mse) if waveform_phase_mse > 0.0 else -math.inf
    )

    epsilon = ppm * 1e-6
    sample_drift_per_s = fs * epsilon
    seconds_per_sample_slip = 1.0 / sample_drift_per_s
    phase_drift_rad_per_s = 2.0 * math.pi * freq * epsilon
    phase_drift_deg_per_s = 360.0 * freq * epsilon
    first_null_s = 0.5 / (freq * epsilon)  # 2 * pi * f * eps * t = pi
    first_null_samples = sample_drift_per_s * first_null_s

    hop_duration_s = hop_i / fs
    drift_samples_per_hop = hop_i * epsilon
    phase_drift_deg_per_hop = 360.0 * freq * hop_duration_s * epsilon
    hops_per_sample_slip = 1.0 / drift_samples_per_hop

    return {
        "fs_hz": fs,
        "sample_period_s": sample_period_s,
        "sample_period_us": sample_period_s * 1e6,
        "freq_hz": freq,
        "ptp_hardware_sync": {
            "node_bound_s": ptp_s,
            "node_bound_us": ptp_s * 1e6,
            "single_node_samples": single_samples,
            "single_node_phase_rad": single_phase_rad,
            "single_node_phase_deg": single_phase_deg,
            "worst_two_node_offset_s": diff_s,
            "worst_two_node_offset_us": diff_s * 1e6,
            "worst_two_node_samples": diff_samples,
            "worst_two_node_phase_rad": diff_phase_rad,
            "worst_two_node_phase_deg": diff_phase_deg,
            "equal_sum_power_gain": coherent_power_gain,
            "equal_sum_attenuation_db": coherent_attenuation_db,
            "target_waveform_phase_mse": waveform_phase_mse,
            "target_waveform_phase_mse_db": waveform_phase_mse_db,
            "accumulates_over_time": False,
        },
        "free_running_sro": {
            "sro_ppm": ppm,
            "sro_ratio_epsilon": epsilon,
            "time_drift_us_per_s": ppm,
            "sample_drift_per_s": sample_drift_per_s,
            "seconds_per_integer_sample_slip": seconds_per_sample_slip,
            "phase_drift_rad_per_s": phase_drift_rad_per_s,
            "phase_drift_deg_per_s": phase_drift_deg_per_s,
            "first_destructive_null_s": first_null_s,
            "first_destructive_null_samples": first_null_samples,
            "wola": {
                "n_fft": n_fft_i,
                "hop": hop_i,
                "hop_duration_ms": hop_duration_s * 1e3,
                "drift_samples_per_hop": drift_samples_per_hop,
                "phase_drift_deg_per_hop": phase_drift_deg_per_hop,
                "hops_per_integer_sample_slip": hops_per_sample_slip,
            },
        },
        "scope": (
            "deterministic narrowband coherent-sum phase budget at f = 2.4 kHz, "
            "contrasting bounded static PTPv2 offset against linear crystal SRO drift"
        ),
    }


def centralized_lcmv(
    r_nn: object,
    constraints: object,
    response: object,
    *,
    condition_limit: float = 1e12,
) -> dict[str, object]:
    """Solve min_w w^H R_nn w subject to C^H w = f."""
    cov = _validate_spd_covariance(r_nn, "r_nn", condition_limit=condition_limit)
    c_mat = _validate_complex_array(constraints, "constraints", ndim=2)
    f_vec = _validate_complex_array(response, "response", ndim=1)
    m_dim, p_dim = c_mat.shape
    if cov.shape[0] != m_dim:
        raise ValueError("constraints row count must match r_nn dimension")
    if f_vec.shape[0] != p_dim:
        raise ValueError("response length must match constraints column count")
    if p_dim > m_dim:
        raise ValueError("number of linear constraints P cannot exceed dimension M")

    k_mat = np.linalg.solve(cov, c_mat)
    gram = 0.5 * (c_mat.conj().T @ k_mat + (c_mat.conj().T @ k_mat).conj().T)
    gram_eigs = np.linalg.eigvalsh(gram)
    if float(gram_eigs[-1]) <= 0.0 or float(gram_eigs[0]) <= 1e-12 * float(gram_eigs[-1]):
        raise ValueError("constraint matrix C is rank-deficient")
    lagrange = np.linalg.solve(gram, f_vec)
    weights = k_mat @ lagrange
    noise_power = float(np.vdot(weights, cov @ weights).real)
    constraint_residual = float(np.max(np.abs(c_mat.conj().T @ weights - f_vec)))
    return {
        "weights": weights,
        "noise_power": noise_power,
        "lagrange_multipliers": lagrange,
        "constraint_max_abs_residual": constraint_residual,
    }


def dasf_lcmv_updates(
    r_nn: object,
    constraints: object,
    response: object,
    *,
    node_sizes: tuple[int, ...] = (2, 2),
    num_filters: int = 1,
    rounds: int = 12,
    initial: list[object] | None = None,
) -> dict[str, object]:
    """Run sequential node-updating DASF-LCMV across a fully connected WASN."""
    cov = _validate_spd_covariance(r_nn, "r_nn")
    c_mat = _validate_complex_array(constraints, "constraints", ndim=2)
    f_vec = _validate_complex_array(response, "response", ndim=1)
    if not isinstance(node_sizes, (tuple, list)) or len(node_sizes) < 2:
        raise ValueError("node_sizes must specify at least 2 nodes")
    sizes = [_validate_positive_int(sz, "node_size", minimum=1) for sz in node_sizes]
    if sum(sizes) != cov.shape[0]:
        raise ValueError("sum(node_sizes) must equal r_nn dimension")
    r_filters = _validate_positive_int(num_filters, "num_filters", minimum=1)
    n_rounds = _validate_positive_int(rounds, "rounds", minimum=1)
    p_dim = c_mat.shape[1]

    offsets = np.cumsum([0] + sizes)
    slices = [slice(int(offsets[u]), int(offsets[u + 1])) for u in range(len(sizes))]

    for u, m_u in enumerate(sizes):
        if r_filters > m_u:
            raise ValueError(f"num_filters ({r_filters}) cannot exceed M_{u} ({m_u})")
        l_u = m_u + (len(sizes) - 1) * r_filters
        if l_u < p_dim:
            raise ValueError(
                f"compressed dimension L_{u}={l_u} is smaller than P={p_dim} constraints"
            )

    compressors: list[np.ndarray] = []
    if initial is None:
        for m_u in sizes:
            compressors.append(np.eye(m_u, r_filters, dtype=np.complex128))
    else:
        if len(initial) != len(sizes):
            raise ValueError("initial compressor list length must match node count")
        for u, (raw_f, m_u) in enumerate(zip(initial, sizes)):
            arr = _validate_complex_array(raw_f, f"initial[{u}]")
            if arr.ndim == 1:
                arr = arr[:, None]
            if arr.shape != (m_u, r_filters):
                raise ValueError(
                    f"initial[{u}] must have shape ({m_u}, {r_filters}), got {arr.shape}"
                )
            compressors.append(_orthonormal_columns(arr))

    central = centralized_lcmv(cov, c_mat, f_vec)
    w_star = central["weights"]
    p_star = float(central["noise_power"])

    history: list[dict[str, object]] = []
    current_weights = np.zeros(cov.shape[0], dtype=np.complex128)
    for rnd in range(n_rounds):
        for q_node, (sl_q, m_q) in enumerate(zip(slices, sizes)):
            l_q = m_q + (len(sizes) - 1) * r_filters
            t_q = np.zeros((cov.shape[0], l_q), dtype=np.complex128)
            t_q[sl_q, :m_q] = np.eye(m_q, dtype=np.complex128)
            col_cursor = m_q
            for v_node, sl_v in enumerate(slices):
                if v_node == q_node:
                    continue
                t_q[sl_v, col_cursor : col_cursor + r_filters] = compressors[v_node]
                col_cursor += r_filters

            r_tilde = 0.5 * (t_q.conj().T @ cov @ t_q + (t_q.conj().T @ cov @ t_q).conj().T)
            c_tilde = t_q.conj().T @ c_mat
            local_sol = centralized_lcmv(r_tilde, c_tilde, f_vec)
            h_q = local_sol["weights"]
            current_weights = t_q @ h_q
            w_local_block = current_weights[sl_q]
            compressors[q_node] = _complete_compressor_columns(
                w_local_block,
                compressors[q_node],
                r_filters,
                ref_norm=float(np.linalg.norm(current_weights)),
            )
            history.append(
                {
                    "round": rnd + 1,
                    "updating_node": q_node + 1,
                    "compressed_dim": l_q,
                    "remaining_dofs": l_q - p_dim,
                    "noise_power": float(local_sol["noise_power"]),
                    "weight_error_l2": float(np.linalg.norm(current_weights - w_star)),
                    "constraint_residual": float(
                        np.max(np.abs(c_mat.conj().T @ current_weights - f_vec))
                    ),
                }
            )

    return {
        "node_sizes": sizes,
        "num_filters": r_filters,
        "num_constraints": p_dim,
        "compressed_dims": [m_u + (len(sizes) - 1) * r_filters for m_u in sizes],
        "remaining_dofs_per_node": [
            m_u + (len(sizes) - 1) * r_filters - p_dim for m_u in sizes
        ],
        "centralized_weights": w_star,
        "centralized_noise_power": p_star,
        "final_weights": current_weights,
        "final_noise_power": float(history[-1]["noise_power"]),
        "final_weight_error_l2": float(history[-1]["weight_error_l2"]),
        "monotone_nonincreasing": all(
            float(history[idx + 1]["noise_power"])
            <= float(history[idx]["noise_power"])
            + 1e-12 * max(1.0, abs(float(history[idx]["noise_power"])))
            for idx in range(len(history) - 1)
        ),
        "history": history,
    }


def dasf_unification_control() -> dict[str, object]:
    """E13-27: MWF/MVDR/GEVD 1D equivalence, 2-constraint DASF-LCMV & P=3 lock-up."""
    # Exact Ch 13 §13.8 2-node 4-mic correlated-noise benchmark:
    a_vec = np.array([1.0, 0.5, 2.0, -0.5], dtype=np.float64)
    r_ss = np.outer(a_vec, a_vec)  # phi_s = 1
    r_nn = np.array(
        [
            [1.0, 0.0, 0.2, 0.8],
            [0.0, 1.0, 0.0, 0.0],
            [0.2, 0.0, 1.0, 0.0],
            [0.8, 0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )
    r_xx = r_ss + r_nn

    # Part 1: Centralized MWF, MVDR, and Principal GEVD share R_nn^(-1) a = [25, 4, 11, -24]^T / 8
    v_opt = np.linalg.solve(r_nn, a_vec)  # [3.125, 0.5, 1.375, -3.0]
    q_snr = float(a_vec @ v_opt)  # 61 / 8 = 7.625
    w_mwf = np.linalg.solve(r_xx, r_ss[:, 0])  # [25, 4, 11, -24]^T / 69
    mvdr_sol = centralized_lcmv(r_nn, a_vec[:, None], np.array([1.0]))
    w_mvdr = mvdr_sol["weights"].real  # [25, 4, 11, -24]^T / 61

    # Principal GEVD of (R_xx, R_nn) via Cholesky whitening L L^T = R_nn, w^T R_nn w = 1
    chol_l = np.linalg.cholesky(r_nn)
    whitened_rx = np.linalg.solve(chol_l, np.linalg.solve(chol_l, r_xx).T)
    whitened_rx = 0.5 * (whitened_rx + whitened_rx.T)
    gevd_vals, gevd_vecs = np.linalg.eigh(whitened_rx)
    lambda_max = float(gevd_vals[-1])  # 69 / 8 = 8.625
    w_gevd = np.linalg.solve(chol_l.T, gevd_vecs[:, -1])
    if w_gevd[0] < 0.0:
        w_gevd = -w_gevd

    # Compressed 3-channel Node-1 model under z_2 = 11 x_3 - 24 x_4
    f2_mvdr = np.array([11.0, -24.0], dtype=np.float64)
    t1_mvdr = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 11.0],
            [0.0, 0.0, -24.0],
        ],
        dtype=np.float64,
    )
    r_ss_tilde = t1_mvdr.T @ r_ss @ t1_mvdr
    r_nn_tilde = t1_mvdr.T @ r_nn @ t1_mvdr
    r_xx_tilde = r_ss_tilde + r_nn_tilde
    a_tilde = t1_mvdr.T @ a_vec  # [1, 0.5, 34]^T

    h_mwf_comp = np.linalg.solve(r_xx_tilde, t1_mvdr.T @ r_ss[:, 0])
    w_mwf_comp = t1_mvdr @ h_mwf_comp
    mvdr_comp_sol = centralized_lcmv(r_nn_tilde, a_tilde[:, None], np.array([1.0]))
    h_mvdr_comp = mvdr_comp_sol["weights"].real
    w_mvdr_comp = t1_mvdr @ h_mvdr_comp

    chol_tilde = np.linalg.cholesky(r_nn_tilde)
    whitened_tilde = np.linalg.solve(
        chol_tilde, np.linalg.solve(chol_tilde, r_xx_tilde).T
    )
    gevd_tilde_vals, gevd_tilde_vecs = np.linalg.eigh(
        0.5 * (whitened_tilde + whitened_tilde.T)
    )
    h_gevd_comp = np.linalg.solve(chol_tilde.T, gevd_tilde_vecs[:, -1])
    w_gevd_comp = t1_mvdr @ h_gevd_comp
    if w_gevd_comp[0] < 0.0:
        w_gevd_comp = -w_gevd_comp

    # Part 2: Two-constraint LCMV C = [a, b], f = [1, 0]^T with interferer b = [0, 1, 1, 0]^T
    b_vec = np.array([0.0, 1.0, 1.0, 0.0], dtype=np.float64)
    c_two = np.column_stack([a_vec, b_vec])
    f_two = np.array([1.0, 0.0], dtype=np.float64)
    lcmv_two_central = centralized_lcmv(r_nn, c_two, f_two)

    # Stale broadcast control: Node 2 keeps [11, -24]^T while Node 1 solves 2-constraint LCMV
    c_two_tilde = t1_mvdr.T @ c_two
    stale_lcmv_sol = centralized_lcmv(r_nn_tilde, c_two_tilde, f_two)
    h_stale = stale_lcmv_sol["weights"].real  # [5210, -2002, 182]^T / 10397
    w_stale = t1_mvdr @ h_stale

    # Sequential DASF-LCMV starting from [1, 0]^T and [11, -24]^T
    dasf_two = dasf_lcmv_updates(
        r_nn,
        c_two,
        f_two,
        node_sizes=(2, 2),
        num_filters=1,
        rounds=1400,
        initial=[np.array([1.0, 0.0]), f2_mvdr],
    )

    # Part 3: P = 3 constraints with R = 1 (L_u = 3 = P -> 0 remaining DOFs -> lock-up)
    # vs R = 2 (L_u = 4 > P = 3 -> centralized optimum recovered)
    c3_vec = np.array([0.0, 0.0, 0.0, 1.0], dtype=np.float64)
    c_three = np.column_stack([a_vec, b_vec, c3_vec])
    f_three = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    dasf_lockup_r1 = dasf_lcmv_updates(
        r_nn,
        c_three,
        f_three,
        node_sizes=(2, 2),
        num_filters=1,
        rounds=6,
        initial=[np.array([1.0, 0.0]), f2_mvdr],
    )
    dasf_unlocked_r2 = dasf_lcmv_updates(
        r_nn,
        c_three,
        f_three,
        node_sizes=(2, 2),
        num_filters=2,
        rounds=2,
    )

    return {
        "rank1_shared_subspace": {
            "r_nn_inv_a": v_opt,
            "q_snr": q_snr,
            "mwf_central_weights": w_mwf,
            "mwf_mse": float(1.0 - a_vec @ w_mwf),
            "mvdr_central_weights": w_mvdr,
            "mvdr_min_noise_power": float(mvdr_sol["noise_power"]),
            "gevd_lambda_max": lambda_max,
            "gevd_central_weights": w_gevd,
            "compressed_node1_mwf_error_l2": float(np.linalg.norm(w_mwf_comp - w_mwf)),
            "compressed_node1_mvdr_error_l2": float(
                np.linalg.norm(w_mvdr_comp - w_mvdr)
            ),
            "compressed_node1_gevd_error_l2": float(
                np.linalg.norm(w_gevd_comp - w_gevd)
            ),
            "compressed_node1_mvdr_h": h_mvdr_comp,
            "compressed_node1_gevd_lambda_max": float(gevd_tilde_vals[-1]),
        },
        "two_constraint_lcmv": {
            "interferer_b": b_vec,
            "central_weights": lcmv_two_central["weights"].real,
            "central_noise_power": float(lcmv_two_central["noise_power"]),
            "stale_mvdr_subspace_h": h_stale,
            "stale_mvdr_subspace_weights": w_stale,
            "stale_mvdr_subspace_noise_power": float(stale_lcmv_sol["noise_power"]),
            "stale_excess_noise_power": float(
                stale_lcmv_sol["noise_power"] - lcmv_two_central["noise_power"]
            ),
            "dasf_r1_final_weights": dasf_two["final_weights"].real,
            "dasf_r1_final_noise_power": dasf_two["final_noise_power"],
            "dasf_r1_weight_error_l2": dasf_two["final_weight_error_l2"],
            "dasf_r1_monotone": dasf_two["monotone_nonincreasing"],
        },
        "three_constraint_dof_boundary": {
            "third_constraint": c3_vec,
            "central_weights": dasf_lockup_r1["centralized_weights"].real,
            "central_noise_power": dasf_lockup_r1["centralized_noise_power"],
            "r1_locked_remaining_dofs": dasf_lockup_r1["remaining_dofs_per_node"],
            "r1_locked_weights": dasf_lockup_r1["final_weights"].real,
            "r1_locked_noise_power": dasf_lockup_r1["final_noise_power"],
            "r2_unlocked_remaining_dofs": dasf_unlocked_r2["remaining_dofs_per_node"],
            "r2_unlocked_weights": dasf_unlocked_r2["final_weights"].real,
            "r2_unlocked_noise_power": dasf_unlocked_r2["final_noise_power"],
            "r2_unlocked_weight_error_l2": dasf_unlocked_r2["final_weight_error_l2"],
        },
        "scope": (
            "known-covariance DASF teaching control: MWF/MVDR/rank-1 GEVD share a 1D "
            "broadcast subspace, 2-constraint LCMV requires subspace adaptation "
            "(34/203 vs stale 2116/10397), and P = L_u = 3 locks single-channel DASF"
        ),
    }
