"""Binaural Signal Matching (BSM-LS, BSM-MagLS) and interaural coherence control.

Under the book's inner-product beamformer convention ``y_e = w_e^H x``, the
synthesized directional ear transfer across Q directions ``H in C^(M x Q)`` is
``p_hat_e = (w_e^H H)^T = H^T w_e^*`` (equivalently ``H^H w_e = p_hat_e^*``).
"""
from __future__ import annotations

import numpy as np

from codes.chapters.ch02.core.conventions import (
    finite_real_array,
    finite_real_scalar,
    hermitian_part,
)


def _validate_complex_array(value, name: str) -> np.ndarray:
    original = np.asarray(value)
    if original.dtype.kind not in "iufc":
        raise ValueError(f"{name} must contain numeric values")
    with np.errstate(over="ignore", invalid="ignore"):
        arr = np.asarray(original, dtype=np.complex128)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def _validate_quadrature_weights(weights, n_dirs: int) -> np.ndarray:
    if weights is None:
        return np.ones(n_dirs, dtype=np.float64)
    w = finite_real_array(weights, "weights")
    if w.shape != (n_dirs,):
        raise ValueError(f"weights must have shape ({n_dirs},)")
    if np.any(w < 0.0) or float(np.sum(w)) <= 0.0:
        raise ValueError("weights must be non-negative with positive sum")
    return w


def _factorize_normal_matrix(H: np.ndarray, w_omega: np.ndarray, reg: float) -> np.ndarray:
    n_mics = H.shape[0]
    with np.errstate(over="raise", invalid="raise"):
        try:
            weighted_h = H * w_omega[np.newaxis, :]
            gram = hermitian_part(weighted_h @ H.conj().T)
            normal = gram + reg * np.eye(n_mics, dtype=np.complex128)
        except FloatingPointError as exc:
            raise ValueError("normal matrix overflowed float64 range") from exc
    try:
        chol = np.linalg.cholesky(normal)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "BSM normal matrix is not positive definite; increase reg or directional span"
        ) from exc
    return chol


def _cholesky_solve(chol: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    y = np.linalg.solve(chol, rhs)
    return np.linalg.solve(chol.conj().T, y)


def bsm_ls_weights(
    H: np.ndarray,
    p_ear: np.ndarray,
    weights: np.ndarray | None = None,
    reg: float = 0.0,
) -> np.ndarray:
    """Solve complex least-squares BSM weights ``w_{e,LS}``.

    Minimizes ``||W_Omega^(1/2) (H^H w_e - p_ear^*)||_2^2 + reg * ||w_e||_2^2``
    under ``p_hat_e = (w_e^H H)^T``.
    """
    h_mat = _validate_complex_array(H, "H")
    if h_mat.ndim != 2 or min(h_mat.shape) < 1:
        raise ValueError("H must be a non-empty 2-D array of shape (M, Q)")
    n_mics, n_dirs = h_mat.shape

    p_vec = _validate_complex_array(p_ear, "p_ear")
    if p_vec.shape != (n_dirs,):
        raise ValueError(f"p_ear must have shape ({n_dirs},)")

    reg_val = finite_real_scalar(reg, "reg")
    if reg_val < 0.0:
        raise ValueError("reg must be non-negative")

    w_omega = _validate_quadrature_weights(weights, n_dirs)
    chol = _factorize_normal_matrix(h_mat, w_omega, reg_val)
    with np.errstate(over="raise", invalid="raise"):
        try:
            rhs = (h_mat * w_omega[np.newaxis, :]) @ p_vec.conj()
            w = _cholesky_solve(chol, rhs)
        except FloatingPointError as exc:
            raise ValueError("BSM-LS right-hand side or solve overflowed float64") from exc
    if w.shape != (n_mics,) or not np.all(np.isfinite(w)):
        raise ValueError("BSM-LS solution is non-finite")
    return w


def bsm_magls_weights(
    H: np.ndarray,
    p_ear_mag: np.ndarray,
    weights: np.ndarray | None = None,
    reg: float = 0.0,
    n_iters: int = 20,
    init_phases: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve Magnitude Least-Squares BSM via alternating phase projection.

    Returns ``(w, objective_history)`` where ``objective_history`` records the
    joint regularized objective
    ``J(w, phi) = sum_q w_q (|w^H h_q| - m_q)^2 + reg * ||w||_2^2``
    after each iteration ``1 .. n_iters``, guaranteed monotonically non-increasing.
    """
    h_mat = _validate_complex_array(H, "H")
    if h_mat.ndim != 2 or min(h_mat.shape) < 1:
        raise ValueError("H must be a non-empty 2-D array of shape (M, Q)")
    n_mics, n_dirs = h_mat.shape

    mag = finite_real_array(p_ear_mag, "p_ear_mag")
    if mag.shape != (n_dirs,) or np.any(mag < 0.0):
        raise ValueError(f"p_ear_mag must be non-negative real of shape ({n_dirs},)")

    reg_val = finite_real_scalar(reg, "reg")
    if reg_val < 0.0:
        raise ValueError("reg must be non-negative")

    if isinstance(n_iters, (bool, np.bool_)) or not isinstance(n_iters, (int, np.integer)) or n_iters < 1:
        raise ValueError("n_iters must be a positive integer")
    n_iters_int = int(n_iters)

    if init_phases is None:
        unit_phasor = np.ones(n_dirs, dtype=np.complex128)
    else:
        phases = finite_real_array(init_phases, "init_phases")
        if phases.shape != (n_dirs,):
            raise ValueError(f"init_phases must have shape ({n_dirs},)")
        unit_phasor = np.exp(1j * phases)

    w_omega = _validate_quadrature_weights(weights, n_dirs)
    chol = _factorize_normal_matrix(h_mat, w_omega, reg_val)
    weighted_h = h_mat * w_omega[np.newaxis, :]

    history = np.empty(n_iters_int, dtype=np.float64)
    w_curr = np.zeros(n_mics, dtype=np.complex128)
    prev_obj = np.inf

    for it in range(n_iters_int):
        with np.errstate(over="raise", invalid="raise"):
            try:
                # Step 1: exact global LS minimization over w given current target phases.
                target_conj = mag * unit_phasor.conj()
                rhs = weighted_h @ target_conj
                w_next = _cholesky_solve(chol, rhs)

                # Synthesized response p_hat = (w^H H)^T = (H^H w)^*
                p_hat = (h_mat.conj().T @ w_next).conj()
                p_hat_abs = np.abs(p_hat)

                obj = float(
                    np.sum(w_omega * ((p_hat_abs - mag) ** 2))
                    + reg_val * np.real(np.vdot(w_next, w_next))
                )
            except FloatingPointError as exc:
                raise ValueError("BSM-MagLS iterate or objective overflowed float64") from exc

        if not (np.all(np.isfinite(w_next)) and np.all(np.isfinite(p_hat_abs)) and np.isfinite(obj)):
            raise ValueError("BSM-MagLS iterate or objective is non-finite")

        # Step 2: pointwise optimal phase projection onto unit circle.
        nonzero = p_hat_abs > 0.0
        next_phasor = unit_phasor.copy()
        next_phasor[nonzero] = p_hat[nonzero] / p_hat_abs[nonzero]

        # Guard against +1e-16 floating-point jitter once already at a fixed point.
        if obj > prev_obj:
            obj = prev_obj
        else:
            w_curr = w_next
            unit_phasor = next_phasor
            prev_obj = obj
        history[it] = obj

    return w_curr, history


def bsm_directional_metrics(
    H: np.ndarray,
    p_ear: np.ndarray,
    w: np.ndarray,
    weights: np.ndarray | None = None,
) -> dict:
    """Evaluate complex LS error and magnitude error of BSM filter ``w``."""
    h_mat = _validate_complex_array(H, "H")
    if h_mat.ndim != 2 or min(h_mat.shape) < 1:
        raise ValueError("H must be 2-D of shape (M, Q)")
    n_mics, n_dirs = h_mat.shape
    p_vec = _validate_complex_array(p_ear, "p_ear")
    w_vec = _validate_complex_array(w, "w")
    if p_vec.shape != (n_dirs,) or w_vec.shape != (n_mics,):
        raise ValueError("incompatible shapes in bsm_directional_metrics")
    w_omega = _validate_quadrature_weights(weights, n_dirs)
    w_norm = w_omega / np.sum(w_omega)

    with np.errstate(over="raise", invalid="raise"):
        try:
            p_hat = (w_vec.conj() @ h_mat).ravel()
            p_hat_abs = np.abs(p_hat)
            complex_err = float(np.sum(w_norm * (np.abs(p_hat - p_vec) ** 2)))
            mag_err = float(np.sum(w_norm * ((p_hat_abs - np.abs(p_vec)) ** 2)))
        except FloatingPointError as exc:
            raise ValueError("BSM directional metrics overflowed float64") from exc
    if not (np.all(np.isfinite(p_hat)) and np.isfinite(complex_err) and np.isfinite(mag_err)):
        raise ValueError("BSM directional metrics are non-finite")
    return {
        "synthesized_response": p_hat,
        "synthesized_magnitude": p_hat_abs,
        "weighted_complex_mse": complex_err,
        "weighted_magnitude_mse": mag_err,
    }


def binaural_coherence_with_mixing(
    R_nn: np.ndarray,
    w_L: np.ndarray,
    w_R: np.ndarray,
    beta: float = 0.0,
    ref_idx: tuple[int, int] = (0, 1),
) -> dict:
    """Compute output binaural noise covariance and interaural coherence under reference mixing.

    Blends ``(1 - beta)`` of beamformer weights ``(w_L, w_R)`` with ``beta`` of the
    raw ear-reference microphones ``ref_idx = (idx_L, idx_R)``:
    ``w_tilde_e = (1 - beta) * w_e + beta * e_{idx_e}``.
    """
    r_raw = _validate_complex_array(R_nn, "R_nn")
    if r_raw.ndim != 2 or r_raw.shape[0] != r_raw.shape[1] or r_raw.shape[0] < 1:
        raise ValueError("R_nn must be a non-empty square matrix")
    n_mics = r_raw.shape[0]

    scale_r = float(np.max(np.abs(r_raw)))
    if scale_r > 0.0:
        r_scaled = r_raw / scale_r
        asym = float(np.linalg.norm(r_scaled - r_scaled.conj().T, ord="fro"))
        norm_r = float(np.linalg.norm(r_scaled, ord="fro"))
        if norm_r > 0.0 and asym > 1e-10 * norm_r:
            raise ValueError("R_nn must be Hermitian")
        r_unit = hermitian_part(r_scaled)
        eigvals = np.linalg.eigvalsh(r_unit)
        scale_eig = max(float(np.max(np.abs(eigvals))), 1.0)
        if float(np.min(eigvals)) < -1e-11 * scale_eig:
            raise ValueError("R_nn must be positive semidefinite")
    else:
        r_unit = np.zeros_like(r_raw)
    r_mat = hermitian_part(r_raw)

    wl = _validate_complex_array(w_L, "w_L")
    wr = _validate_complex_array(w_R, "w_R")
    if wl.shape != (n_mics,) or wr.shape != (n_mics,):
        raise ValueError(f"w_L and w_R must have shape ({n_mics},)")

    beta_val = finite_real_scalar(beta, "beta")
    if not (0.0 <= beta_val <= 1.0):
        raise ValueError("beta must lie in [0, 1]")

    if not isinstance(ref_idx, (tuple, list)) or len(ref_idx) != 2:
        raise ValueError("ref_idx must be a 2-element tuple/list of microphone indices")
    idx_l, idx_r = ref_idx
    for idx in (idx_l, idx_r):
        if isinstance(idx, (bool, np.bool_)) or not isinstance(idx, (int, np.integer)) or not (0 <= idx < n_mics):
            raise ValueError("ref_idx entries must be valid integer microphone indices")

    el = np.zeros(n_mics, dtype=np.complex128)
    er = np.zeros(n_mics, dtype=np.complex128)
    el[int(idx_l)] = 1.0
    er[int(idx_r)] = 1.0

    wl_mix = (1.0 - beta_val) * wl + beta_val * el
    wr_mix = (1.0 - beta_val) * wr + beta_val * er
    if not (np.all(np.isfinite(wl_mix)) and np.all(np.isfinite(wr_mix))):
        raise ValueError("mixed binaural weights are non-finite")

    with np.errstate(over="raise", invalid="raise"):
        try:
            var_l = max(0.0, float(np.real(np.vdot(wl_mix, r_mat @ wl_mix))))
            var_r = max(0.0, float(np.real(np.vdot(wr_mix, r_mat @ wr_mix))))
            cross_lr = complex(np.vdot(wl_mix, r_mat @ wr_mix))
        except FloatingPointError as exc:
            raise ValueError("output binaural covariance overflowed float64") from exc
    if not (np.isfinite(var_l) and np.isfinite(var_r) and np.isfinite(cross_lr)):
        raise ValueError("output binaural covariance is non-finite")

    wl_scale = float(np.max(np.abs(wl_mix)))
    wr_scale = float(np.max(np.abs(wr_mix)))
    if scale_r > 0.0 and wl_scale > 0.0 and wr_scale > 0.0:
        if scale_r > 1e100 or scale_r < 1e-100 or wl_scale > 1e100 or wl_scale < 1e-100 or wr_scale > 1e100 or wr_scale < 1e-100:
            wl_u = wl_mix / wl_scale
            wr_u = wr_mix / wr_scale
            vl_u = max(0.0, float(np.real(np.vdot(wl_u, r_unit @ wl_u))))
            vr_u = max(0.0, float(np.real(np.vdot(wr_u, r_unit @ wr_u))))
            cr_u = complex(np.vdot(wl_u, r_unit @ wr_u))
            sd_l = float(np.sqrt(vl_u))
            sd_r = float(np.sqrt(vr_u))
            coh_complex = (cr_u / sd_l) / sd_r if (sd_l > 0.0 and sd_r > 0.0) else (0.0 + 0.0j)
        else:
            sd_l = float(np.sqrt(var_l))
            sd_r = float(np.sqrt(var_r))
            coh_complex = (cross_lr / sd_l) / sd_r if (sd_l > 0.0 and sd_r > 0.0) else (0.0 + 0.0j)
        ic_out = float(np.clip(np.abs(coh_complex), 0.0, 1.0))
    else:
        coh_complex = 0.0 + 0.0j
        ic_out = 0.0
    if not (np.isfinite(coh_complex) and np.isfinite(ic_out)):
        raise ValueError("output interaural coherence is non-finite")

    return {
        "weights_L": wl_mix,
        "weights_R": wr_mix,
        "var_L": var_l,
        "var_R": var_r,
        "cross_covariance": cross_lr,
        "complex_coherence": coh_complex,
        "ic_out": ic_out,
    }
