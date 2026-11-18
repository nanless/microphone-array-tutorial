"""Kronecker, Khatri-Rao, and column-major vectorization kernels for Appendix A.

Implements the algebraic identities linking coarray DOA (Ch 3), separable planar
beamforming (Ch 5), and covariance matrix fitting / CMF (Ch 12):

    vec(A diag(q) B^T) = (B (*) A) q,
    vec(A diag(q) A^H) = (A^* (*) A) q,

where ``vec(.)`` stacks columns top-to-bottom (column-major / Fortran order) and
``(*)`` denotes the column-wise Khatri-Rao product. Pure NumPy, no SciPy or I/O.
"""
from __future__ import annotations

import math
import numpy as np


def _validate_2d_numeric(value: object, name: str) -> np.ndarray:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be numeric, not bool")
    raw = np.asarray(value)
    if raw.dtype.kind not in "iufc":
        raise ValueError(f"{name} must contain numeric values, not bool or text")
    if isinstance(value, (list, tuple)):
        def _check_no_bool(items: object) -> None:
            if isinstance(items, (list, tuple)):
                for item in items:
                    _check_no_bool(item)
            elif isinstance(items, (bool, np.bool_)):
                raise ValueError(f"{name} must not contain bool")

        _check_no_bool(value)
    target_dtype = np.complex128 if raw.dtype.kind == "c" else np.float64
    with np.errstate(over="ignore", invalid="ignore"):
        arr = np.asarray(raw, dtype=target_dtype)
    if arr.ndim != 2 or min(arr.shape) < 1:
        raise ValueError(f"{name} must be a non-empty 2D matrix")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    return arr


def _complex_record(value: np.ndarray) -> dict[str, object]:
    arr = np.asarray(value, dtype=np.complex128)
    return {"real": arr.real.tolist(), "imag": arr.imag.tolist()}


def column_major_vec(matrix: object) -> np.ndarray:
    """Stack columns of a 2D matrix (M, N) into a 1D vector of length M*N."""
    arr = _validate_2d_numeric(matrix, "matrix")
    return np.reshape(arr, (-1,), order="F").copy()


def khatri_rao_product(left: object, right: object) -> np.ndarray:
    """Compute column-wise Kronecker (Khatri-Rao) product (P*M, K) of (P, K) and (M, K)."""
    left_arr = _validate_2d_numeric(left, "left")
    right_arr = _validate_2d_numeric(right, "right")
    p_rows, k_left = left_arr.shape
    m_rows, k_right = right_arr.shape
    if k_left != k_right:
        raise ValueError(
            f"left and right must have equal column count, got {k_left} and {k_right}"
        )
    out_dtype = (
        np.complex128
        if (np.iscomplexobj(left_arr) or np.iscomplexobj(right_arr))
        else np.float64
    )
    with np.errstate(over="raise", invalid="raise"):
        blocks = left_arr[:, None, :].astype(out_dtype) * right_arr[None, :, :].astype(
            out_dtype
        )
    out = blocks.reshape(p_rows * m_rows, k_left)
    if not np.all(np.isfinite(out)):
        raise ValueError("Khatri-Rao product exceeds float64 range")
    return out


def kronecker_khatri_rao_controls() -> dict[str, object]:
    """E16-21: deterministic JSON-safe verification of Kronecker & Khatri-Rao identities."""
    a_mat = np.array([[1.0, 1.0], [1.0j, -1.0]], dtype=np.complex128)
    q_vec = np.array([2.0, 1.0], dtype=np.float64)

    # 1. Hermitian CSM identity: vec(A diag(q) A^H) = (A^* odot A) q
    r_csm = a_mat @ np.diag(q_vec) @ a_mat.conj().T
    vec_r_col = column_major_vec(r_csm)
    vec_r_row_bug = r_csm.reshape(-1)  # C-order row-major swaps indices (1,0) and (0,1)
    h_dict = khatri_rao_product(a_mat.conj(), a_mat)
    hq_vec = h_dict @ q_vec

    # Missing-conjugate bug: (A odot A) q = vec(A diag(q) A^T), giving -1 at (2, 2)
    h_missing_conj = khatri_rao_product(a_mat, a_mat)
    missing_conj_vec = h_missing_conj @ q_vec

    # 2. General rectangular / non-Hermitian Roth diagonal identity:
    # vec(A diag(q) B^T) = (B odot A) q
    b_mat = np.array([[2.0, 0.0], [1.0, 1.0j]], dtype=np.complex128)
    rect_prod = a_mat @ np.diag(q_vec) @ b_mat.T
    vec_rect = column_major_vec(rect_prod)
    kr_rect = khatri_rao_product(b_mat, a_mat) @ q_vec

    # Full Roth column lemma check: vec(A X B^T) = (B kron A) vec(X) for X = diag(q)
    x_diag = np.diag(q_vec).astype(np.complex128)
    roth_vec = np.kron(b_mat, a_mat) @ column_major_vec(x_diag)

    # 3. Separable 2x2 planar Kronecker beamformer:
    # (w_y kron w_x)^H (a_y kron a_x) = (w_y^H a_y)(w_x^H a_x) = 1
    a_x = np.array([1.0, 1.0j], dtype=np.complex128)
    a_y = np.array([1.0, -1.0], dtype=np.complex128)
    w_x = a_x / 2.0
    w_y = a_y / 2.0
    a_planar = np.kron(a_y, a_x)
    w_planar = np.kron(w_y, w_x)
    response_joint = complex(np.vdot(w_planar, a_planar))
    response_factored = complex(np.vdot(w_y, a_y) * np.vdot(w_x, a_x))
    w_norm_sq = float(np.vdot(w_planar, w_planar).real)
    wng_linear = float(abs(response_joint) ** 2 / w_norm_sq)
    wng_db = 10.0 * math.log10(wng_linear)

    return {
        "exercise_id": "E16-21",
        "A": _complex_record(a_mat),
        "q": q_vec.tolist(),
        "hermitian_csm": {
            "R": _complex_record(r_csm),
            "vec_R_col_major": _complex_record(vec_r_col),
            "vec_R_row_major_bug": _complex_record(vec_r_row_bug),
            "row_vs_col_imaginary_discrepancy_norm": float(
                np.linalg.norm(vec_r_col - vec_r_row_bug)
            ),
            "H_khatri_rao": _complex_record(h_dict),
            "H_times_q": _complex_record(hq_vec),
            "identity_max_abs_error": float(np.max(np.abs(hq_vec - vec_r_col))),
            "missing_conjugate_bug_vector": _complex_record(missing_conj_vec),
            "missing_conjugate_last_entry_real": float(missing_conj_vec[-1].real),
        },
        "rectangular_roth_diagonal": {
            "B": _complex_record(b_mat),
            "A_diag_q_B_transpose": _complex_record(rect_prod),
            "vec_col_major": _complex_record(vec_rect),
            "khatri_rao_B_A_times_q": _complex_record(kr_rect),
            "full_kronecker_roth_vec": _complex_record(roth_vec),
            "khatri_rao_max_abs_error": float(np.max(np.abs(kr_rect - vec_rect))),
            "roth_lemma_max_abs_error": float(np.max(np.abs(roth_vec - vec_rect))),
        },
        "separable_kronecker_beamformer": {
            "a_x": _complex_record(a_x),
            "a_y": _complex_record(a_y),
            "a_planar": _complex_record(a_planar),
            "w_planar": _complex_record(w_planar),
            "distortionless_gain_real": float(response_joint.real),
            "distortionless_gain_imag": float(response_joint.imag),
            "mixed_product_error": float(abs(response_joint - response_factored)),
            "weight_norm_squared": w_norm_sq,
            "wng_linear": wng_linear,
            "wng_db": wng_db,
        },
        "scope": (
            "exact column-major vec / Kronecker / Khatri-Rao algebraic identities; "
            "row-major ravel swaps off-diagonal conjugates and omitting A^* breaks PSD diagonal entries"
        ),
    }
