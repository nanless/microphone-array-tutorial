"""E05-16 & E05-17: first-order differential microphone arrays (DMA) and Kronecker beamforming.

Pure-NumPy deterministic implementations reusing Chapter 5 MVDR/LCMV building blocks.
Follows inner-product output y = w^H x and forward Fourier kernel exp(-j * 2*pi*f*t).
"""

from __future__ import annotations

import math

import numpy as np

from codes.chapters.ch02.core.conventions import SPEED_OF_SOUND, finite_real_scalar
from codes.chapters.ch05.core.beamforming import lcmv_weights as _base_lcmv_weights
from codes.chapters.ch05.core.beamforming import mvdr_weights


def first_order_dma_weights(
    freq_hz: float,
    d_m: float,
    null_angle_deg: float,
    c: float = SPEED_OF_SOUND,
    equalize: bool = True,
) -> np.ndarray:
    """Design 2-mic endfire first-order DMA weights under inner-product convention y = w^H x.

    Endfire angle ``vartheta`` is measured from front target ``0 deg`` to rear ``180 deg``,
    so steering vector is ``a(vartheta) = [1, exp(-j * omega * (d/c) * cos(vartheta))]^T``.
    Null angle must lie in ``[90, 180]`` deg; when ``equalize=True``, enforces ``w^H a(0) == 1``.
    """
    f_hz = finite_real_scalar(freq_hz, "freq_hz")
    spacing = finite_real_scalar(d_m, "d_m")
    null_deg = finite_real_scalar(null_angle_deg, "null_angle_deg")
    sound_speed = finite_real_scalar(c, "c")
    if not isinstance(equalize, (bool, np.bool_)):
        raise ValueError("equalize must be a bool")
    if f_hz <= 0.0 or spacing <= 0.0 or sound_speed <= 0.0:
        raise ValueError("freq_hz, d_m, and c must be strictly positive")
    if not (90.0 <= null_deg <= 180.0):
        raise ValueError("null_angle_deg must lie in rear hemisphere [90.0, 180.0]")

    alpha = -math.cos(math.radians(null_deg))
    propagation_cycles = f_hz * (spacing / sound_speed)
    propagation_phase = 2.0 * math.pi * propagation_cycles
    if not math.isfinite(propagation_phase) or propagation_phase == 0:
        raise ValueError("DMA propagation phase is outside supported float64 range")
    w_raw = np.array([1.0 + 0.0j, -np.exp(1j * (alpha * propagation_phase))], dtype=np.complex128)
    if not equalize:
        return w_raw

    a_target = np.array([1.0 + 0.0j, np.exp(-1j * propagation_phase)], dtype=np.complex128)
    resp_target = complex(np.vdot(w_raw, a_target))
    if abs(resp_target) <= 1e-12:
        raise ValueError("DMA target response vanishes at spatial aliasing frequency; cannot equalize")
    return w_raw / np.conj(resp_target)


def kronecker_mvdr_weights(
    r_x: np.ndarray,
    r_y: np.ndarray,
    a_x: np.ndarray,
    a_y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve separable planar MVDR weights ``w = w_x (x) w_y`` in O(M_x^3 + M_y^3)."""
    w_x = mvdr_weights(r_x, a_x)
    w_y = mvdr_weights(r_y, a_y)
    if w_x.ndim != 1 or w_y.ndim != 1:
        raise ValueError("kronecker_mvdr_weights expects single-frequency 1-D steering vectors a_x and a_y")
    return np.kron(w_x, w_y), w_x, w_y


def lcmv_weights(
    cov: np.ndarray,
    constraints_c: np.ndarray,
    response_f: np.ndarray,
    *,
    relative_diagonal_loading: float = 0.0,
    condition_limit: float = 1e12,
) -> np.ndarray:
    """Thin wrapper around Chapter 5 ``beamforming.lcmv_weights`` solving ``C^H w = f``."""
    return _base_lcmv_weights(
        cov,
        constraints_c,
        response_f,
        relative_diagonal_loading=relative_diagonal_loading,
        condition_limit=condition_limit,
    )


def first_order_dma_example(
    d_m: float = 0.01,
    c: float = SPEED_OF_SOUND,
    gasket_phase_mismatch_deg: float = 5.0,
) -> dict[str, object]:
    """E05-16: classical first-order DMA patterns, +6 dB/oct equalization WNG drop, and 5-deg gasket leak."""
    spacing = finite_real_scalar(d_m, "d_m")
    sound_speed = finite_real_scalar(c, "c")
    mismatch_deg = finite_real_scalar(gasket_phase_mismatch_deg, "gasket_phase_mismatch_deg")
    if spacing <= 0.0 or sound_speed <= 0.0:
        raise ValueError("d_m and c must be strictly positive")

    patterns = []
    for name, alpha in (
        ("dipole", 0.0),
        ("hypercardioid", 1.0 / 3.0),
        ("supercardioid", 1.0 / math.sqrt(3.0)),
        ("cardioid", 1.0),
    ):
        null_deg = math.degrees(math.acos(-alpha))
        q_3d = ((1.0 + alpha) ** 2) / (alpha * alpha + 1.0 / 3.0)
        patterns.append(
            {
                "pattern": name,
                "alpha": alpha,
                "null_angle_deg": null_deg,
                "directivity_factor_3d": q_3d,
                "directivity_index_3d_db": 10.0 * math.log10(q_3d),
            }
        )

    w_raw_200 = first_order_dma_weights(200.0, spacing, 180.0, c=sound_speed, equalize=False)
    w_eq_200 = first_order_dma_weights(200.0, spacing, 180.0, c=sound_speed, equalize=True)
    w_raw_100 = first_order_dma_weights(100.0, spacing, 180.0, c=sound_speed, equalize=False)
    w_eq_100 = first_order_dma_weights(100.0, spacing, 180.0, c=sound_speed, equalize=True)

    a0_200 = np.array([1.0 + 0.0j, np.exp(-1j * 2.0 * math.pi * 200.0 * spacing / sound_speed)])
    a0_100 = np.array([1.0 + 0.0j, np.exp(-1j * 2.0 * math.pi * 100.0 * spacing / sound_speed)])

    raw_mag_200 = abs(complex(np.vdot(w_raw_200, a0_200)))
    raw_mag_100 = abs(complex(np.vdot(w_raw_100, a0_100)))
    wng_200 = 1.0 / float(np.vdot(w_eq_200, w_eq_200).real)
    wng_100 = 1.0 / float(np.vdot(w_eq_100, w_eq_100).real)

    delta_phi = math.radians(mismatch_deg)
    kd_200 = 2.0 * math.pi * 200.0 * spacing / sound_speed
    a_back_ideal_200 = np.array([1.0 + 0.0j, np.exp(+1j * kd_200)], dtype=np.complex128)
    a_back_mismatched_200 = np.array(
        [1.0 + 0.0j, np.exp(1j * delta_phi) * np.exp(+1j * kd_200)],
        dtype=np.complex128,
    )
    ideal_back_leak = abs(complex(np.vdot(w_eq_200, a_back_ideal_200)))
    mismatched_back_leak = abs(complex(np.vdot(w_eq_200, a_back_mismatched_200)))
    mismatched_rear_gain_db = (
        -math.inf if mismatched_back_leak == 0.0 else 20.0 * math.log10(mismatched_back_leak)
    )

    return {
        "spacing_m": spacing,
        "sound_speed_m_s": sound_speed,
        "cardioid_delay_tau0_us": (spacing / sound_speed) * 1e6,
        "classical_patterns": patterns,
        "cardioid_200hz": {
            "kd_rad": kd_200,
            "raw_target_gain": raw_mag_200,
            "raw_target_gain_db": 20.0 * math.log10(raw_mag_200),
            "equalized_wng": wng_200,
            "equalized_wng_db": 10.0 * math.log10(wng_200),
        },
        "cardioid_100hz": {
            "kd_rad": 0.5 * kd_200,
            "raw_target_gain": raw_mag_100,
            "raw_target_gain_db": 20.0 * math.log10(raw_mag_100),
            "equalized_wng": wng_100,
            "equalized_wng_db": 10.0 * math.log10(wng_100),
        },
        "octave_wng_drop_db": 10.0 * math.log10(wng_100 / wng_200),
        "gasket_mismatch_200hz": {
            "phase_mismatch_deg": mismatch_deg,
            "ideal_rear_null_gain": ideal_back_leak,
            "mismatched_rear_gain": mismatched_back_leak,
            "mismatched_rear_gain_db": mismatched_rear_gain_db,
        },
        "limits": (
            "First-order two-mic endfire DMA model for kd << 1; equalized low-frequency WNG falls at "
            "6.02 dB/octave and requires regularizing high-pass roll-off and strict phase matching."
        ),
    }


def kronecker_beamforming_example() -> dict[str, object]:
    """E05-17: separable Kronecker MVDR exactness vs non-separable point-interferer LCMV null coupling."""
    r_x = np.array([[2.0 + 0.0j, -1.0j], [1.0j, 2.0 + 0.0j]], dtype=np.complex128)
    r_y = np.array([[1.0 + 0.0j, 0.0j], [0.0j, 3.0 + 0.0j]], dtype=np.complex128)
    a_x0 = np.ones(2, dtype=np.complex128)
    a_y0 = np.ones(2, dtype=np.complex128)

    w_kron_mvdr, w_x_mvdr, w_y_mvdr = kronecker_mvdr_weights(r_x, r_y, a_x0, a_y0)
    r_joint = np.kron(r_x, r_y)
    a_0 = np.kron(a_x0, a_y0)
    w_joint_mvdr = mvdr_weights(r_joint, a_0)
    mvdr_diff_norm = float(np.linalg.norm(w_joint_mvdr - w_kron_mvdr))
    mvdr_noise_power = float(np.vdot(w_joint_mvdr, r_joint @ w_joint_mvdr).real)

    a_xi = np.array([1.0 + 0.0j, 1.0j], dtype=np.complex128)
    a_yi = np.array([1.0 + 0.0j, 1.0j], dtype=np.complex128)
    a_interferer = np.kron(a_xi, a_yi)
    a_bystander = np.kron(a_xi, a_y0)

    c_1d = np.column_stack((a_x0, a_xi))
    f_1d = np.array([1.0 + 0.0j, 0.0 + 0.0j], dtype=np.complex128)
    w_1d_null = lcmv_weights(np.eye(2, dtype=np.complex128), c_1d, f_1d)
    w_kron_lcmv = np.kron(w_1d_null, w_1d_null)
    w_kron_xonly = np.kron(w_1d_null, a_y0 / 2.0)

    c_joint = np.column_stack((a_0, a_interferer))
    w_joint_lcmv = lcmv_weights(np.eye(4, dtype=np.complex128), c_joint, f_1d)

    wng_kron_both = 1.0 / float(np.vdot(w_kron_lcmv, w_kron_lcmv).real)
    wng_kron_xonly = 1.0 / float(np.vdot(w_kron_xonly, w_kron_xonly).real)
    wng_joint = 1.0 / float(np.vdot(w_joint_lcmv, w_joint_lcmv).real)

    bystander_power_kron_both = abs(complex(np.vdot(w_kron_lcmv, a_bystander))) ** 2
    bystander_power_kron_xonly = abs(complex(np.vdot(w_kron_xonly, a_bystander))) ** 2
    bystander_power_joint = abs(complex(np.vdot(w_joint_lcmv, a_bystander))) ** 2

    return {
        "separable_mvdr": {
            "w_x_mvdr": w_x_mvdr,
            "w_y_mvdr": w_y_mvdr,
            "w_kron_mvdr": w_kron_mvdr,
            "w_joint_mvdr": w_joint_mvdr,
            "l2_difference": mvdr_diff_norm,
            "output_noise_power": mvdr_noise_power,
            "complexity_ratio_8x8": (64**3) / float(8**3 + 8**3),
        },
        "point_interferer_lcmv": {
            "w_kron_doubly_nulled": w_kron_lcmv,
            "w_joint_lcmv": w_joint_lcmv,
            "wng_kron_doubly_nulled": wng_kron_both,
            "wng_kron_x_only": wng_kron_xonly,
            "wng_joint_lcmv": wng_joint,
            "wng_joint_lcmv_db": 10.0 * math.log10(wng_joint),
            "interferer_response_kron_abs": abs(complex(np.vdot(w_kron_lcmv, a_interferer))),
            "interferer_response_joint_abs": abs(complex(np.vdot(w_joint_lcmv, a_interferer))),
            "bystander_power_kron_doubly_nulled": bystander_power_kron_both,
            "bystander_power_kron_x_only": bystander_power_kron_xonly,
            "bystander_power_joint_lcmv": bystander_power_joint,
        },
        "limits": (
            "Rank-1 Kronecker factorization w = w_x (x) w_y is exact for separable R_x (x) R_y, "
            "but forms cross-shaped line nulls rather than point nulls under localized non-separable interference."
        ),
    }
