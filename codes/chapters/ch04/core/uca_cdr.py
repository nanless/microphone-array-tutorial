"""E04-26 & E04-27: UCA phase-mode ESPRIT and blind two-mic coherent-to-diffuse ratio (CDR).

Pure-NumPy deterministic implementations without SciPy or import-time side effects.
References:
- Mathews & Zoltowski, IEEE Trans. Signal Process. 42(9), 1994 (UCA-ESPRIT).
- Schwarz & Habets, Signal Processing 117, 2015, Eq. (25) (DOA-independent CDR).
"""

from __future__ import annotations

import math

import numpy as np

from codes.chapters.ch02.core.conventions import SPEED_OF_SOUND, finite_real_array, finite_real_scalar
from codes.chapters.ch04.core.covariance import _validate_covariance


FIRST_BESSEL_J0_ROOT = 2.404825557695773
FIRST_BESSEL_J1_ROOT = 3.8317059702075125


def bessel_jn(n: int, x: float | np.ndarray, max_terms: int = 28) -> np.ndarray:
    """Evaluate integer-order cylindrical Bessel function J_n(x) in pure NumPy.

    Uses exact alternating power series for ``|x| <= 12`` and Jacobi-Anger Fourier
    integral quadrature for ``|x| > 12`` to avoid alternating-series cancellation.
    """
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)):
        raise ValueError("order n must be an integer, not bool")
    if isinstance(max_terms, (bool, np.bool_)) or not isinstance(max_terms, (int, np.integer)) or max_terms < 10:
        raise ValueError("max_terms must be an integer >= 10")

    order = int(n)
    if order < 0:
        pos = bessel_jn(-order, x, max_terms=max_terms)
        return (-pos) if ((-order) % 2 == 1) else pos

    arr = finite_real_array(x, "x")
    out = np.empty_like(arr, dtype=np.float64)
    abs_x = np.abs(arr)

    mask_series = abs_x <= 12.0
    if np.any(mask_series):
        xs = arr[mask_series]
        half_x = 0.5 * xs
        # Compute (x/2)^n / n! iteratively without integer factorial overflow for large n.
        term = np.ones_like(xs, dtype=np.float64)
        for m in range(1, order + 1):
            term = term * (half_x / float(m))
        acc = term.copy()
        quarter_x2 = half_x * half_x
        terms_limit = max(int(max_terms), order + 20)
        for k in range(1, terms_limit + 1):
            term = -term * quarter_x2 / float(k * (order + k))
            acc = acc + term
            if float(np.max(np.abs(term))) <= 1e-17 * max(1.0, float(np.max(np.abs(acc)))):
                break
        out[mask_series] = acc

    mask_quad = ~mask_series
    if np.any(mask_quad):
        xq = arr[mask_quad]
        # Integral identity J_n(x) = (1/pi) * int_0^pi cos(n*theta - x*sin(theta)) d_theta.
        max_val = float(np.max(np.abs(xq)))
        n_quad = max(128, int(math.ceil(2.0 * max_val + 2 * order + 40)))
        theta = (np.arange(n_quad, dtype=np.float64) + 0.5) * (math.pi / float(n_quad))
        phase = float(order) * theta - xq[..., np.newaxis] * np.sin(theta)
        out[mask_quad] = np.mean(np.cos(phase), axis=-1)

    return out


def uca_steering_vector(num_mics: int, kappa_r: float, azimuth_deg: float) -> np.ndarray:
    """Return horizontal plane-wave UCA steering vector exp(+j * kappa_r * cos(phi - psi_q))."""
    if isinstance(num_mics, (bool, np.bool_)) or not isinstance(num_mics, (int, np.integer)) or num_mics < 3:
        raise ValueError("num_mics must be an integer >= 3")
    kr = finite_real_scalar(kappa_r, "kappa_r")
    az = finite_real_scalar(azimuth_deg, "azimuth_deg")
    if kr <= 0.0:
        raise ValueError("kappa_r must be strictly positive")
    phi = 0.5 * math.pi - math.radians(az)
    psi = 2.0 * np.pi * np.arange(int(num_mics), dtype=np.float64) / float(num_mics)
    return np.exp(1j * kr * np.cos(phi - psi))


def uca_phase_mode_esprit(
    cov: np.ndarray,
    kappa_r: float,
    n_sources: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Estimate horizontal azimuths (deg) via 3-mode (-1, 0, +1) UCA phase-mode ESPRIT.

    Returns ``(azimuth_estimates_deg, mode_covariance_3x3)`` where azimuths follow the
    book convention (0 deg broadside +y, clockwise toward +x, wrapped to (-180, 180]).
    """
    matrix = _validate_covariance(cov)
    m_mics = matrix.shape[0]
    if m_mics < 3:
        raise ValueError("UCA phase-mode ESPRIT requires at least 3 microphones")
    if isinstance(n_sources, (bool, np.bool_)) or not isinstance(n_sources, (int, np.integer)) or not (1 <= n_sources <= 2):
        raise ValueError("3-mode UCA-ESPRIT supports integer n_sources in {1, 2}")
    kr = finite_real_scalar(kappa_r, "kappa_r")
    if kr <= 0.0:
        raise ValueError("kappa_r must be strictly positive")

    j0 = float(bessel_jn(0, kr))
    j1 = float(bessel_jn(1, kr))
    jm1 = float(bessel_jn(-1, kr))
    if abs(j0) < 1e-10 or abs(j1) < 1e-10:
        raise ValueError("Bessel zero singularity at kappa_r")

    psi = 2.0 * np.pi * np.arange(m_mics, dtype=np.float64) / float(m_mics)
    modes = np.array([-1, 0, 1], dtype=np.int64)
    f_h = np.exp(-1j * modes[:, np.newaxis] * psi[np.newaxis, :]) / float(m_mics)
    w_radial = np.array([1.0 / ((-1j) * jm1), 1.0 / j0, 1.0 / (1j * j1)], dtype=np.complex128)
    beamspace = w_radial[:, np.newaxis] * f_h

    r_mode = beamspace @ matrix @ beamspace.conj().T
    r_mode = 0.5 * (r_mode + r_mode.conj().T)
    peak = float(np.max(np.maximum(np.abs(r_mode.real), np.abs(r_mode.imag))))
    if peak == 0.0:
        raise np.linalg.LinAlgError("zero phase-mode covariance has no signal subspace")

    evals, evecs = np.linalg.eigh(r_mode / peak)
    k_src = int(n_sources)
    tol = 100.0 * np.finfo(float).eps * 3.0 * max(1.0, float(np.max(np.abs(evals))))
    if evals[-k_src] <= tol or evals[-k_src] - evals[-k_src - 1] <= tol:
        raise np.linalg.LinAlgError("phase-mode signal subspace is not separated from noise")

    u_s = evecs[:, -k_src:]
    u_upper, u_lower = u_s[:-1], u_s[1:]
    if np.linalg.matrix_rank(u_upper) < k_src or np.linalg.matrix_rank(u_lower) < k_src:
        raise np.linalg.LinAlgError("phase-mode ESPRIT subarray is rank deficient")

    psi_op, *_ = np.linalg.lstsq(u_upper, u_lower, rcond=None)
    eig_mu = np.linalg.eigvals(psi_op)
    if np.any(np.abs(eig_mu) <= 1e-12):
        raise np.linalg.LinAlgError("degenerate zero eigenvalue in phase-mode ESPRIT")

    phi_hat = -np.angle(eig_mu)
    raw_deg = np.rad2deg(0.5 * np.pi - phi_hat)
    wrapped = ((raw_deg + 180.0) % 360.0) - 180.0
    wrapped = np.where((wrapped == -180.0) & (raw_deg > 0.0), 180.0, wrapped)
    return np.sort(wrapped.astype(np.float64)), r_mode


def cdr_nodir(
    gamma_x: complex | np.ndarray,
    gamma_diff: float | np.ndarray,
    eps: float = 1e-12,
) -> np.ndarray:
    """Estimate DOA-independent coherent-to-diffuse ratio (Schwarz & Habets 2015, Eq. 25)."""
    eps_val = finite_real_scalar(eps, "eps")
    if eps_val <= 0.0 or eps_val >= 0.1:
        raise ValueError("eps must lie in (0, 0.1)")

    gx_raw = np.asarray(gamma_x)
    if gx_raw.dtype.kind not in "iufc":
        raise ValueError("gamma_x must be numeric (real or complex), not bool or text")
    with np.errstate(over="ignore", invalid="ignore"):
        gx = np.asarray(gx_raw, dtype=np.complex128)
    if not np.all(np.isfinite(gx)):
        raise ValueError("gamma_x must be finite")

    gd = finite_real_array(gamma_diff, "gamma_diff")
    gx_b, gd_b = np.broadcast_arrays(gx, gd)
    abs_gx2 = gx_b.real * gx_b.real + gx_b.imag * gx_b.imag
    if np.any(abs_gx2 > 1.0 + 1e-9):
        raise ValueError("observed coherence magnitude |gamma_x| cannot exceed 1")
    if np.any(np.abs(gd_b) >= 1.0 - eps_val):
        raise ValueError("diffuse coherence |gamma_diff| must be strictly < 1 (degenerate at DC)")

    abs_gx2_clamped = np.minimum(abs_gx2, 1.0)
    diff = gx_b - gd_b
    a = diff.real * diff.real + diff.imag * diff.imag
    b = abs_gx2_clamped - gd_b * gx_b.real
    c0 = abs_gx2_clamped - 1.0
    disc = np.maximum(b * b - a * c0, 0.0)
    num = b + np.sqrt(disc)
    denom = 1.0 - abs_gx2_clamped

    out = np.empty(gx_b.shape, dtype=np.float64)
    pure_diffuse = a <= eps_val * eps_val
    pure_coherent = (~pure_diffuse) & (denom <= eps_val)
    regular = ~(pure_diffuse | pure_coherent)

    out[pure_diffuse] = 0.0
    out[pure_coherent] = np.inf
    out[regular] = np.maximum(num[regular] / denom[regular], 0.0)
    return out


def uca_phase_mode_example(
    radius_m: float = 0.04,
    kappa_r: float = 1.0,
    true_azimuth_deg: float = 30.0,
    sound_speed: float = SPEED_OF_SOUND,
) -> dict[str, object]:
    """E04-26: UCA circular-harmonics Vandermonde transform, M=8 vs M=4 aliasing, and Bessel zeros."""
    r_m = finite_real_scalar(radius_m, "radius_m")
    kr = finite_real_scalar(kappa_r, "kappa_r")
    az = finite_real_scalar(true_azimuth_deg, "true_azimuth_deg")
    c_m_s = finite_real_scalar(sound_speed, "sound_speed")
    if r_m <= 0.0 or kr <= 0.0 or c_m_s <= 0.0:
        raise ValueError("radius_m, kappa_r, and sound_speed must be strictly positive")

    j0_val = float(bessel_jn(0, kr))
    j1_val = float(bessel_jn(1, kr))
    j3_val = float(bessel_jn(3, kr))
    j7_val = float(bessel_jn(7, kr))
    freq_hz = c_m_s * kr / (2.0 * math.pi * r_m)

    a_m8 = uca_steering_vector(8, kr, az)
    cov_m8 = np.outer(a_m8, a_m8.conj())
    est_m8, _ = uca_phase_mode_esprit(cov_m8, kr, n_sources=1)

    a_m4 = uca_steering_vector(4, kr, az)
    cov_m4 = np.outer(a_m4, a_m4.conj())
    est_m4, _ = uca_phase_mode_esprit(cov_m4, kr, n_sources=1)

    phi = 0.5 * math.pi - math.radians(az)
    psi8 = 2.0 * np.pi * np.arange(8, dtype=np.float64) / 8.0
    modes = np.array([-1, 0, 1], dtype=np.int64)
    w_rad = np.array([1.0 / (1j * j1_val), 1.0 / j0_val, 1.0 / (1j * j1_val)], dtype=np.complex128)
    tilde_a_m8 = w_rad * ((np.exp(-1j * modes[:, np.newaxis] * psi8[np.newaxis, :]) / 8.0) @ a_m8)
    ideal_vandermonde = np.array([np.exp(1j * phi), 1.0 + 0.0j, np.exp(-1j * phi)], dtype=np.complex128)
    vandermonde_err_m8 = np.abs(tilde_a_m8 - ideal_vandermonde)

    f_zero_m0 = c_m_s * FIRST_BESSEL_J0_ROOT / (2.0 * math.pi * r_m)
    f_zero_m1 = c_m_s * FIRST_BESSEL_J1_ROOT / (2.0 * math.pi * r_m)

    try:
        uca_phase_mode_esprit(cov_m8, FIRST_BESSEL_J0_ROOT, n_sources=1)
    except ValueError as err:
        singularity_rejected, singularity_reason = True, str(err)
    else:
        raise AssertionError("Expected Bessel zero singularity rejection at z_{0,1}")

    return {
        "radius_m": r_m,
        "kappa_r": kr,
        "operating_frequency_hz": freq_hz,
        "true_azimuth_deg": az,
        "j0_1": j0_val,
        "j1_1": j1_val,
        "j3_1_alias": j3_val,
        "j7_1_alias": j7_val,
        "radial_equalizer_m0_real": 1.0 / j0_val,
        "radial_equalizer_m1_imag": -1.0 / j1_val,
        "vandermonde_abs_error_m8": vandermonde_err_m8.tolist(),
        "theta_hat_m8_deg": float(est_m8[0]),
        "abs_bias_m8_deg": abs(float(est_m8[0]) - az),
        "theta_hat_m4_deg": float(est_m4[0]),
        "abs_bias_m4_deg": abs(float(est_m4[0]) - az),
        "f_zero_m0_hz": f_zero_m0,
        "f_zero_m1_hz": f_zero_m1,
        "bessel_singularity_rejected": singularity_rejected,
        "bessel_singularity_reason": singularity_reason,
        "limits": (
            "Open-array 3-mode UCA-ESPRIT requires M > 2*kappa*r + 2 to suppress circular aliasing "
            "and fails at open-sphere Bessel zeros where rigid spherical baffles remain non-singular."
        ),
    }


def two_mic_cdr_example() -> dict[str, object]:
    """E04-27: exact rational hand-calculation and boundary regimes for blind two-mic CDR_nodir."""
    gamma_n = 0.5
    true_cdr = 2.0
    gamma_s = 1j
    gamma_x = (true_cdr * gamma_s + gamma_n) / (true_cdr + 1.0)

    cdr_est = float(cdr_nodir(gamma_x, gamma_n))
    tf_weight = cdr_est / (cdr_est + 1.0)

    cdr_pure_diffuse = float(cdr_nodir(gamma_n + 0.0j, gamma_n))
    cdr_pure_anechoic_raw = float(cdr_nodir(1j, gamma_n))
    anechoic_is_infinite = math.isinf(cdr_pure_anechoic_raw)
    cdr_high_coh = float(cdr_nodir(0.999 * 1j, gamma_n))

    try:
        cdr_nodir(1.0 + 0.0j, 1.0)
    except ValueError as err:
        dc_rejected, dc_reason = True, str(err)
    else:
        raise AssertionError("Expected DC diffuse degeneracy rejection at gamma_diff = 1.0")

    return {
        "gamma_diff": float(gamma_n),
        "true_cdr": float(true_cdr),
        "observed_gamma_x_real": float(gamma_x.real),
        "observed_gamma_x_imag": float(gamma_x.imag),
        "observed_abs_gamma_x_sq": float(abs(gamma_x) ** 2),
        "cdr_est": cdr_est,
        "tf_weight": tf_weight,
        "cdr_pure_diffuse": cdr_pure_diffuse,
        "tf_weight_pure_diffuse": 0.0,
        "cdr_pure_anechoic_is_infinite": anechoic_is_infinite,
        "tf_weight_pure_anechoic": 1.0 if anechoic_is_infinite else (cdr_pure_anechoic_raw / (cdr_pure_anechoic_raw + 1.0)),
        "cdr_high_coherence_0p999": cdr_high_coh,
        "dc_degeneracy_rejected": dc_rejected,
        "dc_degeneracy_reason": dc_reason,
        "limits": (
            "Assumes isotropic 3D real diffuse coherence sinc(2*f*d/c) and unit-modulus target coherence; "
            "finite-frame sample coherence variance requires temporal recursive smoothing in practice."
        ),
    }
