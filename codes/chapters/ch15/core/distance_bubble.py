"""Near-field spherical-wave steering and broadside acoustic distance-bubble beamforming.

Conventions:
- Forward Fourier kernel ``exp(-1j * 2*pi*f*t)`` and beamformer output ``y = w^H x``.
- 1-D ULA sensor coordinates ``x_m`` lie along the x-axis; angle ``theta_rad`` is
  measured from the broadside +y normal towards +x in ``[-pi/2, pi/2]``, so a
  point source at range ``r_m`` and angle ``theta_rad`` sits at
  ``(r_m * sin(theta_rad), r_m * cos(theta_rad))``.
"""

from __future__ import annotations

import numpy as np

from codes.chapters.ch02.core.conventions import (
    SPEED_OF_SOUND,
    finite_real_array,
    finite_real_scalar,
    hermitian_part,
)


def fresnel_distance(aperture_m: float, freq_hz: float, c: float = SPEED_OF_SOUND) -> float:
    """Return the Fraunhofer / Fresnel near-field boundary ``r_F = 2 * D^2 / lambda``."""
    d_ap = finite_real_scalar(aperture_m, "aperture_m")
    freq = finite_real_scalar(freq_hz, "freq_hz")
    speed = finite_real_scalar(c, "c")
    if d_ap < 0.0:
        raise ValueError("aperture_m must be non-negative")
    if freq <= 0.0 or speed <= 0.0:
        raise ValueError("freq_hz and c must be strictly positive")
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        wavelength_m = np.float64(speed) / np.float64(freq)
        r_f = 2.0 * (np.float64(d_ap) ** 2) / wavelength_m
    if not (wavelength_m > 0.0 and np.isfinite(wavelength_m) and np.isfinite(r_f)):
        raise ValueError("Fresnel distance exceeds finite float64 range")
    return float(r_f)


def spherical_ula_steering(
    positions_m: np.ndarray,
    r_m: float,
    theta_rad: float,
    freq_hz: float,
    c: float = SPEED_OF_SOUND,
    attenuation: str = "phase_only",
    r_ref_m: float | None = None,
) -> np.ndarray:
    """Evaluate near-field spherical-wave ULA steering vector at ``(r_m, theta_rad)``.

    Parameters
    ----------
    positions_m:
        1-D sensor x-coordinates in meters, shape ``(M,)``.
    r_m:
        Source radial distance from the array origin in meters (``> 0``).
    theta_rad:
        Azimuth angle in radians measured from broadside normal (``0`` = broadside).
    freq_hz:
        Physical frequency in Hz (``>= 0``).
    c:
        Speed of sound in m/s (``> 0``).
    attenuation:
        - ``"phase_only"``: unit-modulus spherical phase relative to ``r_ref_m``;
        - ``"relative_1_over_r"``: ``(r_ref_m / r_sensor) * exp(-1j * k * (r_sensor - r_ref_m))``;
        - ``"absolute_1_over_r"``: ``(1 / r_sensor) * exp(-1j * k * r_sensor)``.
    r_ref_m:
        Reference distance in meters; defaults to ``r_m`` (origin reference).
    """
    pos = finite_real_array(positions_m, "positions_m")
    if pos.ndim != 1 or pos.size < 1:
        raise ValueError("positions_m must be a non-empty 1-D real array")
    r_val = finite_real_scalar(r_m, "r_m")
    theta = finite_real_scalar(theta_rad, "theta_rad")
    freq = finite_real_scalar(freq_hz, "freq_hz")
    speed = finite_real_scalar(c, "c")
    if r_val <= 0.0:
        raise ValueError("r_m must be strictly positive")
    if freq < 0.0 or speed <= 0.0:
        raise ValueError("freq_hz must be non-negative and c must be strictly positive")
    ref_dist = r_val if r_ref_m is None else finite_real_scalar(r_ref_m, "r_ref_m")
    if ref_dist <= 0.0:
        raise ValueError("r_ref_m must be strictly positive")

    src_x = r_val * np.sin(theta)
    src_y = r_val * np.cos(theta)
    sensor_distances = np.hypot(src_x - pos, src_y)
    if np.any(sensor_distances <= 0.0):
        raise ValueError("source position coincides with a sensor coordinate")

    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        wavenumber = 2.0 * np.pi * freq / speed
        if attenuation == "phase_only":
            phase = wavenumber * (sensor_distances - ref_dist)
            if not np.all(np.isfinite(phase)):
                raise ValueError("spherical steering phase exceeds finite float64 range")
            steering = np.exp(-1j * phase)
        elif attenuation == "relative_1_over_r":
            phase = wavenumber * (sensor_distances - ref_dist)
            gain = ref_dist / sensor_distances
            if not (np.all(np.isfinite(phase)) and np.all(np.isfinite(gain))):
                raise ValueError("spherical steering gain or phase exceeds finite float64 range")
            steering = gain * np.exp(-1j * phase)
        elif attenuation == "absolute_1_over_r":
            phase = wavenumber * sensor_distances
            gain = 1.0 / sensor_distances
            if not (np.all(np.isfinite(phase)) and np.all(np.isfinite(gain))):
                raise ValueError("spherical steering gain or phase exceeds finite float64 range")
            steering = gain * np.exp(-1j * phase)
        else:
            raise ValueError(
                "attenuation must be 'phase_only', 'relative_1_over_r', or 'absolute_1_over_r'"
            )
    if not np.all(np.isfinite(steering)):
        raise ValueError("spherical steering vector contains non-finite values")
    return np.asarray(steering, dtype=np.complex128)


def broadside_distance_bubble_weights(
    d_m: float,
    r_target_m: float,
    freq_hz: float,
    c: float = SPEED_OF_SOUND,
    reg: float = 0.0,
    radii_m: np.ndarray | None = None,
    attenuation: str = "phase_only",
) -> dict:
    """Construct 3-mic broadside LCMV distance-bubble weights passing ``r_target_m`` and nulling ``r -> inf``.

    For ``reg == 0.0``, solves ``min ||w||_2^2`` subject to ``w^H a(r_target, 0) = 1``
    and ``w^H [1, 1, 1]^T = 0``. For ``reg > 0.0``, solves the Tikhonov-damped Gram
    system ``w = C (C^H C + reg * I_2)^{-1} [1, 0]^T``.
    """
    spacing = finite_real_scalar(d_m, "d_m")
    r_tar = finite_real_scalar(r_target_m, "r_target_m")
    freq = finite_real_scalar(freq_hz, "freq_hz")
    speed = finite_real_scalar(c, "c")
    reg_val = finite_real_scalar(reg, "reg")
    if spacing <= 0.0 or r_tar <= 0.0 or freq <= 0.0 or speed <= 0.0:
        raise ValueError("d_m, r_target_m, freq_hz, and c must all be strictly positive")
    if reg_val < 0.0:
        raise ValueError("reg must be non-negative")
    if not isinstance(attenuation, str) or attenuation not in ("phase_only", "relative_1_over_r"):
        raise ValueError("broadside_distance_bubble_weights requires unit-origin-normalized steering")

    positions = np.array([-spacing, 0.0, spacing], dtype=np.float64)
    a_near = spherical_ula_steering(positions, r_tar, 0.0, freq, c=speed, attenuation=attenuation)
    a_far = np.ones(3, dtype=np.complex128)
    constraints = np.column_stack((a_near, a_far))
    desired = np.array([1.0 + 0.0j, 0.0 + 0.0j], dtype=np.complex128)

    gram = hermitian_part(constraints.conj().T @ constraints + reg_val * np.eye(2))
    eigvals = np.linalg.eigvalsh(gram)
    if eigvals[0] <= 1e-15 * max(1.0, float(eigvals[-1])):
        raise np.linalg.LinAlgError(
            "near-field and far-field steering vectors are numerically collinear; "
            "increase aperture/frequency or supply positive reg"
        )
    alpha = np.linalg.solve(gram, desired)
    weights = constraints @ alpha
    if not np.all(np.isfinite(weights)):
        raise np.linalg.LinAlgError("distance-bubble weights contain non-finite values")

    norm_w2 = float(np.vdot(weights, weights).real)
    target_resp = complex(np.vdot(weights, a_near))
    far_plane_resp = complex(np.vdot(weights, a_far))
    target_power_gain = float(abs(target_resp) ** 2)
    far_plane_power_gain = float(abs(far_plane_resp) ** 2)
    if norm_w2 <= 0.0 or target_power_gain <= 0.0 or not (np.isfinite(norm_w2) and np.isfinite(target_power_gain)):
        raise np.linalg.LinAlgError("degenerate weight norm or target power gain in distance bubble")
    wng_linear = target_power_gain / norm_w2
    wng_db = float(10.0 * np.log10(wng_linear))
    if not np.isfinite(wng_db):
        raise ValueError("WNG in dB exceeds finite float64 range")

    if radii_m is None:
        grid = np.linspace(0.10, 5.00, 246, dtype=np.float64)
    else:
        grid = finite_real_array(radii_m, "radii_m")
        if grid.ndim != 1 or grid.size < 1 or np.any(grid <= 0.0):
            raise ValueError("radii_m must be a non-empty 1-D array of positive distances")

    radial_power = np.empty(grid.size, dtype=np.float64)
    for idx, radius in enumerate(grid):
        a_r = spherical_ula_steering(positions, float(radius), 0.0, freq, c=speed, attenuation=attenuation)
        radial_power[idx] = float(abs(np.vdot(weights, a_r)) ** 2)
    if not np.all(np.isfinite(radial_power)):
        raise ValueError("radial power response exceeds finite float64 range")

    delta_r_m = float(np.hypot(r_tar, spacing) - r_tar)
    wavelength_m = speed / freq
    outer_phase_rad = float(-2.0 * np.pi * delta_r_m / wavelength_m)
    r_fresnel_m = fresnel_distance(2.0 * spacing, freq, c=speed)

    return {
        "positions_m": positions,
        "aperture_m": float(2.0 * spacing),
        "wavelength_m": float(wavelength_m),
        "fresnel_distance_m": r_fresnel_m,
        "extra_path_m": delta_r_m,
        "outer_phase_rad": outer_phase_rad,
        "outer_phase_deg": float(np.degrees(outer_phase_rad)),
        "weights": weights,
        "weight_norm_squared": norm_w2,
        "target_response": target_resp,
        "target_power_gain": target_power_gain,
        "farfield_plane_power_gain": far_plane_power_gain,
        "wng_linear": wng_linear,
        "wng_db": wng_db,
        "radii_m": grid,
        "radial_power_response": radial_power,
        "gram_condition_number": float(eigvals[-1] / eigvals[0]),
    }


def broadside_two_tap_real_fir(
    weights: np.ndarray,
    freq_hz: float,
    sample_rate_hz: float = 16000.0,
) -> np.ndarray:
    """Realize complex narrowband weights ``y = w^H x`` at interior frequency ``freq_hz`` as 2-tap causal real FIRs.

    Returns array ``taps`` of shape ``(M, 2)`` such that
    ``taps[m, 0] + taps[m, 1] * exp(-1j * 2*pi*freq_hz/sample_rate_hz) == conj(weights[m])``.
    """
    raw = np.asarray(weights)
    if raw.dtype.kind not in "iufc" or raw.ndim != 1 or raw.size < 1:
        raise ValueError("weights must be a non-empty 1-D numeric array")
    w = np.asarray(raw, dtype=np.complex128)
    if not np.all(np.isfinite(w)):
        raise ValueError("weights must be finite")
    freq = finite_real_scalar(freq_hz, "freq_hz")
    fs = finite_real_scalar(sample_rate_hz, "sample_rate_hz")
    if fs <= 0.0 or not (0.0 < freq < 0.5 * fs):
        raise ValueError("freq_hz must lie strictly inside (0, sample_rate_hz / 2)")

    omega = 2.0 * np.pi * freq / fs
    sin_om = float(np.sin(omega))
    cos_om = float(np.cos(omega))
    if abs(sin_om) <= 1e-12:
        raise ValueError("freq_hz is too close to 0 or Nyquist for 2-tap real FIR realization")
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        h1 = w.imag / sin_om
        h0 = w.real - h1 * cos_om
        taps = np.column_stack((h0, h1)).astype(np.float64)
    if not np.all(np.isfinite(taps)):
        raise ValueError("2-tap FIR coefficients exceed finite float64 range")
    return taps
