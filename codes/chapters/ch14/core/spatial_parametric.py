"""FOA encoding, 3-D rotation, active acoustic intensity, and diffuseness.

Conventions follow ``codes/chapters/ch02/core/conventions.py``:
- Forward Fourier transform kernel ``exp(-1j * 2*pi*f*t)`` (time harmonic ``exp(+1j*w*t)``).
- Standard 3-D Ambisonics right-handed Cartesian axes: ``+x`` front, ``+y`` left,
  ``+z`` up. Spherical angles ``(azimuth, elevation)`` satisfy
  ``n_src = [cos(el)*cos(az), cos(el)*sin(az), sin(el)]`` pointing towards the source.
- Under Euler's momentum equation ``j*w*rho0*u = -grad(p)``, a plane wave arriving
  from ``n_src`` propagates along ``d_prop = -n_src`` and has physical particle
  velocity ``u = -(p / (rho0 * c)) * n_src = +(p / (rho0 * c)) * d_prop``.
"""
from __future__ import annotations

import numpy as np

from codes.chapters.ch02.core.conventions import (
    finite_real_array,
    finite_real_scalar,
)

# Regular tetrahedron capsule look directions: FLU, FRD, BLD, BRU.
TETRA_DIRECTIONS = np.array(
    [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ],
    dtype=np.float64,
) / np.sqrt(3.0)


def _validate_complex_or_real(value, name: str = "input") -> np.ndarray:
    original = np.asarray(value)
    if original.dtype.kind not in "iufc":
        raise ValueError(f"{name} must contain numeric values")
    target_dtype = np.complex128 if original.dtype.kind == "c" else np.float64
    with np.errstate(over="ignore", invalid="ignore"):
        array = np.asarray(original, dtype=target_dtype)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _validate_order(order: str) -> str:
    if not isinstance(order, str) or order not in ("WXYZ", "ACN"):
        raise ValueError("order must be 'WXYZ' or 'ACN'")
    return order


def reorder_foa(b_format: np.ndarray, *, src_order: str = "WXYZ", dst_order: str = "ACN") -> np.ndarray:
    """Convert 4-channel FOA between ``[W, X, Y, Z]`` and ACN ``[W, Y, Z, X]``."""
    _validate_order(src_order)
    _validate_order(dst_order)
    b = _validate_complex_or_real(b_format, "b_format")
    if b.ndim not in (1, 2) or b.shape[0] != 4:
        raise ValueError("b_format must have shape (4,) or (4, N)")
    if src_order == dst_order:
        return b.copy()
    if src_order == "WXYZ" and dst_order == "ACN":
        return b[[0, 2, 3, 1]].copy()
    return b[[0, 3, 1, 2]].copy()


def tetrahedral_a_to_b_matrix(*, cardioid_alpha: float = 0.5, order: str = "WXYZ") -> np.ndarray:
    """Return exact 4x4 SN3D FOA encoding matrix for coincident tetrahedral capsules.

    Each capsule q in (FLU, FRD, BLD, BRU) has response
    ``a_q = ((1 - alpha) + alpha * (u_q . n_src)) * p``.
    For standard cardioid capsules (``cardioid_alpha = 0.5``), ``W = p`` and
    ``[X, Y, Z] = p * [n_x, n_y, n_z]`` under SN3D normalization.
    """
    _validate_order(order)
    alpha = finite_real_scalar(cardioid_alpha, "cardioid_alpha")
    if not (0.0 < alpha < 1.0):
        raise ValueError("cardioid_alpha must lie strictly in (0, 1)")
    g_w = 1.0 / (4.0 * (1.0 - alpha))
    g_xyz = np.sqrt(3.0) / (4.0 * alpha)
    mat_wxyz = np.array(
        [
            [g_w, g_w, g_w, g_w],
            [g_xyz, g_xyz, -g_xyz, -g_xyz],
            [g_xyz, -g_xyz, g_xyz, -g_xyz],
            [g_xyz, -g_xyz, -g_xyz, g_xyz],
        ],
        dtype=np.float64,
    )
    if order == "ACN":
        return mat_wxyz[[0, 2, 3, 1]]
    return mat_wxyz


def encode_tetrahedral_a_to_b(
    a_format: np.ndarray,
    *,
    cardioid_alpha: float = 0.5,
    order: str = "WXYZ",
) -> np.ndarray:
    """Encode 4-channel tetrahedral A-format (FLU, FRD, BLD, BRU) to SN3D FOA."""
    a = _validate_complex_or_real(a_format, "a_format")
    if a.ndim not in (1, 2) or a.shape[0] != 4:
        raise ValueError("a_format must have shape (4,) or (4, N)")
    if a.ndim == 2 and a.shape[1] < 1:
        raise ValueError("a_format snapshot dimension must be non-empty")
    encoder = tetrahedral_a_to_b_matrix(cardioid_alpha=cardioid_alpha, order=order)
    with np.errstate(over="raise", invalid="raise"):
        try:
            out = encoder @ a
        except FloatingPointError as exc:
            raise ValueError("FOA encoding overflowed float64 range") from exc
    if not np.all(np.isfinite(out)):
        raise ValueError("FOA encoding produced non-finite values")
    return out


def synthesize_tetrahedral_a_format(
    azimuth_rad: float,
    elevation_rad: float = 0.0,
    amplitude: complex = 1.0,
    *,
    cardioid_alpha: float = 0.5,
) -> np.ndarray:
    """Synthesize coincident tetrahedral A-format response for a single plane wave."""
    az = finite_real_scalar(azimuth_rad, "azimuth_rad")
    el = finite_real_scalar(elevation_rad, "elevation_rad")
    alpha = finite_real_scalar(cardioid_alpha, "cardioid_alpha")
    if not (0.0 < alpha < 1.0):
        raise ValueError("cardioid_alpha must lie strictly in (0, 1)")
    amp_arr = _validate_complex_or_real(amplitude, "amplitude")
    if amp_arr.ndim != 0:
        raise ValueError("amplitude must be a scalar")
    n_src = np.array(
        [np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)],
        dtype=np.float64,
    )
    gains = (1.0 - alpha) + alpha * (TETRA_DIRECTIONS @ n_src)
    with np.errstate(over="raise", invalid="raise"):
        try:
            out = gains * amp_arr.item()
        except FloatingPointError as exc:
            raise ValueError("A-format synthesis overflowed float64") from exc
    if not np.all(np.isfinite(out)):
        raise ValueError("A-format synthesis produced non-finite values")
    return out


def open_sphere_foa_radial_weights(kr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return open-sphere FOA modal weights ``(b0(kr), b1(kr))`` for ``kr >= 0``.

    For omnidirectional capsules on a sphere of radius R, the 0th-order mode scales
    as ``b0(kr) = j0(kr) = sin(kr)/kr`` (vanishing at ``kr = pi``!), and the 1st-order
    dipole mode scales as ``b1(kr) = j * j1(kr) = j * (sin(kr) - kr*cos(kr)) / kr^2``.
    """
    x = np.atleast_1d(finite_real_array(kr, "kr"))
    if x.ndim != 1 or x.size < 1 or np.any(x < 0.0):
        raise ValueError("kr must be a non-empty 1-D non-negative real array")
    b0 = np.sinc(x / np.pi)
    b1 = np.empty(x.shape, dtype=np.complex128)
    small = x < 1e-4
    xs = x[small]
    # Taylor series: j1(x) = x/3 - x^3/30 + x^5/840
    b1[small] = 1j * (xs / 3.0 - (xs**3) / 30.0 + (xs**5) / 840.0)
    mid = (~small) & (x <= 1e100)
    xm = x[mid]
    b1[mid] = 1j * (np.sin(xm) - xm * np.cos(xm)) / (xm**2)
    huge = x > 1e100
    xh = x[huge]
    b1[huge] = 1j * (np.sin(xh) / xh - np.cos(xh)) / xh
    if not (np.all(np.isfinite(b0)) and np.all(np.isfinite(b1))):
        raise ValueError("radial modal weights must be finite")
    return b0, b1


def foa_rotation_matrix(
    yaw_rad: float = 0.0,
    pitch_rad: float = 0.0,
    roll_rad: float = 0.0,
    *,
    listener_head_rotation: bool = False,
    order: str = "WXYZ",
) -> np.ndarray:
    """Return 4x4 orthogonal FOA rotation matrix for yaw (Z), pitch (Y), roll (X).

    When ``listener_head_rotation=False`` (default), applies active scene rotation
    ``R_scene = Rz(yaw) @ Ry(pitch) @ Rx(roll)`` to source directions.
    When ``listener_head_rotation=True``, transforms world FOA into a listener-fixed
    frame whose head rotated by ``(yaw, pitch, roll)``, applying ``R_scene.T``.
    """
    if type(listener_head_rotation) is not bool:
        raise ValueError("listener_head_rotation must be bool")
    _validate_order(order)
    yaw = finite_real_scalar(yaw_rad, "yaw_rad")
    pitch = finite_real_scalar(pitch_rad, "pitch_rad")
    roll = finite_real_scalar(roll_rad, "roll_rad")

    cy, sy = np.cos(yaw), np.sin(yaw)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cr, sr = np.cos(roll), np.sin(roll)

    rz = np.array([[cy, -sy, 0.0], [sy, cy, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)
    ry = np.array([[cp, 0.0, sp], [0.0, 1.0, 0.0], [-sp, 0.0, cp]], dtype=np.float64)
    rx = np.array([[1.0, 0.0, 0.0], [0.0, cr, -sr], [0.0, sr, cr]], dtype=np.float64)

    r3 = rz @ ry @ rx
    if listener_head_rotation:
        r3 = r3.T

    rot_wxyz = np.eye(4, dtype=np.float64)
    rot_wxyz[1:, 1:] = r3
    if order == "ACN":
        perm = [0, 2, 3, 1]
        return rot_wxyz[np.ix_(perm, perm)]
    return rot_wxyz


def rotate_foa(
    b_format: np.ndarray,
    yaw_rad: float = 0.0,
    pitch_rad: float = 0.0,
    roll_rad: float = 0.0,
    *,
    listener_head_rotation: bool = False,
    order: str = "WXYZ",
) -> np.ndarray:
    """Apply 3-D yaw-pitch-roll rotation to FOA ``[W, X, Y, Z]`` (or ACN)."""
    b = _validate_complex_or_real(b_format, "b_format")
    if b.ndim not in (1, 2) or b.shape[0] != 4:
        raise ValueError("b_format must have shape (4,) or (4, N)")
    if b.ndim == 2 and b.shape[1] < 1:
        raise ValueError("b_format snapshot dimension must be non-empty")
    rot = foa_rotation_matrix(
        yaw_rad,
        pitch_rad,
        roll_rad,
        listener_head_rotation=listener_head_rotation,
        order=order,
    )
    with np.errstate(over="raise", invalid="raise"):
        try:
            out = rot @ b
        except FloatingPointError as exc:
            raise ValueError("FOA rotation overflowed float64 range") from exc
    if not np.all(np.isfinite(out)):
        raise ValueError("FOA rotation produced non-finite values")
    return out


def foa_to_pressure_velocity(
    b_format: np.ndarray,
    *,
    rho0: float = 1.0,
    c: float = 1.0,
    order: str = "WXYZ",
) -> tuple[np.ndarray, np.ndarray]:
    """Convert SN3D FOA to acoustic pressure ``p`` and physical velocity ``u``.

    Since ``[X, Y, Z] = p * n_src`` points towards the source while a plane wave
    propagates along ``d_prop = -n_src``, physical particle velocity is
    ``u = -[X, Y, Z] / (rho0 * c)``.
    """
    rho0_val = finite_real_scalar(rho0, "rho0")
    c_val = finite_real_scalar(c, "c")
    if rho0_val <= 0.0 or c_val <= 0.0:
        raise ValueError("rho0 and c must be strictly positive")
    z0 = rho0_val * c_val
    if not np.isfinite(z0) or z0 <= 0.0:
        raise ValueError("acoustic impedance rho0 * c must be finite and positive")
    b_wxyz = reorder_foa(b_format, src_order=order, dst_order="WXYZ")
    p = b_wxyz[0]
    with np.errstate(over="raise", invalid="raise"):
        try:
            u = -b_wxyz[1:4] / z0
        except FloatingPointError as exc:
            raise ValueError("particle velocity overflowed float64") from exc
    if not np.all(np.isfinite(u)):
        raise ValueError("particle velocity must be finite")
    return p, u


def active_intensity_and_diffuseness(
    p: np.ndarray | complex | float,
    u: np.ndarray,
    rho0: float = 1.0,
    c: float = 1.0,
) -> dict:
    """Compute active acoustic intensity, total energy density, diffuseness, and DOA.

    Parameters
    ----------
    p : complex scalar or 1-D array of shape ``(N,)``
        Complex acoustic pressure phasor(s).
    u : complex/real array of shape ``(3,)`` or ``(3, N)``
        3-D particle velocity phasor(s) matching ``p``.
    rho0 : float, default 1.0
        Ambient medium density (kg/m^3, or normalized 1.0).
    c : float, default 1.0
        Speed of sound (m/s, or normalized 1.0).

    Returns
    -------
    dict
        Contains ``active_intensity`` (``E[I_a]``), ``reactive_intensity`` (``E[I_r]``),
        ``energy_density`` (``E[E]``), ``diffuseness`` in ``[0, 1]``,
        ``intensity_unit`` (propagation unit vector ``+E[I_a]/||E[I_a]||``),
        ``doa_unit`` (source unit vector ``-E[I_a]/||E[I_a]||``), and ``valid_direction``.
    """
    rho0_val = finite_real_scalar(rho0, "rho0")
    c_val = finite_real_scalar(c, "c")
    if rho0_val <= 0.0 or c_val <= 0.0:
        raise ValueError("rho0 and c must be strictly positive")

    p_arr = _validate_complex_or_real(p, "p")
    u_arr = _validate_complex_or_real(u, "u")

    if p_arr.ndim == 0:
        if u_arr.shape != (3,):
            raise ValueError("for scalar p, u must have shape (3,)")
        p_snap = p_arr.reshape(1)
        u_snap = u_arr.reshape(3, 1)
    elif p_arr.ndim == 1:
        if p_arr.size < 1:
            raise ValueError("p must be non-empty")
        if u_arr.shape != (3, p_arr.size):
            raise ValueError("for 1-D p of length N, u must have shape (3, N)")
        p_snap = p_arr
        u_snap = u_arr
    else:
        raise ValueError("p must be a scalar or 1-D snapshot array")

    # Scale jointly to prevent intermediate square overflow/underflow on extreme finite inputs.
    scale_p = float(np.max(np.abs(p_snap)))
    scale_u = float(np.max(np.abs(u_snap)))
    z0 = rho0_val * c_val
    if not np.isfinite(z0) or z0 <= 0.0:
        raise ValueError("acoustic impedance rho0 * c must be finite and positive")
    scale = max(scale_p, scale_u * z0)

    if scale == 0.0:
        zeros3 = np.zeros(3, dtype=np.float64)
        return {
            "active_intensity": zeros3,
            "reactive_intensity": zeros3.copy(),
            "energy_density": 0.0,
            "diffuseness": 1.0,
            "intensity_unit": zeros3.copy(),
            "doa_unit": zeros3.copy(),
            "valid_direction": False,
        }

    p_norm = p_snap / scale
    u_norm = u_snap * (z0 / scale)

    cross_norm = np.conj(p_norm)[np.newaxis, :] * u_norm
    ia_norm = 0.5 * np.mean(cross_norm.real, axis=1)
    ir_norm = 0.5 * np.mean(cross_norm.imag, axis=1)
    e_pot_norm = 0.25 * float(np.mean(np.abs(p_norm) ** 2))
    e_kin_norm = 0.25 * float(np.mean(np.sum(np.abs(u_norm) ** 2, axis=0)))
    e_tot_norm = e_pot_norm + e_kin_norm

    ia_mag_norm = float(np.linalg.norm(ia_norm))
    # Dimensionless ratio ||E[I_a]|| / (c * E[E]) is invariant to scale!
    transport_ratio = ia_mag_norm / e_tot_norm if e_tot_norm > 0.0 else 0.0
    # Clamp tiny floating-point roundoff around 0 and 1 (e.g., 1 + 1e-16 -> 1.0).
    diffuseness = float(np.clip(1.0 - transport_ratio, 0.0, 1.0))

    with np.errstate(over="raise", invalid="raise"):
        try:
            factor_i = (scale**2) / z0
            factor_e = (scale**2) / (rho0_val * (c_val**2))
            active_i = ia_norm * factor_i
            reactive_i = ir_norm * factor_i
            energy_e = e_tot_norm * factor_e
        except FloatingPointError as exc:
            raise ValueError("intensity or energy density overflowed float64") from exc
    if not (np.all(np.isfinite(active_i)) and np.all(np.isfinite(reactive_i)) and np.isfinite(energy_e)):
        raise ValueError("intensity or energy density overflowed float64")

    if ia_mag_norm > 0.0:
        intensity_unit = ia_norm / ia_mag_norm
        doa_unit = -intensity_unit
        valid_dir = True
    else:
        intensity_unit = np.zeros(3, dtype=np.float64)
        doa_unit = np.zeros(3, dtype=np.float64)
        valid_dir = False

    return {
        "active_intensity": active_i,
        "reactive_intensity": reactive_i,
        "energy_density": float(energy_e),
        "diffuseness": diffuseness,
        "intensity_unit": intensity_unit,
        "doa_unit": doa_unit,
        "valid_direction": valid_dir,
    }
