"""FOA loudspeaker decoding (Basic / max-r_E / in-phase), SOFA geometry audit, and head-tracking latency.

Numeric validation helpers follow ``codes/chapters/ch02/core/conventions.py``, while 3-D
spatial axes follow AES69-2022 (SOFA 2.1) and ACN/SN3D Ambisonics:
- Right-handed listener-local 3-D Cartesian axes: ``+x`` front, ``+y`` left, ``+z`` up
  (distinct from the clockwise-from-``+y`` horizontal ULA convention in Chapters 1-5).
- SOFA local spherical coordinates ``(azimuth_deg, elevation_deg, radius_m)`` measure
  azimuth counter-clockwise from ``+x`` (front) toward ``+y`` (left) in the horizontal
  plane and elevation in ``[-90, 90]`` degrees from the ``xy``-plane toward ``+z`` (up),
  with radius ``r >= 0``:
  ``[x, y, z] = r * [cos(el)*cos(az), cos(el)*sin(az), sin(el)]``.
- On an isotropic ``D``-dimensional ``L``-loudspeaker layout (``D = 2`` regular polygon
  ring with ``L >= 4``, or ``D = 3`` spherical 3-design such as the ``L = 6``
  regular octahedron), First-Order Ambisonics (``N = 1``) decoding with order-1 modal
  weight ``g1 >= 0`` feeds speaker ``l`` with gain
  ``d_l(n_s) = (1 + D * g1 * (u_l . n_s)) / L``.
"""
from __future__ import annotations

import numpy as np

from codes.chapters.ch02.core.conventions import (
    finite_real_array,
    finite_real_scalar,
)

# Regular 2-D horizontal hexagon (L = 6): azimuths 0, 60, 120, 180, 240, 300 deg.
REGULAR_HEXAGON_2D_DIRS = np.array(
    [
        [np.cos(np.pi * k / 3.0), np.sin(np.pi * k / 3.0), 0.0]
        for k in range(6)
    ],
    dtype=np.float64,
)

# Regular 3-D octahedron (L = 6): [+x, -x, +y, -y, +z, -z] (exact spherical 3-design).
REGULAR_OCTAHEDRON_3D_DIRS = np.array(
    [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ],
    dtype=np.float64,
)


def _validate_spatial_dim(spatial_dim: int) -> int:
    if isinstance(spatial_dim, (bool, np.bool_)) or not isinstance(spatial_dim, (int, np.integer)):
        raise ValueError("spatial_dim must be an integer (2 or 3)")
    dim = int(spatial_dim)
    if dim not in (2, 3):
        raise ValueError("spatial_dim must be 2 or 3")
    return dim


def foa_max_re_order_weight(spatial_dim: int = 3, mode: str = "max_re") -> float:
    """Return exact closed-form FOA (N=1) dipole weight ``g1`` for 2-D or 3-D decoding.

    Parameters
    ----------
    spatial_dim : {2, 3}, default 3
        ``2`` for horizontal regular-ring decoding (Chebyshev ``T_2(g1) = 0``),
        ``3`` for full-sphere periphonic decoding (Legendre ``P_2(g1) = 0``).
    mode : {"basic", "max_re", "in_phase"}, default "max_re"
        - ``"basic"``: ``g1 = 1.0`` (unit low-frequency Gerzon velocity vector ``||r_V|| = 1``).
        - ``"max_re"``: maximizes high-frequency Gerzon energy vector ``||r_E||``
          (``1 / sqrt(2)`` in 2-D, ``1 / sqrt(3)`` in 3-D).
        - ``"in_phase"``: eliminates rear anti-phase sidelobes so ``min_l d_l >= 0``
          (``1 / 2`` in 2-D, ``1 / 3`` in 3-D).
    """
    dim = _validate_spatial_dim(spatial_dim)
    if not isinstance(mode, str) or mode not in ("basic", "max_re", "in_phase"):
        raise ValueError("mode must be 'basic', 'max_re', or 'in_phase'")
    if mode == "basic":
        return 1.0
    if mode == "max_re":
        return float(1.0 / np.sqrt(float(dim)))
    return float(1.0 / float(dim))


def foa_loudspeaker_decode_metrics(
    spk_dirs: np.ndarray,
    src_dir: np.ndarray,
    *,
    spatial_dim: int = 3,
    g1: float = 1.0,
    atol: float = 1e-9,
) -> dict:
    """Evaluate FOA loudspeaker panning gains and Gerzon ``(r_V, r_E)`` localization vectors.

    Parameters
    ----------
    spk_dirs : array_like of shape ``(L, 3)``
        Unit direction vectors ``u_l`` of ``L`` equidistant loudspeakers.
    src_dir : array_like of shape ``(3,)``
        Unit direction vector ``n_s`` of the virtual plane-wave source.
    spatial_dim : {2, 3}, default 3
        Decoding manifold dimension ``D`` in ``d_l = (1 + D * g1 * (u_l . n_s)) / L``.
    g1 : float, default 1.0
        Non-negative 1st-order dipole modal weight relative to 0th-order ``g0 = 1``.
    atol : float, default 1e-9
        Tolerance for negative-lobe and layout-isotropy checks.
    """
    dim = _validate_spatial_dim(spatial_dim)
    g1_val = finite_real_scalar(g1, "g1")
    tol = finite_real_scalar(atol, "atol")
    if g1_val < 0.0:
        raise ValueError("g1 must be non-negative")
    if tol < 0.0:
        raise ValueError("atol must be non-negative")

    u_mat = finite_real_array(spk_dirs, "spk_dirs")
    if u_mat.ndim != 2 or u_mat.shape[1] != 3 or u_mat.shape[0] < 2:
        raise ValueError("spk_dirs must have shape (L, 3) with L >= 2")
    u_norms = np.linalg.norm(u_mat, axis=1)
    if not np.all(np.isfinite(u_norms)) or np.any(np.abs(u_norms - 1.0) > 1e-6):
        raise ValueError("every row of spk_dirs must be a unit vector")

    n_src = finite_real_array(src_dir, "src_dir")
    if n_src.ndim != 1 or n_src.shape != (3,):
        raise ValueError("src_dir must have shape (3,)")
    src_norm = float(np.linalg.norm(n_src))
    if not np.isfinite(src_norm) or abs(src_norm - 1.0) > 1e-6:
        raise ValueError("src_dir must be a unit vector")

    if dim == 2:
        if np.any(np.abs(u_mat[:, 2]) > 1e-6) or abs(float(n_src[2])) > 1e-6:
            raise ValueError("for spatial_dim=2, spk_dirs and src_dir must lie in the horizontal z=0 plane")

    n_spk = int(u_mat.shape[0])
    cosines = u_mat @ n_src
    with np.errstate(over="raise", invalid="raise"):
        try:
            gains = (1.0 + float(dim) * g1_val * cosines) / float(n_spk)
            gains_sq = gains * gains
        except FloatingPointError as exc:
            raise ValueError("FOA loudspeaker decoding overflowed float64") from exc
    if not (np.all(np.isfinite(gains)) and np.all(np.isfinite(gains_sq))):
        raise ValueError("FOA loudspeaker gains must be finite")

    amp_sum = float(np.sum(gains))
    energy_sum = float(np.sum(gains_sq))
    if abs(amp_sum) <= 1e-15 or energy_sum <= 1e-30:
        raise ValueError("degenerate loudspeaker gains with vanishing amplitude or energy sum")

    r_v = (gains @ u_mat) / amp_sum
    r_v_norm = float(np.linalg.norm(r_v))
    r_e = (gains_sq @ u_mat) / energy_sum
    r_e_norm = float(np.linalg.norm(r_e))
    r_e_acos_deg = float(np.rad2deg(np.arccos(np.clip(r_e_norm, -1.0, 1.0))))

    centroid = np.sum(u_mat, axis=0)
    gram = u_mat.T @ u_mat
    if dim == 3:
        target_gram = (float(n_spk) / 3.0) * np.eye(3, dtype=np.float64)
    else:
        target_gram = np.diag([float(n_spk) / 2.0, float(n_spk) / 2.0, 0.0])
    third_moment = np.einsum("li,lj,lk->ijk", u_mat, u_mat, u_mat)
    is_isotropic = bool(
        np.allclose(centroid, 0.0, atol=1e-6)
        and np.allclose(gram, target_gram, atol=1e-6)
        and np.allclose(third_moment, 0.0, atol=1e-6)
    )
    analytic_re = float((2.0 * g1_val) / (1.0 + float(dim) * (g1_val ** 2)))

    min_gain = float(np.min(gains))
    max_gain = float(np.max(gains))
    return {
        "gains": gains,
        "amp_sum": amp_sum,
        "energy_sum": energy_sum,
        "r_V": r_v,
        "r_V_norm": r_v_norm,
        "r_E": r_e,
        "r_E_norm": r_e_norm,
        "r_E_acos_deg": r_e_acos_deg,
        "energy_spread_2acos_deg": float(2.0 * r_e_acos_deg),
        "analytic_isotropic_r_E_norm": analytic_re,
        "is_isotropic_layout": is_isotropic,
        "min_gain": min_gain,
        "max_gain": max_gain,
        "has_negative_lobe": bool(min_gain < -tol),
    }


def sofa_spherical_to_cartesian(coords: np.ndarray) -> np.ndarray:
    """Convert AES69-2022 SOFA spherical ``[az_deg, el_deg, r_m]`` to Cartesian ``[x, y, z]`` m.

    Azimuth ``az_deg`` is measured counter-clockwise from ``+x`` (front) toward ``+y`` (left),
    elevation ``el_deg`` lies in ``[-90, 90]`` degrees toward ``+z`` (up), and ``r_m >= 0``.
    """
    arr = finite_real_array(coords, "coords")
    if arr.ndim < 1 or arr.shape[-1] != 3 or arr.size < 3:
        raise ValueError("coords must have trailing dimension 3 as [az_deg, el_deg, r_m]")
    az_deg = arr[..., 0]
    el_deg = arr[..., 1]
    radius = arr[..., 2]
    if np.any(radius < 0.0):
        raise ValueError("SOFA spherical radius r_m must be non-negative")
    if np.any(np.abs(el_deg) > 90.0 + 1e-12):
        raise ValueError("SOFA spherical elevation el_deg must lie within [-90, 90] degrees")
    el_clamped = np.clip(el_deg, -90.0, 90.0)
    az_rad = np.deg2rad(az_deg)
    el_rad = np.deg2rad(el_clamped)
    cos_el = np.cos(el_rad)
    x = radius * cos_el * np.cos(az_rad)
    y = radius * cos_el * np.sin(az_rad)
    z = radius * np.sin(el_rad)
    out = np.stack([x, y, z], axis=-1)
    if not np.all(np.isfinite(out)):
        raise ValueError("SOFA spherical-to-Cartesian conversion overflowed float64")
    return out


def sofa_cartesian_to_spherical(coords: np.ndarray) -> np.ndarray:
    """Convert AES69-2022 SOFA Cartesian ``[x, y, z]`` m to spherical ``[az_deg, el_deg, r_m]``.

    Returns counter-clockwise azimuth ``az_deg`` wrapped to ``(-180, 180]`` degrees (``0.0`` on
    the ``z``-axis), elevation ``el_deg`` in ``[-90, 90]`` degrees (``0.0`` at the origin), and
    Euclidean radius ``r_m >= 0`` in the ``(+x front, +y left, +z up)`` listener frame.
    """
    arr = finite_real_array(coords, "coords")
    if arr.ndim < 1 or arr.shape[-1] != 3 or arr.size < 3:
        raise ValueError("coords must have trailing dimension 3 as [x, y, z]")
    x = arr[..., 0]
    y = arr[..., 1]
    z = arr[..., 2]
    r_xy = np.hypot(x, y)
    radius = np.hypot(r_xy, z)
    if not np.all(np.isfinite(radius)):
        raise ValueError("SOFA Cartesian radius overflowed float64")
    az_deg = np.where(r_xy > 0.0, np.rad2deg(np.arctan2(y, x)), 0.0)
    az_deg = np.where(az_deg <= -180.0, az_deg + 360.0, az_deg)
    el_deg = np.where(radius > 0.0, np.rad2deg(np.arctan2(z, r_xy)), 0.0)
    return np.stack([az_deg, el_deg, radius], axis=-1)


def audit_sofa_simple_freefield_geometry(
    listener_view: np.ndarray,
    listener_up: np.ndarray,
    receiver_pos: np.ndarray,
    *,
    coord_type: str = "cartesian",
    min_head_radius_m: float = 0.03,
    max_head_radius_m: float = 0.20,
    atol: float = 1e-6,
    strict: bool = False,
) -> dict:
    """Check listener geometry and a declared symmetric-head integration policy.

    AES69-2022 ReceiverPosition is listener-local; ListenerView/ListenerUp are
    global directions. Radius/symmetry checks below are application heuristics,
    not normative AES69 conformance requirements. ``strict`` enforces this
    entire teaching policy, not a complete SOFA file validator.

    Detects four recurrent real-world SOFA integration bugs:
    1. Non-unit or non-orthogonal ``ListenerView`` / ``ListenerUp`` vectors.
    2. Left/right ear channel swap (Left ear ``R=0`` failing ``y_left > 0`` along ``up x view``).
    3. Missing ``ReceiverPosition:Type = "spherical"`` attribute causing ``[90, 0, 0.09]``
       to be parsed as ``x = 90 m`` in Cartesian metres instead of ``az = +90 deg``.
    4. Asymmetric ear placement across the listener head origin.
    """
    if type(strict) is not bool:
        raise ValueError("strict must be bool")
    if not isinstance(coord_type, str) or coord_type not in ("cartesian", "spherical"):
        raise ValueError("coord_type must be 'cartesian' or 'spherical'")
    tol = finite_real_scalar(atol, "atol")
    r_min = finite_real_scalar(min_head_radius_m, "min_head_radius_m")
    r_max = finite_real_scalar(max_head_radius_m, "max_head_radius_m")
    if tol <= 0.0:
        raise ValueError("atol must be strictly positive")
    if r_min <= 0.0 or r_max <= r_min:
        raise ValueError("require 0 < min_head_radius_m < max_head_radius_m")

    view = finite_real_array(listener_view, "listener_view")
    up = finite_real_array(listener_up, "listener_up")
    if view.shape != (3,) or up.shape != (3,):
        raise ValueError("listener_view and listener_up must each have shape (3,)")
    view_norm = float(np.linalg.norm(view))
    up_norm = float(np.linalg.norm(up))
    if view_norm <= 1e-15 or up_norm <= 1e-15:
        raise ValueError("listener_view and listener_up must be non-zero 3-D vectors")

    recv_raw = finite_real_array(receiver_pos, "receiver_pos")
    if recv_raw.shape != (2, 3):
        raise ValueError("receiver_pos for binaural SimpleFreeFieldHRIR must have shape (2, 3)")

    if coord_type == "spherical":
        recv_cart = sofa_spherical_to_cartesian(recv_raw)
        recv_sph = recv_raw.copy()
    else:
        recv_cart = recv_raw.copy()
        recv_sph = sofa_cartesian_to_spherical(recv_cart)

    view_unit_ok = abs(view_norm - 1.0) <= tol
    up_unit_ok = abs(up_norm - 1.0) <= tol
    view_hat = view / view_norm
    up_hat = up / up_norm
    ortho_dot = float(np.dot(view_hat, up_hat))
    ortho_ok = abs(ortho_dot) <= tol
    axes_orthonormal = bool(view_unit_ok and up_unit_ok and ortho_ok)

    left_raw = np.cross(up_hat, view_hat)
    left_raw_norm = float(np.linalg.norm(left_raw))
    left_axis = left_raw / left_raw_norm if left_raw_norm > 1e-12 else np.zeros(3, dtype=np.float64)

    # ReceiverPosition already uses local +x front, +y left, +z up. Projecting
    # it onto a global left axis would falsely reject a rotated listener.
    ear_left_proj = recv_cart[:, 1]
    ear_radii_m = np.linalg.norm(recv_cart, axis=1)
    ear_order_valid = bool(
        axes_orthonormal and ear_left_proj[0] > tol and ear_left_proj[1] < -tol
    )
    ears_symmetric = bool(np.allclose(recv_cart[0] + recv_cart[1], 0.0, atol=tol))
    plausible_head_radius = bool(np.all((ear_radii_m >= r_min) & (ear_radii_m <= r_max)))
    likely_unconverted_spherical = bool(coord_type == "cartesian" and np.any(ear_radii_m > 1.0))

    issues: list[str] = []
    if not (view_unit_ok and up_unit_ok):
        issues.append("non_unit_listener_axes")
    if not ortho_ok:
        issues.append("non_orthogonal_listener_axes")
    if not ear_order_valid:
        issues.append("ear_sign_swap_or_nonpositive_left_ear")
    if not plausible_head_radius:
        issues.append("implausible_head_radius_m")
    if likely_unconverted_spherical:
        issues.append("likely_unconverted_spherical_degrees")
    if not ears_symmetric:
        issues.append("asymmetric_receiver_pair")

    passed = len(issues) == 0
    if strict and not passed:
        raise ValueError("SOFA SimpleFreeFieldHRIR geometry audit failed: " + ", ".join(issues))

    return {
        "passed": passed,
        "coord_type": coord_type,
        "axes_orthonormal": axes_orthonormal,
        "view_up_dot": ortho_dot,
        "left_axis": left_axis,
        "left_axis_coordinate_system": "global",
        "receiver_coordinate_system": "listener_local",
        "policy": "symmetric binaural teaching geometry; not full AES69 conformance",
        "receiver_cartesian_m": recv_cart,
        "receiver_spherical": recv_sph,
        "ear_left_projection_m": ear_left_proj,
        "ear_radii_m": ear_radii_m,
        "ear_order_valid": ear_order_valid,
        "ears_symmetric": ears_symmetric,
        "plausible_head_radius": plausible_head_radius,
        "likely_unconverted_spherical_degrees": likely_unconverted_spherical,
        "issues": issues,
    }


def head_tracking_latency_audit(
    yaw_rate_deg_s: float,
    latency_s: np.ndarray | float,
    *,
    elevation_rad: float = 0.0,
    amplitude: float = 1.0,
    maa_deg: float = 2.0,
    prediction_horizon_s: float = 0.0,
) -> dict:
    """Compute motion-to-sound angular lag and FOA chord error during head rotation.

    Under head yaw velocity ``yaw_rate_deg_s`` and pipeline latency ``latency_s``,
    a world-locked plane wave at elevation ``phi`` and pressure magnitude ``|p|``
    suffers angular lag ``Delta_theta = |yaw_rate| * latency_s`` and exact FOA WXYZ
    Euclidean error ``||b_lag - b_ideal||_2 = 2 * |p| * cos(phi) * |sin(Delta_theta / 2)|``.
    """
    yaw_rate = finite_real_scalar(yaw_rate_deg_s, "yaw_rate_deg_s")
    el = finite_real_scalar(elevation_rad, "elevation_rad")
    amp = finite_real_scalar(amplitude, "amplitude")
    maa = finite_real_scalar(maa_deg, "maa_deg")
    pred_s = finite_real_scalar(prediction_horizon_s, "prediction_horizon_s")

    if abs(el) > 0.5 * np.pi + 1e-12:
        raise ValueError("elevation_rad must lie in [-pi/2, pi/2]")
    if amp < 0.0:
        raise ValueError("amplitude must be non-negative")
    if maa <= 0.0:
        raise ValueError("maa_deg must be strictly positive")
    if pred_s < 0.0:
        raise ValueError("prediction_horizon_s must be non-negative")

    lat_arr = finite_real_array(latency_s, "latency_s")
    is_scalar = lat_arr.ndim == 0
    if lat_arr.ndim > 1 or lat_arr.size < 1:
        raise ValueError("latency_s must be a scalar or non-empty 1-D array")
    if np.any(lat_arr < 0.0):
        raise ValueError("latency_s must be non-negative")

    el_clamped = float(np.clip(el, -0.5 * np.pi, 0.5 * np.pi))
    cos_el = float(np.cos(el_clamped))
    abs_rate = abs(yaw_rate)

    with np.errstate(over="raise", invalid="raise"):
        try:
            lag_deg = abs_rate * lat_arr
            lag_rad = np.deg2rad(lag_deg)
            chord_err = 2.0 * amp * cos_el * np.abs(np.sin(0.5 * lag_rad))
            maa_ratio = lag_deg / maa
            exceeds = lag_deg > maa

            res_lat_s = np.abs(lat_arr - pred_s)
            comp_lag_deg = abs_rate * res_lat_s
            comp_lag_rad = np.deg2rad(comp_lag_deg)
            comp_chord_err = 2.0 * amp * cos_el * np.abs(np.sin(0.5 * comp_lag_rad))
            overshoot_deg = float(abs_rate * pred_s)
        except FloatingPointError as exc:
            raise ValueError("head-tracking latency audit overflowed float64") from exc

    if not (
        np.all(np.isfinite(lag_deg))
        and np.all(np.isfinite(chord_err))
        and np.all(np.isfinite(maa_ratio))
        and np.isfinite(overshoot_deg)
    ):
        raise ValueError("head-tracking latency metrics must be finite")

    if is_scalar:
        return {
            "yaw_rate_deg_s": yaw_rate,
            "latency_s": float(lat_arr),
            "angular_lag_deg": float(lag_deg),
            "angular_lag_rad": float(lag_rad),
            "foa_chord_error": float(chord_err),
            "lag_to_maa_ratio": float(maa_ratio),
            "exceeds_maa": bool(exceeds),
            "prediction_horizon_s": pred_s,
            "compensated_lag_deg": float(comp_lag_deg),
            "compensated_foa_chord_error": float(comp_chord_err),
            "abrupt_stop_overshoot_deg": overshoot_deg,
        }

    return {
        "yaw_rate_deg_s": yaw_rate,
        "latency_s": lat_arr,
        "angular_lag_deg": lag_deg,
        "angular_lag_rad": lag_rad,
        "foa_chord_error": chord_err,
        "lag_to_maa_ratio": maa_ratio,
        "exceeds_maa": exceeds,
        "prediction_horizon_s": pred_s,
        "compensated_lag_deg": comp_lag_deg,
        "compensated_foa_chord_error": comp_chord_err,
        "abrupt_stop_overshoot_deg": overshoot_deg,
    }
