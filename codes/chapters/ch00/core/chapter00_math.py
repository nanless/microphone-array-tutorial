"""Pure-NumPy mathematical foundations kernel for Chapter 0 (E00-01..E00-16).

Implements deterministic geometric, spectral, linear-algebraic, Wirtinger-calculus,
probabilistic, Bayesian, combinatorial, and special-function primitives used across
Chapters 1..15 and Appendices A/B without SciPy or external dependencies.
"""
from __future__ import annotations

import math
import numpy as np

DEFAULT_SOUND_SPEED_MPS = 343.0


# ---------------------------------------------------------------------------
# Strict input validation helpers
# ---------------------------------------------------------------------------

def _check_bool(value, name: str) -> bool:
    if type(value) is not bool:
        raise ValueError(f"{name} must be a Python bool")
    return value


def _check_int(
    value,
    name: str,
    *,
    min_value: int | None = None,
    max_value: int | None = None,
) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer (bool rejected)")
    out = int(value)
    if min_value is not None and out < min_value:
        raise ValueError(f"{name} must be >= {min_value}, got {out}")
    if max_value is not None and out > max_value:
        raise ValueError(f"{name} must be <= {max_value}, got {out}")
    return out


def _check_real_scalar(
    value,
    name: str,
    *,
    min_value: float | None = None,
    strict_min: bool = False,
    max_value: float | None = None,
    strict_max: bool = False,
) -> float:
    if isinstance(value, (bool, np.bool_, complex, np.complexfloating)):
        raise ValueError(f"{name} must be a finite real number (bool/complex rejected)")
    if not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError(f"{name} must be a finite real number")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    if min_value is not None:
        if strict_min and not (out > min_value):
            raise ValueError(f"{name} must be > {min_value}, got {out}")
        if not strict_min and out < min_value:
            raise ValueError(f"{name} must be >= {min_value}, got {out}")
    if max_value is not None:
        if strict_max and not (out < max_value):
            raise ValueError(f"{name} must be < {max_value}, got {out}")
        if not strict_max and out > max_value:
            raise ValueError(f"{name} must be <= {max_value}, got {out}")
    return out


def _reject_bool_in_raw_array(raw: np.ndarray, name: str, *, allow_complex: bool) -> None:
    if raw.dtype == np.bool_:
        raise ValueError(f"{name} must not contain booleans")
    if not allow_complex and np.iscomplexobj(raw):
        raise ValueError(f"{name} must be real-valued (complex rejected)")
    if raw.dtype == object:
        for item in raw.flat:
            if isinstance(item, (bool, np.bool_)):
                raise ValueError(f"{name} must not contain booleans")
            if not allow_complex and isinstance(item, (complex, np.complexfloating)):
                raise ValueError(f"{name} must be real-valued (complex rejected)")


def _check_real_array(
    arr,
    name: str,
    *,
    ndim: int | None = None,
    shape: tuple[int, ...] | None = None,
    nonempty: bool = True,
) -> np.ndarray:
    raw = np.asarray(arr)
    _reject_bool_in_raw_array(raw, name, allow_complex=False)
    try:
        out = np.asarray(raw, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} could not be converted to float64") from exc
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real numbers")
    if ndim is not None and out.ndim != ndim:
        raise ValueError(f"{name} must have ndim={ndim}, got {out.ndim}")
    if shape is not None and out.shape != shape:
        raise ValueError(f"{name} must have shape={shape}, got {out.shape}")
    if nonempty and out.size == 0:
        raise ValueError(f"{name} must be non-empty")
    return out


def _check_complex_array(
    arr,
    name: str,
    *,
    ndim: int | None = None,
    shape: tuple[int, ...] | None = None,
    nonempty: bool = True,
) -> np.ndarray:
    raw = np.asarray(arr)
    _reject_bool_in_raw_array(raw, name, allow_complex=True)
    try:
        out = np.asarray(raw, dtype=np.complex128)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} could not be converted to complex128") from exc
    if not (np.all(np.isfinite(out.real)) and np.all(np.isfinite(out.imag))):
        raise ValueError(f"{name} must contain only finite complex numbers")
    if ndim is not None and out.ndim != ndim:
        raise ValueError(f"{name} must have ndim={ndim}, got {out.ndim}")
    if shape is not None and out.shape != shape:
        raise ValueError(f"{name} must have shape={shape}, got {out.shape}")
    if nonempty and out.size == 0:
        raise ValueError(f"{name} must be non-empty")
    return out


def _check_hermitian(
    mat,
    name: str,
    *,
    psd: bool = False,
    strictly_pd: bool = False,
    tol: float = 1e-11,
) -> np.ndarray:
    arr = _check_complex_array(mat, name, ndim=2)
    if arr.shape[0] != arr.shape[1]:
        raise ValueError(f"{name} must be square, got {arr.shape}")
    max_abs = float(np.max(np.abs(arr)))
    if max_abs == 0.0:
        if strictly_pd:
            raise ValueError(f"{name} must be strictly positive definite (zero matrix)")
        return np.zeros_like(arr)
    scale = max_abs
    asym = float(np.max(np.abs(arr - arr.conj().T)))
    if asym > tol * scale:
        raise ValueError(f"{name} is not Hermitian (max asymmetry {asym:.3e})")
    herm = 0.5 * (arr + arr.conj().T)
    if psd or strictly_pd:
        eigs = np.linalg.eigvalsh(herm)
        min_eig = float(np.min(eigs))
        if strictly_pd and not (min_eig > tol * scale):
            raise ValueError(f"{name} must be strictly positive definite (min eig {min_eig:.3e})")
        if psd and min_eig < -tol * scale:
            raise ValueError(f"{name} must be positive semidefinite (min eig {min_eig:.3e})")
    return herm


# ---------------------------------------------------------------------------
# §0.1 Geometry, coordinate conventions, Fresnel range & SO(3) rotations
# ---------------------------------------------------------------------------

def wrap_angle_pm180_deg(angle_deg: float) -> float:
    ang = _check_real_scalar(angle_deg, "angle_deg")
    wrapped = (ang + 180.0) % 360.0 - 180.0
    if wrapped <= -180.0:
        wrapped += 360.0
    return float(wrapped)


def coordinate_convention_bridge(
    theta_book_deg: float,
    phi_el_deg: float,
    radius_m: float = 1.0,
) -> dict:
    """Convert planar-array broadside convention (+y=0°, CW to +x) to Cartesian,
    spherical-harmonic (alpha, colatitude), and AES69/AmbiX spherical coordinates."""
    th_deg = _check_real_scalar(theta_book_deg, "theta_book_deg", min_value=-360.0, max_value=360.0)
    el_deg = _check_real_scalar(phi_el_deg, "phi_el_deg", min_value=-90.0, max_value=90.0)
    r_m = _check_real_scalar(radius_m, "radius_m", min_value=0.0, strict_min=True)

    th = math.radians(th_deg)
    el = math.radians(el_deg)
    u_vec = np.array(
        [math.cos(el) * math.sin(th), math.cos(el) * math.cos(th), math.sin(el)],
        dtype=np.float64,
    )
    pos_m = r_m * u_vec

    alpha_deg = wrap_angle_pm180_deg(90.0 - th_deg)
    colatitude_deg = 90.0 - el_deg
    alpha_rad = math.radians(alpha_deg)
    colat_rad = math.radians(colatitude_deg)
    u_from_sh = np.array(
        [
            math.sin(colat_rad) * math.cos(alpha_rad),
            math.sin(colat_rad) * math.sin(alpha_rad),
            math.cos(colat_rad),
        ],
        dtype=np.float64,
    )
    r_rec = float(np.linalg.norm(pos_m))
    el_rec_deg = math.degrees(math.asin(np.clip(pos_m[2] / r_rec, -1.0, 1.0)))
    th_rec_deg = math.degrees(math.atan2(pos_m[0], pos_m[1]))

    # AES69 / SOFA & AmbiX listener-local coordinates:
    # +x_A = front (+y_book), +y_A = left (-x_book), +z_A = up (+z_book);
    # phi_SOFA = azimuth CCW from front (+x_A -> +y_A) = -theta_book, delta_SOFA = elevation = phi_el.
    sofa_az_deg = wrap_angle_pm180_deg(-th_deg)
    sofa_el_deg = el_deg
    u_ambix_sofa = np.array([u_vec[1], -u_vec[0], u_vec[2]], dtype=np.float64)
    pos_ambix_sofa_m = r_m * u_ambix_sofa
    sofa_az_rad = math.radians(sofa_az_deg)
    sofa_el_rad = math.radians(sofa_el_deg)
    u_from_sofa = np.array(
        [
            math.cos(sofa_el_rad) * math.cos(sofa_az_rad),
            math.cos(sofa_el_rad) * math.sin(sofa_az_rad),
            math.sin(sofa_el_rad),
        ],
        dtype=np.float64,
    )
    return {
        "theta_book_deg": th_deg,
        "phi_el_deg": el_deg,
        "radius_m": r_m,
        "alpha_ccw_from_px_deg": alpha_deg,
        "colatitude_deg": colatitude_deg,
        "sofa_azimuth_deg": sofa_az_deg,
        "sofa_elevation_deg": sofa_el_deg,
        "unit_direction": u_vec,
        "position_cartesian_m": pos_m,
        "ambix_sofa_unit_direction": u_ambix_sofa,
        "ambix_sofa_position_cartesian_m": pos_ambix_sofa_m,
        "sh_consistency_error": float(np.linalg.norm(u_vec - u_from_sh)),
        "ambix_sofa_consistency_error": float(np.linalg.norm(u_ambix_sofa - u_from_sofa)),
        "roundtrip_theta_book_deg": th_rec_deg,
        "roundtrip_phi_el_deg": el_rec_deg,
        "roundtrip_radius_m": r_rec,
    }


def fresnel_spherical_wave_approx(
    source_range_m: float,
    sensor_x_m: float,
    theta_book_deg: float,
    freq_hz: float,
    aperture_m: float,
    sound_speed_mps: float = DEFAULT_SOUND_SPEED_MPS,
) -> dict:
    """Compare exact spherical-wave path difference against 1st-order plane-wave
    and 2nd-order Fresnel Taylor approximations."""
    r = _check_real_scalar(source_range_m, "source_range_m", min_value=0.0, strict_min=True)
    x_m = _check_real_scalar(sensor_x_m, "sensor_x_m")
    th_deg = _check_real_scalar(theta_book_deg, "theta_book_deg", min_value=-90.0, max_value=90.0)
    f_hz = _check_real_scalar(freq_hz, "freq_hz", min_value=0.0, strict_min=True)
    ap_m = _check_real_scalar(aperture_m, "aperture_m", min_value=0.0, strict_min=True)
    c = _check_real_scalar(sound_speed_mps, "sound_speed_mps", min_value=0.0, strict_min=True)
    if abs(x_m) >= r:
        raise ValueError("sensor_x_m magnitude must be strictly smaller than source_range_m")

    th = math.radians(th_deg)
    wavelength_m = c / f_hz
    wavenumber = 2.0 * math.pi / wavelength_m
    fraunhofer_m = 2.0 * (ap_m ** 2) / wavelength_m

    r_sensor = math.sqrt(r * r - 2.0 * r * x_m * math.sin(th) + x_m * x_m)
    delta_exact = r_sensor - r
    delta_farfield = -x_m * math.sin(th)
    delta_fresnel = -x_m * math.sin(th) + (x_m * x_m * (math.cos(th) ** 2)) / (2.0 * r)

    err_ff_m = abs(delta_exact - delta_farfield)
    err_fr_m = abs(delta_exact - delta_fresnel)
    return {
        "wavelength_m": wavelength_m,
        "fraunhofer_distance_m": fraunhofer_m,
        "in_nearfield_fresnel_zone": bool(r < fraunhofer_m),
        "sensor_distance_m": r_sensor,
        "delta_r_exact_m": delta_exact,
        "delta_r_plane_wave_m": delta_farfield,
        "delta_r_fresnel_m": delta_fresnel,
        "plane_wave_path_error_m": err_ff_m,
        "fresnel_path_error_m": err_fr_m,
        "plane_wave_phase_error_deg": math.degrees(wavenumber * err_ff_m),
        "fresnel_phase_error_deg": math.degrees(wavenumber * err_fr_m),
    }


def spherical_and_fresnel_geometry_example(
    theta_book_deg: float = 30.0,
    phi_el_deg: float = 45.0,
    radius_m: float = 1.5,
    fresnel_range_m: float = 0.50,
    sensor_x_m: float = 0.10,
    fresnel_theta_deg: float = 30.0,
    freq_hz: float = 3430.0,
    aperture_m: float = 0.20,
    sound_speed_mps: float = DEFAULT_SOUND_SPEED_MPS,
) -> dict:
    """E00-01: Spherical/AES69 coordinate bridge, spherical law of cosines, and 2nd-order Fresnel range."""
    coord = coordinate_convention_bridge(theta_book_deg, phi_el_deg, radius_m)
    peer = coordinate_convention_bridge(-60.0, 0.0, 1.0)
    cos_gamma = float(np.dot(coord["unit_direction"], peer["unit_direction"]))
    th1, el1 = math.radians(coord["theta_book_deg"]), math.radians(coord["phi_el_deg"])
    th2, el2 = math.radians(-60.0), 0.0
    cos_gamma_law = (
        math.sin(el1) * math.sin(el2)
        + math.cos(el1) * math.cos(el2) * math.cos(th1 - th2)
    )
    fresnel = fresnel_spherical_wave_approx(
        fresnel_range_m,
        sensor_x_m,
        fresnel_theta_deg,
        freq_hz,
        aperture_m,
        sound_speed_mps,
    )
    return {
        "coordinate_bridge": coord,
        "spherical_law_of_cosines": {
            "dot_product_cos": cos_gamma,
            "trig_law_cos": cos_gamma_law,
            "central_angle_deg": math.degrees(math.acos(np.clip(cos_gamma, -1.0, 1.0))),
            "abs_discrepancy": abs(cos_gamma - cos_gamma_law),
        },
        "fresnel_approximation": fresnel,
    }


def elementary_rotation_3d(axis: str, angle_rad: float) -> np.ndarray:
    ang = _check_real_scalar(angle_rad, "angle_rad")
    c, s = math.cos(ang), math.sin(ang)
    if axis == "x":
        return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=np.float64)
    if axis == "y":
        return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]], dtype=np.float64)
    if axis == "z":
        return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)
    raise ValueError(f"axis must be 'x', 'y', or 'z', got {axis!r}")


def so3_rotation_matrix_composition_example(
    yaw_deg: float = 90.0,
    pitch_deg: float = 90.0,
    roll_deg: float = 0.0,
    vector: tuple[float, float, float] | list[float] | np.ndarray = (1.0, 2.0, 2.0),
) -> dict:
    """E00-02: SO(3) rotation composition, non-commutativity, norm preservation, and passive transpose inverse."""
    yaw = math.radians(_check_real_scalar(yaw_deg, "yaw_deg"))
    pitch = math.radians(_check_real_scalar(pitch_deg, "pitch_deg"))
    roll = math.radians(_check_real_scalar(roll_deg, "roll_deg"))
    vec = _check_real_array(vector, "vector", shape=(3,))

    rz = elementary_rotation_3d("z", yaw)
    ry = elementary_rotation_3d("y", pitch)
    rx = elementary_rotation_3d("x", roll)

    r_zyx = rz @ ry @ rx
    r_yzx = ry @ rz @ rx
    e_x = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    rotated_vec = r_zyx @ vec
    recovered_vec = r_zyx.T @ rotated_vec

    return {
        "r_zyx": r_zyx,
        "r_yzx": r_yzx,
        "orthogonality_fro_error": float(np.linalg.norm(r_zyx.T @ r_zyx - np.eye(3), ord="fro")),
        "determinant": float(np.linalg.det(r_zyx)),
        "ex_under_rz_then_ry": r_zyx @ e_x,
        "ex_under_ry_then_rz": r_yzx @ e_x,
        "commutator_frobenius_norm": float(np.linalg.norm(r_zyx - r_yzx, ord="fro")),
        "input_vector": vec,
        "input_norm": float(np.linalg.norm(vec)),
        "rotated_vector": rotated_vec,
        "rotated_norm": float(np.linalg.norm(rotated_vec)),
        "passive_inverse_recovered": recovered_vec,
        "inverse_recovery_error": float(np.linalg.norm(recovered_vec - vec)),
    }


# ---------------------------------------------------------------------------
# §0.2 Complex phasors, decibels, DFT Parseval, linear conv & STFT COLA
# ---------------------------------------------------------------------------

def phasor_delay_and_decibel_budget_example(
    freq_hz: float = 1000.0,
    quadrature_delay_s: float = 2.5e-4,
    destructive_delay_s: float = 5.0e-4,
    slope_delay_s: float = 1.25e-4,
    num_mics: int = 4,
) -> dict:
    """E00-03: Complex phasor addition, linear phase delay slope, and 20log10 vs 10log10 dB budgets."""
    f0 = _check_real_scalar(freq_hz, "freq_hz", min_value=0.0, strict_min=True)
    tau_q = _check_real_scalar(quadrature_delay_s, "quadrature_delay_s", min_value=0.0)
    tau_d = _check_real_scalar(destructive_delay_s, "destructive_delay_s", min_value=0.0)
    tau_g = _check_real_scalar(slope_delay_s, "slope_delay_s", min_value=0.0)
    m_cnt = _check_int(num_mics, "num_mics", min_value=1)

    z_quad = 1.0 + np.exp(-1j * 2.0 * math.pi * f0 * tau_q)
    z_dest = 1.0 + np.exp(-1j * 2.0 * math.pi * f0 * tau_d)

    f_grid = np.array([500.0, 1500.0, 2500.0], dtype=np.float64)
    phases = np.unwrap(np.angle(np.exp(-1j * 2.0 * np.pi * f_grid * tau_g)))
    slope_rad_per_hz = float((phases[-1] - phases[0]) / (f_grid[-1] - f_grid[0]))
    recovered_group_delay_s = -slope_rad_per_hz / (2.0 * math.pi)

    coherent_amp_db = 20.0 * math.log10(m_cnt)
    coherent_pow_db = 10.0 * math.log10(m_cnt * m_cnt)
    uncorr_noise_pow_db = 10.0 * math.log10(m_cnt)
    net_snr_gain_db = coherent_pow_db - uncorr_noise_pow_db

    # Two-block Jensen inequality illustration: power ratios 1.0 (0 dB) and 99.0 (19.9564 dB)
    blk_ratios = np.array([1.0, 99.0], dtype=np.float64)
    mean_of_db = float(np.mean(10.0 * np.log10(blk_ratios)))
    db_of_mean_power = float(10.0 * math.log10(np.mean(blk_ratios)))

    return {
        "quadrature_phasor": complex(z_quad),
        "quadrature_magnitude": float(abs(z_quad)),
        "quadrature_power": float(abs(z_quad) ** 2),
        "quadrature_phase_deg": math.degrees(math.atan2(z_quad.imag, z_quad.real)),
        "destructive_magnitude": float(abs(z_dest)),
        "frequency_grid_hz": f_grid,
        "unwrapped_phases_rad": phases,
        "recovered_group_delay_us": recovered_group_delay_s * 1e6,
        "coherent_sum_db": coherent_amp_db,
        "coherent_power_db_check": coherent_pow_db,
        "uncorrelated_noise_sum_db": uncorr_noise_pow_db,
        "array_white_noise_gain_db": net_snr_gain_db,
        "block_average_of_db": mean_of_db,
        "pooled_energy_db": db_of_mean_power,
    }


def dft_parseval_and_zero_padded_conv_example(
    x: list[float] | np.ndarray = (1.0, 2.0, -1.0, 0.0),
    h: list[float] | np.ndarray = (1.0, 0.5, -0.5),
    n_pad: int = 8,
) -> dict:
    """E00-04: DFT Parseval energy identity, Hermitian symmetry, circular aliasing vs zero-padded linear conv."""
    x_arr = _check_real_array(x, "x", ndim=1)
    h_arr = _check_real_array(h, "h", ndim=1)
    lx, lh = x_arr.size, h_arr.size
    linear_len = lx + lh - 1
    n_fft_pad = _check_int(n_pad, "n_pad", min_value=linear_len)

    x_dft = np.fft.fft(x_arr)
    energy_time = float(np.sum(x_arr ** 2))
    energy_freq = float(np.sum(np.abs(x_dft) ** 2) / lx)

    # Conjugate symmetry X[N-k] == conj(X[k]) for k = 1..N-1
    k_int = np.arange(1, lx)
    conj_sym_err = float(np.max(np.abs(x_dft[lx - k_int] - np.conj(x_dft[k_int]))))

    y_linear = np.convolve(x_arr, h_arr, mode="full")
    n_circ = max(lx, lh)
    y_circ = np.fft.ifft(np.fft.fft(x_arr, n=n_circ) * np.fft.fft(h_arr, n=n_circ)).real
    y_zp_full = np.fft.ifft(np.fft.fft(x_arr, n=n_fft_pad) * np.fft.fft(h_arr, n=n_fft_pad)).real
    y_zp = y_zp_full[:linear_len]

    return {
        "x": x_arr,
        "h": h_arr,
        "dft_x": x_dft,
        "time_energy": energy_time,
        "parseval_freq_energy": energy_freq,
        "parseval_abs_error": abs(energy_time - energy_freq),
        "conjugate_symmetry_error": conj_sym_err,
        "linear_convolution": y_linear,
        "circular_convolution_no_pad": y_circ,
        "circular_aliasing_l2_error": float(np.linalg.norm(y_circ - y_linear[:n_circ])),
        "zero_padded_fft_convolution": y_zp,
        "zero_padded_max_abs_error": float(np.max(np.abs(y_zp - y_linear))),
    }


def periodic_hann_window(win_length: int) -> np.ndarray:
    n_win = _check_int(win_length, "win_length", min_value=2)
    idx = np.arange(n_win, dtype=np.float64)
    return 0.5 * (1.0 - np.cos(2.0 * math.pi * idx / n_win))


def stft_wola_roundtrip(
    signal: np.ndarray,
    analysis_win: np.ndarray,
    synthesis_win: np.ndarray,
    hop_size: int,
    *,
    normalize_denominator: bool = False,
) -> dict:
    x = _check_real_array(signal, "signal", ndim=1)
    w_a = _check_real_array(analysis_win, "analysis_win", ndim=1)
    w_s = _check_real_array(synthesis_win, "synthesis_win", ndim=1, shape=w_a.shape)
    n_win = w_a.size
    hop = _check_int(hop_size, "hop_size", min_value=1, max_value=n_win)
    norm_flag = _check_bool(normalize_denominator, "normalize_denominator")
    if x.size < 2 * n_win:
        raise ValueError("signal length must be at least 2 * win_length to form a steady interior")

    num_frames = 1 + (x.size - n_win) // hop
    recon = np.zeros(x.size, dtype=np.float64)
    cola_sum = np.zeros(x.size, dtype=np.float64)
    for m in range(num_frames):
        start = m * hop
        frame = x[start : start + n_win] * w_a
        spec = np.fft.rfft(frame, n=n_win)
        time_frame = np.fft.irfft(spec, n=n_win) * w_s
        recon[start : start + n_win] += time_frame
        cola_sum[start : start + n_win] += w_a * w_s

    last_covered = (num_frames - 1) * hop
    steady_start = n_win - hop
    steady_end = last_covered + hop
    if norm_flag:
        min_den = float(np.min(cola_sum[steady_start:steady_end]))
        if min_den <= 1e-12:
            raise ValueError("WOLA denominator vanishes in steady interior")
        recon_eval = recon.copy()
        recon_eval[steady_start:steady_end] /= cola_sum[steady_start:steady_end]
    else:
        recon_eval = recon

    steady_den = cola_sum[steady_start:steady_end]
    steady_err = float(np.max(np.abs(recon_eval[steady_start:steady_end] - x[steady_start:steady_end])))
    return {
        "reconstructed": recon_eval,
        "cola_denominator": cola_sum,
        "steady_slice": [int(steady_start), int(steady_end)],
        "steady_cola_min": float(np.min(steady_den)),
        "steady_cola_max": float(np.max(steady_den)),
        "steady_cola_ripple": float(np.max(steady_den) - np.min(steady_den)),
        "steady_max_abs_error": steady_err,
    }


def stft_cola_window_reconstruction_example(
    win_length: int = 8,
    signal_length: int = 32,
) -> dict:
    """E00-05: Periodic Hann and root-Hann WOLA COLA exact reconstruction vs unnormalized squared-Hann ripple."""
    n_win = _check_int(win_length, "win_length", min_value=4)
    if n_win % 4 != 0:
        raise ValueError("win_length must be divisible by 4 for 50% and 75% overlap checks")
    sig_len = _check_int(signal_length, "signal_length", min_value=2 * n_win)

    idx = np.arange(sig_len, dtype=np.float64)
    sig = np.cos(2.0 * math.pi * 3.0 * idx / sig_len) + 0.5 * np.sin(2.0 * math.pi * 7.0 * idx / sig_len)

    w_hann = periodic_hann_window(n_win)
    w_rect = np.ones(n_win, dtype=np.float64)
    w_sqrt = np.sqrt(w_hann)

    ola_hann_50 = stft_wola_roundtrip(sig, w_hann, w_rect, n_win // 2, normalize_denominator=False)
    wola_sqrt_50 = stft_wola_roundtrip(sig, w_sqrt, w_sqrt, n_win // 2, normalize_denominator=False)
    wrong_sq_50 = stft_wola_roundtrip(sig, w_hann, w_hann, n_win // 2, normalize_denominator=False)
    wola_sq_75_raw = stft_wola_roundtrip(sig, w_hann, w_hann, n_win // 4, normalize_denominator=False)
    wola_sq_75_norm = stft_wola_roundtrip(sig, w_hann, w_hann, n_win // 4, normalize_denominator=True)

    return {
        "win_length": n_win,
        "periodic_hann_taps": w_hann,
        "hann_ola_50pct": {
            "steady_cola_min": ola_hann_50["steady_cola_min"],
            "steady_cola_max": ola_hann_50["steady_cola_max"],
            "steady_max_abs_error": ola_hann_50["steady_max_abs_error"],
        },
        "sqrt_hann_wola_50pct": {
            "steady_cola_min": wola_sqrt_50["steady_cola_min"],
            "steady_cola_max": wola_sqrt_50["steady_cola_max"],
            "steady_max_abs_error": wola_sqrt_50["steady_max_abs_error"],
        },
        "unnormalized_squared_hann_50pct_counterexample": {
            "steady_cola_min": wrong_sq_50["steady_cola_min"],
            "steady_cola_max": wrong_sq_50["steady_cola_max"],
            "steady_cola_ripple": wrong_sq_50["steady_cola_ripple"],
            "ripple_ratio_db": 20.0 * math.log10(wrong_sq_50["steady_cola_max"] / wrong_sq_50["steady_cola_min"]),
            "steady_max_abs_error": wrong_sq_50["steady_max_abs_error"],
        },
        "squared_hann_75pct": {
            "constant_cola_gain": wola_sq_75_raw["steady_cola_min"],
            "normalized_max_abs_error": wola_sq_75_norm["steady_max_abs_error"],
        },
    }


# ---------------------------------------------------------------------------
# §0.3 Complex subspaces, orthogonal projections, SVD/Tikhonov & Woodbury/Schur
# ---------------------------------------------------------------------------

def complex_column_space_projection(a_mat: np.ndarray, b_vec: np.ndarray) -> dict:
    a = _check_complex_array(a_mat, "a_mat", ndim=2)
    m, n = a.shape
    b = _check_complex_array(b_vec, "b_vec", shape=(m,))
    if m < n:
        raise ValueError("a_mat must be tall or square (m >= n) with full column rank")
    gram = _check_hermitian(a.conj().T @ a, "A^H A", strictly_pd=True)
    gram_inv = np.linalg.inv(gram)
    z_ls = gram_inv @ (a.conj().T @ b)
    p_col = a @ gram_inv @ a.conj().T
    p_perp = np.eye(m, dtype=np.complex128) - p_col
    b_par = p_col @ b
    b_perp = p_perp @ b
    return {
        "gram_matrix": gram,
        "gram_det": float(np.linalg.det(gram).real),
        "ls_coordinates": z_ls,
        "projector_col": p_col,
        "projector_left_null": p_perp,
        "b_parallel": b_par,
        "b_orthogonal": b_perp,
        "idempotency_fro_error": float(np.linalg.norm(p_col @ p_col - p_col, ord="fro")),
        "hermitian_fro_error": float(np.linalg.norm(p_col.conj().T - p_col, ord="fro")),
        "orthogonality_inner_abs": float(abs(np.vdot(b_par, b_perp))),
        "normal_equation_residual_norm": float(np.linalg.norm(a.conj().T @ b_perp)),
        "norm_sq_b": float(np.vdot(b, b).real),
        "norm_sq_parallel": float(np.vdot(b_par, b_par).real),
        "norm_sq_orthogonal": float(np.vdot(b_perp, b_perp).real),
        "trace_projector_col": float(np.trace(p_col).real),
        "trace_projector_left_null": float(np.trace(p_perp).real),
    }


def orthogonal_projection_and_four_subspaces_example() -> dict:
    """E00-06: Complex column-space projector P_R(A), Hermitian idempotency, and four-subspace Pythagoras split."""
    a_mat = np.array([[1.0, 0.0], [1.0j, 1.0], [0.0, -1.0j]], dtype=np.complex128)
    b_vec = np.array([1.0, 0.0, 1.0j], dtype=np.complex128)
    res = complex_column_space_projection(a_mat, b_vec)
    u_left_null = np.array([1.0, -1.0j, 1.0], dtype=np.complex128) / math.sqrt(3.0)
    return {
        "a_mat": a_mat,
        "b_vec": b_vec,
        "left_null_unit_basis": u_left_null,
        "left_null_check_norm": float(np.linalg.norm(a_mat.conj().T @ u_left_null)),
        **res,
    }


def tikhonov_svd_solve(
    g_mat: np.ndarray,
    y_obs: np.ndarray,
    loading_eps: float,
) -> dict:
    g = _check_complex_array(g_mat, "g_mat", ndim=2)
    y = _check_complex_array(y_obs, "y_obs", shape=(g.shape[0],))
    eps = _check_real_scalar(loading_eps, "loading_eps", min_value=0.0)

    u, s, vh = np.linalg.svd(g, full_matrices=False)
    if eps == 0.0 and (
        g.shape[0] < g.shape[1]
        or float(s[0]) == 0.0
        or float(s[-1]) <= 1e-12 * max(1.0, float(s[0]))
    ):
        raise ValueError("unregularized inverse requires m >= n and full column rank")
    filt = s / (s * s + eps)
    x_svd = vh.conj().T @ (filt * (u.conj().T @ y))
    n = g.shape[1]
    normal_mat = g.conj().T @ g + eps * np.eye(n, dtype=np.complex128)
    x_normal = np.linalg.solve(normal_mat, g.conj().T @ y)
    cond_2 = (
        float(s[0] / s[-1])
        if (g.shape[0] >= g.shape[1] and float(s[-1]) > 0.0)
        else float("inf")
    )
    return {
        "singular_values": s,
        "condition_number_2": cond_2,
        "filter_factors": filt,
        "solution": x_svd,
        "normal_eq_consistency_error": float(np.linalg.norm(x_svd - x_normal)),
    }


def svd_condition_and_tikhonov_loading_example(
    sigma_max: float = 10.0,
    sigma_min: float = 0.01,
    noise_norm: float = 0.10,
    selected_loading: float = 0.10,
) -> dict:
    """E00-07: SVD 2-norm condition number, worst-case noise amplification, and Tikhonov spectral filtering."""
    s1 = _check_real_scalar(sigma_max, "sigma_max", min_value=0.0, strict_min=True)
    s2 = _check_real_scalar(sigma_min, "sigma_min", min_value=0.0, strict_min=True)
    if s2 > s1:
        raise ValueError("sigma_min must not exceed sigma_max")
    eta = _check_real_scalar(noise_norm, "noise_norm", min_value=0.0)
    eps_sel = _check_real_scalar(selected_loading, "selected_loading", min_value=0.0, strict_min=True)

    inv_rt2 = 1.0 / math.sqrt(2.0)
    u_mat = inv_rt2 * np.array([[1.0, 1.0], [1.0, -1.0]], dtype=np.float64)
    g_mat = u_mat @ np.diag([s1, s2])
    x_true = np.array([1.0, 0.0], dtype=np.float64)
    y_clean = g_mat @ x_true
    delta_y = eta * u_mat[:, 1]
    y_noisy = y_clean + delta_y

    unreg = tikhonov_svd_solve(g_mat, y_noisy, 0.0)
    reg_sel = tikhonov_svd_solve(g_mat, y_noisy, eps_sel)

    err_unreg = float(np.linalg.norm(unreg["solution"] - x_true))
    err_reg = float(np.linalg.norm(reg_sel["solution"] - x_true))
    rel_obs_noise = float(np.linalg.norm(delta_y) / np.linalg.norm(y_clean))
    rel_sol_err_unreg = float(err_unreg / np.linalg.norm(x_true))

    sweep = []
    for eps in (0.0, 1e-4, 1e-2, 1e-1, 1.0):
        sol = tikhonov_svd_solve(g_mat, y_noisy, eps)["solution"]
        sweep.append(
            {
                "loading_eps": eps,
                "solution_real": sol.real.tolist(),
                "error_norm": float(np.linalg.norm(sol - x_true)),
            }
        )

    return {
        "g_matrix": g_mat,
        "singular_values": unreg["singular_values"],
        "condition_number_2": unreg["condition_number_2"],
        "relative_observation_perturbation": rel_obs_noise,
        "unregularized_solution": unreg["solution"].real,
        "unregularized_error_norm": err_unreg,
        "relative_unregularized_error": rel_sol_err_unreg,
        "bound_tightness_ratio": rel_sol_err_unreg / (unreg["condition_number_2"] * rel_obs_noise)
        if eta > 0.0
        else 1.0,
        "tikhonov_loading_eps": eps_sel,
        "tikhonov_solution": reg_sel["solution"].real,
        "tikhonov_error_norm": err_reg,
        "error_reduction_factor": err_unreg / err_reg if err_reg > 0.0 else float("inf"),
        "loading_sweep": sweep,
    }


def _require_nonsingular(mat: np.ndarray, name: str, rcond: float = 1e-12) -> None:
    sv = np.linalg.svd(mat, compute_uv=False)
    if sv.size == 0 or not (float(sv[0]) > 0.0 and float(sv[-1]) > rcond * float(sv[0])):
        raise ValueError(f"{name} must be non-singular")


def woodbury_inverse_update(
    a_inv: np.ndarray,
    u_mat: np.ndarray,
    c_inv: np.ndarray,
    v_herm: np.ndarray,
) -> np.ndarray:
    """Compute (A + U C V)^{-1} = A^{-1} - A^{-1} U (C^{-1} + V A^{-1} U)^{-1} V A^{-1}."""
    ai = _check_complex_array(a_inv, "a_inv", ndim=2)
    u = _check_complex_array(u_mat, "u_mat", ndim=2)
    ci = _check_complex_array(c_inv, "c_inv", ndim=2)
    vh = _check_complex_array(v_herm, "v_herm", ndim=2)
    m = ai.shape[0]
    r = u.shape[1]
    if ai.shape != (m, m) or u.shape != (m, r) or ci.shape != (r, r) or vh.shape != (r, m):
        raise ValueError("incompatible shapes in woodbury_inverse_update")
    _require_nonsingular(ai, "a_inv")
    _require_nonsingular(ci, "c_inv")
    capacitance = ci + vh @ ai @ u
    _require_nonsingular(capacitance, "capacitance")
    return ai - (ai @ u) @ np.linalg.solve(capacitance, vh @ ai)


def woodbury_and_schur_complement_example(
    noise_var: float = 1.0,
    interferer_var: float = 2.0,
) -> dict:
    """E00-08: Sherman-Morrison-Woodbury rank-1 covariance inverse and partitioned Schur determinant."""
    s_n2 = _check_real_scalar(noise_var, "noise_var", min_value=0.0, strict_min=True)
    s_i2 = _check_real_scalar(interferer_var, "interferer_var", min_value=0.0)

    u_vec = np.array([1.0, 1.0j, -1.0], dtype=np.complex128)
    r_mat = s_n2 * np.eye(3, dtype=np.complex128) + s_i2 * np.outer(u_vec, u_vec.conj())
    r0_inv = (1.0 / s_n2) * np.eye(3, dtype=np.complex128)

    # Sherman-Morrison closed form & general Woodbury helper
    denom = 1.0 + (s_i2 / s_n2) * float(np.vdot(u_vec, u_vec).real)
    r_inv_sm = r0_inv - ((s_i2 / (s_n2 * s_n2)) / denom) * np.outer(u_vec, u_vec.conj())
    r_inv_wb = woodbury_inverse_update(
        r0_inv,
        u_vec[:, None],
        np.array([[1.0 / s_i2]], dtype=np.complex128) if s_i2 > 0 else np.eye(1, dtype=np.complex128),
        u_vec.conj()[None, :] if s_i2 > 0 else np.zeros((1, 3), dtype=np.complex128),
    )
    r_inv_direct = np.linalg.inv(r_mat)

    # Partition R into 2x2 top-left A11 and scalar bottom-right a22
    a11 = r_mat[:2, :2]
    a12 = r_mat[:2, 2]
    a22 = float(r_mat[2, 2].real)
    schur_s = float((a22 - np.vdot(a12, np.linalg.solve(a11, a12))).real)
    det_a11 = float(np.linalg.det(a11).real)
    det_via_schur = det_a11 * schur_s
    det_sylvester = (s_n2 ** 3) * denom
    det_direct = float(np.linalg.det(r_mat).real)

    return {
        "r_matrix": r_mat,
        "sherman_morrison_denominator": denom,
        "r_inv_sherman_morrison": r_inv_sm,
        "woodbury_vs_direct_fro_error": float(np.linalg.norm(r_inv_wb - r_inv_direct, ord="fro")),
        "sm_vs_direct_fro_error": float(np.linalg.norm(r_inv_sm - r_inv_direct, ord="fro")),
        "interferer_attenuation_ratio": float(np.linalg.norm(r_inv_sm @ u_vec) / np.linalg.norm(r0_inv @ u_vec)),
        "partition_det_a11": det_a11,
        "schur_complement_s": schur_s,
        "det_via_schur": det_via_schur,
        "det_sylvester": det_sylvester,
        "det_abs_error": abs(det_via_schur - det_direct),
    }


# ---------------------------------------------------------------------------
# §0.4 Rayleigh quotient, GEV joint diagonalization, vec/Kronecker/Khatri-Rao
# ---------------------------------------------------------------------------

def solve_hermitian_gev(r_a: np.ndarray, r_b: np.ndarray) -> dict:
    """Solve R_a w = lambda R_b w for Hermitian R_a and Hermitian positive-definite R_b
    via Cholesky whitening R_b = L L^H, returning ascending generalized eigenvalues and
    B-orthonormal eigenvectors Phi satisfying Phi^H R_b Phi = I and Phi^H R_a Phi = diag(lambda)."""
    a = _check_hermitian(r_a, "r_a")
    b = _check_hermitian(r_b, "r_b", strictly_pd=True)
    if a.shape != b.shape:
        raise ValueError("r_a and r_b must have identical dimensions")

    chol_l = np.linalg.cholesky(b)
    chol_inv = np.linalg.inv(chol_l)
    whitened = 0.5 * (chol_inv @ a @ chol_inv.conj().T + (chol_inv @ a @ chol_inv.conj().T).conj().T)
    eigvals, q_unitary = np.linalg.eigh(whitened)
    phi = np.linalg.solve(chol_l.conj().T, q_unitary)
    return {
        "cholesky_l": chol_l,
        "whitened_matrix": whitened,
        "gev_eigenvalues": eigvals,
        "joint_diagonalizer_phi": phi,
        "b_whiten_fro_error": float(np.linalg.norm(phi.conj().T @ b @ phi - np.eye(a.shape[0]), ord="fro")),
        "a_diag_fro_error": float(np.linalg.norm(phi.conj().T @ a @ phi - np.diag(eigvals), ord="fro")),
    }


def rayleigh_quotient_and_gev_diagonalization_example() -> dict:
    """E00-09: Hermitian Rayleigh bounds, GEV whitening, and correlated-noise SINR gain over unwhitened matched filter."""
    r_nn = np.array([[2.0, 1.0], [1.0, 2.0]], dtype=np.complex128)
    a_tar = np.array([1.0, 1.0j], dtype=np.complex128)
    r_ss = np.outer(a_tar, a_tar.conj())

    nn_eigs = np.linalg.eigvalsh(r_nn)
    gev = solve_hermitian_gev(r_ss, r_nn)
    w_gev = gev["joint_diagonalizer_phi"][:, -1]

    sinr_gev = float((np.vdot(w_gev, r_ss @ w_gev).real) / (np.vdot(w_gev, r_nn @ w_gev).real))
    sinr_unwhitened = float((np.vdot(a_tar, r_ss @ a_tar).real) / (np.vdot(a_tar, r_nn @ a_tar).real))
    closed_form_max_sinr = float(np.vdot(a_tar, np.linalg.solve(r_nn, a_tar)).real)

    return {
        "r_nn": r_nn,
        "r_ss": r_ss,
        "r_nn_eigenvalues": nn_eigs,
        "gev_eigenvalues": gev["gev_eigenvalues"],
        "principal_gev_weight": w_gev,
        "max_gev_sinr": sinr_gev,
        "closed_form_aH_Rinv_a": closed_form_max_sinr,
        "unwhitened_matched_sinr": sinr_unwhitened,
        "whitening_sinr_gain_db": 10.0 * math.log10(sinr_gev / sinr_unwhitened),
        "b_whiten_fro_error": gev["b_whiten_fro_error"],
        "a_diag_fro_error": gev["a_diag_fro_error"],
    }


def vec_col_major(mat: np.ndarray) -> np.ndarray:
    arr = _check_complex_array(mat, "mat", ndim=2)
    return arr.reshape(-1, order="F")


def khatri_rao_columnwise(b_mat: np.ndarray, a_mat: np.ndarray) -> np.ndarray:
    """Column-wise Kronecker (Khatri-Rao) product B ⊙ A of shape (P*M, K)."""
    b = _check_complex_array(b_mat, "b_mat", ndim=2)
    a = _check_complex_array(a_mat, "a_mat", ndim=2)
    if b.shape[1] != a.shape[1]:
        raise ValueError(f"column counts must match for Khatri-Rao product: {b.shape} vs {a.shape}")
    return np.einsum("pk,mk->pmk", b, a).reshape(b.shape[0] * a.shape[0], b.shape[1])


def vec_kronecker_and_khatri_rao_example() -> dict:
    """E00-10: Column-major vec(AXB)=(B^T ⊗ A)vec(X), Kronecker spectrum, and Khatri-Rao coarray expansion."""
    a_mat = np.array([[1.0, 2.0j], [-1.0j, 0.5]], dtype=np.complex128)
    x_mat = np.array([[2.0, -1.0 + 1.0j], [0.5j, 3.0]], dtype=np.complex128)
    b_mat = np.array([[1.0 + 1.0j, -0.5], [2.0, 1.0 - 2.0j]], dtype=np.complex128)

    lhs_vec = vec_col_major(a_mat @ x_mat @ b_mat)
    rhs_vec = np.kron(b_mat.T, a_mat) @ vec_col_major(x_mat)
    wrong_herm_vec = np.kron(b_mat.conj().T, a_mat) @ vec_col_major(x_mat)

    p_mat = np.diag([2.0, 5.0])
    q_mat = np.array([[3.0, 1.0], [1.0, 3.0]], dtype=np.float64)
    kron_pq = np.kron(p_mat, q_mat)
    kron_eigs = np.linalg.eigvalsh(kron_pq)
    outer_eigs = np.sort(np.outer(np.linalg.eigvalsh(p_mat), np.linalg.eigvalsh(q_mat)).ravel())

    # 3-sensor sparse linear array at integer lags {0, 1, 3} -> difference coarray {-3..+3} (7 unique lags)
    sensor_pos = np.array([0, 1, 3], dtype=np.int64)
    diff_lags = sorted({int(m - n) for m in sensor_pos for n in sensor_pos})
    phases = np.array([0.35 * math.pi, -0.55 * math.pi], dtype=np.float64)
    steer = np.exp(1j * sensor_pos[:, None] * phases[None, :])
    src_powers = np.array([2.0, 1.0], dtype=np.float64)
    noise_var = 0.25

    r_cov = steer @ np.diag(src_powers) @ steer.conj().T + noise_var * np.eye(3, dtype=np.complex128)
    kr_virtual = khatri_rao_columnwise(steer.conj(), steer)
    vec_r_model = kr_virtual @ src_powers + noise_var * vec_col_major(np.eye(3, dtype=np.complex128))

    # Schur product theorem check: Hadamard product of two PSD matrices is PSD
    hadamard_psd = (steer @ steer.conj().T) * (steer @ steer.conj().T).conj()
    hadamard_min_eig = float(np.min(np.linalg.eigvalsh(hadamard_psd)))

    return {
        "vec_axb_max_abs_error": float(np.max(np.abs(lhs_vec - rhs_vec))),
        "wrong_hermitian_transpose_error": float(np.linalg.norm(lhs_vec - wrong_herm_vec)),
        "kronecker_trace": float(np.trace(kron_pq)),
        "product_of_traces": float(np.trace(p_mat) * np.trace(q_mat)),
        "kronecker_eigenvalues": kron_eigs,
        "expected_product_eigenvalues": outer_eigs,
        "sparse_sensor_positions": sensor_pos,
        "difference_coarray_unique_lags": diff_lags,
        "num_unique_virtual_lags": len(diff_lags),
        "khatri_rao_vec_max_abs_error": float(np.max(np.abs(vec_col_major(r_cov) - vec_r_model))),
        "hadamard_psd_min_eigenvalue": hadamard_min_eig,
    }


# ---------------------------------------------------------------------------
# §0.5 Wirtinger calculus, constrained Lagrange (MVDR/LCMV), NNLS KKT & MM
# ---------------------------------------------------------------------------

def solve_lcmv_weights(
    r_cov: np.ndarray,
    c_mat: np.ndarray,
    f_vec: np.ndarray,
) -> np.ndarray:
    """Closed-form complex Lagrange LCMV weights w = R^{-1} C (C^H R^{-1} C)^{-1} f
    for C^H w = f."""
    r = _check_hermitian(r_cov, "r_cov", strictly_pd=True)
    c = _check_complex_array(c_mat, "c_mat", ndim=2)
    m, k = c.shape
    if r.shape[0] != m:
        raise ValueError("c_mat row count must match r_cov dimension")
    if k > m:
        raise ValueError("number of linear constraints cannot exceed sensor count")
    f = _check_complex_array(f_vec, "f_vec", shape=(k,))
    r_inv_c = np.linalg.solve(r, c)
    gram_c = _check_hermitian(c.conj().T @ r_inv_c, "C^H R^{-1} C", strictly_pd=True)
    lagrange_dual = np.linalg.solve(gram_c, f)
    return r_inv_c @ lagrange_dual


def wirtinger_lagrange_mvdr_and_lcmv_example(fd_eps: float = 1e-6) -> dict:
    """E00-11: Wirtinger conjugate gradient verification and closed-form MVDR vs null-constrained LCMV."""
    eps = _check_real_scalar(fd_eps, "fd_eps", min_value=1e-10, max_value=1e-2)
    a_tar = np.array([1.0, 1.0, 1.0], dtype=np.complex128)
    a_int = np.array([1.0, 1.0j, -1.0], dtype=np.complex128)
    r_cov = np.eye(3, dtype=np.complex128) + 2.0 * np.outer(a_int, a_int.conj())

    w_test = np.array([0.4 - 0.2j, -0.3 + 0.5j, 0.2 + 0.1j], dtype=np.complex128)
    delta_w = np.array([0.3 + 0.4j, -0.2 + 0.1j, 0.5 - 0.3j], dtype=np.complex128)

    def quad_cost(w: np.ndarray) -> float:
        return float(np.vdot(w, r_cov @ w).real)

    # Wirtinger gradient nabla_{w^*} (w^H R w) = R w; real directional derivative = 2 Re{(R w)^H delta_w}
    grad_w_conj = r_cov @ w_test
    analytic_dir_deriv = 2.0 * float(np.vdot(grad_w_conj, delta_w).real)
    central_fd = (quad_cost(w_test + eps * delta_w) - quad_cost(w_test - eps * delta_w)) / (2.0 * eps)

    # Numeric partial check: 0.5 * (dJ/dw_R + j dJ/dw_I) == R w
    fd_wirtinger = np.zeros(3, dtype=np.complex128)
    for idx in range(3):
        e_k = np.zeros(3, dtype=np.complex128)
        e_k[idx] = 1.0
        dj_dr = (quad_cost(w_test + eps * e_k) - quad_cost(w_test - eps * e_k)) / (2.0 * eps)
        dj_di = (quad_cost(w_test + 1j * eps * e_k) - quad_cost(w_test - 1j * eps * e_k)) / (2.0 * eps)
        fd_wirtinger[idx] = 0.5 * (dj_dr + 1j * dj_di)

    w_mvdr = solve_lcmv_weights(r_cov, a_tar[:, None], np.array([1.0 + 0.0j]))
    c_lcmv = np.column_stack([a_tar, a_int])
    f_lcmv = np.array([1.0 + 0.0j, 0.0 + 0.0j])
    w_lcmv = solve_lcmv_weights(r_cov, c_lcmv, f_lcmv)

    p_mvdr = quad_cost(w_mvdr)
    p_lcmv = quad_cost(w_lcmv)

    return {
        "wirtinger_gradient_rw": grad_w_conj,
        "wirtinger_fd_max_abs_error": float(np.max(np.abs(grad_w_conj - fd_wirtinger))),
        "directional_derivative_analytic": analytic_dir_deriv,
        "directional_derivative_fd_error": abs(analytic_dir_deriv - central_fd),
        "w_mvdr": w_mvdr,
        "mvdr_target_response": complex(np.vdot(w_mvdr, a_tar)),
        "mvdr_interferer_response": complex(np.vdot(w_mvdr, a_int)),
        "mvdr_output_power": p_mvdr,
        "w_lcmv": w_lcmv,
        "lcmv_target_response": complex(np.vdot(w_lcmv, a_tar)),
        "lcmv_interferer_response_abs": float(abs(np.vdot(w_lcmv, a_int))),
        "lcmv_output_power": p_lcmv,
        "lcmv_minus_mvdr_power": p_lcmv - p_mvdr,
    }


def solve_nnls_projected_gradient(
    p_mat: np.ndarray,
    b_vec: np.ndarray,
    *,
    max_iters: int = 500,
    tol: float = 1e-13,
) -> dict:
    """Minimize 0.5 * ||P q - b||_2^2 subject to q >= 0 via exact-Lipschitz projected gradient
    with active-set KKT polish and full KKT verification."""
    p = _check_real_array(p_mat, "p_mat", ndim=2)
    b = _check_real_array(b_vec, "b_vec", shape=(p.shape[0],))
    n_it = _check_int(max_iters, "max_iters", min_value=1)
    eps_tol = _check_real_scalar(tol, "tol", min_value=1e-15, max_value=1e-4)

    gram = p.T @ p
    rhs = p.T @ b
    lip = float(np.max(np.linalg.eigvalsh(gram)))
    if lip <= 0.0:
        raise ValueError("p_mat must not be identically zero")
    step = 1.0 / lip
    q = np.zeros(p.shape[1], dtype=np.float64)
    for _ in range(n_it):
        grad = gram @ q - rhs
        q_next = np.maximum(0.0, q - step * grad)
        if float(np.linalg.norm(q_next - q)) <= eps_tol:
            q = q_next
            break
        q = q_next

    # Active-set subspace polish on positive support
    support = q > 1e-10
    if np.any(support):
        p_sup = p[:, support]
        q_sup, *_ = np.linalg.lstsq(p_sup, b, rcond=None)
        if np.all(q_sup >= 0.0):
            q_cand = np.zeros_like(q)
            q_cand[support] = q_sup
            grad_cand = gram @ q_cand - rhs
            if np.all(grad_cand >= -1e-12):
                q = q_cand

    grad_final = gram @ q - rhs
    comp_slack = float(np.max(np.abs(q * grad_final)))
    return {
        "q_opt": q,
        "dual_multiplier_mu": grad_final,
        "objective_value": 0.5 * float(np.sum((p @ q - b) ** 2)),
        "primal_min_entry": float(np.min(q)),
        "dual_min_entry": float(np.min(grad_final)),
        "complementary_slackness_max_abs": comp_slack,
    }


def mm_smoothed_l1_solve(
    a_mat: np.ndarray,
    y_vec: np.ndarray,
    reg_lambda: float,
    smooth_delta: float = 1e-3,
    n_iters: int = 20,
) -> dict:
    """Minimize J(x) = 0.5*||A x - y||_2^2 + reg_lambda * sum_i sqrt(|x_i|^2 + delta^2)
    via quadratic Majorization-Minimization (IRLS upper-bound minimization)."""
    a = _check_complex_array(a_mat, "a_mat", ndim=2)
    y = _check_complex_array(y_vec, "y_vec", shape=(a.shape[0],))
    lam = _check_real_scalar(reg_lambda, "reg_lambda", min_value=0.0, strict_min=True)
    delta = _check_real_scalar(smooth_delta, "smooth_delta", min_value=0.0, strict_min=True)
    iters = _check_int(n_iters, "n_iters", min_value=1)

    m, n = a.shape
    gram = a.conj().T @ a
    ah_y = a.conj().T @ y
    x = np.zeros(n, dtype=np.complex128)

    def obj(vec: np.ndarray) -> float:
        res = a @ vec - y
        pen = np.sum(np.sqrt(np.abs(vec) ** 2 + delta * delta))
        return 0.5 * float(np.vdot(res, res).real) + lam * float(pen)

    history = [obj(x)]
    surrogate_touch_gaps = []
    for _ in range(iters):
        weights = 1.0 / np.sqrt(np.abs(x) ** 2 + delta * delta)
        # Surrogate Q(z | x) touches J(x) at z = x and upper-bounds J(z) everywhere
        sys_mat = gram + lam * np.diag(weights)
        x_next = np.linalg.solve(sys_mat, ah_y)
        res_next = a @ x_next - y
        q_val = (
            0.5 * float(np.vdot(res_next, res_next).real)
            + lam * float(np.sum(np.sqrt(np.abs(x) ** 2 + delta * delta)
                                 + 0.5 * weights * (np.abs(x_next) ** 2 - np.abs(x) ** 2)))
        )
        j_next = obj(x_next)
        surrogate_touch_gaps.append(float(q_val - j_next))
        x = x_next
        history.append(j_next)

    diffs = np.diff(history)
    return {
        "solution": x,
        "objective_history": np.asarray(history, dtype=np.float64),
        "monotonically_non_increasing": bool(np.all(diffs <= 1e-14)),
        "surrogate_always_upper_bounds": bool(min(surrogate_touch_gaps) >= -1e-14),
        "initial_objective": history[0],
        "final_objective": history[-1],
    }


def nnls_kkt_and_mm_reweighted_l1_example() -> dict:
    """E00-12: NNLS KKT complementary slackness vs naive negative clipping, and MM surrogate monotonic descent."""
    p_mat = np.array([[1.0, 0.5], [0.5, 1.0]], dtype=np.float64)
    b_vec = np.array([1.0, -0.2], dtype=np.float64)

    q_ls = np.linalg.solve(p_mat, b_vec)
    q_naive_clip = np.maximum(0.0, q_ls)
    naive_clip_obj = 0.5 * float(np.sum((p_mat @ q_naive_clip - b_vec) ** 2))
    naive_clip_grad = p_mat.T @ (p_mat @ q_naive_clip - b_vec)

    nnls = solve_nnls_projected_gradient(p_mat, b_vec)

    a_mm = np.array([[1.0, 0.6, -0.3], [0.2, 1.0, 0.5]], dtype=np.complex128)
    x_sparse = np.array([1.5 + 0.0j, 0.0 + 0.0j, -0.8 + 0.0j], dtype=np.complex128)
    y_mm = a_mm @ x_sparse
    mm_res = mm_smoothed_l1_solve(a_mm, y_mm, reg_lambda=0.05, smooth_delta=1e-2, n_iters=25)

    return {
        "unconstrained_ls": q_ls,
        "naive_clipped_ls": q_naive_clip,
        "naive_clipped_objective": naive_clip_obj,
        "naive_clipped_kkt_violation_grad0": float(abs(naive_clip_grad[0])),
        "nnls_kkt": nnls,
        "mm_reweighted_l1": {
            "solution_real": mm_res["solution"].real,
            "initial_objective": mm_res["initial_objective"],
            "final_objective": mm_res["final_objective"],
            "monotonically_non_increasing": mm_res["monotonically_non_increasing"],
            "surrogate_always_upper_bounds": mm_res["surrogate_always_upper_bounds"],
        },
    }


# ---------------------------------------------------------------------------
# §0.6 Complex Gaussian propriety, DFT real bins, CRLB, Wiener & Kalman
# ---------------------------------------------------------------------------

def complex_second_order_moments(
    var_real: float,
    var_imag: float,
    corr_ri: float = 0.0,
) -> dict:
    """Analytical covariance E[|z|^2] and pseudo-covariance E[z^2] for zero-mean z = u + j v."""
    su2 = _check_real_scalar(var_real, "var_real", min_value=0.0)
    sv2 = _check_real_scalar(var_imag, "var_imag", min_value=0.0)
    rho = _check_real_scalar(corr_ri, "corr_ri", min_value=-1.0, max_value=1.0)
    cov_r = su2 + sv2
    pseudo_c = complex(su2 - sv2, 2.0 * rho * math.sqrt(su2 * sv2))
    return {
        "variance_r": cov_r,
        "pseudo_variance_c": pseudo_c,
        "pseudo_variance_abs": abs(pseudo_c),
        "is_proper_second_order": bool(abs(pseudo_c) <= 1e-15),
    }


def dft_real_white_noise_pseudo_covariance(n_fft: int = 8, sigma2: float = 1.0) -> dict:
    """Exact analytical covariance E[|X[k]|^2] and pseudo-covariance E[X[k]^2] of the unitary
    N-point DFT of real zero-mean white noise x ~ N(0, sigma2 I_N)."""
    n = _check_int(n_fft, "n_fft", min_value=2)
    s2 = _check_real_scalar(sigma2, "sigma2", min_value=0.0, strict_min=True)
    idx = np.arange(n, dtype=np.float64)
    f_unitary = np.exp(-1j * 2.0 * math.pi * np.outer(idx, idx) / n) / math.sqrt(n)
    # Covariance E[X X^H] = s2 * F F^H = s2 * I; Pseudo-covariance E[X X^T] = s2 * F F^T
    cov_diag = np.diag(s2 * (f_unitary @ f_unitary.conj().T)).real
    pseudo_diag = np.diag(s2 * (f_unitary @ f_unitary.T)).real
    return {
        "n_fft": n,
        "variance_per_bin": cov_diag,
        "pseudo_variance_per_bin": pseudo_diag,
        "dc_pseudo_variance": float(pseudo_diag[0]),
        "nyquist_pseudo_variance": float(pseudo_diag[n // 2]) if n % 2 == 0 else None,
        "interior_max_abs_pseudo_variance": float(
            np.max(np.abs([pseudo_diag[k] for k in range(n) if k != 0 and (n % 2 != 0 or k != n // 2)]))
        )
        if n > 2
        else 0.0,
    }


def narrowband_tdoa_crlb(
    freq_hz: float,
    snr_linear: float,
    num_snapshots: int = 1,
) -> dict:
    """Fisher information and CRLB for 2-microphone narrowband TDOA estimation x = s_0 [1, e^{-j2pi f tau}]^T + v
    comparing known source phase vs unknown deterministic source phase (Schur reduction)."""
    f = _check_real_scalar(freq_hz, "freq_hz", min_value=0.0, strict_min=True)
    snr = _check_real_scalar(snr_linear, "snr_linear", min_value=0.0, strict_min=True)
    snaps = _check_int(num_snapshots, "num_snapshots", min_value=1)

    omega = 2.0 * math.pi * f
    fim_tau_known_phase = 2.0 * snaps * (omega ** 2) * snr
    # Joint 2x2 FIM wrt [tau, varphi]^T: [[2 L w^2 SNR, -2 L w SNR], [-2 L w SNR, 4 L SNR]]
    fim_joint = np.array(
        [
            [2.0 * snaps * (omega ** 2) * snr, -2.0 * snaps * omega * snr],
            [-2.0 * snaps * omega * snr, 4.0 * snaps * snr],
        ],
        dtype=np.float64,
    )
    fim_eff_unknown_phase = float(fim_joint[0, 0] - (fim_joint[0, 1] ** 2) / fim_joint[1, 1])
    crlb_var_known = 1.0 / fim_tau_known_phase
    crlb_var_unknown = 1.0 / fim_eff_unknown_phase
    return {
        "freq_hz": f,
        "snr_linear": snr,
        "num_snapshots": snaps,
        "fim_known_phase": fim_tau_known_phase,
        "fim_unknown_phase_schur": fim_eff_unknown_phase,
        "crlb_std_known_phase_us": math.sqrt(crlb_var_known) * 1e6,
        "crlb_std_unknown_phase_us": math.sqrt(crlb_var_unknown) * 1e6,
        "unknown_vs_known_variance_ratio": crlb_var_unknown / crlb_var_known,
    }


def proper_complex_gaussian_and_crlb_example() -> dict:
    """E00-13: Proper vs Improper complex Gaussian pseudo-covariance, DFT DC/Nyquist exception, and TDOA CRLB."""
    proper = complex_second_order_moments(0.5, 0.5, 0.0)
    improper = complex_second_order_moments(0.8, 0.2, 0.5)
    dft_check = dft_real_white_noise_pseudo_covariance(n_fft=8, sigma2=1.0)
    crlb = narrowband_tdoa_crlb(freq_hz=1000.0, snr_linear=10.0, num_snapshots=25)
    return {
        "proper_circular_case": proper,
        "improper_elliptical_case": improper,
        "real_dft_8pt_propriety": dft_check,
        "narrowband_tdoa_crlb": crlb,
    }


def wiener_lmmse_filter(
    r_xx: np.ndarray,
    r_xs: np.ndarray,
    var_s: float,
) -> dict:
    r = _check_hermitian(r_xx, "r_xx", strictly_pd=True)
    p = _check_complex_array(r_xs, "r_xs", shape=(r.shape[0],))
    phi_s = _check_real_scalar(var_s, "var_s", min_value=0.0, strict_min=True)

    w_opt = np.linalg.solve(r, p)
    mmse = float(phi_s - np.vdot(w_opt, p).real)
    # Cross-covariance between estimation error e = s - w^H x and observation x^H: E[e x^H] = r_xs^H - w^H R_xx
    orth_residual = p.conj() - (w_opt.conj() @ r)
    return {
        "wiener_weights": w_opt,
        "mmse": mmse,
        "orthogonality_residual_norm": float(np.linalg.norm(orth_residual)),
    }


def kalman_joseph_step(
    f_mat: np.ndarray,
    h_mat: np.ndarray,
    q_cov: np.ndarray,
    r_cov: np.ndarray,
    s_prev: np.ndarray,
    p_prev: np.ndarray,
    z_obs: np.ndarray,
) -> dict:
    f = _check_real_array(f_mat, "f_mat", ndim=2)
    n_dim = f.shape[0]
    if f.shape != (n_dim, n_dim):
        raise ValueError("f_mat must be square")
    h = _check_real_array(h_mat, "h_mat", ndim=2)
    m_dim = h.shape[0]
    if h.shape[1] != n_dim:
        raise ValueError("h_mat column count must match state dimension")
    q = _check_hermitian(q_cov, "q_cov", psd=True).real
    r = _check_hermitian(r_cov, "r_cov", strictly_pd=True).real
    s0 = _check_real_array(s_prev, "s_prev", shape=(n_dim,))
    p0 = _check_hermitian(p_prev, "p_prev", psd=True).real
    z = _check_real_array(z_obs, "z_obs", shape=(m_dim,))
    if q.shape != (n_dim, n_dim) or p0.shape != (n_dim, n_dim) or r.shape != (m_dim, m_dim):
        raise ValueError("covariance dimensions mismatch state/measurement sizes")

    s_pred = f @ s0
    p_pred = f @ p0 @ f.T + q
    innovation = z - h @ s_pred
    s_innov = h @ p_pred @ h.T + r
    k_gain = np.linalg.solve(s_innov.T, (p_pred @ h.T).T).T
    s_post = s_pred + k_gain @ innovation

    i_kh = np.eye(n_dim, dtype=np.float64) - k_gain @ h
    p_joseph = i_kh @ p_pred @ i_kh.T + k_gain @ r @ k_gain.T
    p_standard = i_kh @ p_pred
    return {
        "predicted_state": s_pred,
        "predicted_cov": p_pred,
        "innovation": innovation,
        "innovation_cov": s_innov,
        "kalman_gain": k_gain,
        "posterior_state": s_post,
        "posterior_cov_joseph": p_joseph,
        "joseph_vs_standard_fro_diff": float(np.linalg.norm(p_joseph - p_standard, ord="fro")),
        "posterior_min_eigenvalue": float(np.min(np.linalg.eigvalsh(0.5 * (p_joseph + p_joseph.T)))),
    }


def wiener_orthogonality_and_kalman_joseph_example() -> dict:
    """E00-14: LMMSE Wiener orthogonality principle, MVDR+postfilter identity, and Kalman Joseph covariance update."""
    phi_s = 2.0
    a_vec = np.array([1.0, 1.0], dtype=np.complex128)
    r_vv = np.array([[1.0, 0.5], [0.5, 1.0]], dtype=np.complex128)
    r_xx = phi_s * np.outer(a_vec, a_vec.conj()) + r_vv
    r_xs = phi_s * a_vec

    wiener = wiener_lmmse_filter(r_xx, r_xs, phi_s)
    w_mvdr = solve_lcmv_weights(r_vv, a_vec[:, None], np.array([1.0 + 0.0j]))
    sigma2_mvdr = float(np.vdot(w_mvdr, r_vv @ w_mvdr).real)
    sinr_mvdr = phi_s / sigma2_mvdr
    scalar_post = sinr_mvdr / (1.0 + sinr_mvdr)
    mvdr_post_diff = float(np.linalg.norm(wiener["wiener_weights"] - scalar_post * w_mvdr))

    f_mat = np.array([[1.0, 1.0], [0.0, 1.0]], dtype=np.float64)
    h_mat = np.array([[1.0, 0.0]], dtype=np.float64)
    q_cov = np.diag([0.0, 1.0])
    r_cov = np.array([[3.0]], dtype=np.float64)
    s0 = np.array([10.0, 2.0], dtype=np.float64)
    p0 = np.diag([4.0, 1.0])
    z1 = np.array([13.5], dtype=np.float64)

    kf = kalman_joseph_step(f_mat, h_mat, q_cov, r_cov, s0, p0, z1)

    # Perturbed/quantized gain check: true MSE covariance equals Joseph formula, NOT (I - K_tilde H) P_pred
    k_perturbed = np.array([[0.75], [0.20]], dtype=np.float64)
    i_kt = np.eye(2) - k_perturbed @ h_mat
    true_mse_cov = i_kt @ kf["predicted_cov"] @ i_kt.T + k_perturbed @ r_cov @ k_perturbed.T
    naive_cov = i_kt @ kf["predicted_cov"]

    return {
        "wiener": {
            **wiener,
            "mvdr_weights": w_mvdr,
            "mvdr_sinr": sinr_mvdr,
            "scalar_wiener_postfilter": scalar_post,
            "mvdr_postfilter_equivalence_error": mvdr_post_diff,
        },
        "kalman": {
            **kf,
            "perturbed_gain_naive_asymmetry": float(np.linalg.norm(naive_cov - naive_cov.T, ord="fro")),
            "perturbed_gain_naive_cov_error": float(np.linalg.norm(true_mse_cov - naive_cov, ord="fro")),
        },
    }


def hungarian_linear_sum_assignment(cost_matrix: np.ndarray) -> dict:
    """Exact O(K^3) Kuhn-Munkres / Jonker-Volgenant shortest-augmenting-path solver
    for a square K x K real cost matrix, returning optimal assignment and dual potentials."""
    cost = _check_real_array(cost_matrix, "cost_matrix", ndim=2)
    k = cost.shape[0]
    if cost.shape[1] != k:
        raise ValueError(f"cost_matrix must be square, got {cost.shape}")

    # 1-based indexing classic shortest augmenting path algorithm
    u = np.zeros(k + 1, dtype=np.float64)
    v = np.zeros(k + 1, dtype=np.float64)
    p = np.zeros(k + 1, dtype=np.int64)
    way = np.zeros(k + 1, dtype=np.int64)

    for i in range(1, k + 1):
        p[0] = i
        j0 = 0
        minv = np.full(k + 1, np.inf, dtype=np.float64)
        used = np.zeros(k + 1, dtype=bool)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = np.inf
            j1 = 0
            for j in range(1, k + 1):
                if not used[j]:
                    cur = float(cost[i0 - 1, j - 1]) - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            for j in range(k + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break

    col_for_row = np.zeros(k, dtype=np.int64)
    for j in range(1, k + 1):
        col_for_row[p[j] - 1] = j - 1

    dual_u = u[1:].copy()
    dual_v = v[1:].copy()
    total_cost = float(sum(cost[r, col_for_row[r]] for r in range(k)))
    dual_sum = float(np.sum(dual_u) + np.sum(dual_v))
    slack_min = float(np.min(cost - (dual_u[:, None] + dual_v[None, :])))
    return {
        "row_ind": np.arange(k, dtype=np.int64),
        "col_ind": col_for_row,
        "total_cost": total_cost,
        "dual_u": dual_u,
        "dual_v": dual_v,
        "dual_objective_sum": dual_sum,
        "min_reduced_cost": slack_min,
        "duality_gap": abs(total_cost - dual_sum),
    }


def greedy_bipartite_assignment(cost_matrix: np.ndarray) -> dict:
    cost = _check_real_array(cost_matrix, "cost_matrix", ndim=2)
    k = cost.shape[0]
    if cost.shape[1] != k:
        raise ValueError("cost_matrix must be square")
    used_rows: set[int] = set()
    used_cols: set[int] = set()
    col_for_row = np.full(k, -1, dtype=np.int64)
    for _ in range(k):
        best_val = np.inf
        best_pair = (-1, -1)
        for r in range(k):
            if r in used_rows:
                continue
            for c in range(k):
                if c in used_cols:
                    continue
                if cost[r, c] < best_val:
                    best_val = float(cost[r, c])
                    best_pair = (r, c)
        used_rows.add(best_pair[0])
        used_cols.add(best_pair[1])
        col_for_row[best_pair[0]] = best_pair[1]
    return {
        "col_ind": col_for_row,
        "total_cost": float(sum(cost[r, col_for_row[r]] for r in range(k))),
    }


def hungarian_bipartite_assignment_example() -> dict:
    """E00-15: Permutation-invariant cost matrix, greedy assignment suboptimality, and Hungarian strong duality."""
    cost = np.array(
        [
            [1.0, 2.0, 9.0],
            [2.0, 9.0, 9.0],
            [9.0, 3.0, 4.0],
        ],
        dtype=np.float64,
    )
    greedy = greedy_bipartite_assignment(cost)
    opt = hungarian_linear_sum_assignment(cost)
    perm = np.array([2, 0, 1], dtype=np.int64)
    opt_permuted = hungarian_linear_sum_assignment(cost[perm][:, perm])
    return {
        "cost_matrix": cost,
        "greedy_assignment": greedy["col_ind"],
        "greedy_total_cost": greedy["total_cost"],
        "hungarian_assignment": opt["col_ind"],
        "hungarian_total_cost": opt["total_cost"],
        "cost_saving_over_greedy": greedy["total_cost"] - opt["total_cost"],
        "dual_u": opt["dual_u"],
        "dual_v": opt["dual_v"],
        "min_reduced_cost": opt["min_reduced_cost"],
        "duality_gap": opt["duality_gap"],
        "permuted_optimal_cost": opt_permuted["total_cost"],
        "permutation_invariance_error": abs(opt["total_cost"] - opt_permuted["total_cost"]),
    }


# ---------------------------------------------------------------------------
# §0.7 Acoustic special functions: Bessel, spherical Bessel, Legendre & Y_n^m
# ---------------------------------------------------------------------------

def cylindrical_bessel_j(order: int, x) -> np.ndarray | float:
    """Pure-NumPy cylindrical Bessel function J_n(x) for integer order n."""
    n = _check_int(order, "order")
    scalar_input = np.ndim(x) == 0
    arr = _check_real_array([x] if scalar_input else x, "x")
    if np.max(np.abs(arr)) > 30.0:
        raise ValueError("cylindrical_bessel_j series/quadrature domain requires |x| <= 30")

    m_abs = abs(n)
    if float(np.max(np.abs(arr))) <= 6.0 and m_abs <= 40:
        half_x = 0.5 * arr
        total = np.zeros_like(arr, dtype=np.float64)
        # Power series J_m(x) = sum_{k=0}^{55} (-1)^k / (k! (m+k)!) (x/2)^{2k+m}
        for k in range(56):
            denom = float(math.factorial(k) * math.factorial(m_abs + k))
            term = ((-1.0) ** k / denom) * (half_x ** (2 * k + m_abs))
            total += term
            if float(np.max(np.abs(term))) <= 1e-17 * max(1.0, float(np.max(np.abs(total)))):
                break
    else:
        # Periodic integral J_m(x) = (1/pi) \int_0^pi cos(m alpha - x sin alpha) d alpha
        n_quad = max(128, 2 * m_abs + 64)
        alpha = (np.arange(n_quad, dtype=np.float64) + 0.5) * (math.pi / n_quad)
        total = np.cos(m_abs * alpha - arr[..., None] * np.sin(alpha)).mean(axis=-1)
    if n < 0 and (m_abs % 2 == 1):
        total = -total
    return float(total[0]) if scalar_input else total


def spherical_bessel_jy(order: int, x) -> tuple[np.ndarray | float, np.ndarray | float]:
    """Pure-NumPy spherical Bessel functions j_n(x) and y_n(x) for non-negative integer n and x > 0."""
    n = _check_int(order, "order", min_value=0, max_value=12)
    scalar_input = np.ndim(x) == 0
    arr = _check_real_array([x] if scalar_input else x, "x")
    if np.any(arr <= 0.0):
        raise ValueError("spherical_bessel_jy requires strictly positive x > 0")

    y0 = -np.cos(arr) / arr
    if n == 0:
        j_out = np.sin(arr) / arr
        y_out = y0
    else:
        y1 = -np.cos(arr) / (arr * arr) - np.sin(arr) / arr
        if n == 1:
            y_out = y1
        else:
            ym1, ym0 = y0, y1
            for l_idx in range(1, n):
                yp1 = ((2.0 * l_idx + 1.0) / arr) * ym0 - ym1
                ym1, ym0 = ym0, yp1
            y_out = ym0

        # Power series for small/moderate x avoids 1/x^n subtraction cancellation;
        # forward-stable upward recurrence for x > max(3.5, n) avoids alternating-series cancellation.
        cutoff = max(3.5, float(n))
        small_mask = arr <= cutoff
        j_out = np.zeros_like(arr, dtype=np.float64)
        if np.all(small_mask):
            for k in range(56):
                double_fact = math.factorial(2 * n + 2 * k + 1) // (
                    (2 ** (n + k)) * math.factorial(n + k)
                )
                denom = float((2 ** k) * math.factorial(k) * double_fact)
                term = ((-1.0) ** k / denom) * (arr ** (2 * k + n))
                j_out += term
                if float(np.max(np.abs(term))) <= 1e-17 * max(1.0, float(np.max(np.abs(j_out)))):
                    break
        else:
            if np.any(small_mask):
                arr_s = arr[small_mask]
                tot_s = np.zeros_like(arr_s, dtype=np.float64)
                for k in range(56):
                    double_fact = math.factorial(2 * n + 2 * k + 1) // (
                        (2 ** (n + k)) * math.factorial(n + k)
                    )
                    denom = float((2 ** k) * math.factorial(k) * double_fact)
                    term = ((-1.0) ** k / denom) * (arr_s ** (2 * k + n))
                    tot_s += term
                    if float(np.max(np.abs(term))) <= 1e-17 * max(1.0, float(np.max(np.abs(tot_s)))):
                        break
                j_out[small_mask] = tot_s
            large_mask = ~small_mask
            arr_l = arr[large_mask]
            j0_l = np.sin(arr_l) / arr_l
            j1_l = np.sin(arr_l) / (arr_l * arr_l) - np.cos(arr_l) / arr_l
            if n == 1:
                j_out[large_mask] = j1_l
            else:
                jm1, jm0 = j0_l, j1_l
                for l_idx in range(1, n):
                    jp1 = ((2.0 * l_idx + 1.0) / arr_l) * jm0 - jm1
                    jm1, jm0 = jm0, jp1
                j_out[large_mask] = jm0

    if scalar_input:
        return float(j_out[0]), float(y_out[0])
    return j_out, y_out


def spherical_bessel_derivative(order: int, x) -> tuple[np.ndarray | float, np.ndarray | float]:
    """Derivatives j_n'(x) and y_n'(x) via exact recurrence f_n'(x) = (n/x) f_n(x) - f_{n+1}(x)."""
    n = _check_int(order, "order", min_value=0, max_value=11)
    scalar_input = np.ndim(x) == 0
    arr = _check_real_array([x] if scalar_input else x, "x")
    jn, yn = spherical_bessel_jy(n, arr)
    jnp1, ynp1 = spherical_bessel_jy(n + 1, arr)
    dj = (n / arr) * jn - jnp1
    dy = (n / arr) * yn - ynp1
    if scalar_input:
        return float(dj[0]), float(dy[0])
    return dj, dy


def legendre_polynomial(order: int, x) -> np.ndarray | float:
    """Bonnet three-term recurrence for Legendre polynomial P_n(x) on [-1, 1]."""
    n = _check_int(order, "order", min_value=0)
    scalar_input = np.ndim(x) == 0
    arr = _check_real_array([x] if scalar_input else x, "x")
    if np.any(np.abs(arr) > 1.0 + 1e-15):
        raise ValueError("Legendre polynomial argument x must lie in [-1, 1]")
    arr = np.clip(arr, -1.0, 1.0)
    p0 = np.ones_like(arr)
    if n == 0:
        return float(p0[0]) if scalar_input else p0
    p1 = arr.copy()
    if n == 1:
        return float(p1[0]) if scalar_input else p1
    pm1, pm0 = p0, p1
    for k in range(1, n):
        pp1 = ((2 * k + 1) * arr * pm0 - k * pm1) / (k + 1)
        pm1, pm0 = pm0, pp1
    return float(pm0[0]) if scalar_input else pm0


def associated_legendre_pnm(order: int, degree: int, x) -> np.ndarray | float:
    """Associated Legendre function P_n^m(x) for 0 <= m <= n on [-1, 1] with Condon-Shortley phase (-1)^m."""
    n = _check_int(order, "order", min_value=0)
    m = _check_int(degree, "degree", min_value=0, max_value=n)
    scalar_input = np.ndim(x) == 0
    arr = _check_real_array([x] if scalar_input else x, "x")
    if np.any(np.abs(arr) > 1.0 + 1e-15):
        raise ValueError("associated_legendre_pnm requires x in [-1, 1]")
    arr = np.clip(arr, -1.0, 1.0)

    # Sectoral seed P_m^m(x) = (-1)^m (2m-1)!! (1 - x^2)^{m/2}
    double_fact = 1.0
    for k in range(1, m + 1):
        double_fact *= (2 * k - 1)
    pmm = ((-1.0) ** m) * double_fact * ((1.0 - arr * arr) ** (0.5 * m))
    if n == m:
        return float(pmm[0]) if scalar_input else pmm
    pmmp1 = arr * (2 * m + 1) * pmm
    if n == m + 1:
        return float(pmmp1[0]) if scalar_input else pmmp1
    prev2, prev1 = pmm, pmmp1
    for l_idx in range(m + 2, n + 1):
        curr = ((2 * l_idx - 1) * arr * prev1 - (l_idx + m - 1) * prev2) / (l_idx - m)
        prev2, prev1 = prev1, curr
    return float(prev1[0]) if scalar_input else prev1


def complex_spherical_harmonic_ynm(
    order: int,
    degree: int,
    colatitude_rad: float,
    azimuth_rad: float,
) -> complex:
    """Orthonormal complex spherical harmonic Y_n^m(vartheta, alpha) for -n <= m <= n."""
    n = _check_int(order, "order", min_value=0)
    m = _check_int(degree, "degree", min_value=-n, max_value=n)
    vartheta = _check_real_scalar(colatitude_rad, "colatitude_rad", min_value=0.0, max_value=math.pi)
    alpha = _check_real_scalar(azimuth_rad, "azimuth_rad")

    m_abs = abs(m)
    norm = math.sqrt(((2 * n + 1) / (4.0 * math.pi)) * (math.factorial(n - m_abs) / math.factorial(n + m_abs)))
    pnm = float(associated_legendre_pnm(n, m_abs, math.cos(vartheta)))
    y_pos = norm * pnm * complex(math.cos(m_abs * alpha), math.sin(m_abs * alpha))
    if m >= 0:
        return y_pos
    return ((-1) ** m_abs) * y_pos.conjugate()


def bessel_legendre_and_spherical_harmonics_example() -> dict:
    """E00-16: Cylindrical/spherical Bessel recurrences & rigid-sphere Wronskian, Legendre orthogonality, and Y_n^m addition theorem."""
    x_val = 2.0
    j0 = cylindrical_bessel_j(0, x_val)
    j1 = cylindrical_bessel_j(1, x_val)
    j2 = cylindrical_bessel_j(2, x_val)
    recurrence_err = abs((j0 + j2) - (2.0 * 1.0 / x_val) * j1)

    # Jacobi-Anger expansion check: exp(j * kr * cos(phi)) = sum_{m=-M}^{M} j^m J_m(kr) exp(j m phi)
    kr = 1.8
    phi_ja = 0.65
    ja_sum = sum(
        (1j ** m) * cylindrical_bessel_j(m, kr) * complex(math.cos(m * phi_ja), math.sin(m * phi_ja))
        for m in range(-18, 19)
    )
    ja_exact = complex(math.cos(kr * math.cos(phi_ja)), math.sin(kr * math.cos(phi_ja)))

    # Open vs rigid sphere radial weights at kr = pi (first open-sphere null!)
    kr_null = math.pi
    sj0, sy0 = spherical_bessel_jy(0, kr_null)
    dj0, dy0 = spherical_bessel_derivative(0, kr_null)
    wronskian_val = sj0 * dy0 - dj0 * sy0
    wronskian_exact = 1.0 / (kr_null * kr_null)
    h0_2 = complex(sj0, -sy0)
    dh0_2 = complex(dj0, -dy0)
    b0_open_pi = 4.0 * math.pi * sj0
    b0_rigid_pi = 4.0 * math.pi * (sj0 - (dj0 / dh0_2) * h0_2)

    # Exact 5-point Gauss-Legendre quadrature on [-1, 1] integrates polynomials up to degree 9
    s_inner = math.sqrt(5.0 - 2.0 * math.sqrt(10.0 / 7.0)) / 3.0
    s_outer = math.sqrt(5.0 + 2.0 * math.sqrt(10.0 / 7.0)) / 3.0
    w_inner = (322.0 + 13.0 * math.sqrt(70.0)) / 900.0
    w_outer = (322.0 - 13.0 * math.sqrt(70.0)) / 900.0
    gl_nodes = np.array([-s_outer, -s_inner, 0.0, s_inner, s_outer], dtype=np.float64)
    gl_weights = np.array([w_outer, w_inner, 128.0 / 225.0, w_inner, w_outer], dtype=np.float64)

    gram_leg = np.zeros((4, 4), dtype=np.float64)
    for m_idx in range(4):
        pm = legendre_polynomial(m_idx, gl_nodes)
        for n_idx in range(4):
            pn = legendre_polynomial(n_idx, gl_nodes)
            gram_leg[m_idx, n_idx] = float(np.sum(gl_weights * pm * pn))
    expected_leg_diag = np.diag([2.0 / (2 * n_idx + 1) for n_idx in range(4)])
    odd_parity_err = abs(legendre_polynomial(3, -0.6) + legendre_polynomial(3, 0.6))

    # Spherical Harmonic Addition Theorem check for orders n = 0..3
    vt1, az1 = 0.70, 0.45
    vt2, az2 = 1.25, -0.85
    u1 = np.array(
        [math.sin(vt1) * math.cos(az1), math.sin(vt1) * math.sin(az1), math.cos(vt1)],
        dtype=np.float64,
    )
    u2 = np.array(
        [math.sin(vt2) * math.cos(az2), math.sin(vt2) * math.sin(az2), math.cos(vt2)],
        dtype=np.float64,
    )
    cos_gamma = float(np.dot(u1, u2))
    addition_errors = []
    for n_ord in range(4):
        lhs = sum(
            complex_spherical_harmonic_ynm(n_ord, m_deg, vt1, az1)
            * complex_spherical_harmonic_ynm(n_ord, m_deg, vt2, az2).conjugate()
            for m_deg in range(-n_ord, n_ord + 1)
        )
        rhs = ((2 * n_ord + 1) / (4.0 * math.pi)) * legendre_polynomial(n_ord, cos_gamma)
        addition_errors.append(abs(lhs - rhs))

    return {
        "bessel_j0_j1_j2_at_2": [j0, j1, j2],
        "bessel_three_term_recurrence_error": recurrence_err,
        "jacobi_anger_abs_error": abs(ja_sum - ja_exact),
        "wronskian_at_pi": wronskian_val,
        "wronskian_identity_abs_error": abs(wronskian_val - wronskian_exact),
        "open_sphere_b0_abs_at_pi": abs(b0_open_pi),
        "rigid_sphere_b0_abs_at_pi": abs(b0_rigid_pi),
        "rigid_sphere_b0_closed_form_abs_at_pi": 4.0 * math.pi / ((kr_null ** 2) * abs(dh0_2)),
        "legendre_gram_matrix_0_to_3": gram_leg,
        "legendre_orthogonality_fro_error": float(np.linalg.norm(gram_leg - expected_leg_diag, ord="fro")),
        "legendre_p3_odd_parity_error": odd_parity_err,
        "spherical_harmonic_addition_max_error": max(addition_errors),
    }


# ---------------------------------------------------------------------------
# Additional Ch00 foundation primitives (Dirichlet, R^{1/2}, matrix Wirtinger,
# Gauss-Markov GLS, Wilson CI / zero-failure bound, real FOA spherical harmonics)
# ---------------------------------------------------------------------------

def dirichlet_ula_array_factor(num_mics: int, psi_rad: float) -> dict:
    """Compare direct M-element geometric sum sum_{m=0}^{M-1} exp(-j m psi)
    against the closed-form Dirichlet kernel exp(-j (M-1) psi / 2) * sin(M psi / 2) / sin(psi / 2)."""
    m_cnt = _check_int(num_mics, "num_mics", min_value=1)
    psi = _check_real_scalar(psi_rad, "psi_rad")

    idx = np.arange(m_cnt, dtype=np.float64)
    direct = complex(np.sum(np.exp(-1j * idx * psi)))
    half_psi = 0.5 * psi
    sin_half = math.sin(half_psi)
    if abs(sin_half) <= 1e-12:
        # At psi = 2 k pi, every term exp(-j m 2 k pi) = 1, so both sum and limit equal +M
        closed = complex(float(m_cnt), 0.0)
    else:
        phase = complex(math.cos(-0.5 * (m_cnt - 1) * psi), math.sin(-0.5 * (m_cnt - 1) * psi))
        closed = phase * (math.sin(0.5 * m_cnt * psi) / sin_half)
    return {
        "num_mics": m_cnt,
        "psi_rad": psi,
        "direct_sum": direct,
        "dirichlet_closed_form": closed,
        "abs_error": abs(direct - closed),
        "normalized_magnitude": abs(direct) / m_cnt,
    }


def hermitian_matrix_sqrt(mat: np.ndarray) -> dict:
    """Compute Hermitian positive-definite square root R^{1/2} = U Lambda^{1/2} U^H
    and its inverse R^{-1/2} = U Lambda^{-1/2} U^H."""
    r = _check_hermitian(mat, "mat", strictly_pd=True)
    eigvals, eigvecs = np.linalg.eigh(r)
    sqrt_eigs = np.sqrt(eigvals)
    r_sqrt = (eigvecs * sqrt_eigs[None, :]) @ eigvecs.conj().T
    r_inv_sqrt = (eigvecs * (1.0 / sqrt_eigs)[None, :]) @ eigvecs.conj().T
    return {
        "r_sqrt": r_sqrt,
        "r_inv_sqrt": r_inv_sqrt,
        "eigenvalues": eigvals,
        "square_fro_error": float(np.linalg.norm(r_sqrt @ r_sqrt - r, ord="fro")),
        "whitening_fro_error": float(
            np.linalg.norm(r_inv_sqrt @ r @ r_inv_sqrt - np.eye(r.shape[0]), ord="fro")
        ),
    }


def matrix_wirtinger_logdet_check(w_mat: np.ndarray, fd_eps: float = 1e-6) -> dict:
    """Verify matrix Wirtinger log-determinant gradient:
    d ln|det(W)| / d W^* = 0.5 * W^{-H} against entry-wise real/imag central finite differences."""
    w = _check_complex_array(w_mat, "w_mat", ndim=2)
    m = w.shape[0]
    if w.shape != (m, m):
        raise ValueError(f"w_mat must be square, got {w.shape}")
    _require_nonsingular(w, "w_mat")
    eps = _check_real_scalar(fd_eps, "fd_eps", min_value=1e-10, max_value=1e-2)

    def f_logabsdet(mat: np.ndarray) -> float:
        return float(math.log(abs(np.linalg.det(mat))))

    grad_analytic = 0.5 * np.linalg.inv(w).conj().T
    grad_fd = np.zeros_like(w)
    for row in range(m):
        for col in range(m):
            e_rc = np.zeros_like(w)
            e_rc[row, col] = 1.0
            df_dr = (f_logabsdet(w + eps * e_rc) - f_logabsdet(w - eps * e_rc)) / (2.0 * eps)
            df_di = (f_logabsdet(w + 1j * eps * e_rc) - f_logabsdet(w - 1j * eps * e_rc)) / (2.0 * eps)
            grad_fd[row, col] = 0.5 * (df_dr + 1j * df_di)

    return {
        "analytic_gradient": grad_analytic,
        "fd_gradient": grad_fd,
        "max_abs_error": float(np.max(np.abs(grad_analytic - grad_fd))),
    }


def gauss_markov_gls_solve(
    a_mat: np.ndarray,
    b_vec: np.ndarray,
    noise_cov: np.ndarray,
) -> dict:
    """Generalized / Weighted Least Squares (Gauss-Markov BLUE) via Cholesky whitening C = L L^H."""
    a = _check_complex_array(a_mat, "a_mat", ndim=2)
    m, p = a.shape
    if m < p:
        raise ValueError("a_mat must have m >= p with full column rank")
    b = _check_complex_array(b_vec, "b_vec", shape=(m,))
    c_nn = _check_hermitian(noise_cov, "noise_cov", strictly_pd=True)
    if c_nn.shape != (m, m):
        raise ValueError("noise_cov dimension must match a_mat row count")

    chol_l = np.linalg.cholesky(c_nn)
    a_tilde = np.linalg.solve(chol_l, a)
    b_tilde = np.linalg.solve(chol_l, b)

    gram_gls = _check_hermitian(a_tilde.conj().T @ a_tilde, "A^H C^{-1} A", strictly_pd=True)
    x_gls = np.linalg.solve(gram_gls, a_tilde.conj().T @ b_tilde)
    cov_gls = np.linalg.inv(gram_gls)

    # Direct metric orthogonality: A^H C^{-1} (A x_gls - b) == 0
    resid_gls = a @ x_gls - b
    weighted_orth = a.conj().T @ np.linalg.solve(c_nn, resid_gls)

    # Compare against OLS covariance (A^H A)^{-1} A^H C A (A^H A)^{-1}
    gram_ols = _check_hermitian(a.conj().T @ a, "A^H A", strictly_pd=True)
    p_ols = np.linalg.solve(gram_ols, a.conj().T)
    x_ols = p_ols @ b
    cov_ols = p_ols @ c_nn @ p_ols.conj().T
    cov_diff = 0.5 * ((cov_ols - cov_gls) + (cov_ols - cov_gls).conj().T)

    return {
        "cholesky_l": chol_l,
        "x_gls": x_gls,
        "x_ols": x_ols,
        "cov_gls": cov_gls,
        "cov_ols": cov_ols,
        "weighted_orthogonality_norm": float(np.linalg.norm(weighted_orth)),
        "cov_ols_minus_gls_min_eig": float(np.min(np.linalg.eigvalsh(cov_diff))),
    }


def binomial_zero_failure_upper_bound(n_trials: int, alpha: float = 0.05) -> float:
    """Exact one-sided (1 - alpha) upper bound p_U = 1 - alpha^(1/n) when 0 failures occur in n trials."""
    n = _check_int(n_trials, "n_trials", min_value=1)
    a = _check_real_scalar(alpha, "alpha", min_value=0.0, strict_min=True, max_value=1.0, strict_max=True)
    return float(1.0 - (a ** (1.0 / n)))


def wilson_score_interval(
    k_success: int,
    n_trials: int,
    z_value: float = 1.959963984540054,
) -> dict:
    """Wilson score confidence interval for binomial proportion k / n."""
    n = _check_int(n_trials, "n_trials", min_value=1)
    k = _check_int(k_success, "k_success", min_value=0, max_value=n)
    z = _check_real_scalar(z_value, "z_value", min_value=0.0, strict_min=True)

    p_hat = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (p_hat + z2 / (2.0 * n)) / denom
    half_width = (z / denom) * math.sqrt((p_hat * (1.0 - p_hat) / n) + (z2 / (4.0 * n * n)))
    return {
        "p_hat": p_hat,
        "lower": max(0.0, center - half_width),
        "upper": min(1.0, center + half_width),
        "center": center,
        "half_width": half_width,
    }


def real_spherical_harmonic_foa(colatitude_rad: float, azimuth_rad: float) -> dict:
    """Real first-order spherical harmonics (ACN 0..3: W, Y, Z, X) in SN3D, N3D, and orthonormal conventions.

    Angles ``(colatitude_rad, azimuth_rad)`` are in the listener-local AmbiX/SOFA
    frame (``+x`` front, ``+y`` left, ``+z`` up), with colatitude measured from ``+z``
    and azimuth measured counterclockwise from ``+x`` toward ``+y``.
    """
    vt = _check_real_scalar(colatitude_rad, "colatitude_rad", min_value=0.0, max_value=math.pi)
    az = _check_real_scalar(azimuth_rad, "azimuth_rad")
    nx = math.sin(vt) * math.cos(az)
    ny = math.sin(vt) * math.sin(az)
    nz = math.cos(vt)
    # ACN ordering: 0=(0,0) W, 1=(1,-1) Y, 2=(1,0) Z, 3=(1,+1) X
    sn3d = np.array([1.0, ny, nz, nx], dtype=np.float64)
    n3d = np.array([1.0, math.sqrt(3.0) * ny, math.sqrt(3.0) * nz, math.sqrt(3.0) * nx], dtype=np.float64)
    ortho = n3d / math.sqrt(4.0 * math.pi)
    return {
        "unit_vec_xyz": np.array([nx, ny, nz], dtype=np.float64),
        "acn_sn3d": sn3d,
        "acn_n3d": n3d,
        "acn_orthonormal": ortho,
    }
