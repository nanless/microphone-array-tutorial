"""Conditional target speech extraction (TSE) modulation, dual-state loss, and zero-overlap bypass.

Conventions:
- Real waveforms are 1-D arrays of shape ``(samples,)`` in float64.
- Feature maps for FiLM conditioning have shape ``(channels, frames)`` or
  ``(channels, bins, frames)`` in float64, modulated along axis 0.
- Power and energy ratios in decibels follow ``10 * log10(E_num / E_den)``.
"""

from __future__ import annotations

import math

import numpy as np

from codes.chapters.ch02.core.conventions import finite_real_array, finite_real_scalar


def _finite_nonbool_scalar(value: object, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite real scalar, not bool")
    scalar = finite_real_scalar(value, name)
    if not math.isfinite(scalar):
        raise ValueError(f"{name} must be finite")
    return scalar


def _safe_energy(signal: np.ndarray) -> float:
    """Compute sum of squares with peak scaling to prevent intermediate overflow."""
    peak = float(np.max(np.abs(signal)))
    if peak == 0.0:
        return 0.0
    scaled = signal / peak
    energy = float(np.dot(scaled, scaled)) * (peak * peak)
    if not np.isfinite(energy):
        raise ValueError("signal energy exceeds finite float64 range")
    return energy


def film_conditioning(
    features: np.ndarray,
    gamma: np.ndarray,
    beta: np.ndarray,
) -> np.ndarray:
    """Apply Feature-wise Linear Modulation ``Z = gamma * H + beta`` along channel axis 0.

    Parameters
    ----------
    features:
        Real feature array ``H`` of shape ``(C, T)`` or ``(C, F, T)``.
    gamma:
        Channel scaling vector of shape ``(C,)``.
    beta:
        Channel shift vector of shape ``(C,)``.
    """
    h = finite_real_array(features, "features")
    g = finite_real_array(gamma, "gamma")
    b = finite_real_array(beta, "beta")
    if h.ndim not in (2, 3) or min(h.shape) < 1:
        raise ValueError("features must have shape (C, T) or (C, F, T) with non-empty axes")
    channels = h.shape[0]
    if g.shape != (channels,) or b.shape != (channels,):
        raise ValueError(f"gamma and beta must both have shape ({channels},)")
    expand = (channels,) + (1,) * (h.ndim - 1)
    with np.errstate(over="ignore", invalid="ignore"):
        modulated = g.reshape(expand) * h + b.reshape(expand)
    if not np.all(np.isfinite(modulated)):
        raise ValueError("modulated features exceed finite float64 range")
    return modulated


def inspect_standard_sisdr(s_hat: np.ndarray, s_ref: np.ndarray) -> dict:
    """Diagnose standard SI-SDR and explicitly show why ``s_ref = 0`` is undefined."""
    est = finite_real_array(s_hat, "s_hat")
    ref = finite_real_array(s_ref, "s_ref")
    if est.ndim != 1 or ref.ndim != 1 or est.size < 1 or est.shape != ref.shape:
        raise ValueError("s_hat and s_ref must be 1-D real arrays of equal positive length")
    ref_peak = float(np.max(np.abs(ref)))
    est_peak = float(np.max(np.abs(est)))
    ref_energy = _safe_energy(ref)
    est_energy = _safe_energy(est)
    if ref_peak == 0.0:
        return {
            "ref_energy": 0.0,
            "est_energy": est_energy,
            "projection_undefined": True,
            "standard_sisdr_db": None,
            "reason": "0/0 singularity in <s_hat, s_ref> / ||s_ref||_2^2 when target is absent",
        }
    if est_peak == 0.0:
        return {
            "ref_energy": ref_energy,
            "est_energy": 0.0,
            "projection_undefined": False,
            "standard_sisdr_db": -np.inf,
            "reason": "zero target projection makes log10(0 / ||e_noise||_2^2) = -inf",
        }
    ref_u = ref / ref_peak
    est_u = est / est_peak
    ref_u_energy = float(np.dot(ref_u, ref_u))
    alpha_u = float(np.dot(est_u, ref_u) / ref_u_energy)
    s_proj_u = alpha_u * ref_u
    e_noise_u = est_u - s_proj_u
    proj_u_energy = float(np.dot(s_proj_u, s_proj_u))
    noise_u_energy = float(np.dot(e_noise_u, e_noise_u))
    if proj_u_energy == 0.0:
        return {
            "ref_energy": ref_energy,
            "est_energy": est_energy,
            "projection_undefined": False,
            "standard_sisdr_db": -np.inf,
            "reason": "zero target projection makes log10(0 / ||e_noise||_2^2) = -inf",
        }
    if noise_u_energy == 0.0:
        return {
            "ref_energy": ref_energy,
            "est_energy": est_energy,
            "projection_undefined": False,
            "standard_sisdr_db": np.inf,
            "reason": "exact collinearity makes residual noise energy zero (+inf dB)",
        }
    return {
        "ref_energy": ref_energy,
        "est_energy": est_energy,
        "projection_undefined": False,
        "standard_sisdr_db": float(10.0 * np.log10(proj_u_energy / noise_u_energy)),
        "reason": "well-defined active target reference",
    }


def tse_dual_state_loss(
    s_hat: np.ndarray,
    s_ref: np.ndarray,
    x_mix: np.ndarray,
    eps_s: float = 1e-8,
    tau_floor: float = 1e-3,
    ratio_eps: float = 0.0,
) -> dict:
    """Compute dual-state TSE loss handling both target-present and target-absent segments.

    - Target-present (``||s_ref||_2^2 > eps_s``): ``-SI-SDR(s_hat, s_ref)``.
    - Target-absent (``||s_ref||_2^2 <= eps_s``): mixture-referenced energy-floor penalty
      ``10 * log10((||s_hat||_2^2 + tau_floor * ||x_mix||_2^2) / (tau_floor * ||x_mix||_2^2))``.
    """
    est = finite_real_array(s_hat, "s_hat")
    ref = finite_real_array(s_ref, "s_ref")
    mix = finite_real_array(x_mix, "x_mix")
    eps_s_val = _finite_nonbool_scalar(eps_s, "eps_s")
    tau_val = _finite_nonbool_scalar(tau_floor, "tau_floor")
    ratio_eps_val = _finite_nonbool_scalar(ratio_eps, "ratio_eps")
    if est.ndim != 1 or ref.ndim != 1 or mix.ndim != 1 or est.size < 1:
        raise ValueError("s_hat, s_ref, and x_mix must be non-empty 1-D real arrays")
    if not (est.shape == ref.shape == mix.shape):
        raise ValueError("s_hat, s_ref, and x_mix must have identical shapes")
    if eps_s_val < 0.0:
        raise ValueError("eps_s must be non-negative")
    if tau_val <= 0.0:
        raise ValueError("tau_floor must be strictly positive")
    if ratio_eps_val < 0.0:
        raise ValueError("ratio_eps must be non-negative")

    ref_energy = _safe_energy(ref)
    est_energy = _safe_energy(est)
    mix_energy = _safe_energy(mix)
    sisdr_diag = inspect_standard_sisdr(est, ref)

    if ref_energy > eps_s_val:
        ref_peak = float(np.max(np.abs(ref)))
        ref_u = ref / ref_peak
        ref_u_energy = float(np.dot(ref_u, ref_u))
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            alpha = float(np.dot(est, ref_u) / (ref_peak * ref_u_energy))
            s_proj = alpha * ref
            e_noise = est - s_proj
        if not (np.isfinite(alpha) and np.all(np.isfinite(s_proj)) and np.all(np.isfinite(e_noise))):
            raise ValueError("SI-SDR projection exceeds finite float64 range")
        proj_energy = _safe_energy(s_proj)
        noise_energy = _safe_energy(e_noise)
        num = proj_energy + ratio_eps_val
        den = noise_energy + ratio_eps_val
        if num <= 0.0 or den <= 0.0 or not (np.isfinite(num) and np.isfinite(den)):
            raise ValueError(
                "target-present SI-SDR hits zero projection or zero residual; "
                "pass positive ratio_eps to stabilize degenerate inputs"
            )
        si_sdr_db = float(10.0 * np.log10(num / den))
        if not np.isfinite(si_sdr_db):
            raise ValueError("SI-SDR exceeds finite float64 range")
        return {
            "target_present": True,
            "loss_db": -si_sdr_db,
            "si_sdr_db": si_sdr_db,
            "standard_sisdr_undefined": bool(sisdr_diag["projection_undefined"]),
            "ref_energy": ref_energy,
            "est_energy": est_energy,
            "mix_energy": mix_energy,
            "proj_energy": proj_energy,
            "noise_energy": noise_energy,
            "floor_energy": None,
            "scale_alpha": alpha,
        }

    if mix_energy <= 0.0:
        raise ValueError("target-absent energy-floor penalty requires positive mixture energy")
    with np.errstate(over="ignore", invalid="ignore"):
        floor_energy = tau_val * mix_energy
    if floor_energy <= 0.0 or not np.isfinite(floor_energy):
        raise ValueError("mixture floor energy underflowed or overflowed")
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        ratio_over_floor = est_energy / floor_energy
        loss_db = float(
            10.0 * (np.log10(ratio_over_floor) if ratio_over_floor > 1e100 else np.log10(1.0 + ratio_over_floor))
        )
    if not np.isfinite(loss_db):
        raise ValueError("target-absent loss exceeds finite float64 range")
    return {
        "target_present": False,
        "loss_db": loss_db,
        "si_sdr_db": None,
        "standard_sisdr_undefined": bool(sisdr_diag["projection_undefined"]),
        "ref_energy": ref_energy,
        "est_energy": est_energy,
        "mix_energy": mix_energy,
        "proj_energy": None,
        "noise_energy": None,
        "floor_energy": floor_energy,
        "scale_alpha": None,
    }


def zero_overlap_bypass_decision(
    y_primary: np.ndarray,
    y_secondary: np.ndarray,
    threshold_db: float = -18.0,
) -> dict:
    """Protect single-speaker segments against consonant over-splitting via slot energy ratio.

    When ``10 * log10(||y_secondary||_2^2 / ||y_primary||_2^2) < threshold_db``,
    the secondary separator slot is classified as leaked unvoiced/fricative energy
    rather than an independent concurrent speaker, and ``y_primary + y_secondary``
    is recombined to preserve speech naturalness.
    """
    pri = finite_real_array(y_primary, "y_primary")
    sec = finite_real_array(y_secondary, "y_secondary")
    thresh = _finite_nonbool_scalar(threshold_db, "threshold_db")
    if pri.ndim != 1 or sec.ndim != 1 or pri.size < 1 or pri.shape != sec.shape:
        raise ValueError("y_primary and y_secondary must be 1-D real arrays of equal positive length")

    peak_pri = float(np.max(np.abs(pri)))
    peak_sec = float(np.max(np.abs(sec)))
    if peak_pri <= 0.0:
        raise ValueError("y_primary must have strictly positive energy")

    e_pri = _safe_energy(pri)
    e_sec = _safe_energy(sec)
    if peak_sec == 0.0:
        ratio_db = -np.inf
        bypass = True
    else:
        pri_u = pri / peak_pri
        sec_u = sec / peak_sec
        norm_pri_sq = float(np.dot(pri_u, pri_u))
        norm_sec_sq = float(np.dot(sec_u, sec_u))
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            peak_ratio = peak_sec / peak_pri
        log10_peak_ratio = (
            float(np.log10(peak_ratio))
            if (peak_ratio > 0.0 and np.isfinite(peak_ratio))
            else float(np.log10(peak_sec) - np.log10(peak_pri))
        )
        ratio_db = float(20.0 * log10_peak_ratio + 10.0 * np.log10(norm_sec_sq / norm_pri_sq))
        if not np.isfinite(ratio_db):
            raise ValueError("slot energy ratio exceeds finite float64 range")
        bypass = bool(ratio_db < thresh)

    with np.errstate(over="ignore", invalid="ignore"):
        output = (pri + sec) if bypass else pri.copy()
    if not np.all(np.isfinite(output)):
        raise ValueError("recombined bypass output exceeds finite float64 range")
    return {
        "bypass_separation": bypass,
        "secondary_active": not bypass,
        "energy_ratio_db": ratio_db,
        "threshold_db": thresh,
        "primary_energy": e_pri,
        "secondary_energy": e_sec,
        "output": output,
    }
