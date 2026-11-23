"""Pure-NumPy teaching kernels for convolutive CI-SDR and streaming separation gates.

Covers Chapter 08 exercises E08-25 (full-tail causal Toeplitz CI-SDR vs scalar
SI-SDR and a separately named truncated projection) and E08-26 (narrowband-crossband
vs Mel-subband streaming MAC/s accounting and the -18 dB single-speaker
zero-overlap fricative-splitting bypass gate). No SciPy imports; float64.
"""
from __future__ import annotations

import math
import numbers
import numpy as np


def _finite_real_array(value, name: str = "input") -> np.ndarray:
    original = np.asarray(value)
    if original.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numeric values, not complex, bool or text")
    with np.errstate(over="ignore", invalid="ignore"):
        array = np.asarray(original, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite and representable as float64")
    return array


def _finite_real_scalar(value, name: str = "input") -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a real scalar, not bool")
    array = _finite_real_array(value, name)
    if array.ndim != 0:
        raise ValueError(f"{name} must be a scalar")
    return float(array)


def _integer(value, name: str, *, minimum: int | None = None) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Integral):
        raise ValueError(f"{name} must be an integer")
    result = int(value)
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return result


def causal_fir_toeplitz(
    reference: np.ndarray, filter_length: int, *, full_convolution: bool = True,
) -> np.ndarray:
    """Build ``(N+L-1, L)`` full causal convolution, or an explicit ``(N,L)`` crop.

    ``full_convolution=False`` is only the truncated comparison in E08-25;
    CI-SDR always retains the full convolution tail.
    """
    ref = _finite_real_array(reference, "reference")
    if ref.ndim != 1 or ref.size == 0:
        raise ValueError("reference must be a non-empty 1-D real array")
    taps = _integer(filter_length, "filter_length", minimum=1)
    if taps > ref.size:
        raise ValueError("filter_length must satisfy 1 <= filter_length <= len(reference)")
    if not isinstance(full_convolution, (bool, np.bool_)):
        raise ValueError("full_convolution must be boolean")
    n_samples = ref.size
    rows = n_samples + taps - 1 if full_convolution else n_samples
    toeplitz = np.zeros((rows, taps), dtype=np.float64)
    for lag in range(taps):
        count = min(n_samples, rows - lag)
        toeplitz[lag:lag + count, lag] = ref[:count]
    return toeplitz


def ci_sdr_toeplitz(
    estimate: np.ndarray,
    reference: np.ndarray,
    filter_length: int = 512,
    *,
    zero_mean: bool = False,
    epsilon: float = 1e-12,
) -> dict:
    """Full-tail causal L-tap least-squares projection CI-SDR teaching kernel.

    Matches Chapter 08 ``si_sdr`` relative-energy flooring at ``filter_length=1``
    (caps at ``+/- 10*log10(1/epsilon)`` dB). Validates 1-D finite real arrays,
    ``1 <= filter_length <= len(reference)``, and non-silent signals. The
    estimate is zero-padded by L-1 samples before fitting and scoring. Returned
    projection/residual arrays therefore have N+L-1 samples. Absolute report
    energy fields outside float64 support cause an explicit ``ValueError``.
    An unrepresentable raw Gram matrix is None with a representation status;
    scaled least-squares projection and scoring remain available.
    """
    if not isinstance(zero_mean, (bool, np.bool_)):
        raise ValueError("zero_mean must be boolean")
    est = _finite_real_array(estimate, "estimate")
    ref = _finite_real_array(reference, "reference")
    if est.ndim != 1 or est.shape != ref.shape or est.size == 0:
        raise ValueError("estimate and reference must be equal-length non-empty 1-D arrays")
    taps = _integer(filter_length, "filter_length", minimum=1)
    if taps > ref.size:
        raise ValueError("filter_length must satisfy 1 <= filter_length <= len(reference)")
    eps = _finite_real_scalar(epsilon, "epsilon")
    if not 0.0 < eps < 1.0:
        raise ValueError("epsilon must lie strictly between 0 and 1")

    est_peak = float(np.max(np.abs(est)))
    ref_peak = float(np.max(np.abs(ref)))
    if est_peak == 0.0 or ref_peak == 0.0:
        raise ValueError("CI-SDR is undefined for a silent or constant-after-centering signal")

    est_unit = est / est_peak
    ref_unit = ref / ref_peak
    if zero_mean:
        est_unit = est_unit - np.mean(est_unit)
        ref_unit = ref_unit - np.mean(ref_unit)
    est_unit_energy = float(est_unit @ est_unit)
    if est_unit_energy == 0.0 or float(ref_unit @ ref_unit) == 0.0:
        raise ValueError("CI-SDR is undefined for a zero-energy signal")

    design_unit = causal_fir_toeplitz(ref_unit, taps)
    est_padded_unit = np.pad(est_unit, (0, taps - 1))
    h_unit, _, rank, _ = np.linalg.lstsq(design_unit, est_padded_unit, rcond=None)

    proj_unit = design_unit @ h_unit
    res_unit = est_padded_unit - proj_unit
    target_fraction = max(float(proj_unit @ proj_unit) / est_unit_energy, eps)
    residual_fraction = max(float(res_unit @ res_unit) / est_unit_energy, eps)
    score_db = float(10.0 * (math.log10(target_fraction) - math.log10(residual_fraction)))

    with np.errstate(over="raise", invalid="raise"):
        try:
            scale_ratio = est_peak / ref_peak
            filter_taps = h_unit * scale_ratio
            target_projection = proj_unit * est_peak
            est_eval = est_padded_unit * est_peak
            ref_eval = ref_unit * ref_peak
            residual = est_eval - target_projection
            design_raw = causal_fir_toeplitz(ref_eval, taps)
            target_energy = float(target_projection @ target_projection)
            residual_energy = float(residual @ residual)
            estimate_energy = float(est_eval @ est_eval)
        except FloatingPointError as error:
            raise ValueError("unscaled CI-SDR energies exceed float64 range") from error
    if not all(np.all(np.isfinite(value)) for value in
               (filter_taps, target_projection, residual,
                target_energy, residual_energy, estimate_energy)):
        raise ValueError("unscaled CI-SDR report exceeds float64 range")

    for vector, energy in ((target_projection, target_energy),
                           (residual, residual_energy), (est_eval, estimate_energy)):
        if energy == 0 and np.any(vector != 0):
            raise ValueError("positive unscaled CI-SDR energy is below float64 support")

    # The full convolution has a positive column energy for every nonzero
    # reference. Do not let raw-unit Gram under/overflow erase that fact.
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        gram_raw = design_raw.T @ design_raw
    if not np.all(np.isfinite(gram_raw)):
        gram_raw, gram_representation = None, "overflow"
    elif np.any(np.diag(gram_raw) == 0):
        gram_raw, gram_representation = None, "positive_diagonal_underflow"
    else:
        gram_representation = "finite"

    return {
        "ci_sdr_db": score_db,
        "filter_length": taps,
        "filter_taps": filter_taps,
        "target_projection": target_projection,
        "residual": residual,
        "target_energy": target_energy,
        "residual_energy": residual_energy,
        "estimate_energy": estimate_energy,
        "gram_matrix": gram_raw,
        "gram_matrix_representation": gram_representation,
        "rank": int(rank),
        "evaluation_samples": int(est_eval.size),
        "full_convolution": True,
        "filter_as_long_as_reference": bool(taps >= ref.size),
        "overparameterized": bool(taps > design_unit.shape[0]),
    }


def streaming_subband_mac_budget(
    sample_rate_hz: int,
    hop_samples: int,
    num_subbands: int,
    num_blocks: int,
    hidden_dim: int,
    projections_per_block: int = 4,
) -> dict:
    """Exact MAC/s accounting for narrowband-crossband streaming separators."""
    fs = _integer(sample_rate_hz, "sample_rate_hz", minimum=1)
    hop = _integer(hop_samples, "hop_samples", minimum=1)
    subbands = _integer(num_subbands, "num_subbands", minimum=1)
    blocks = _integer(num_blocks, "num_blocks", minimum=1)
    dim = _integer(hidden_dim, "hidden_dim", minimum=1)
    projs = _integer(projections_per_block, "projections_per_block", minimum=1)

    frame_rate_hz = fs / hop
    macs_per_frame_per_subband = blocks * projs * (dim ** 2)
    macs_per_frame = subbands * macs_per_frame_per_subband
    # Exact integer MAC/s when hop divides sample_rate_hz.
    if fs % hop == 0:
        macs_per_second: int | float = (fs // hop) * macs_per_frame
    else:
        macs_per_second = float(frame_rate_hz * macs_per_frame)
    return {
        "sample_rate_hz": fs,
        "hop_samples": hop,
        "frame_hop_ms": float(1000.0 * hop / fs),
        "frame_rate_hz": float(frame_rate_hz),
        "num_subbands": subbands,
        "num_blocks": blocks,
        "hidden_dim": dim,
        "projections_per_block": projs,
        "macs_per_frame": int(macs_per_frame),
        "macs_per_second": macs_per_second,
        "gmac_per_second": float(macs_per_second / 1e9),
    }


def single_speaker_bypass_gate(
    slot_waveforms: np.ndarray,
    *,
    threshold_db: float = -18.0,
    min_primary_energy: float = 1e-12,
) -> dict:
    """Evaluate a secondary/primary energy gate on real shape ``(2, T)``.

    Ratios use a common scaled domain. An unrepresentable absolute energy is
    returned as ``None`` with its ``energy_representation`` status; logarithmic
    energies retain positive values whose absolute report would under/overflow.
    A merged output outside float64 range is explicitly rejected.
    """
    slots = _finite_real_array(slot_waveforms, "slot_waveforms")
    if slots.ndim != 2 or slots.shape[0] != 2 or slots.shape[1] == 0:
        raise ValueError("slot_waveforms must have shape (2, T) with T >= 1")
    tau_db = _finite_real_scalar(threshold_db, "threshold_db")
    min_energy = _finite_real_scalar(min_primary_energy, "min_primary_energy")
    if tau_db >= 0.0:
        raise ValueError("threshold_db must be strictly negative")
    if min_energy <= 0.0:
        raise ValueError("min_primary_energy must be strictly positive")

    peaks = np.max(np.abs(slots), axis=1)
    common_peak = float(np.max(peaks))
    if common_peak == 0:
        raise ValueError("primary slot energy is below min_primary_energy (silent input)")
    scaled = slots / common_peak
    scaled_energies = np.sum(scaled * scaled, axis=1)
    order = np.argsort(-scaled_energies, kind="mergesort")
    primary_idx = int(order[0])
    secondary_idx = int(order[1])
    energies, representations, log_energies = [], [], []
    for index, peak in enumerate(peaks):
        if peak == 0:
            energies.append(0.0)
            representations.append("zero")
            log_energies.append(None)
            continue
        unit = slots[index] / peak
        unit_energy = float(unit @ unit)
        log_energies.append(2.0 * math.log(float(peak)) + math.log(unit_energy))
        root_energy = math.sqrt(unit_energy)
        if peak > math.sqrt(np.finfo(float).max) / root_energy:
            energies.append(None)
            representations.append("overflow")
        else:
            norm = float(peak) * root_energy
            energy = norm * norm
            finite_positive = math.isfinite(energy) and energy > 0
            energies.append(energy if finite_positive else None)
            representations.append("finite" if finite_positive else
                                   "overflow" if math.isinf(energy) else "underflow")
    e_primary = energies[primary_idx]
    e_secondary = energies[secondary_idx]
    if representations[primary_idx] == "underflow" or (e_primary is not None and e_primary < min_energy):
        raise ValueError("primary slot energy is below min_primary_energy (silent input)")

    linear_ratio = float(scaled_energies[secondary_idx] / scaled_energies[primary_idx])
    if peaks[secondary_idx] == 0:
        ratio_db = -math.inf
    elif linear_ratio > 0:
        ratio_db = 10.0 * math.log10(linear_ratio)
    else:
        log_ratio = log_energies[secondary_idx] - log_energies[primary_idx]
        linear_ratio = math.exp(log_ratio)
        ratio_db = 10.0 / math.log(10.0) * log_ratio
    bypass = bool(ratio_db <= tau_db)
    if bypass:
        output_slots = np.zeros_like(slots)
        with np.errstate(over="raise", invalid="raise"):
            try:
                output_slots[primary_idx] = slots[0] + slots[1]
            except FloatingPointError as error:
                raise ValueError("merged slot output exceeds float64 range") from error
    else:
        output_slots = slots.copy()

    return {
        "slot_energies": np.asarray(energies, dtype=object if None in energies else float),
        "energy_representation": representations,
        "slot_log_energies": log_energies,
        "primary_slot_index": primary_idx,
        "secondary_slot_index": secondary_idx,
        "primary_energy": e_primary,
        "secondary_energy": e_secondary,
        "secondary_to_primary_ratio": linear_ratio,
        "ratio_underflow": bool(linear_ratio == 0 and peaks[secondary_idx] > 0),
        "secondary_to_primary_db": ratio_db if math.isfinite(ratio_db) else None,
        "threshold_db": tau_db,
        "threshold_linear": float(10.0 ** (tau_db / 10.0)),
        "bypass": bypass,
        "output_waveforms": output_slots,
    }


def ci_sdr_vs_si_sdr_experiment() -> dict:
    """E08-25: Exact 4-sample hand-calculable comparison of SI-SDR (L=1) vs CI-SDR (L=2, L=4)."""
    s_ref = np.array([1.0, 0.0, -1.0, 0.0], dtype=np.float64)
    s_delay1 = np.array([0.0, 1.0, 0.0, -1.0], dtype=np.float64)
    w_orth = 0.1 * np.array([1.0, -1.0, 1.0, -1.0], dtype=np.float64)
    s_room = s_ref + s_delay1 + w_orth
    s_bad = 10.0 * w_orth

    def _serialize(res: dict) -> dict:
        return {
            "ci_sdr_db": res["ci_sdr_db"],
            "filter_length": res["filter_length"],
            "filter_taps": res["filter_taps"].tolist(),
            "target_projection": res["target_projection"].tolist(),
            "residual": res["residual"].tolist(),
            "target_energy": res["target_energy"],
            "residual_energy": res["residual_energy"],
            "estimate_energy": res["estimate_energy"],
            "gram_matrix": None if res["gram_matrix"] is None else res["gram_matrix"].tolist(),
            "gram_matrix_representation": res["gram_matrix_representation"],
            "rank": res["rank"],
            "evaluation_samples": res["evaluation_samples"],
            "full_convolution": res["full_convolution"],
            "filter_as_long_as_reference": res["filter_as_long_as_reference"],
            "overparameterized": res["overparameterized"],
        }

    t4 = causal_fir_toeplitz(s_ref, 4)
    truncated_t4 = causal_fir_toeplitz(s_ref, 4, full_convolution=False)
    truncated_h = np.linalg.solve(truncated_t4, s_bad)
    return {
        "reference": s_ref.tolist(),
        "delayed_reference": s_delay1.tolist(),
        "orthogonal_distortion": w_orth.tolist(),
        "cases": {
            "pure_one_sample_delay": {
                "estimate": s_delay1.tolist(),
                "si_sdr_L1": _serialize(ci_sdr_toeplitz(s_delay1, s_ref, filter_length=1)),
                "ci_sdr_L2": _serialize(ci_sdr_toeplitz(s_delay1, s_ref, filter_length=2)),
            },
            "two_tap_room_plus_distortion": {
                "estimate": s_room.tolist(),
                "si_sdr_L1": _serialize(ci_sdr_toeplitz(s_room, s_ref, filter_length=1)),
                "ci_sdr_L2": _serialize(ci_sdr_toeplitz(s_room, s_ref, filter_length=2)),
            },
            "full_tail_vs_truncated_L_equals_N": {
                "estimate": s_bad.tolist(),
                "toeplitz_L4": t4.tolist(),
                "gram_L4_determinant": float(np.linalg.det(t4.T @ t4)),
                "truncated_toeplitz_L4": truncated_t4.tolist(),
                "truncated_L4_determinant": float(np.linalg.det(truncated_t4)),
                "truncated_L4_filter_taps": truncated_h.tolist(),
                "truncated_L4_residual_energy": float(np.sum((s_bad - truncated_t4 @ truncated_h)**2)),
                "si_sdr_L1": _serialize(ci_sdr_toeplitz(s_bad, s_ref, filter_length=1)),
                "ci_sdr_L2": _serialize(ci_sdr_toeplitz(s_bad, s_ref, filter_length=2)),
                "ci_sdr_L4": _serialize(ci_sdr_toeplitz(s_bad, s_ref, filter_length=4)),
            },
        },
        "scope": (
            "four-sample full-tail Toeplitz projection; short-FIR delay/equalization "
            "examples and a separately labeled cropped projection that overfits at L=N"
        ),
    }


def streaming_mac_and_bypass_experiment() -> dict:
    """E08-26: Narrowband-crossband MAC/s budget & -18 dB single-speaker fricative bypass gate."""
    fullband_budget = streaming_subband_mac_budget(
        sample_rate_hz=16000,
        hop_samples=128,
        num_subbands=256,
        num_blocks=4,
        hidden_dim=64,
        projections_per_block=4,
    )
    mel_budget = streaming_subband_mac_budget(
        sample_rate_hz=16000,
        hop_samples=128,
        num_subbands=64,
        num_blocks=4,
        hidden_dim=64,
        projections_per_block=4,
    )

    # Disjoint-support 4-sample orthogonal representation of vowel (energy 99) + fricative /s/ (energy 1).
    s_vowel = np.array([math.sqrt(99.0), 0.0, 0.0, 0.0], dtype=np.float64)
    s_fric = np.array([0.0, 1.0, 0.0, 0.0], dtype=np.float64)
    s_clean = s_vowel + s_fric
    split_slots = np.vstack((s_vowel, s_fric))
    slot1_eval = ci_sdr_toeplitz(split_slots[0], s_clean, filter_length=1, zero_mean=False)
    gate_single = single_speaker_bypass_gate(split_slots, threshold_db=-18.0)
    bypassed_error_energy = float(
        np.sum((gate_single["output_waveforms"][0] - s_clean) ** 2)
    )

    # True 2-speaker overlap at -10 dB SIR (energies 99 and 9.9).
    overlap_slots = np.vstack(
        (
            np.array([math.sqrt(99.0), 0.0, 0.0, 0.0], dtype=np.float64),
            np.array([0.0, 0.0, math.sqrt(9.9), 0.0], dtype=np.float64),
        )
    )
    gate_overlap = single_speaker_bypass_gate(overlap_slots, threshold_db=-18.0)

    def _serialize_gate(g: dict) -> dict:
        return {
            "slot_energies": g["slot_energies"].tolist(),
            "primary_slot_index": g["primary_slot_index"],
            "secondary_slot_index": g["secondary_slot_index"],
            "primary_energy": g["primary_energy"],
            "secondary_energy": g["secondary_energy"],
            "secondary_to_primary_ratio": g["secondary_to_primary_ratio"],
            "secondary_to_primary_db": g["secondary_to_primary_db"],
            "threshold_db": g["threshold_db"],
            "threshold_linear": g["threshold_linear"],
            "bypass": g["bypass"],
            "output_waveforms": g["output_waveforms"].tolist(),
        }

    return {
        "mac_budget": {
            "fullband_256": fullband_budget,
            "mel_band_64": mel_budget,
            "compression_factor": float(
                fullband_budget["macs_per_second"] / mel_budget["macs_per_second"]
            ),
        },
        "fricative_splitting_case": {
            "vowel_energy": 99.0,
            "fricative_energy": 1.0,
            "total_clean_energy": 100.0,
            "fricative_to_vowel_db": float(10.0 * math.log10(1.0 / 99.0)),
            "slot1_uncentered_si_sdr_db": slot1_eval["ci_sdr_db"],
            "reference_word_count": 4,
            "unbypassed_missing_fricative_substitutions": 1,
            "unbypassed_wer": 0.25,
            "bypass_gate": _serialize_gate(gate_single),
            "bypassed_reconstruction_error_energy": bypassed_error_energy,
            "bypassed_wer": 0.0,
        },
        "two_speaker_overlap_case": {
            "sir_db": -10.0,
            "bypass_gate": _serialize_gate(gate_overlap),
        },
        "scope": (
            "linear projection MAC/s accounting and deterministic 2-slot energy-ratio "
            "bypass control; word-error counts refer to the declared 4-word toy transcript"
        ),
    }
