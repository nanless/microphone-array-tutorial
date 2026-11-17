"""Deterministic PDM CIC decimation and ITU-T P.340/P.808 evaluation helpers.

All calculations are closed-form mathematical controls from stated inputs:
Hogenauer register bit growth, exact discrete CIC transfer droop, linear-phase
group delay, stereo rising/falling-edge half-clock skew, ITU-T P.340 double-talk
attenuation categories, and ITU-T P.808 trap/gold worker screening. No hardware
capture, neural MOS inference, or subjective listening test is executed.
"""

from __future__ import annotations

import cmath
from collections.abc import Sequence
from fractions import Fraction
import math
from numbers import Integral

import numpy as np

from codes.chapters.ch02.core.conventions import SPEED_OF_SOUND, finite_real_scalar

P340_TYPES = ("Type 1", "Type 2a", "Type 2b", "Type 2c", "Type 3")


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer, not bool or float")
    result = int(value)
    if result < 1:
        raise ValueError(f"{name} must be >= 1")
    return result


def _exact_decimal_fraction(value: object, name: str) -> Fraction:
    """Convert finite int/float/Fraction to exact decimal Fraction without bools."""
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must not be a boolean")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value))
    if isinstance(value, (float, np.floating)):
        if isinstance(value, np.floating) and value.dtype.itemsize > 8:
            raise ValueError(f"{name} extended-precision floats are not supported")
        scalar = float(value)
        if not math.isfinite(scalar):
            raise ValueError(f"{name} must be finite")
        return Fraction(str(scalar))
    raise ValueError(f"{name} must be a finite int, float, or Fraction")


def cic_decimation_metrics(
    pdm_clock_hz: float,
    decimation: int,
    order: int,
    eval_freq_hz: float,
    *,
    diff_delay: int = 1,
    input_bits: int = 1,
    mic_spacing_m: float = 0.012,
    sound_speed: float = SPEED_OF_SOUND,
    hal_freq_hz: float = 1000.0,
    hal_gains: tuple[float, float] = (1.0, 0.8),
) -> dict[str, float | int | bool]:
    """Compute exact CIC register width, passband droop, group delay, and stereo skew."""
    f_pdm = finite_real_scalar(pdm_clock_hz, "pdm_clock_hz")
    f_eval = finite_real_scalar(eval_freq_hz, "eval_freq_hz")
    d_m = finite_real_scalar(mic_spacing_m, "mic_spacing_m")
    c_m_s = finite_real_scalar(sound_speed, "sound_speed")
    f_hal = finite_real_scalar(hal_freq_hz, "hal_freq_hz")
    r = _positive_int(decimation, "decimation")
    n_stages = _positive_int(order, "order")
    m_delay = _positive_int(diff_delay, "diff_delay")
    b_in = _positive_int(input_bits, "input_bits")

    if f_pdm <= 0.0:
        raise ValueError("pdm_clock_hz must be positive")
    if d_m <= 0.0 or c_m_s <= 0.0:
        raise ValueError("mic_spacing_m and sound_speed must be positive")
    if f_hal <= 0.0:
        raise ValueError("hal_freq_hz must be positive")

    rm = r * m_delay
    nyquist_after_cic = f_pdm / (2.0 * rm)
    if not (0.0 < f_eval < nyquist_after_cic):
        raise ValueError("eval_freq_hz must lie strictly inside (0, pdm_clock_hz / (2 * decimation * diff_delay))")

    if not isinstance(hal_gains, Sequence) or isinstance(hal_gains, (str, bytes)) or len(hal_gains) != 2:
        raise ValueError("hal_gains must be a two-element numeric sequence")
    g1 = finite_real_scalar(hal_gains[0], "hal_gains[0]")
    g2 = finite_real_scalar(hal_gains[1], "hal_gains[1]")
    if g1 <= 0.0 or g2 <= 0.0 or g1 == g2:
        raise ValueError("hal_gains must be distinct positive channel gains")

    dc_gain = rm ** n_stages
    bit_growth = (dc_gain - 1).bit_length()
    max_register_bits = b_in + bit_growth

    arg_small = math.pi * f_eval / f_pdm
    arg_large = arg_small * rm
    single_stage_gain = abs(math.sin(arg_large) / (rm * math.sin(arg_small)))
    passband_gain = single_stage_gain ** n_stages
    passband_droop_db = 20.0 * math.log10(passband_gain)
    cont_sinc_stage = abs(math.sin(arg_large) / arg_large)
    continuous_sinc_droop_db = 20.0 * n_stages * math.log10(cont_sinc_stage)

    group_delay_s = n_stages * (rm - 1) / (2.0 * f_pdm)
    half_clock_skew_s = 1.0 / (2.0 * f_pdm)
    skew_path_m = c_m_s * half_clock_skew_s
    if skew_path_m > 2.0 * d_m:
        raise ValueError("half-clock skew path exceeds 2 * mic_spacing_m; no real endfire null exists")

    cos_null = -(1.0 - skew_path_m / d_m)
    cos_null = max(-1.0, min(1.0, cos_null))
    endfire_null_deg = math.degrees(math.acos(cos_null))
    endfire_null_shift_deg = 180.0 - endfire_null_deg

    omega_tau0 = 2.0 * math.pi * f_hal * (d_m / c_m_s)
    front_mag = abs(g1 - g2 * cmath.exp(-1j * 2.0 * omega_tau0))
    if front_mag <= 0.0:
        raise ValueError("hal_freq_hz hits an endfire spatial alias null; front response is zero")
    rear_mag = abs(g1 - g2)
    hal_rel_linear = rear_mag / front_mag
    hal_rel_db = 20.0 * math.log10(hal_rel_linear)

    return {
        "pdm_clock_hz": f_pdm,
        "output_rate_hz": f_pdm / r,
        "decimation": r,
        "diff_delay": m_delay,
        "order": n_stages,
        "input_bits": b_in,
        "eval_freq_hz": f_eval,
        "dc_gain": dc_gain,
        "bit_growth": bit_growth,
        "max_register_bits": max_register_bits,
        "fits_signed_24bit": max_register_bits <= 24,
        "single_stage_gain": single_stage_gain,
        "passband_gain": passband_gain,
        "passband_droop_db": passband_droop_db,
        "continuous_sinc_droop_db": continuous_sinc_droop_db,
        "group_delay_s": group_delay_s,
        "group_delay_us": group_delay_s * 1e6,
        "half_clock_skew_s": half_clock_skew_s,
        "half_clock_skew_ns": half_clock_skew_s * 1e9,
        "skew_path_mm": skew_path_m * 1e3,
        "skew_phase_error_deg": 360.0 * f_eval * half_clock_skew_s,
        "endfire_null_deg": endfire_null_deg,
        "endfire_null_shift_deg": endfire_null_shift_deg,
        "hal_gain_mismatch_180_linear": hal_rel_linear,
        "hal_gain_mismatch_180_db": hal_rel_db,
    }


def pdm_cic_teaching_case() -> dict[str, object]:
    """Return the E10-35 stereo PDM CIC decimation and skew control."""
    stage1_r16 = cic_decimation_metrics(3_072_000.0, 16, 4, 7000.0)
    stage1_r64 = cic_decimation_metrics(3_072_000.0, 64, 4, 7000.0)
    stage1_r64_20k = cic_decimation_metrics(3_072_000.0, 64, 4, 20000.0)
    single_16k_r64 = cic_decimation_metrics(1_024_000.0, 64, 4, 7000.0)
    return {
        "stage1_R16": stage1_r16,
        "stage1_R64": stage1_r64,
        "stage1_R64_20khz": stage1_r64_20k,
        "single_16khz_R64_7khz": single_16k_r64,
        "half_clock_skew_ns": stage1_r16["half_clock_skew_ns"],
        "skew_path_mm": stage1_r16["skew_path_mm"],
        "phase_error_7khz_deg": stage1_r16["skew_phase_error_deg"],
        "endfire_null_deg": stage1_r16["endfire_null_deg"],
        "endfire_null_shift_deg": stage1_r16["endfire_null_shift_deg"],
        "hal_unprocessed_180_linear": 0.0,
        "hal_voice_comm_180_rel_linear": stage1_r16["hal_gain_mismatch_180_linear"],
        "hal_voice_comm_180_rel_db": stage1_r16["hal_gain_mismatch_180_db"],
        "scope": "closed-form Hogenauer CIC, half-clock PDM edge skew, and 2-mic cardioid gain mismatch",
    }


def _classify_p340_sending(attenuation_db: float) -> str:
    if attenuation_db <= 3.0:
        return "Type 1"
    if attenuation_db <= 6.0:
        return "Type 2a"
    if attenuation_db <= 9.0:
        return "Type 2b"
    if attenuation_db <= 12.0:
        return "Type 2c"
    return "Type 3"


def _classify_p340_receiving(attenuation_db: float) -> str:
    if attenuation_db <= 3.0:
        return "Type 1"
    if attenuation_db <= 5.0:
        return "Type 2a"
    if attenuation_db <= 8.0:
        return "Type 2b"
    if attenuation_db <= 10.0:
        return "Type 2c"
    return "Type 3"


def p340_and_p808_evaluation(
    single_talk_rms: tuple[float, float],
    candidates_dt_rms: dict[str, tuple[float, float]],
    worker_submissions: list[dict[str, object]],
    *,
    trap_expected: int = 1,
    gold_targets: tuple[float, float] = (4.8, 1.4),
    gold_tol: float = 1.0,
) -> dict[str, object]:
    """Classify ITU-T P.340 duplex types and filter ITU-T P.808 crowd ratings."""
    if (not isinstance(single_talk_rms, Sequence)
            or isinstance(single_talk_rms, (str, bytes))
            or len(single_talk_rms) != 2):
        raise ValueError("single_talk_rms must be (u_send_st, u_recv_st)")
    u_s_st = finite_real_scalar(single_talk_rms[0], "single_talk_rms[0]")
    u_r_st = finite_real_scalar(single_talk_rms[1], "single_talk_rms[1]")
    if u_s_st <= 0.0 or u_r_st <= 0.0:
        raise ValueError("single_talk_rms values must be strictly positive")

    if not isinstance(candidates_dt_rms, dict) or not candidates_dt_rms:
        raise ValueError("candidates_dt_rms must be a non-empty dict")

    p340_results: dict[str, dict[str, float | str]] = {}
    for name, pair in candidates_dt_rms.items():
        if not isinstance(name, str) or not name:
            raise ValueError("candidate names must be non-empty strings")
        if not isinstance(pair, Sequence) or isinstance(pair, (str, bytes)) or len(pair) != 2:
            raise ValueError(f"candidate {name} double-talk RMS must be a 2-tuple")
        u_s_dt = finite_real_scalar(pair[0], f"{name}.u_send_dt")
        u_r_dt = finite_real_scalar(pair[1], f"{name}.u_recv_dt")
        if u_s_dt <= 0.0 or u_r_dt <= 0.0:
            raise ValueError(f"candidate {name} double-talk RMS values must be positive")
        ah_s_db = 20.0 * math.log10(u_s_st / u_s_dt)
        ah_r_db = 20.0 * math.log10(u_r_st / u_r_dt)
        send_type = _classify_p340_sending(ah_s_db)
        recv_type = _classify_p340_receiving(ah_r_db)
        overall_type = P340_TYPES[max(P340_TYPES.index(send_type), P340_TYPES.index(recv_type))]
        p340_results[name] = {
            "u_send_dt": u_s_dt,
            "u_recv_dt": u_r_dt,
            "sending_dt_attenuation_db": ah_s_db,
            "receiving_dt_attenuation_db": ah_r_db,
            "sending_type": send_type,
            "receiving_type": recv_type,
            "overall_type": overall_type,
        }

    trap_target = _positive_int(trap_expected, "trap_expected")
    if trap_target > 5:
        raise ValueError("trap_expected ACR score must be in 1..5")
    if not isinstance(gold_targets, Sequence) or isinstance(gold_targets, (str, bytes)) or len(gold_targets) != 2:
        raise ValueError("gold_targets must be (clean_mos, degraded_mos)")
    g_clean = _exact_decimal_fraction(gold_targets[0], "gold_targets[0]")
    g_deg = _exact_decimal_fraction(gold_targets[1], "gold_targets[1]")
    g_tol = _exact_decimal_fraction(gold_tol, "gold_tol")
    if g_tol < 0 or g_clean <= g_deg:
        raise ValueError("gold_tol must be non-negative and clean gold target must exceed degraded target")

    if not isinstance(worker_submissions, Sequence) or isinstance(worker_submissions, (str, bytes)) or not worker_submissions:
        raise ValueError("worker_submissions must be a non-empty sequence")

    seen_workers: set[str] = set()
    clip_count: int | None = None
    raw_votes: list[Fraction] = []
    screened_votes: list[Fraction] = []
    worker_rows: list[dict[str, object]] = []
    accepted_ids: list[str] = []
    rejected_ids: list[str] = []

    for entry in worker_submissions:
        if not isinstance(entry, dict):
            raise ValueError("each worker submission must be a dict")
        wid = entry.get("worker_id")
        if not isinstance(wid, str) or not wid or wid in seen_workers:
            raise ValueError("worker_id must be a unique non-empty string")
        seen_workers.add(wid)

        trap_val = _positive_int(entry.get("trap"), f"{wid}.trap")
        if trap_val > 5:
            raise ValueError(f"{wid}.trap must be in 1..5")

        gold_pair = entry.get("gold")
        if not isinstance(gold_pair, Sequence) or isinstance(gold_pair, (str, bytes)) or len(gold_pair) != 2:
            raise ValueError(f"{wid}.gold must be a 2-element sequence")
        r_clean = _exact_decimal_fraction(gold_pair[0], f"{wid}.gold[0]")
        r_deg = _exact_decimal_fraction(gold_pair[1], f"{wid}.gold[1]")
        if not (Fraction(1) <= r_clean <= Fraction(5) and Fraction(1) <= r_deg <= Fraction(5)):
            raise ValueError(f"{wid}.gold ratings must lie in [1, 5]")

        ratings_seq = entry.get("clips")
        if not isinstance(ratings_seq, Sequence) or isinstance(ratings_seq, (str, bytes)) or not ratings_seq:
            raise ValueError(f"{wid}.clips must be a non-empty sequence")
        if clip_count is None:
            clip_count = len(ratings_seq)
        elif len(ratings_seq) != clip_count:
            raise ValueError("all workers must rate the same number of clips")

        clip_fracs = [_exact_decimal_fraction(v, f"{wid}.clips") for v in ratings_seq]
        if any(not (Fraction(1) <= v <= Fraction(5)) for v in clip_fracs):
            raise ValueError(f"{wid}.clips ACR ratings must lie in [1, 5]")

        trap_ok = trap_val == trap_target
        gold_ok = (
            abs(r_clean - g_clean) <= g_tol
            and abs(r_deg - g_deg) <= g_tol
            and r_clean > r_deg
        )
        passed = trap_ok and gold_ok
        if passed:
            reason = "pass"
            accepted_ids.append(wid)
            screened_votes.extend(clip_fracs)
        elif not trap_ok and not gold_ok:
            reason = "trap_and_gold"
            rejected_ids.append(wid)
        elif not trap_ok:
            reason = "trap_mismatch"
            rejected_ids.append(wid)
        else:
            reason = "gold_mismatch"
            rejected_ids.append(wid)

        raw_votes.extend(clip_fracs)
        worker_rows.append({
            "worker_id": wid,
            "passed": passed,
            "reason": reason,
            "clip_sum": float(sum(clip_fracs)),
        })

    if not screened_votes:
        raise ValueError("all workers failed P.808 screening; screened MOS is undefined")

    raw_mean = sum(raw_votes, Fraction(0)) / len(raw_votes)
    screened_mean = sum(screened_votes, Fraction(0)) / len(screened_votes)
    debias = screened_mean - raw_mean
    retention = Fraction(len(accepted_ids), len(worker_rows))

    return {
        "single_talk_rms": {"send": u_s_st, "recv": u_r_st},
        "p340": p340_results,
        "p808": {
            "worker_rows": worker_rows,
            "accepted_workers": accepted_ids,
            "rejected_workers": rejected_ids,
            "worker_retention_fraction": str(retention),
            "worker_retention_rate": float(retention),
            "raw_vote_count": len(raw_votes),
            "screened_vote_count": len(screened_votes),
            "raw_mos_fraction": str(raw_mean),
            "raw_mos": float(raw_mean),
            "screened_mos_fraction": str(screened_mean),
            "screened_mos": float(screened_mean),
            "mos_debias_fraction": str(debias),
            "mos_debias_delta": float(debias),
        },
    }


def p340_p808_teaching_case() -> dict[str, object]:
    """Return the E10-36 P.340 double-talk and P.808 crowd-screening control."""
    result = p340_and_p808_evaluation(
        single_talk_rms=(0.20, 0.20),
        candidates_dt_rms={
            "A": (0.04, 0.20),
            "B": (0.08, 0.10),
            "C": (0.16, 0.18),
        },
        worker_submissions=[
            {"worker_id": "w1", "trap": 1, "gold": (5, 1), "clips": [4, 5, 4, 5]},
            {"worker_id": "w2", "trap": 1, "gold": (4, 2), "clips": [4, 4, 5, 4]},
            {"worker_id": "w3", "trap": 1, "gold": (5, 2), "clips": [5, 4, 4, 4]},
            {"worker_id": "w4", "trap": 1, "gold": (5, 1), "clips": [4, 4, 4, 5]},
            {"worker_id": "w5", "trap": 3, "gold": (3, 3), "clips": [3, 3, 3, 3]},
            {"worker_id": "w6", "trap": 1, "gold": (2, 4), "clips": [1, 2, 1, 2]},
        ],
    )
    result["single_talk_erle_db"] = {"A": 48.0, "B": 35.0, "C": 32.0}
    result["evaluator_blind_spots"] = {
        "POLQA_P863": "intrusive auto-alignment masks fixed 250 ms buffer latency",
        "ViSQOL_v3": "patch spectro-temporal similarity tolerates constant lag and mild phase warp",
        "DNSMOS_P835": "non-intrusive neural MOS can reward generative vocoder hallucinations without WER check",
        "NISQA": "blind dimension predictor is sensitive to out-of-domain codec/far-field priors without ASR",
    }
    result["scope"] = "constructed RMS levels and ACR ballot table; no ITU software, neural model, or human study run"
    return result
