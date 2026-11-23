"""Audio-visual (AV) frame synchronization, camera-lag onset loss, and occlusion fallback.

Conventions:
- Audio STFT hop index ``l`` (``0 <= l < n_audio_frames``) has physical hop-center
  timestamp ``t_a[l] = (l + 0.5) * audio_hop_s`` seconds.
- Video frame ``m >= 0`` at frame rate ``video_fps`` spans physical exposure interval
  ``[m * T_v, (m + 1) * T_v)`` with center ``(m + 0.5) * T_v`` (``T_v = 1 / video_fps``).
- Positive ``av_offset_s > 0`` denotes camera/ISP pipeline lag relative to audio:
  a video frame exposed at physical time ``t_lip`` arrives at host timestamp
  ``t_lip + av_offset_s``, so uncompensated arrival-time matching queries past lip
  state at ``t_a[l] - av_offset_s``.
"""

from __future__ import annotations

import numpy as np

from codes.chapters.ch02.core.conventions import finite_real_array, finite_real_scalar


def map_video_to_audio_frames(
    n_audio_frames: int,
    audio_hop_s: float,
    video_fps: float,
    av_offset_s: float = 0.0,
    onset_frame: int = 10,
    consonant_frames: int = 8,
    *,
    alignment_mode: str = "online_last_arrived",
) -> dict:
    """Map video to audio using explicit online arrival or offline association.

    Online mode selects the last frame whose center timestamp plus ISP lag
    has arrived. This is a center-arrival teaching model, not an assertion
    that a real camera delivers a full exposure at its center. Offline mode
    associates the nearest exposure center, allowing a future-arriving frame.
    Audio decisions in this finite model occur at the stated hop-center times;
    a real STFT/inference availability delay is not measured here.
    An unavailable online frame has index -1 and zero timestamp placeholders;
    consumers must use ``valid_video_mask`` before indexing or plotting.

    Parameters
    ----------
    n_audio_frames:
        Total number of audio STFT frames ``L >= 1``.
    audio_hop_s:
        Audio STFT hop duration ``H_a`` in seconds (``> 0``).
    video_fps:
        Camera frame rate ``f_v`` in frames per second (``> 0``).
    av_offset_s:
        Camera pipeline lag ``Delta t_av`` in seconds (``> 0`` means video lags audio).
    onset_frame:
        Reference word-initial consonant start STFT frame index for onset-loss audit.
    consonant_frames:
        Duration of the word-initial plosive/fricative segment in STFT frames.
    """
    if isinstance(n_audio_frames, (bool, np.bool_)) or not isinstance(n_audio_frames, (int, np.integer)) or n_audio_frames < 1:
        raise ValueError("n_audio_frames must be a positive integer")
    if isinstance(onset_frame, (bool, np.bool_)) or not isinstance(onset_frame, (int, np.integer)) or onset_frame < 0:
        raise ValueError("onset_frame must be a non-negative integer")
    if isinstance(consonant_frames, (bool, np.bool_)) or not isinstance(consonant_frames, (int, np.integer)) or consonant_frames < 1:
        raise ValueError("consonant_frames must be a positive integer")

    hop = finite_real_scalar(audio_hop_s, "audio_hop_s")
    fps = finite_real_scalar(video_fps, "video_fps")
    offset = finite_real_scalar(av_offset_s, "av_offset_s")
    if not isinstance(alignment_mode, str) or alignment_mode not in ("online_last_arrived", "offline_nearest_exposure"):
        raise ValueError("unknown alignment_mode")
    if hop <= 0.0 or fps <= 0.0:
        raise ValueError("audio_hop_s and video_fps must be strictly positive")

    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        video_period_s = 1.0 / fps
        hops_per_vid = float(video_period_s / hop)
        audio_indices = np.arange(int(n_audio_frames), dtype=np.int64)
        audio_times_s = (audio_indices + 0.5) * hop
        queried_lip_times_s = audio_times_s - offset
        scaled_vid_query = np.round(queried_lip_times_s / video_period_s, 12)
        t_onset_s = float(onset_frame * hop)
        cons_dur_s = float(consonant_frames * hop)
        cons_indices = np.arange(int(onset_frame), int(onset_frame + consonant_frames), dtype=np.int64)
        cons_times_s = (cons_indices + 0.5) * hop
        cons_queried_lip_s = cons_times_s - offset
        scaled_cons_query = np.round(cons_queried_lip_s / video_period_s, 12)

    if (
        video_period_s <= 0.0
        or not np.isfinite(video_period_s)
        or not np.isfinite(hops_per_vid)
        or not np.isfinite(t_onset_s)
        or not np.isfinite(cons_dur_s)
        or not np.all(np.isfinite(audio_times_s))
        or not np.all(np.isfinite(queried_lip_times_s))
        or not np.all(np.isfinite(scaled_vid_query))
        or not np.all(np.isfinite(scaled_cons_query))
        or np.any(np.abs(scaled_vid_query) > 9.0e15)
        or np.any(np.abs(scaled_cons_query) > 9.0e15)
    ):
        raise ValueError("audio-visual timestamp mapping exceeds finite float64/int64 support")

    def associate(scaled_query: np.ndarray, audio_times: np.ndarray):
        online = alignment_mode == "online_last_arrived"
        indices = np.floor(scaled_query - (0.5 if online else 0.0)).astype(np.int64)
        if online:
            # Rounding a quotient must never admit even a slightly future arrival.
            arrival = (indices + 0.5) * video_period_s + offset
            indices = indices - (arrival > audio_times)
            valid = indices >= 0
            indices = np.where(valid, indices, -1)
        else:
            indices = np.maximum(0, indices)
            valid = np.ones(indices.shape, dtype=bool)
        centers = np.where(valid, (indices + 0.5) * video_period_s, 0.0)
        arrivals = np.where(valid, centers + offset, 0.0)
        return indices, valid, centers, arrivals

    video_indices, valid_video, video_exposure_centers_s, video_arrivals_s = associate(
        scaled_vid_query, audio_times_s,
    )
    effective_lip_lag_s = np.where(valid_video, audio_times_s - video_exposure_centers_s, 0.0)

    # Quantify word-initial consonant clipping when talker mouth opens at t_onset = onset_frame * hop
    _, cons_valid, cons_vid_centers_s, _ = associate(scaled_cons_query, cons_times_s)
    if not all(np.all(np.isfinite(values)) for values in (
        video_exposure_centers_s, cons_vid_centers_s, video_arrivals_s, effective_lip_lag_s,
    )):
        raise ValueError("mapped video timestamps exceed finite float64 range")
    # A video exposure centered strictly before t_onset_s captures pre-speech closed lips -> gates onset off
    clipped_mask = (~cons_valid) | (cons_vid_centers_s < (t_onset_s - 1e-12))
    clipped_consonant_frames = int(np.sum(clipped_mask))
    consonant_loss_fraction = float(clipped_consonant_frames / int(consonant_frames))

    return {
        "n_audio_frames": int(n_audio_frames),
        "audio_hop_s": hop,
        "video_fps": fps,
        "video_period_s": video_period_s,
        "hops_per_video_frame": hops_per_vid,
        "av_offset_s": offset,
        "alignment_mode": alignment_mode,
        "audio_center_times_s": audio_times_s,
        "audio_decision_times_s": audio_times_s.copy(),
        "queried_lip_times_s": queried_lip_times_s,
        "video_frame_indices": video_indices,
        "valid_video_mask": valid_video,
        "unavailable_video_mask": ~valid_video,
        "mapped_video_centers_s": video_exposure_centers_s,
        "mapped_video_arrival_times_s": video_arrivals_s,
        "effective_lip_lag_s": effective_lip_lag_s,
        "onset_audit": {
            "onset_frame": int(onset_frame),
            "onset_time_s": t_onset_s,
            "consonant_frames": int(consonant_frames),
            "consonant_duration_s": float(consonant_frames * hop),
            "clipped_consonant_frames": clipped_consonant_frames,
            "clipped_duration_s": float(clipped_consonant_frames * hop),
            "consonant_loss_fraction": consonant_loss_fraction,
        },
    }


def av_confidence_gated_mask(
    mask_audio: np.ndarray,
    mask_av: np.ndarray,
    alpha_v: np.ndarray | float,
) -> np.ndarray:
    """Blend audio-only and audio-visual TF masks via visual confidence ``alpha_v in [0, 1]``.

    Implements ``M_fused = (1 - alpha_v) * mask_audio + alpha_v * mask_av``, smoothly
    falling back to ``mask_audio`` under lip occlusion, head turn, or camera dropout.
    """
    m_aud = finite_real_array(mask_audio, "mask_audio")
    m_av = finite_real_array(mask_av, "mask_av")
    alpha = finite_real_array(alpha_v, "alpha_v")
    if m_aud.ndim not in (1, 2) or m_aud.size < 1 or m_aud.shape != m_av.shape:
        raise ValueError("mask_audio and mask_av must have identical non-empty 1-D or 2-D shapes")
    if np.any(m_aud < 0.0) or np.any(m_aud > 1.0) or np.any(m_av < 0.0) or np.any(m_av > 1.0):
        raise ValueError("masks must lie within [0, 1]")
    if np.any(alpha < 0.0) or np.any(alpha > 1.0):
        raise ValueError("visual confidence alpha_v must lie within [0, 1]")

    if alpha.ndim == 0:
        w_v = float(alpha)
        return (1.0 - w_v) * m_aud + w_v * m_av
    if m_aud.ndim == 1 and alpha.shape == m_aud.shape:
        return (1.0 - alpha) * m_aud + alpha * m_av
    if m_aud.ndim == 2 and alpha.shape == (m_aud.shape[1],):
        return (1.0 - alpha[None, :]) * m_aud + alpha[None, :] * m_av
    if m_aud.ndim == 2 and alpha.shape == m_aud.shape:
        return (1.0 - alpha) * m_aud + alpha * m_av
    raise ValueError("alpha_v must be scalar, per-frame (T,), or match mask shape")
