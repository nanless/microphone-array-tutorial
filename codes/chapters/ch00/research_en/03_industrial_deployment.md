# Industrial Audio Implementation: From Capture and State Management to Deployment and Scoring

🌐 [简体中文](../research/03_industrial_deployment.md) ｜ **English**

Static re-verification date for DeepFilterNet host callback counting and `faster-enhancer.c` lifecycle/pacing: 2026-10-05; static re-verification date for AEC playout-volume notifications and XMOS `ALT` hold counting: 2026-10-04; initial verification date: 2026-09-22; re-verification date for AEC interfaces I03/I05 and scorer I20: 2026-09-23; verification date for I01 timestamps, I18 I/O binding, and I28 source interfaces: 2026-09-26; re-verification date for Chapter 10 industrial boundaries and I30–I32: 2026-09-28; FastEnhancer conference acceptance status separately verified on 2026-09-29. On 2026-10-01, direct Chapter 10 interfaces were re-checked and the five bounded original-source contracts in I32 were executed, while historical reports preserve their original dates and environments. Corresponding main text: Chapters 10, 11, and 13, plus Appendix B. This guide examines the engineering layers that must be supplied once textbook algorithms enter a continuous audio product: device drivers, reference routing, cross-block state persistence, fixed-point DSP kernels, neural inference runtimes, and objective evaluation harnesses. Exact upstream commit hashes are governed by [`SOURCES.lock.json`](../SOURCES.lock.json); web documentation versions indicate our reading references rather than implying local installation of matching system packages.

Except for local interface probes that explicitly report inputs, host environments, and measured numerical outputs, the test procedures below specify recommended hardware acceptance steps rather than already-measured product benchmarks. Source acquisition, dependency installation, compilation, runtime execution, and acoustic chamber measurement represent distinct verification states. Without target hardware, checkpoints, or licensed speech corpora, one can complete static source verification, yet must never mark an item as having passed device acceptance.

When multiple networked nodes collaborate on waveform enhancement, start from the [observation and update models in Extended Topic II](../../../../chapters_en/13_distributed-enhancement.md), proceed to the [distributed network timeline](#distributed-network-deployment) in §8 below, and consult [pinned distributed source contracts](04_source_reproduction.md#distributed-reproduction) to separate actual helper execution from static inspection.

## 1. Audio Capture, Clocking, and Continuous Streams

### I01: PortAudio Callbacks and Hardware Timestamps

PortAudio provides a cross-platform C abstraction over host audio backends. This repository pins the commit corresponding to release `19.7.0`. Begin with the callback signature in `include/portaudio.h`, then inspect the capture and playback examples under `examples/`; although PortAudio itself uses the MIT License, host OS backends carry their own platform dependencies. Official [callback guidelines](https://portaudio.com/docs/v19-doxydocs/writing_a_callback.html) strictly forbid operations with unbounded or unpredictable latency—such as dynamic memory allocation, file I/O, mutex locking, or console printing—inside the real-time audio callback.

Every hardware test log must record the active host API, device ID, sample rate, sample format, channel count, requested vs. negotiated latency, and callback block length. Note that `frameCount` denotes time samples *per channel*—never multiply `frameCount` by channel count and mistake the scalar sample total for a frame count. Log input ADC timestamps, output DAC timestamps, and callback status flags; never treat timestamps across two distinct physical devices as sharing a common clock without explicit calibration.

End-to-end latency acceptance can inject an acoustic or electrical loopback impulse while stressing disk I/O and CPU load, correlating callback status flags, inter-callback intervals, ring-buffer occupancy, and audible dropouts. A conservative architecture copies callback frames into pre-allocated lock-free ring buffers, though bounded real-time DSP with verified worst-case execution time may run in-callback; blocking Python inference or synchronous logging must never sit on the callback thread. After hot-unplugging and reconnecting a USB array, re-verify physical channel mapping rather than merely re-opening the previous numeric device index.

Pinned [`portaudio.h`](https://github.com/PortAudio/portaudio/blob/3f7bee79a65327d2e0965e8a74299723ed6f072d/include/portaudio.h) permits `paFramesPerBufferUnspecified`, under which `frameCount` varies dynamically from callback to callback. Real-time factor (RTF) must divide processing time by the newly arrived audio duration `frameCount / sampleRate`; similarly, overlapped STFT pipelines must normalize compute budget against the *hop duration* of newly arrived samples rather than inflating the denominator to the full analysis window length (which reuses historical samples). Hardware callback blocks, algorithm hop sizes, and playback blocks can be decoupled cleanly via bounded FIFO queues without forcing all stages onto a giant least-common-multiple buffer.

Crucially, `paInputOverflow` signals only that input samples were dropped—it does *not* provide a portable, exact dropped-sample count. Under fixed-size callback buffers, upstream documentation notes that dropped samples may fall at one or multiple interior positions within the delivered buffer, whereas variable-size buffers place drops prior to the first sample. Status flags, timestamps, hardware sample counters, and stream restarts must all be logged separately: record an exact dropped-sample count only when corroborated by independent hardware counters, and otherwise log `"discontinuity_detected, dropped_samples=unknown"`—never substituting `0` and pretending AEC or resampler states remain phase-continuous.

**Log timestamp tap points alongside their exact sample-index correspondence.** Official [`PaStreamCallbackTimeInfo`](https://files.portaudio.com/docs/v19-doxydocs/structPaStreamCallbackTimeInfo.html) exposes three distinct instants in seconds on the stream's clock: `inputBufferAdcTime` (ADC capture instant of the first sample in `inputBuffer`), `currentTime` (instant the callback was invoked), and `outputBufferDacTime` (DAC playout instant of the first sample in `outputBuffer`). Subtracting `outputBufferDacTime - inputBufferAdcTime` yields true end-to-end sample latency *only* when both belong to the same clock domain and you track which input sample maps to which output sample. If a worker thread defers output emission by one callback block, the sample-index offset shifts by one block as well—one cannot keep subtracting the two first-sample timestamps of the *same* callback invocation.

[Exercise E10-17](../../../../chapters_en/10_engineering-practice.md#sec-10-11) works out by hand a `15` ms input age plus `30` ms output lead yielding `45` ms corresponding-sample latency—rising to `55` ms when queued one extra `10` ms block—and shows why `3` ms of in-window compute sitting between those two hardware instants must NOT be double-counted. PortAudio's [Buffering, Latency and Timing Implementation Guidelines](https://github.com/PortAudio/portaudio/wiki/BufferingLatencyAndTimingImplementationGuidelines) further caution that some host backends report only partial/nominal driver buffering, so struct fields never substitute for physical loopback impulse verification (no soundcard was opened in E10-17).

### I02: ALSA PCM State Machine and XRUN Recovery

On Linux, device-level behavior can be audited directly in `alsa-project/alsa-lib` via `src/pcm/pcm.c` and [`test/pcm.c`](https://github.com/alsa-project/alsa-lib/blob/f84cd4ced7b36fddb8e4ee24404cf7c091d27020/test/pcm.c). Root `COPYING` specifies LGPL-2.1, with per-file headers governing exact version clauses and exceptions. ALSA's PCM (Pulse-Code Modulation) API cleanly separates hardware ring-buffer parameters (`hw_params`), software start/wakeup thresholds (`sw_params`), and runtime stream states.

The [official ALSA PCM documentation](https://www.alsa-project.org/alsa-doc/alsa-lib/pcm.html) distinguishes error codes by root cause: `-EPIPE` signals playback underrun or capture overrun (XRUN); `-ESTRPIPE` indicates system suspend; `-ENODEV` indicates physical device removal; and hardware echo-reference loopback subdevices may return `-ENODATA` when the paired playback stream has not yet started. Calling a single blind recovery wrapper on every negative return code—and assuming sample continuity afterward—breaks downstream stateful DSP.

Deployment specs must pin `period_size`, `buffer_size`, `avail_min`, `start_threshold`, interleaved vs. non-interleaved access layout, and monotonic timestamping mode. Test short reads/writes, `SIGSTOP`/`SIGCONT` process stalls, OS suspend/resume, USB hot-unplugging, and opening the echo-reference capture stream *before* starting far-end playback. After any XRUN recovery, log a stream discontinuity event, report an exact lost-sample count only if hardware frame counters or monotonic timestamps uniquely determine the gap (recording the derivation basis), and notify AEC, resamplers, and trackers to reset or re-align—calling `snd_pcm_prepare()` restarts DMA without magically restoring lost acoustic history. No Linux soundcard or hot-plug hardware tests were executed here.

### I03: PipeWire's Four AEC Streams

PipeWire's canonical repository resides on freedesktop GitLab, mirrored officially at `PipeWire/pipewire` on GitHub. Key entry points are [`src/modules/module-echo-cancel.c`](https://github.com/PipeWire/pipewire/blob/62316ae8cf659ca8e5f1ee866ae743d2438d47a4/src/modules/module-echo-cancel.c) and [`spa/plugins/aec/aec-webrtc.cpp`](https://github.com/PipeWire/pipewire/blob/62316ae8cf659ca8e5f1ee866ae743d2438d47a4/spa/plugins/aec/aec-webrtc.cpp). Root `COPYING` is MIT, while optional WebRTC audio-processing backends carry their own licenses.

[PipeWire 1.6.9 `module-echo-cancel` documentation](https://docs.pipewire.org/page_module_echo_cancel.html) separates four distinct graph streams: physical microphone `capture`, virtual echo-cancelled `source` exposed to applications, virtual `sink` receiving application audio, and physical loudspeaker `playback`. Any system audio intended for cancellation must traverse the designated reference path; `monitor.mode` alters how the reference tap is acquired, so `audio.rate`, channel map, target node bindings, and `node.latency` must all be archived alongside a dump of the active PipeWire graph.

In pinned `meson.build`, PipeWire selects dynamically among `webrtc-audio-processing-2`, `-1`, or legacy system libraries at build time—never assume PipeWire runs the exact same WebRTC AEC3 commit pinned elsewhere in this book merely because both say "WebRTC." Moreover, in this pinned commit, [`aec-webrtc.cpp`](https://github.com/PipeWire/pipewire/blob/62316ae8cf659ca8e5f1ee866ae743d2438d47a4/spa/plugins/aec/aec-webrtc.cpp) enables high-pass filtering and noise suppression by default, disables AGC, and passes `(num_blocks - 1) * 10` ms as the nominal stream delay to WebRTC APM rather than querying hardware DAC/ADC timestamps! Consequently, any benchmark comparing a PipeWire virtual mic against direct in-app WebRTC APM must record the distro library version, submodule toggles, `num_blocks`, reference routing, and exact output tap rather than attributing every difference to the linear AEC core.

Verify that VoIP speech, OS notification chimes, and media playback all enter the reference sink, and re-verify after switching default output devices (e.g., speakers to HDMI/Bluetooth). Inject reference dropouts and playout-delay steps to observe recovery. Crucially, when the OS sound server already runs AEC, enabling a second in-app AEC alters residual statistics and compounds near-end speech distortion—always compare single-AEC vs. cascaded-AEC configurations rather than merely checking that a virtual microphone node appears.

When an application embeds pinned WebRTC Audio Processing Module (APM) directly instead of relying on PipeWire, follow [`api/audio/audio_processing.h`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/api/audio/audio_processing.h): feed far-end render frames destined for the loudspeaker into `ProcessReverseStream()`, call `set_stream_delay_ms()` with the current hardware buffer delay immediately before each capture frame, and then invoke `ProcessStream()`. Externally, APM operates on `10` ms linear PCM frames—channel-interleaved for `int16_t` APIs and deinterleaved per-channel pointers for `float` APIs.

At pinned commit `0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e`, the `int16_t` overload restricts sample rates strictly to `8, 16, 32, 48` kHz (`NativeRate`), requiring identical rates across capture input, capture output, and reverse render streams plus matching capture I/O channel layouts; the `float` overload accepts arbitrary valid rates and distinct I/O channel counts. Thus `44.1` kHz integer PCM cannot enter the `int16_t` overload directly without resampling or switching to the `float` interface; see [`Initialize()` constraints and `NativeRate` enum](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/api/audio/audio_processing.h#509).

Note that `set_stream_delay_ms()` represents the sum of render buffering *after* `ProcessReverseStream()` until physical DAC output plus capture buffering *after* physical ADC input until `ProcessStream()`—never room acoustic time-of-flight or WAV header metadata.

In pinned [`audio_processing_impl.cc`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/audio_processing_impl.cc), delay arguments outside `[0, 500]` ms are clamped to `[0, 500]` ms with a warning return code; always log the requested delay, return status, and read-back `stream_delay_ms()`. To compare linear adaptive cancellation against post-nonlinear suppression, set `echo_canceller.export_linear_aec_output = true` and call `GetLinearAecOutput()` right after `ProcessStream()`: because that exported linear tap is fixed at `16` kHz (`160` samples per `10` ms block), never splice or diff it naively at `48` kHz device rate! See [Enhancement Guide A06](02_aec_wpe_separation.md#sec-u-dd6f9a17f8) for full AEC3 internals.

**Audit playout-volume notifications separately from reference waveforms.** In pinned WebRTC, `CreatePlayoutVolumeChange()` posts an unnormalized integer volume notification; on the capture side, detecting a change propagates a boolean `gain_change` flag rather than applying sample-wise gain compensation to the reference buffer. Inside [`refined_filter_update_gain.cc` line 52](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/aec3/refined_filter_update_gain.cc#52), the `gain_change` branch remains an explicit `TODO` alongside ERLE-estimator state resets—so it must never be described as automatic filter rescaling! Always record whether volume scaling occurs upstream or downstream of the AEC reference tap alongside exact sample indices of volume steps; see [Enhancement Guide A06](02_aec_wpe_separation.md#sec-u-dd6f9a17f8) (static source re-verification only; zero new hardware or C++ runs).

### I04: Stateful Sample-Rate Conversion in `libsamplerate`

`libsndfile/libsamplerate` identifies version `0.2.2` in its README under the `BSD-2-Clause` license. Core state lifecycle functions `src_new`, `src_process`, and `src_reset` reside in [`src/samplerate.c`](https://github.com/libsndfile/libsamplerate/blob/0844c208f683527c08ea8a80acc13b398aa9c8bf/src/samplerate.c), backed by polyphase sinc/interpolator kernels under `src/`.

The [Full API specification](https://libsndfile.github.io/libsamplerate/api_full.html) defines `src_ratio = output_sample_rate / input_sample_rate` and returns per-call `input_frames_used` and `output_frames_gen`. Callers must neither assume 100% consumption of every input buffer nor instantiate a fresh `SRC_STATE` per chunk. On the final input buffer, set `end_of_input = 1` and continue invoking `src_process` on the same state with zero new input until `output_frames_gen == 0` to drain internal FIR history; merely rewriting a WAV header's sample rate or halting as soon as input samples run out violates the streaming protocol.

When `SRC_DATA.src_ratio` changes between consecutive `src_process` calls, `libsamplerate` linearly interpolates the conversion ratio across the output block; conversely, calling `src_set_ratio()` forces an immediate step to the new ratio on the next block, bypassing linear slewing. Because the two APIs serve distinct control objectives, updating a ratio field never substitutes for blind clock-drift estimation or artifact-free phase tracking; see the [`Set Ratio` section of the Full API](https://libsndfile.github.io/libsamplerate/api_full.html).

Under this book's SRO sign convention, a capture clock running `+100` ppm fast (`16001.6` Hz vs. `16000` Hz) requires resampling ratio `1 / 1.0001 ≈ 0.99990001` to align onto the nominal reference clock. Validate whole-utterance vs. chunked streaming residuals after compensating fixed FIR group delay, measure stopband alias rejection near Nyquist, and test abrupt steps, slewed ratio updates, and multi-hour FIFO occupancy—a matching file duration over a few seconds never proves long-term phase coherence.

To verify the C streaming contract concretely, this repository executes a [re-runnable `libsamplerate` SRO experiment](../../ch10/examples/run_libsamplerate_sro.py) on a 10 s mono synthetic impulse train: relative to a `16` kHz reference clock, the simulated capture clock runs `+100` ppm fast (`16001.6` samples/s) during `0–5` s and `+150` ppm fast (`16002.4` samples/s) during `5–10` s, with impulses placed at exact reference instants `1, 3, 4, 6, 8, 9` s. Integrating piecewise clock rates and rounding yields `160020` captured frames, with an analytical uncompensated drift increment of `16` frames between the `1` s and `9` s impulses. Crucially, the piecewise clock rates here are **oracle ground truth supplied directly to the resampler**, NOT estimated blindly from audio.

Compiling pinned `libsamplerate` `0.2.2` (`0844c208f683527c08ea8a80acc13b398aa9c8bf`) as an isolated static library in a temporary directory on macOS arm64 (AppleClang 21.0.0, CMake 4.3.2), the harness runs `SRC_SINC_MEDIUM_QUALITY` in streaming mono float32 with `257`-frame input chunks, `64`-frame output buffers, and persistent cross-chunk state. The uncompensated control holds `src_ratio = 1.0`, whereas the compensated run passes `16000 / 16001.6 ≈ 0.99990001` during the first 5 s and `16000 / 16002.4 ≈ 0.99985002` during the second 5 s, setting `end_of_input = 1` on the final chunk and draining until zero output frames remain (exercising `SRC_DATA.src_ratio` linear slewing rather than `src_set_ratio` step jumps):

| Condition | Total Output Frames | Peak Sample-Index Offsets of the 6 Impulses vs. Reference Clock (frames) | Drift Increment Over `1–9` s (frames) | Partial-Consumption Calls | Tail Frames Drained After Input Exhausted |
|---|---:|---|---:|---:|---:|
| Uncompensated (`ratio = 1`) | 160 020 | `2, 5, 6, 10, 15, 18` | 16 | 1 875 | 84 |
| Oracle Clock-Ratio Compensation | 159 999 | `0, 0, 0, 0, 0, 0` | 0 | 1 874 | 127 |

Here "partial-consumption calls" counts invocations where `input_frames_used < input_frames`, with cumulative consumed input verified at exactly `160020` frames in both runs alongside per-call ratio and drain-order assertions. Notice that the compensated stream emits `159999` frames (`1` frame shy of `10 s × 16 kHz = 160000`)—never claim output length equals `160000` merely because all six interior impulse peaks land at `0` offset! Only isolated impulse peak indices are measured here (without post-hoc delay/gain fitting, inter-impulse waveform scoring, stopband rejection checks, or listening tests).

Re-run via `.venv/bin/python codes/chapters/ch10/examples/run_libsamplerate_sro.py`, which verifies upstream commit/worktree cleanliness, builds in `/tmp`, emits JSON, and leaves zero binaries or trace artifacts inside the repository; [standalone failure-mode unit tests](../../../../tests/test_codes_libsamplerate_sro.py) verify detection of inverted ratios, dropped unconsumed tails, premature termination, and omitted EOF draining. Real hardware still requires testing blind SRO estimation noise, ratio-update artifacts, anti-aliasing stopbands, timestamp jitter loops, and multi-hour queue stability.

A companion pure-NumPy [estimate–compensate–verify teaching loop](../../ch10/sro_closed_loop_demo.py) runs inside `.venv` without external C builds or soundcards, delegating `fit_delay_ppm` directly to the [single canonical `estimate_sro_ppm` kernel](../../ch10/core/engineering.py) (preserving the exact numerical results originally logged on 2026-09-28). Synthesizing `12` s of analytical waveforms with **noiseless per-sample synthetic timestamps**—reference rate `2000` Hz, device rate `2000.3` Hz (`+150` ppm), initial arrival offset `+2` ms after reference $t=0$, and a single dropped sample at physical device index `12000`—nominal reference index $n$ sits at $n/2000$ s while pre-drop device index $n$ sits at $0.002 + n/2000.3$ s. Over nine anchors in the first `2` s, linear regression recovers intercept `-2` ms and first-order slope `149.9775` ppm, which inverts via $s/(1-s)$ to exact `150` ppm (`+2` ms initial offset)—after first detecting integer-multiple timestamp gaps to restore physical sample indices and restricting regression strictly to **pre-drop** anchors so a dropped-sample step is never misfit as clock drift.

Calibrating on device indices `0..4000`, the streaming linear interpolator preserves cross-block sample history and fractional phase to emit reference samples starting at index `4500`; under `257`-sample chunks, the single-sample gap marks reference indices `12002` and `12003` invalid rather than interpolating across missing data. Against the analytical reference waveform, dimensionless squared-error MSE compares raw index-by-index alignment (confounded by `+2` ms offset, `+150` ppm drift, and the post-drop `1`-sample slip) against compensated valid samples:

| Reference Sample Interval | Valid Compensated Samples | Uncompensated MSE | Compensated MSE |
|---|---:|---:|---:|
| `5 000–10 999` (pre-drop) | 6 000 | 0.02869 | 0.000002773 |
| `13 000–21 999` (post-drop) | 9 000 | 0.02131 | 0.000002477 |

Run `.venv/bin/python codes/chapters/ch10/sro_closed_loop_demo.py` to reproduce all parameters, estimates, gap masks, and segment MSEs, backed by [unit tests](../../../../tests/test_codes_sro_closed_loop.py) verifying exact slopes, `2` ms intercept, single-sample drop recovery, and `257`- vs. `509`-sample chunk equivalence. Because timestamps are noiseless, drift is constant, and residual MSE includes 2-tap linear interpolation error, neither this synthetic closed loop nor the oracle-ratio `libsamplerate` run substitutes for real hardware clock-recovery acceptance.

On 2026-09-23, a read-only host audio enumeration check (`system_profiler SPAudioDataType -json` and CoreAudio property queries, without opening microphones or recording audio) found a single nominal `48` kHz mono `"MacBook Air Microphone"` and `"MacBook Air Speakers"` sharing an identical non-zero `kAudioDevicePropertyClockDomain` (noting also that restricted sandboxed subprocesses returned an empty device list, illustrating why empty sandboxed queries never prove hardware absence). Per [Apple CoreAudio `clockDomain` documentation](https://developer.apple.com/documentation/coreaudio/audiohardwareclock/clockdomain), devices sharing a non-zero clock domain are hardware-synchronized and cannot serve as two independent asynchronous crystal clocks; lacking a second independent capture device on that host, zero real-hardware ppm figures were fabricated. Real two-device SRO measurement must inspect `clockDomain`, disable hidden OS/Aggregate-Device resampling, log raw frame counters and host timestamps on both streams, and separate initial offset, steady slope, and XRUN steps across multiple windows.

Note finally that while pinned WebRTC AEC3's [`render_delay_controller.cc`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/aec3/render_delay_controller.cc) exposes `HasClockdrift()` to `block_processor.cc`, detecting buffer-delay drift or resetting adaptive echo paths is *not* equivalent to resampling two asynchronous clocks onto a common sample grid; cross-clock-domain AEC products must still validate true ASRC closed-loop stability.

### I05: SpeexDSP Streaming Resampler and Echo Ring Buffers

In `SpeexDSP` (`BSD-3-Clause`, pinned in `SOURCES.lock.json`), `libspeexdsp/resample.c` and `include/speex/speex_resampler.h` implement an independent polyphase resampler, while `libspeexdsp/mdf.c` implements partitioned-block frequency-domain AUMDF echo cancellation (contrasted against Chapter 6's time-domain NLMS). Crucially, `SpeexResamplerState` and `SpeexEchoState` are completely separate structs—linking SpeexDSP never automatically compensates playback/capture SRO inside the AEC!

Resampler setup pins input/output rates, channels, quality `[0..10]`, `int16` vs. `float` API, and per-call buffer lengths; AEC setup separately pins `frame_size`, `filter_length` (tail length in samples), sample rate, and render reference buffering. Per [pinned `speex_echo.h`](https://gitlab.xiph.org/xiph/speexdsp/-/blob/8e29a256ef0235ebbe7fcb8417b5ac7731eb8307/include/speex/speex_echo.h) and the [Speex Manual](https://www.speex.org/docs/manual/speex-manual/node7.html), synchronous `speex_echo_cancellation()` expects the caller to supply already-aligned capture and playback frames without adding internal playback delay, whereas asynchronous helpers `speex_echo_playback()` / `speex_echo_capture()` insert a fixed internal two-frame ring buffer—requiring that every far-end frame reach the AEC *strictly before* its acoustic echo arrives at the microphone while avoiding excessive causal bulk delay that wastes FIR tail taps.

Watch out for a classic silent pitfall in pinned `mdf.c`: `speex_echo_state_init()` initializes `sampling_rate = 8000` Hz! Feeding `16` kHz WAV audio without explicitly calling `speex_echo_ctl(st, SPEEX_ECHO_SET_SAMPLING_RATE, &rate)` (and verifying via `SPEEX_ECHO_GET_SAMPLING_RATE`) misconfigures internal notch/adaptation time constants; see [upstream `mdf.c`](https://gitlab.xiph.org/xiph/speexdsp/-/blob/8e29a256ef0235ebbe7fcb8417b5ac7731eb8307/libspeexdsp/mdf.c). Likewise, upstream `testecho.c` is a minimal API demo that ignores `fread` short reads at EOF; rigorous evaluation requires exact PCM16 frame accounting, as implemented in our [macOS arm64 synchronous real-pair Speex experiment](02_aec_wpe_separation.md#aec) (which tests synchronous 1-mic/1-speaker execution only, not asynchronous or hardware loops).

For multi-channel setups (`speex_echo_state_init_mc(frame_size, filter_length, nb_mic, nb_speakers)`), static audit of [lines 687–733 of pinned `mdf.c`](https://gitlab.xiph.org/xiph/speexdsp/-/blob/8e29a256ef0235ebbe7fcb8417b5ac7731eb8307/libspeexdsp/mdf.c#L687-733) uncovers a severe upstream bug in the asynchronous wrappers: while synchronous `speex_echo_cancellation()` correctly indexes interleaved speakers and mics via `i*K + speak` and `i*C + chan`, asynchronous `speex_echo_playback()` copies only `frame_size` scalar samples (omitting `* nb_speakers`), and `speex_echo_capture()` copies only `frame_size` scalar samples on ring-buffer underflow (omitting `* nb_mic`)! Consequently, **pinned asynchronous `speex_echo_playback`/`capture` MUST NOT be fed multi-channel interleaved buffers**; applications must maintain their own full-frame interleaved reference ring buffer and call synchronous `speex_echo_cancellation()`—verifying channel ordering with impulses at *distinct sample indices per channel* (identical simultaneous impulses mask channel swaps!), and noting that `SPEEX_ECHO_GET_IMPULSE_RESPONSE` is marked unimplemented for multi-channel in this commit.

Even in mono, asynchronous `speex_echo_playback()` silently drops playback frames arriving prior to the first `speex_echo_capture()` call, and `speex_echo_capture()` copies raw microphone input straight to output with a warning on reference underflow; see [asynchronous entry points in `mdf.c`](https://gitlab.xiph.org/xiph/speexdsp/-/blob/8e29a256ef0235ebbe7fcb8417b5ac7731eb8307/libspeexdsp/mdf.c) and [Speex debugging notes](https://www.speex.org/docs/manual/speex-manual/node7.html). Comparing synchronous vs. asynchronous runs on identical PCM while logging playback indices, ring occupancy, warnings, and exact `out == in` pass-through frames immediately separates caller queue bugs from AUMDF convergence speed.

Finally, Speex AUMDF regulates step size continuously via estimated residual-to-error leakage rather than a binary double-talk detector (DTD); enabling nonlinear residual echo suppression requires binding `SpeexEchoState` into `SpeexPreprocessState` via `SPEEX_PREPROCESS_SET_ECHO_STATE` so `preprocess.c` queries `speex_echo_get_residual()`. Always tap and archive *both* linear MDF output and post-preprocessor output—never inferring linear tap convergence from aggressive nonlinear suppression! For same-input three-AEC synthetic ground-truth comparisons, see [Enhancement Guide A07](02_aec_wpe_separation.md#aec).

## 2. Noise Estimation, Voice Activity, and Gain Control

Chapter 10's [power spectral subtraction hand derivation and Exercise E10-13](../../../../chapters_en/10_engineering-practice.md#sec-10-3-2) fix initial noise-only frames, per-bin mean noise power, and an observation-relative spectral floor inside [this book's baseline `noise_suppression.py`](../../ch10/core/noise_suppression.py). Industrial modules below—WebRTC VAD, AGC2, Silero VAD, and DeepFilterNet—rely on fundamentally different state machines and data priors whose perceptual quality or CPU load can never be extrapolated from stationary spectral subtraction on synthetic tones.

### From Stationary Noise Floors to Controlled Online Noise Tracking

When background noise rises mid-call, freezing an initial noise estimate $D_f$ leaves unattenuated residual noise, whereas blindly averaging every active speech frame into $D_f$ cancels target speech. **MCRA** (*Minima Controlled Recursive Averaging*) solves this by modulating the **update rate of the noise power estimator** itself: reading page 13, Sections II–III, Eqs. (3)–(14) of Cohen and Berdugo's [2002 IEEE SPL paper](https://israelcohen.com/wp-content/uploads/2018/05/SPL_Jan2002.pdf) (verified 2026-10-01) traces conditional noise updates, smoothed speech-presence probability, time-frequency local averaging, two-window running minima tracking, and ratio thresholding—far richer than multiplying a final Wiener gain by a binary VAD bit.

Given speech-presence probability $p_t \in [0, 1]$ and baseline noise smoothing factor $\alpha_d$, MCRA forms the time-varying smoothing factor $\widetilde{\alpha}_t = \alpha_d + (1 - \alpha_d)p_t$. For $\alpha_d = 0.8$, previous noise power `2`, and current noisy power `10` (in identical linear power units at a single bin), $p_t = 0$ updates noise power to `3.6`, whereas $p_t = 0.9$ slows the update to `2.16`. That single-step conditional recursion is exercised over multi-frame sequences in [Exercise E10-28](../../../../chapters_en/10_engineering-practice.md#e10-28), while full minima-tracking and two-iteration IMCRA bias compensation are analyzed in [Spatial Processing §§50–51](01_spatial_and_tracking.md#sec-u-2973a53a28).

Not every C library implements soft-probability MCRA: in pinned SpeexDSP (`preprocess.c`), `preproc_update_prob` sets a hard binary flag per bin whenever $0.4 S > S_{\min}$, and `preproc_noise_update` updates the noise estimate only when that flag is `0` **or when current power falls below the previous noise estimate**. Our bounded native helper test in I32 confirms downward noise tracking even while flagged, without mislabeling Speex's heuristic as full Cohen MCRA/IMCRA. (Visiting [Israel Cohen's software page](https://webee.technion.ac.il/Sites/People/IsraelCohen/Info/Software.html) yielded no license-verified redistributable archive during our check, so we cite the published equations without redistributing author MATLAB scripts.)

### I06: WebRTC Legacy GMM VAD

Located at [`common_audio/vad/include/webrtc_vad.h`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/common_audio/vad/include/webrtc_vad.h) alongside `webrtc_vad.c` and `vad_core.c` (`BSD-3-Clause` with third-party notices), this classic 6-band Gaussian Mixture Model VAD is distinct from newer RNN-based voice probability estimators inside WebRTC APM—always specify exact header/symbol names rather than saying "WebRTC VAD" generically.

The C API accepts mono signed 16-bit PCM strictly at rates and frame lengths validated by `WebRtcVad_ValidRateAndFrameLength` (`10, 20, 30` ms—i.e., `160, 320, 480` samples at `16` kHz), with aggressiveness modes `0..3` adjusting likelihood-ratio thresholds rather than dBFS energy cutoffs. Acceptance testing must sweep false-alarm vs. miss trade-offs across modes `0..3`, onset/offset word clipping, and stationary noise triggers—taking care to scale, round, and saturate normalized `[-1, 1]` floats to `[-32768, 32767]` (naive `(int16_t)x` C casting zeros everything in `(-1, 1)`!) and documenting whether VAD taps raw mic 0 or beamformed output.

### I07: WebRTC AGC2

Defined in [`modules/audio_processing/gain_controller2.h`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/gain_controller2.h) and `modules/audio_processing/agc2/`, WebRTC's second-generation automatic gain controller combines fixed digital gain, adaptive speech-level estimation, noise-floor gating, and lookahead limiting—going well beyond Chapter 10's peak-limiter baseline.

Deployments must log fixed digital gain, adaptive digital gain settings, whether analog/OS input-volume control is enabled (remembering that analog PGA changes ahead of the ADC alter the echo path seen by AEC!), sample rate, and channel count, re-initializing state whenever rate or device topology changes. Evaluate level steps, quiet utterances after prolonged silence, impulsive clicks, pre-clipped inputs, and double-talk—comparing noise floors at *matched output speech loudness* rather than mistaking expected linear gain makeup of background noise for NR failure.

### I08: Silero VAD Chunk Sizes and Recurrent State

`snakers4/silero-vad` (`MIT` code license; verify pinned ONNX/JIT model assets separately) separates low-level `OnnxWrapper`, streaming `VADIterator`, and offline `get_speech_timestamps` inside [`src/silero_vad/utils_vad.py`](https://github.com/snakers4/silero_vad/blob/60b7ffa243625ebdc1070275a29f18c8786a/src/silero_vad/utils_vad.py).

At the pinned revision, `OnnxWrapper` enforces exact `32` ms chunks—`512` samples at `16` kHz or `256` samples at `8` kHz—while carrying internal LSTM state and a context tail across consecutive chunks: never reset state per chunk within a stream, never share one wrapper instance across concurrent streams without state multiplexing, and document onset/offset hysteresis thresholds, minimum speech/silence durations, and pre-/post-roll padding separately from raw frame probabilities. Note also that upstream's convenience integer-multiple decimation shortcut (`sr // 16000` strided slicing) lacks an anti-aliasing low-pass filter, so hardware `48` kHz streams should pass through proper bandlimited SRC first, buffering `10` ms (`160`-sample) driver callbacks into `512`-sample model chunks with explicit timestamp bookkeeping.

### I09: DeepFilterNet End-to-End Execution Paths

`Rikorose/DeepFilterNet` (`Apache-2.0 OR MIT` dual-licensed code; model weights and corpora governed separately) cleanly separates Python training/offline scripts (`DeepFilterNet/df/`), the Rust `libDF` inference engine ([`libDF/src/`](https://github.com/Rikorose/DeepFilterNet/tree/d375b2d8309e0935d165700c91da9de862a99c31/libDF/src)), and real-time Linux LADSPA plugin wrappers ([`ladspa/`](https://github.com/Rikorose/DeepFilterNet/tree/d375b2d8309e0935d165700c91da9de862a99c31/ladspa)).

Targeting `48` kHz full-band speech enhancement, deployments must record FFT length, hop size, ERB/DF lookahead frames, deep-filtering tap order, and delay-compensation flags—noting that offline `--compensate-delay` file slicing eliminates zero latency in live streaming! In pinned `ladspa/README.md`, even zero-lookahead causal configurations incur a `20` ms minimum STFT framing delay *before* adding PipeWire/OS quantum buffering (`ladspa/filter-chain-configs/deepfilter-mono-source.conf`, built via `cargo build --release -p deep-filter-ladspa` when Rust toolchains are installed). Benchmark stationary/nonstationary noise, competing talkers, music, and silence against downstream WER, subjective quality, physical loopback latency, and release-build Rust CPU load rather than Python offline script wall-time; see [upstream repository](https://github.com/Rikorose/DeepFilterNet).

#### Blocking Sleeps and Latency-Counter Quirks in the Pinned LADSPA Plugin

Static inspection of [`ladspa/src/lib.rs`](https://github.com/Rikorose/DeepFilterNet/blob/d375b2d8309e0935d165700c91da9de862a99c31/ladspa/src/lib.rs#L213-L224) exposes crucial differences from hard-real-time audio callback rules: at initialization, each channel's output queue is pre-filled with one `hop` of zeros; inside `run()`, incoming samples are pushed onto a shared mutex-protected input queue while a background worker runs inference, and `run()` loops locking and polling the output queue—calling `thread::sleep` whenever output samples are not yet ready! Because mutex locking and `sleep` provide zero bounded-wait guarantee, wiring this LADSPA plugin into PipeWire does not automatically satisfy I01 hard-real-time callback rules. Moreover, lines `434–459` inject an extra `hop` of silence buffering whenever per-call `rtf >= 1.0`, triggering an explicit `panic!` if `proc_delay` reaches `1` second of samples (static source findings; Rust was not compiled or run here).

Even more subtly, lines `461–479` contain a concrete channel-vs.-sample axis bug in the **latency-reduction** branch: `o_q` is a `VecDeque` *per channel*, so `o_q.iter_mut().take(self.frame_size)` iterates at most `n_channels` queues and pops only **1 sample** (`pop_front()`) per channel, yet immediately decrements `proc_delay` by an entire `self.frame_size` (`480` samples)! On a 2-channel stream (`frame_size = 480`) starting with `960` queued samples per channel, our independent Python queue trace confirms `959` samples remain per channel (next sample index `1`) while `proc_delay` drops from `960` to `480` (whereas dropping a full hop would leave `480` samples starting at index `480`). Preserved untouched without claiming neural-core bugs or hardware measurements, this static audit is bound by SHA-256 in [I32](#industrial-upstream-audit).

#### Why Host Callback Counts Cannot Stand In for Model Hop Counts

Inspecting [lines 461–488 of `ladspa/src/lib.rs`](https://github.com/Rikorose/DeepFilterNet/blob/d375b2d8309e0935d165700c91da9de862a99c31/ladspa/src/lib.rs#L461-L488) reveals a second timing assumption in that same latency-reduction branch: it checks `t_proc_change > 10 * sr / frame_size` while incrementing `t_proc_change += 1` **once per host `run(sample_count)` callback** rather than accumulating `sample_count` or elapsed wall-clock seconds!

At `sr = 48000` Hz and model hop `frame_size = 480` samples, the threshold RHS equals `1000`. Starting from `t_proc_change = 0` under sustained `rtf < 0.5` with zero resets, callback `1001` enters the check with `t_proc_change == 1000` (failing strict `>`), so callback `1002` (`t_proc_change == 1001`) is the earliest call that triggers latency reduction—corresponding to wildly different elapsed audio durations depending on the host PipeWire/JACK quantum `sample_count`:

| Fixed Host Callback Block (`sample_count`) | Audio Advanced Over 1000 Callbacks (s) | Audio Advanced from Initial `0` Through Callback 1002 (s) |
|---:|---:|---:|
| 128 | 2.666667 | 2.672000 |
| 480 | 10.000000 | 10.020000 |
| 1024 | 21.333333 | 21.376000 |

(Note also that when an underrun resets `t_proc_change = 0` mid-function before the trailing `+= 1`, `t_proc_change` exits that call as `1`, reaching `1001` on the `1001`st subsequent stable callback; variable host quanta require summing `sample_count` per call.) Furthermore, upstream's local `rtf` timer wraps queue locking and `sleep` waits rather than pure inference compute, and its `"Underrun detected"` log string reflects worker lag rather than hardware ALSA `xrun_count`.

## 3. Embedded DSP Firmware and Inter-Module Control Coupling

### I10: XMOS `lib_voice` AEC and Delay Control

`xmos/fwk_voice` officially migrated to [`xmos/lib_voice`](https://github.com/xmos/lib_voice) (`README` v1.1.0, specifying XTC Tools `15.3.1`, Python `3.11`, and `lib_xcore_math`). Core AEC and Automatic Delay Estimation & Control (ADEC) modules reside under `lib_voice/api/aec/`, `lib_voice/src/aec/`, and `adec/`.

Configurations must pin microphone and reference-speaker channel counts, frame advance, main/shadow filter tail partitions, BFP memory pools, xCORE hardware thread mapping, and reference bulk delay—testing ADEC re-alignment triggers, shadow-to-main swaps, playback-delay steps, acoustic path changes, silent references, and double-talk rather than swapping in textbook NLMS while keeping vendor control constants. Note that [`LICENSE.rst`](https://github.com/xmos/lib_voice/blob/c9f1a9bf95cd88c7950adf4bf631c217f900ad25/LICENSE.rst) is the **XMOS Public Licence v1**, restricting commercial deployment to XMOS silicon with special-use exclusions—valuable for vendor firmware architecture study, but never classifiable as hardware-agnostic MIT/BSD open source.

**Post-silence AEC hold counting applies strictly to the non-default `ALT` architecture branch.** In pinned [`stage1.h`](https://github.com/xmos/lib_voice/blob/c9f1a9bf95cd88c7950adf4bf631c217f900ad25/lib_voice/api/stage1/stage1.h), `ALT_ARCH_MODE` defaults to `0`; inside [`stage1.c` line 49](https://github.com/xmos/lib_voice/blob/c9f1a9bf95cd88c7950adf4bf631c217f900ad25/lib_voice/src/stage1/stage1.c#L49), only when `ALT_ARCH_MODE` is enabled does control flow pass reference frames through the delay buffer before testing far-end activity against a hold limit of $16000 \times 3 / 240 = 200$ frames (`15` ms per `240`-sample advance). Because upstream evaluates strict `count > limit` *before* `count++`, starting from active AEC (`count = 0`) keeps AEC active—and forces downstream reference-activity flags high—across the first **201** consecutive inactive reference frames, bypassing AEC only on frame **202**! Thus `"3 s"` is a nominal configuration constant rather than an exact 200-frame cutoff or universal default across `ALT_ARCH_MODE=0`, and governs AEC/IC routing hold rather than physical RT60 or double-talk detection (static source audit only; zero xCORE execution, and unreadable external XVF3510 PDFs are not cited as code evidence).

### I11: XMOS Interference Canceller (IC) and VNR Control Coupling

Under [`lib_voice/src/`](https://github.com/xmos/lib_voice/tree/c9f1a9bf95cd88c7950adf4bf631c217f900ad25/lib_voice/src) (`ic/` and `vnr/`, same version and XMOS license as I10), adaptive Interference Cancellation (IC) relies on Voice-to-Noise Ratio (VNR) estimates to gate filter adaptation—remembering that VNR is a noisy statistical ratio rather than an oracle VAD bit, so target speech leaking into the blocking/reference channel can cause target cancellation unless adaptation gating is swept across target angles, distances, interferer directions, and mic gain mismatches.

Notice in pinned `lib_voice/api/ic/ic_defines.h` that `IC_FRAME_ADVANCE = 240` (`15` ms of new audio at `16` kHz) whereas `IC_FRAME_LENGTH = 512` (`32` ms analysis window)—never confuse analysis window length with per-call hop latency! Archiving auto-adapted, adaptation-frozen, and bypassed outputs side-by-side cleanly isolates linear filter distortion from VNR gating errors, including re-convergence lag after near-end speech ends.

### I12: XMOS NS, AGC, and Inter-Stage Metadata Plumbing

Inspecting `lib_voice/api/ns/`, `api/agc/`, and `api/stage1/` ([API tree](https://github.com/xmos/lib_voice/tree/c9f1a9bf95cd88c7950adf4bf631c217f900ad25/lib_voice/api)) shows that embedded voice pipelines pass structured metadata—VAD/VNR flags, AEC residual/convergence states, and gain states—alongside PCM buffers. Porting a standalone NS or AGC function without reconstructing its upstream metadata contract silently alters its operating mode; evaluate NS-only, AGC-only, and NS+AGC combinations against equal-gain unprocessed capture on noise pumping, word-tail clipping, level steps, and double-talk.

### I13: Sound Open Firmware (SOF) Fixed Time-Domain Beamforming (TDFB)

Sound Open Firmware (`thesofproject/sof`, predominantly `BSD-3-Clause` with mixed per-file notices documented in `LICENCE`) implements embedded time-domain filter-and-sum beamforming under [`src/audio/tdfb/`](https://github.com/thesofproject/sof/tree/b6c6a05d52536313fe8e8752b1c4e069b1cc4002/src/audio/tdfb) and [`tdfb_generic.c`](https://github.com/thesofproject/sof/blob/b6c6a05d52536313fe8e8752b1c4e069b1cc4002/src/audio/tdfb/tdfb_generic.c).

Unlike online covariance-adaptive MVDR, SOF TDFB applies offline-designed FIR filter banks selected by topology blobs, array geometry, and steering-angle controls. Verify each mic-to-beam FIR impulse response directly, sweep acoustic incidence angles, inject channel swaps and gain mismatches, and cross-check `tdfb_generic.c` against HiFi/SIMD DSP intrinsics within fixed-point tolerances.

Crucially, audit the offline MATLAB design scripts before trusting generated blobs! In Chapter 5, our [static audit script](../../ch05/examples/audit_sof_tdfb_design.py) and [initial report](../../ch05/reports/sof_tdfb_design_audit.json) identified two mathematical discrepancies in pinned `sof_bf_design.m`—an unnormalized vs. normalized `sinc` argument convention and a WNG denominator variable mismatch—and our 2026-10-04 [current SOF contract report](../../ch05/reports/sof_tdfb_design_current.json) further documented an optional azimuth-sweep simulation loop that re-initializes `nmi` inside the loop body, retaining only the final azimuth's result (all verified statically with scalar Python math checks without running MATLAB/Octave, exporting FIR blobs, or flashing DSP hardware).

Between ideal frequency-domain weights and physical DSP output lie per-mic causal delay alignment, finite FIR tap windowing, int16 coefficient quantization, multi-branch accumulation, and format-specific rounding/saturation (`S16`/`S24` rounding+saturation vs. `S32` shift+saturation and output channel bitmasks in `tdfb_generic.c`)—see [Spatial Processing §60](01_spatial_and_tracking.md#per-channel-mask-evidence).

### I14: SOF Firmware Sample-Rate Conversion (`src`) and Topology Routing

Under [`src/audio/src/`](https://github.com/thesofproject/sof/tree/b6c6a05d52536313fe8e8752b1c4e069b1cc4002/src/audio/src), `src` stands for **Sample Rate Conversion** (not generic "source code"!), implementing fixed-ratio firmware polyphase conversion—distinct from asynchronous clock-drift recovery or hardware DAI synchronization. Adding an `src` topology widget alters group delay, SRAM footprint, and passband/stopband edges; verify supported vs. rejected rate pairs, block boundaries, stop/restart resets, near-Nyquist aliasing, and multi-channel phase lock without post-hoc per-channel time alignment.

### I15: ARM CMSIS-DSP `Q15` Fixed-Point FIR

[`Source/FilteringFunctions/arm_fir_q15.c`](https://github.com/ARM-software/CMSIS-DSP/blob/83a2d7bc98c81b4bbe4a6f48b1f2ecf179868a0b/Source/FilteringFunctions/arm_fir_q15.c) in `ARM-software/CMSIS-DSP` (`Apache-2.0`) provides an industry-standard fixed-point baseline—multiplying `1.15` operands into a `34.30` (`q63_t`) accumulator before right-shifting by `15` (`>> 15` floor truncation) and saturating to `[-32768, 32767]`. Because this textbook's teaching `q15_dot` uses round-to-nearest-even, bit-exact parity requires pinning rounding mode, tap order, state-buffer length, and SIMD macro branches explicitly rather than treating `"Q15"` as self-specifying.

Specifically, [`arm_fir_init_q15.c`](https://github.com/ARM-software/CMSIS-DSP/blob/83a2d7bc98c81b4bbe4a6f48b1f2ecf179868a0b/Source/FilteringFunctions/arm_fir_init_q15.c) requires `numTaps` to be an **even integer $\ge 4$** with coefficients stored in **time-reversed order**! Moreover, required `pState` buffer length varies by build configuration: `numTaps + blockSize - 1` under scalar C (no `ARM_MATH_LOOPUNROLL`), `numTaps + blockSize` when `ARM_MATH_LOOPUNROLL` is defined, and `8`-tap zero-padded multiples under Helium MVE—never size `pState` from one `#ifdef` branch and deploy on another!

In I32's native contract test, compiling unmodified `arm_fir_q15.c` and `arm_fir_init_q15.c` under scalar C11 with a 4-tap causal response $[0.5, 0.5, 0, 0]$ (stored time-reversed) on `Q15` integer input $[1, 0, -1, 0, 1]$ produces $[0, 0, -1, -1, 0]$: arithmetic right-shift truncates $+0.5\text{ LSB} \to 0$ while flooring $-0.5\text{ LSB} \to -1$, contrasting sharply with round-to-even, alongside verified `-32768 * -32768` positive full-scale saturation (`32767`) and whole-buffer vs. chunked state equivalence (scalar host check only; zero Cortex-M/MVE cycle benchmarks).

## 4. Neural Inference Runtimes and Memory Budgets

### I16: ARM CMSIS-NN Quantized Kernels

`ARM-software/CMSIS-NN` (`Apache-2.0`; `Include/arm_nnfunctions.h`, `Source/`, `Tests/UnitTest/`) implements `int8`/`int16` neural operators targeting Cortex-M DSP/Helium extensions, bit-aligned against TensorFlow Lite Micro reference quantization ([official repository](https://github.com/ARM-software/CMSIS-NN)). Deployment audits must pin per-tensor vs. per-channel scale/zero-point parameters, activation clamping, scratch-buffer byte requirements, MVE/DSP compiler flags, and *exact front-end STFT/mel/log fixed-point scaling*—since bit-exact integer convolutions still degrade accuracy if front-end feature quantization drifts—and fail explicitly on unsupported ops or insufficient scratch RAM.

### I17: TensorFlow Lite Micro `tensor_arena` and Audio Front Ends

[`tensorflow/lite/micro/examples/micro_speech/`](https://github.com/tensorflow/tflite-micro/tree/9f638f18154dff868e7053572f089be845b2fbdf/tensorflow/lite/micro/examples/micro_speech) (`Apache-2.0` code; checkpoints/datasets licensed separately) illustrates bare-metal keyword-spotting execution without dynamic heap allocation: callers pre-allocate a contiguous byte array (`tensor_arena`) partitioned by TFLM's [memory planner](https://github.com/tensorflow/tflite-micro/blob/9f638f18154dff868e7053572f089be845b2fbdf/tensorflow/lite/micro/docs/memory_management.md) into reusable non-persistent activation `head` buffers, temporary scratch space, and persistent `tail` metadata/state—meaning `.tflite` flatbuffer size (e.g., an illustrative `100` KiB file) never bounds runtime arena SRAM, nor does the arena automatically include DMA audio ring buffers, FFT twiddle tables, or RTOS thread stacks. Profile peak high-water marks via `RecordingMicroInterpreter`, verify clean error codes on undersized arenas, and test continuous streaming across pause/resume and ADC overruns.

### I18: ONNX Runtime Threading, Execution Providers, and I/O Binding

Pinned `microsoft/onnxruntime` (`MIT`; C API in [`onnxruntime_c_api.h`](https://github.com/microsoft/onnxruntime/blob/48795e0281bcaa6b3b63b28af5e078b4c41e995f/include/onnxruntime/core/session/onnxruntime_c_api.h)) requires recording execution provider (EP), opset version, dynamic vs. static axes, feature normalization, and recurrent cache tensors alongside [intra-op vs. inter-op thread-pool settings](https://onnxruntime.ai/docs/performance/tune-performance/threading.html)—measuring session creation, cold first inference, warmed steady-state latency, CPU fallback copies, and multi-session thread contention rather than copying single-session offline thread counts blindly.

**Which copy overhead does I/O Binding actually eliminate?** Standard `Run()` calls with a GPU/NPU EP may copy host inputs to device memory and copy outputs back to host RAM on every `10–16` ms frame. Per [ONNX Runtime I/O Binding docs](https://onnxruntime.ai/docs/performance/tune-performance/iobinding.html), binding device-resident `OrtValue` buffers keeps recurrent streaming caches (`cache_in_*` / `cache_out_*`) on-device across frames when supported by the model signature and EP. Auditing pinned commit `48795e0281bcaa6b3b63b28af5e078b4c41e995f` highlights five C API contracts:

| C API Entry Point | Bound Objects to Verify | Mandatory Preconditions |
|---|---|---|
| `CreateIoBinding`, `BindInput` | Input tensor name, `OrtValue`, `OrtMemoryInfo` | Exact `dtype`, shape, and device match; caller buffer must outlive `RunWithBinding` |
| `BindOutput` | Pre-allocated fixed-shape output `OrtValue` | Output shape must fit pre-allocated capacity; never mistake stale buffer contents for fresh output on error |
| `BindOutputToDevice` | `OrtMemoryInfo` for dynamic-shape outputs | Delegates per-call device allocation to ORT—so binding to device does NOT guarantee zero per-frame allocations! |
| `RunWithBinding` | Bound inputs/outputs and `OrtRunOptions` | Log active EP, silent CPU-fallback nodes, and cross-frame state rollover |
| `SynchronizeBoundInputs`, `SynchronizeBoundOutputs` | Stream/queue synchronization barriers | Pinned header comments note backend-dependent behavior (no-op on some EPs)—never treat them as universal GPU timers! |

Keeping recurrent caches on-GPU still leaves DMA/host-to-device transfers for raw microphone PCM/STFT inputs and enhanced output frames, nor does `IoBinding` permit aliasing input and output pointers in-place without explicit operator guarantees. Note that this section is a pinned C API audit and three-way experimental design (default `Run` vs. persistent device `IoBinding` vs. forced per-frame host round-trip); zero GPU benchmarks were executed locally.

## 5. Evaluation Datasets and Scorer Reproducibility

### I19: Local `DNSMOS` Prediction

[`DNSMOS/dnsmos_local.py`](https://github.com/microsoft/DNS-Challenge/tree/591184a9fcb2cbdec02520fed81a32bbbf9d73ff/DNSMOS) in `microsoft/DNS-Challenge` (`MIT` under `LICENSE-CODE`; datasets governed separately) runs non-intrusive neural quality estimation locally. Always distinguish baseline `DNSMOS` (P.808), three-score `DNSMOS P.835` (`SIG, BAK, OVRL`), and personalized `pDNSMOS` (which penalizes un-enrolled competing talkers as background interference!), pinning ONNX file SHA-256, resampling rate, sliding-window framing, short-clip handling, and full-corpus failure accounting—never confusing neural MOS predictions with human listening panels or spatial localization fidelity.

### I20: `AECMOS` Input Alignment and Challenge Scoring Windows

[`AECMOS/AECMOS_local/`](https://github.com/microsoft/AEC-Challenge/tree/6c633d0a9d2a143a0e364899b91b06f127315b18/AECMOS/AECMOS_local) (`MIT` code; ONNX checkpoints and challenge recordings governed separately; web API deprecated per [official README](https://github.com/microsoft/AEC-Challenge/tree/6c633d0a9d2a143a0e364899b91b06f127315b18/AECMOS)) predicts separate echo degradation and near-end/other degradation scores from `(lpb, mic, enh)` triplets.

Inspecting [`aecmos.py` lines 60–81](https://github.com/microsoft/AEC-Challenge/blob/6c633d0a9d2a143a0e364899b91b06f127315b18/AECMOS/AECMOS_local/aecmos.py#L60-L81) at pinned commit `6c633d0a9d2a143a0e364899b91b06f127315b18` reveals three crucial silent behaviors: (1) all three WAVs are loaded via `librosa.load(..., sr=model_sr)` with default `mono=True` ([`librosa.load` 0.9.1 docs](https://librosa.org/doc-playground/0.9.1/generated/librosa.load.html)), silently downmixing multi-channel files to mono and resampling; (2) all three signals are truncated to their common minimum length; and (3) `run()` slices and scores **strictly the first `20` seconds** of inputs $\ge 20$ s! Moreover, scoring Interspeech 2021 or ICASSP 2022 blind test sets requires **pre-slicing all three streams synchronously on the original timeline** per [official cropping rules](https://github.com/microsoft/AEC-Challenge/blob/6c633d0a9d2a143a0e364899b91b06f127315b18/AECMOS/README.md)—retaining the second half for far-end single-talk, the trailing `(duration_s - 15) / 2` seconds for double-talk, and full duration for near-end single-talk—before checking whether any cropped clip still exceeds `20` s (static source audit only; zero `AECMOS` ONNX inference was run here).

### I21: CHiME-8 Data Preparation, Text Normalization, and Multi-Talker Meeting Scoring

`chimechallenge/chime-utils` (`MIT`, pinned at `152882404f572d40769ef02bf91c5a9a9cfc9c78`; [`Scoring` README](https://github.com/chimechallenge/chime-utils/blob/152882404f572d40769ef02bf91c5a9a9cfc9c78/README.md#scoring)) and underlying `fgnt/meeteval` (`MIT`, pinned at `6e3dc81284f2d6928f7ef9e620fd3b6906daa429`, verified 2026-09-28; meeting speech corpora licensed separately) standardize `SegLST` JSON formatting (`session_id`, `speaker`, `start_time`, `end_time`, `words`) and multi-speaker WER evaluation.

#### Output Stream Slots vs. Speaker Identities Govern Which WER Metric Applies

Section 3.2 Eq. (3), Sections 3.3–3.5, and Figure 1 of the [MeetEval paper (v3, 2024-01-25)](https://arxiv.org/pdf/2307.11394) separate speaker-attributed diarization streams from unattributed continuous speech separation (CSS) slots: **cpWER** optimizes a single global bijection between reference speakers and hypothesis streams across an entire session; **ORC-WER** (*Optimal Reference Combination WER*) assigns each intact reference utterance to an output stream while preserving intra-speaker/reference utterance order; and **MIMO-WER** generalizes assignment across multiple input/output streams while preserving each reference speaker's internal utterance order.

Because a 2-output CSS front end reuses its two output slots across $>2$ meeting participants over time ([CSS ICASSP 2020 paper §§1–2](https://www.microsoft.com/en-us/research/uploads/prod/2020/04/ICASSP2020__Continuous_speech_separation__dataset_and_analysis.pdf)), suppose speakers A and B overlap saying `"a"` and `"b"`, followed later by speaker C saying `"c"`, and CSS cleanly outputs slot 0 = `"a c"`, slot 1 = `"b"`. Every word is preserved without crosstalk, so **ORC-WER = 0**; yet feeding those two CSS slots directly into **cpWER** (padded with an empty third hypothesis speaker) incurs `1` deletion and `1` insertion on `"c"`, yielding $\text{cpWER} = 2/3$! Thus cpWER penalizes *lack of speaker diarization* rather than acoustic separation failure—and conversely, ORC-WER forbids splitting reference utterance `"b c"` post-hoc into `"b"` and `"c"` to artificially lower errors against `"a b"` / `"c"` (whose true ORC error count is `2`).

Nor does **DI-cpWER** measure speaker-identification accuracy directly: in pinned [`di_cp.py`](https://github.com/fgnt/meeteval/blob/6e3dc81284f2d6928f7ef9e620fd3b6906daa429/meeteval/wer/wer/di_cp.py), `greedy_di_cp_word_error_rate` swaps reference and hypothesis roles inside greedy ORC before swapping insertion/deletion counts back—yielding a greedy upper bound rather than an exact optimum (and note that the pedagogical `di_cp_lev` snippet inside upstream's doc/simplified module contains a parameter/variable typo raising `NameError` when invoked directly, an isolated doc-function bug kept distinct from the Cython-backed production entry point).

#### Utterance Granularity, Word Timestamps, and Diagnostic Metric Boundaries

Reading Sections IV-F, V-A Eq. (10), VI-A, and VII Algorithm 1 of the authors' follow-up manuscript [*Word Error Rate Definitions and Algorithms for Long-Form Multi-talker Speech Recognition* (2025, v1)](https://arxiv.org/html/2508.02112v1) (supplemental verification 2026-10-02) reinforces three rules: (1) reference utterance boundaries are fixed evaluation inputs ([Exercise E11-23](../../../../chapters_en/11_selection-guide.md#e11-23)); (2) DI-cpWER is a diagnostic tool rather than a primary leaderboard metric because re-segmenting hypotheses alters assignment degrees of freedom; and (3) standard Levenshtein ignores timestamps completely—matching reference word `"spring"` on `[0, 1]` s against hypothesis `"spring"` at `[10, 11]` s with `0` errors, whereas **time-constrained Levenshtein** rejects disjoint intervals (`1` deletion + `1` insertion = `2` errors) while matching overlapping `[0, 1]` s and `[0.5, 1.5]` s with `0` errors ([Exercise E11-24](../../../../chapters_en/11_selection-guide.md#e11-24), operating on pre-tokenized whole words).

Pay close attention to boundary-touching intervals in pinned [`levenshtein.h::overlaps` (line 94)](https://github.com/fgnt/meeteval/blob/6e3dc81284f2d6928f7ef9e620fd3b6906daa429/meeteval/wer/matching/levenshtein.h#L94): upstream tests two **strict** inequalities (`a_begin < b_end && b_begin < a_end`). Consequently, abutting positive-length intervals `[0, 1]` and `[1, 2]` do **not** overlap (`2` errors), whereas a zero-duration point hypothesis `[1.5, 1.5]` strictly inside `[0, 3]` **does** satisfy `0 < 1.5 && 1.5 < 3` (`0` errors), while point hypotheses at endpoints `0.0` or `3.0` fail (`2` errors)! Indeed, pinned [`time_constrained.py` lines 75–102](https://github.com/fgnt/meeteval/blob/6e3dc81284f2d6928f7ef9e620fd3b6906daa429/meeteval/wer/wer/time_constrained.py#L75-L102) subdivide unaligned utterance spans proportionally to character length for references (`"abc b"` over `[0, 4]` s $\to$ `[0, 3]` and `[3, 4]`) and Collapse hypothesis words to **interval midpoints** (`[1.5, 1.5]` and `[3.5, 3.5]`) so hypotheses cannot game overlap by inflating word spans ([2025 preprint v1 §V-C](https://arxiv.org/html/2508.02112v1); note also that the `5` s collar used in [Neumann et al. CHiME-7 2023 §4.2, p. 30](https://www.isca-archive.org/chime_2023/neumann23_chime.pdf) reflects corpus annotation jitter rather than a universal law).

#### Missing Sessions Must Never Vanish Silently from Corpus Denominators

Invoking pinned [`apply_multi_file` (`seglst.py` line 565)](https://github.com/fgnt/meeteval/blob/6e3dc81284f2d6928f7ef9e620fd3b6906daa429/meeteval/io/seglst.py#L565) directly with a lightweight segment-counting callback confirms upstream's missing-session guardrails:

| Session-Set Boundary Condition | Pinned `apply_multi_file` Behavior |
|---|---|
| 10 reference sessions, 1 hypothesis session missing (`10%`) | Evaluates all 10 sessions, substituting an empty `SegLST` for the missing hypothesis with a warning |
| 2 reference sessions, 1 hypothesis session missing (`50% > 0.1`) | Raises `ValueError` exceeding default `missing_hypothesis` ratio threshold `0.1` |
| Hypothesis contains extra session ID absent from reference | Raises exception by default |
| Same extra hypothesis session with explicit `partial=True` | Evaluates intersection of reference sessions and logs a warning listing dropped IDs |
| Empty reference session set | Raises exception |

Hence silent system outputs on hard sessions must emit explicit empty hypothesis entries rather than dropping files or passing `partial=True` to shrink the denominator.

#### Control-Flow Quirks in the Pinned CHiME Scoring Wrapper

Static and mocked-I/O inspection of [`chime_utils/scoring/meeteval.py`](https://github.com/chimechallenge/chime-utils/blob/152882404f572d40769ef02bf91c5a9a9cfc9c78/chime_utils/scoring/meeteval.py) highlights three wrapper specifics: (1) `_wer` computes scenario-level WER first (`chime6`, `dipco`, `mixer6`) and **macro-averages** across scenarios rather than pooling word counts micro-wise; (2) when `ignore_missing=False`, `_load_and_prepare` logs an error without returning—raising `UnboundLocalError` if scenario 0 (`chime6`) is missing a file, or **silently reusing `chime6`'s reference/hypothesis dict** if scenario 1 (`mixer6`) is missing a file (!); and (3) its iterative text-normalization fixed-point loop checks `if words == words:` (comparing the post-normalization string to *itself* rather than to the previous iteration!), exiting after a single pass (confirmed by passing `"aaa"` to a stateful mock normalizer that strips one leading character per call, returning `"aa"` instead of `""`). Always pre-validate file completeness and normalizer idempotence before running challenge wrappers.

<a id="meeting-scoring-audit"></a>

#### Execution Tiers and Re-runnable Entry Points for the Meeting-Scoring Audit

The [meeting-scoring interface audit](../../ch11/examples/audit_meeting_scoring_interfaces.py) and its [current report](../../ch11/reports/meeting_scoring_interfaces_current.json) separate direct `apply_multi_file` calls, AST-extracted doc functions, pseudo-timestamp helpers, mocked-I/O `chime-utils` control-flow probes, and the expected `ModuleNotFoundError` (`cy_levenshtein` Cython extension uncompiled) when calling production `cp_word_error_rate`—preserving the [2026-09-28 historical report](../../ch11/reports/meeting_scoring_interfaces.json) untouched and verifying expected values via [independent offline DP/permutation tests](../../../../tests/test_codes_meeting_scoring_interfaces.py) inside a clean subprocess (logging 32 imported upstream `.py` modules alongside `source_selection_mismatch` status without installing Cython extensions or running ASR).

<a id="meeting-kernel-contracts"></a>

#### Native Execution Contract for the Original C++ Word-Edit Alignment Kernels

Complementing the Python audit above, our [native C++ meeting-kernel audit](../../ch11/examples/audit_meeting_kernel_contracts.py) `#include`s unmodified pinned `levenshtein.h` and natively executes `levenshtein_distance_` and `time_constrained_levenshtein_distance_v2_` across 13 deterministic cases (the original 9 cases plus 4 zero-duration point-timestamp controls, unit edit costs `1/1/1`, verified against [independent 2D DP unit tests](../../../../tests/test_codes_meeting_kernel_contracts.py)):

| Pre-Tokenized Input & Interval Relationship | Standard Levenshtein Errors | Time-Constrained Kernel Errors |
|---|---:|---:|
| Identical word, disjoint intervals | 0 | 2 |
| Identical word, abutting endpoints (`[0, 1]` vs. `[1, 2]`) | 0 | 2 |
| Identical word, strictly positive overlap | 0 | 0 |
| Distinct words, strictly positive overlap | 1 | 1 |
| Distinct words, disjoint intervals | 1 | 2 |
| Empty reference, 2-word hypothesis | 2 | 2 |
| 2-word reference, empty hypothesis | 2 | 2 |
| Both empty | 0 | 0 |
| Identical 2 words, both intervals disjoint | 0 | 4 |
| Identical word, ref `[0, 3]`, interior point hyp `[1.5, 1.5]` | 0 | 0 |
| Identical word, ref `[0, 3]`, left-endpoint point hyp `[0, 0]` | 0 | 2 |
| Identical word, ref `[0, 3]`, right-endpoint point hyp `[3, 3]` | 0 | 2 |
| Identical 2 words, ref `[0, 3], [3, 4]`, midpoint hyps `1.5, 3.5` | 0 | 0 |

Archived in the [current kernel report](../../ch11/reports/meeting_kernel_contracts_current.json) while preserving the [9-case historical report](../../ch11/reports/meeting_kernel_contracts.json), the runner sanitizes compiler include-path env vars, records real `.d` compiler dependencies across the 7 verified upstream files, compiles in `/tmp`, and writes inside the repo solely under explicit `--report`:

```bash
.venv/bin/python -B -m codes.chapters.ch11.examples.audit_meeting_kernel_contracts
.venv/bin/python -B -m codes.chapters.ch11.examples.audit_meeting_kernel_contracts --report codes/chapters/ch11/reports/meeting_kernel_contracts_current.json
```

## 6. Turning Upstream Repositories into Executable Pipeline Experiments

Always isolate **one concrete end-to-end chain** at a time—for instance, on a Linux conferencing endpoint: `ALSA/PortAudio capture -> PipeWire or in-app AEC -> fixed TDFB -> single-channel NS -> VAD/AGC -> ASR`, deciding explicitly up front whether PipeWire or the application owns AEC and reference capture.

Every experimental log must archive at minimum: upstream commit SHA, checkpoint digest, dependency versions and compiler flags, hardware and firmware revisions, room acoustics setup, channel map and clock domain, dataset split, callback block size, state initialization/reset rules, evaluation metrics, per-file scores, and CPU/RAM/power telemetry. Stateful regression suites must cover full-utterance, chunked streaming, pause/resume, and fault injection—unit tests on stateless math functions never substitute for them.

Mismatched block lengths require FIFO adaptation, yet never require waiting for the least common multiple (LCM) of all stage block sizes! At `16` kHz, for instance, a `10` ms capture callback is `160` samples, a `15` ms processing hop is `240` samples, and a `32` ms VAD window is `512` samples; their LCM is `7680` samples (`480` ms!). Waiting `480` ms to fill an LCM super-frame before launching downstream modules injects nearly half a second of completely avoidable latency.

Instead, streaming ring queues invoke each consumer immediately once its own required input hop has accumulated, retaining leftover samples in the queue and associating VAD control decisions via explicit sample timestamps. Acceptance tests must verify zero duplicate or dropped samples per consumer and exact sample-interval attribution for every control decision (note that this FIFO scheduling example does not imply chaining incompatible vendor firmware binaries together).

In this textbook, [Engineering Exercises E10-07 and E10-08](../../../../chapters_en/10_engineering-practice.md#sec-10-11) work out exact byte budgets and consumer sample intervals: running `.venv/bin/python -m codes.chapters.ch00.cross_chapter.exercises_engineering`, `events` in E10-08 logs half-open `[start, end)` sample slices and callback ticks while `leftovers` tracks buffered remainders—showing the two consumers emitting their first valid blocks at `20` ms and `40` ms respectively (assuming instant zero-compute delivery rather than claiming a physical soundcard latency measurement).

In that same runner, E10-14 links six-frame serial execution, in-flight buffer capacity, concurrent resident RAM, and periodic duty-cycle power estimation ([full hand derivation in Chapter 10](../../../../chapters_en/10_engineering-practice.md#sec-10-4-3)), while Chapter 11 Exercises E11-05 through E11-07 screen candidate designs against circular-array inter-mic aliasing spacing, single-channel conferencing output/latency hard constraints, and cross-device clock drift—all verifying deterministic arithmetic constraints on textbook specifications rather than certifying hardware acoustics.

Likewise, [state and timestamp exercises (`tracking_time_exercises.py`)](../cross_chapter/tracking_time_exercises.py) isolate two classic edge cases: E10-15 traces frame-by-frame VAD pre-roll, onset trigger, and hangover hold to emit the exact sample slice `[160, 1120)` without duplicating the trigger frame, while E11-08 derives the exact one-sided binomial risk upper bound from zero observed failures—separating `"zero failures observed in N trials"` from `"true failure probability proven negligible"` alongside E09-08's tracker state timestamps (without treating correlated consecutive audio frames as i.i.d. Bernoulli trials!).

Complementing those checks, [engineering boundary exercises (`engineering_boundary_exercises.py`)](../cross_chapter/engineering_boundary_exercises.py) verify exact boundary counterexamples: E10-16 confirms zero dropped frames when service time equals period $T$ exactly vs. deadline misses when slightly exceeding $T$; E10-17 computes corresponding-sample loopback latency; and E11-09 preserves six paired sessions including one tie, enumerating all $2^5 = 32$ sign permutations across the five non-tied sessions to compute the exact one-sided sign-test $p$-value—demonstrating why pooled micro-WER magnitude and per-session win/loss consistency answer distinct statistical questions.

Real-time callbacks must always distinguish [`PaStreamCallbackTimeInfo`](https://files.portaudio.com/docs/v19-doxydocs/structPaStreamCallbackTimeInfo.html)'s `inputBufferAdcTime`, `currentTime`, and `outputBufferDacTime` (in seconds on the stream's clock, verified 2026-09-22), never subtracting uncalibrated clocks across distinct devices, and re-running quantization, timestamp, and FIFO consumption checks whenever models or runtimes are upgraded.

Everything cataloged above separates static source study and experimental design from actual execution: consult `codes/chapters/ch00/upstream/` status files for local checkout state, and rely strictly on logged reports for build, runtime, and hardware results.

## 7. Independent Comparisons of Resampling, Loudness, Scoring, File I/O, and Block Floating-Point Interfaces

The six entries below cover sample-rate conversion, loudness metering, intrusive perceptual scoring, audio file I/O, and block floating-point arithmetic. They neither form a mandatory serial chain nor belong exclusively to multi-microphone array processing; entries I22, I23, and I26 include executed host C interface measurements bounded strictly to their logged stimuli, versions, and metrics without claiming hardware or perceptual speech-quality certification.

### I22: `libsoxr` Ratio Direction, Slewing, and Output Delay

`libsoxr` is the standalone SoX resampling library. This book audited [official SourceForge commit `945b592b70470e29f917f4de89b4281fbbd540c0` (v0.1.3)](https://sourceforge.net/p/soxr/code/ci/945b592b70470e29f917f4de89b4281fbbd540c0/tree/), identical to the `chirlu/soxr` GitHub mirror. Read `src/soxr.h` first, then `examples/5-variable-rate.c`, and finally filter implementations under `src/`. Upstream [`LICENCE`](https://sourceforge.net/p/soxr/code/ci/945b592b70470e29f917f4de89b4281fbbd540c0/tree/LICENCE) specifies `LGPL-2.1-or-later` while pointing additionally to per-file notices in `pffft.c`—permitting local source study while imposing LGPL/notice obligations on binary redistribution or static linking.

When comparing `libsoxr` against `libsamplerate` (I04), check **ratio direction** first! Whereas `libsamplerate` defines `src_ratio = out_rate / in_rate`, `libsoxr`'s `soxr_set_io_ratio` defines `io_ratio = in_rate / out_rate`—the exact reciprocal! Compensating a `+100` ppm fast capture clock onto nominal rate uses `≈ 0.99990001` in `libsamplerate` versus `1.0001` in `libsoxr` (whereas `soxr_create` takes separate `input_rate` and `output_rate` arguments directly).

Variable-rate operation requires enabling `SOXR_VR` at creation time: per [lines 37–66 of `examples/5-variable-rate.c`](https://github.com/chirlu/soxr/blob/945b592b70470e29f917f4de89b4281fbbd540c0/examples/5-variable-rate.c), initialize `soxr_create` with the maximum anticipated `io_ratio`, then call `soxr_set_io_ratio(soxr, io_ratio, slew_len)`, where `slew_len` specifies the number of **output samples** over which to ramp smoothly to the new ratio (`0` = immediate step; `160` output samples at `16` kHz = `10` ms transition).

During streaming calls ([`src/soxr.h`](https://github.com/chirlu/soxr/blob/945b592b70470e29f917f4de89b4281fbbd540c0/src/soxr.h)), accumulate consumed input frames `idone` and emitted output frames `odone`, carry unconsumed input forward, and signal end-of-stream strictly by passing `in = NULL` (passing `ilen = 0` with a non-`NULL` pointer does NOT trigger EOF draining!) while looping until `odone == 0`. Remember that `soxr_delay()` returns pending delay denominated in **output samples** rather than milliseconds (excluding driver/queue latency), and that `soxr_create` defaults to `SOXR_HQ` whereas `soxr_oneshot` defaults to `SOXR_LQ`—so whole-utterance vs. chunked comparisons must equalize quality recipes, phase responses, data types, threading, and integer dither explicitly.

On 2026-09-22, this book built and executed pinned `libsoxr`, `libebur128`, and `libsndfile` (I22, I23, I26) on macOS arm64 (AppleClang 21.0.0, CMake 4.3.2) as isolated `Release` static libraries without installing system packages—disabling OpenMP and SIMD in `soxr`, disabling external compressed codecs in `libsndfile`, and passing `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` on the command line without touching upstream files (verified clean before and after build). Driven by [`run_industrial_interfaces.py`](../../ch10/examples/run_industrial_interfaces.py) and [`industrial_interfaces.c`](../../ch10/examples/industrial_interfaces.c), all SHA-256 digests and raw outputs live in [`industrial_interfaces.json`](../../ch10/reports/industrial_interfaces.json).

On a 1 s, `48` kHz stereo float64 stimulus (left channel: `0.5`-amplitude impulse at 0-based frame `12000`; right channel: `0.1`-amplitude, `1` kHz sine), downsampling to `16` kHz under explicit `SOXR_HQ`, `flags=0`, single-threaded interleaved `float64`, and `113`-frame output buffers across whole-utterance, `127`-frame chunks, and `509`-frame chunks yields:

| Input Feeding Mode | Output Frames Before EOF Drain | Frames Emitted During `in=NULL` Drain | Total Output Frames | Max Absolute Diff vs. Whole-Buffer |
|---|---:|---:|---:|---:|
| Whole buffer | 15471 | 529 | 16000 | 0 |
| 127-frame chunks | 15659 | 341 | 16000 | 0 |
| 509-frame chunks | 15622 | 378 | 16000 | 0 |

Post-drain `soxr_delay()` reads `0` across all three runs (pre-drain differences reflect `113`-frame output buffer staging rather than distinct algorithmic latencies). The left-channel impulse peak lands at exact output frame `4000` ($12000 \times 16000 / 48000 = 4000$), right-channel interior `[320, 15680)` matches analytical `1` kHz sine within $6.84 \times 10^{-8}$ without gain/phase fitting, and whole-buffer feeding exercises `141` partial-consumption iterations. Conversely, calling `soxr_clear()` mid-stream at input frame `24000` without draining loses `434` output frames (`15566` total) and produces an unaligned index-wise discrepancy of `0.07669`—flagging protocol misuse rather than perceptual distortion.

#### Executed: Variable-Ratio Resampling and Slewing Under Two-Stage Known Clock Drift

A second native probe on pinned `libsoxr` tests `SOXR_HQ | SOXR_VR` on known synthetic clock drift: nominal `48` kHz input runs `+100` ppm fast ($48004.8$ Hz, `240024` frames) over `0–5` s and `+150` ppm fast ($48007.2$ Hz, `240036` frames) over `5–10` s (`480060` input frames total), targeting `16` kHz stereo `float64` (left: `0.5`-amplitude impulses at exact true seconds `1..9` s; right: `0.1`-amplitude, `3` kHz sine evaluated at exact per-sample physical instants).

Initializing `soxr_create` with max `io_ratio = 48007.2 / 16000 = 3.00045` and initial `io_ratio = 48004.8 / 16000 = 3.0003`, output buffers are capped right at cumulative output frame `80000` (`5.0` s) so the second `soxr_set_io_ratio` fires at the exact `80000`-frame boundary across four modes—holding stale `3.0003`, stepping immediately (`slew_len=0`), slewing over `160` output frames (`10` ms), and slewing over `320` output frames (`20` ms):

| Ratio-Update Mode at Output Frame 80000 | Total Drained Output Frames | `1` s Impulse Index | `5` s Impulse Index | `9` s Impulse Index |
|---|---:|---:|---:|---:|
| Hold initial `3.0003` | 160004 | 16002 | 80002 | 144005 |
| Immediate step (`slew_len = 0`) | 160000 | 16002 | 80002 | 144002 |
| `160`-frame slew (`10` ms) | 160000 | 16002 | 80002 | 144002 |
| `320`-frame slew (`20` ms) | 160000 | 16002 | 80002 | 144002 |

All four modes share an initial `2`-frame filter/peak offset (`16002` at `1` s); comparing *relative drift* from `1` s to `9` s, holding stale `3.0003` drifts by `+3` frames (`144005 - 16002 - 128000 = 3`, consistent with `50` ppm over `4–5` s quantized onto integer impulse indices and `160004` total output frames), whereas all three updated modes eliminate that drift (`160000` total frames, `9` s peak at `144002`). Whole-buffer and `127`-frame chunking under `160`-frame slewing agree to `0` max difference; compared index-by-index against `slew_len = 0`, `160`-frame and `320`-frame slews differ by at most `0.000624623` and `0.001245994` respectively—consistent with analytical linear-slew cumulative input-phase offsets $(3.00045 - 3.0003) \times 160 / 2 = 0.012$ and $0.024$ input frames. Notably, after the final `in = NULL` call returns `odone == 0` under `SOXR_VR`, `soxr_delay()` still reads `100` output samples (unlike `0` in fixed-rate mode)—logged verbatim in [`industrial_interfaces.json`](../../ch10/reports/industrial_interfaces.json) without misinterpreting `100` as unwritten audio or hardware delay.

### I23: `libebur128` Integrated Loudness and True Peak

Peak-limiting AGCs bound sample amplitude, yet do not regulate audio volume according to human loudness perception. `jiixyj/libebur128` provides loudness and peak metering via [`ebur128/ebur128.h`](https://github.com/jiixyj/libebur128/blob/67b33abe1558160ed76ada1322329b0e9e058b02/ebur128/ebur128.h) and `ebur128/ebur128.c`. The version macro in the pinned snapshot is `1.2.6`, and root [`COPYING`](https://github.com/jiixyj/libebur128/blob/67b33abe1558160ed76ada1322329b0e9e058b02/COPYING) specifies the MIT License; bind implementations and metering modes strictly to this commit rather than claiming full standards certification across every configuration merely because the library name mentions R128.

Momentary loudness, short-term loudness, and integrated loudness operate over distinct statistical time windows; Loudness Range (`LRA`) characterizes temporal variability, while sample peak and true peak inspect waveform amplitude. In the C API, integrated loudness is returned in `LUFS` (*Loudness Units relative to Full Scale*), whereas `ebur128_true_peak` returns linear full-scale amplitude—converting to `dBTP` (*decibels True Peak*) requires taking $20\log_{10}$ relative to digital full scale. Neither quantity substitutes for the other, and neither yields acoustic `dB SPL` without physical microphone sensitivity calibration.

The underlying measurement algorithm follows [ITU-R BS.1770-5](https://www.itu.int/rec/R-REC-BS.1770-5-202311-I) (approved 2023-11-22, verified in force on the official ITU page on 2026-09-28). A metering standard does not automatically dictate target loudness for every product: items `g–m` on pages 3–4 of [EBU R 128-2023](https://tech.ebu.ch/docs/r/r128.pdf) address broadcast programme delivery—specifying whole-programme integration, a `-23.0` LUFS target, and a `-1.0` dBTP ceiling in linear PCM production alongside delivery tolerances. Never copy broadcast mastering targets blindly as default AGC setpoints for ASR, VoIP telephony, or wake-word front ends (the `≈ -23` LUFS sine probe below simply reflects the amplitude of `0.1` chosen for that test stimulus).

Initialize the desired metering modes via `ebur128_init`, then stream interleaved audio continuously. Frame counts passed to `ebur128_add_frames_*` denote time frames *per channel*, not total scalar samples across channels. Channel role assignments directly weight loudness summation: because the header's default channel map assumes standard surround/programme speaker layouts, never feed four raw microphones of an array directly as four surround programme channels! Measure a single enhanced output stream as mono `CENTER`, and maintain separate single-channel meter instances when comparing individual microphones.

Two-stage gating determines which time blocks enter integrated loudness, so digital silence or sub-threshold signals can return `-inf`. Callers must record meter state, active duration, and undefined/failure reasons explicitly—never clamping `-inf` to `0` LUFS and averaging it into corpus means! Likewise, true-peak oversampling estimates reconstructed inter-sample waveform crests that can exceed the maximum stored discrete sample magnitude, so unclipped sample peaks alone do not guarantee D/A headroom.

Recommended protocol check: feed a sufficiently long, steady-amplitude non-silent signal well above the gating threshold, then scale the entire waveform by `0.5`. Provided the gated active region remains identical, integrated loudness drops analytically by $20\log_{10} 2 \approx 6.0206$ LU while both sample peak and linear true peak halve. Test digital silence, short clips, per-chunk state resets, and misassigned channel maps separately (verifying metering contracts rather than subjective listening quality).

Core `libebur128` builds cleanly with C and CMake without neural weights; root `CMakeLists.txt` points optional unit tests toward `test/`, whose file-reading dependencies and standard EBU test vectors are audited separately from compiling the static library itself.

Our executed synthetic check feeds a `48` kHz, `10` s, `1` kHz mono sine explicitly configured as `CENTER` under `EBUR128_MODE_I | EBUR128_MODE_TRUE_PEAK`. Amplitudes `0.1` and `0.05` yield integrated loudness `-23.00360` LUFS and `-29.02420` LUFS—differing by `-6.020599913` LU in exact agreement with $20\log_{10}(0.5)$ because both sit far above the `-70` LUFS absolute gate (note that halving real speech containing near-threshold passages can alter absolute-gate retention and shift integrated loudness by other than `6.02` LU).

Streaming `127`-frame chunks while preserving state matches whole-buffer integrated loudness to `0` LU difference, with sample peaks `0.1` and `0.05` and estimated true peaks `≈ 0.1000000015` and `0.0500000007` (linear ratio `0.5`). Digital silence returns `-inf`, serialized as `null` alongside reason `negative_infinity_no_gated_energy` and `0` peaks in [`industrial_interfaces.json`](../../ch10/reports/industrial_interfaces.json).

#### Executed: Constructed Signal with Inter-Sample True Peaks

Because `48` kHz samples a `1` kHz sine very densely near its crest, a second stimulus constructs a clear inter-sample peak gap: a `48` kHz, `10` s, **`12` kHz** mono sine with initial phase $\pi/4$ and amplitude `0.95`, tapered by `100` ms raised-cosine fade-in/fade-out windows at both ends to prevent step-onset overshoot in the interpolation FIR. Across the steady interior, discrete sample phases advance by $\pi/2$ per step, giving every sample exact magnitude $0.95/\sqrt{2} = 0.6717514421$ (`-3.45583` dBFS) while the continuous sinusoid reaches `0.95` (`-0.44553` dBFS) midway between samples—an analytical `3.01030` dB gap unaffected in the interior by boundary fades.

In pinned `libebur128` `1.2.6`, the [`ebur128.c` polyphase interpolator](https://github.com/jiixyj/libebur128/blob/67b33abe1558160ed76ada1322329b0e9e058b02/ebur128/ebur128.c#L119-L174) applies a 49-tap, 4× oversampling FIR at `48` kHz—a finite-order estimator that does not claim exact sinc reconstruction on every bandlimited waveform. Native C execution returns sample peak `0.671751442206` (`-3.45583` dBFS) and true-peak estimate `0.946371197701` (`-0.47877` dBTP), exceeding sample peak by `2.97706` dB and matching bit-for-bit between whole-buffer and `127`-frame streaming (differing by `0.00363` from continuous `0.95`, so the two figures are kept distinct).

Omitting the `100` ms fade-in produces a discontinuous step onset: independent convolution with those exact 49 taps shows the initial transient reaching `≈ 0.96139 > 0.95`, explaining why smooth tapers are essential to isolate steady-state inter-sample peaks; see [`industrial_interfaces.json`](../../ch10/reports/industrial_interfaces.json) and [`industrial_interfaces.c`](../../ch10/examples/industrial_interfaces.c).

### I24: `pystoi` STOI, ESTOI, and Invalid-Score Handling

Chapter 10's Short-Time Objective Intelligibility (`STOI`) metric can be inspected through `mpariente/pystoi`, which implements both standard `STOI` and Extended STOI (`ESTOI`) in [`pystoi/stoi.py`](https://github.com/mpariente/pystoi/blob/74872b000753a7a42ff51aa0868af8c82c7f9053/pystoi/stoi.py). Maintained as an independent Python package whose [`README.md`](https://github.com/mpariente/pystoi/blob/74872b000753a7a42ff51aa0868af8c82c7f9053/README.md) credits Cees Taal's MATLAB reference for equivalence tests, the two implementations and runtimes must not be conflated.

Marked `0.4.1` in `setup.py` under root `MIT` `LICENSE` and depending on NumPy and SciPy, its Octave/MATLAB cross-check suite requires separate external assets not executed automatically upon installing Python `pystoi`.

Inputs are 1D equal-length clean reference and processed speech waveforms sharing an identical sample rate and time alignment. Internally, `pystoi` resamples to `10` kHz, strips silent frames relative to the clean reference envelope, and computes short-time third-octave correlations; because `extended=False` and `extended=True` apply distinct normalizations, always report `STOI` and `ESTOI` in separate columns.

Always verify whether a returned float reflects a valid evaluation: whenever fewer than `30` STFT frames remain after silence removal, pinned `pystoi` issues a `RuntimeWarning` and returns `1e-5`—an implementation failure sentinel rather than a genuine `"near-zero intelligibility"` measurement! Scoring harnesses must trap warnings, log clip IDs/lengths/reasons, and retain failed clips in the corpus denominator rather than averaging `1e-5` into the mean.

#### Executed on Original Python API: Exceptions, Sentinels, and Finite Zero-Input Scores

On 2026-09-28 inside a pre-existing isolated environment (Python `3.13.12`, NumPy `2.3.5`, SciPy `1.17.1`, without installing SciPy into shared `.venv`), 12 calls to unmodified `pystoi.stoi` on identical `10` kHz `float64` 1D pairs (bypassing resampling) revealed three distinct zero-input regimes:

| All-Zero Input Length | Standard `STOI` | `ESTOI` |
|---|---|---|
| `256` samples | `AxisError`, no score returned | Identical `AxisError` |
| `1024` samples | `RuntimeWarning` + `1e-5` sentinel | Identical `RuntimeWarning` + `1e-5` sentinel |
| `10000` samples | Finite `0.0`, **no warning** | RNG-seed-dependent finite float, **no warning** |

At `256` samples, `utils.remove_silent_frames` raises `AxisError` right after frame formation before reaching the `30`-frame check. Conversely, a `10000`-sample all-zero reference has zero maximum frame energy, so relative silence gating removes zero frames! While `STOI` returns `0.0` silently, `ESTOI`'s `utils.row_col_normalize` adds machine-epsilon `np.random.randn` jitter independently to both arguments—yielding `≈ 0.00856519` after `np.random.seed(0)` and `≈ −0.000209994` after `np.random.seed(1)` on all-zero self-comparison! Flagged explicitly as `finite_zero_energy_invalid` in our report, these numbers carry zero intelligibility meaning.

On non-degenerate random input `np.random.default_rng(4).normal(size=10000)`, self-comparison returns exact `1.0` in both `STOI` and `ESTOI`, matching unit-vector correlation identity; see [I32 script and report](#industrial-upstream-audit) for per-call digests, seeds, exceptions, and warnings.

For full speech evaluation, prepare licensed clean natural speech paired with self-identity, known-SNR noisy, all-zero muted, and time-misaligned versions, logging mode, SciPy version, alignment rule, and short-clip rejection—remembering that `STOI` is neither word recognition accuracy nor meaningful on synthetic pure tones.

### I25: ViSQOL Operating Modes, Model Assets, and Bazel Dependencies

Google's **ViSQOL** (*Virtual Speech Quality Objective Listener*) is a full-reference perceptual quality estimator that compares clean and degraded spectro-temporal neurograms to predict `MOS-LQO` (*Mean Opinion Score—Listening Quality Objective*). Requiring a clean reference (unlike non-intrusive `DNSMOS` in I19), it neither implements ITU-T POLQA nor replaces human listening panels; see [`README.md` Guidelines, License, and FAQ](https://github.com/google/visqol/blob/38d0b0163e441047d4429bf07ad09e5b9031d02c/README.md).

Upstream `audio` mode operates at `48` kHz while `speech` mode operates at `16` kHz, downmixing multi-channel inputs to mono prior to comparison—so ViSQOL scores never certify array directivity, binaural ILD/ITD cues, or spatial image preservation! In `speech` mode, record VAD gating, mapping model file, and `scaled` vs. `unscaled` MOS mapping explicitly without mixing modes.

Crucially, similarity-to-MOS mapping requires model weights: inspecting `src/`, `src/proto/`, and Python bindings suffices to study API plumbing, yet `model/` contains both `.tflite` neural models and `.txt` LIBSVM support-vector files that are genuine numerical model weights (never mistake `.txt` SVR files for prose documentation!). When sparse checkouts omit `model/` and `testdata/`, document that exclusion explicitly rather than crashing at runtime and blaming the algorithm.

Inspecting root `WORKSPACE` and `BUILD` at the pinned commit ([`WORKSPACE`](https://github.com/google/visqol/blob/38d0b0163e441047d4429bf07ad09e5b9031d02c/WORKSPACE)) shows Bazel pulling TensorFlow `2.11.0`, `protobuf`, `pybind11`, `Abseil`, `LIBSVM`, `Armadillo`, and `PFFFT`—record source size, external fetch footprint, Bazel cache size, and runtime RAM separately under root `Apache-2.0` code licensing (which does not relicense external dependencies or speech corpora).

Upstream's FAQ notes that ViSQOL was designed around codec and VoIP network artifacts and behaves unevenly on neural noise suppression or generative speech artifacts, requiring corroboration from WER and human listening tests.

Furthermore, the [ViSQOL v3 QoMEX 2020 paper](https://arxiv.org/pdf/2004.09584) establishes sharper technical boundaries: §II.C recommends scoring `3–10` s active speech/audio segments; §II.B notes limited sensitivity to certain subtle music-codec bitrates and cautions against using scores as blind CI regression gates; §III.C performs fine-grained patch-wise temporal realignment internally (meaning ViSQOL scores cannot detect constant latency shifts or subtle chunk-boundary slips!); and §III.B defines conformance-version tracking bound to model and code hashes (static reading only; ViSQOL was not built or run here).

### I26: `libsndfile` Frame Counts, Amplitude Normalization, and Short Reads

Arrays entering DSP routines depend directly on audio file decoding contracts. We study `libsndfile/libsndfile` (`include/sndfile.h`, `src/sndfile.c`, `src/pcm.c`) strictly as an I/O protocol baseline rather than an array algorithm: the [official API](https://libsndfile.github.com/libsndfile/api.html) separates scalar item-counted (`sf_read_*`), time-frame-counted (`sf_readf_*`), and raw byte interfaces.

Reading `160` stereo frames via `sf_readf_float(..., 160)` requires a `320`-float buffer, corresponding to `sf_read_float(..., 320)`! At EOF, returned counts may fall strictly below the requested count while trailing buffer slots remain untouched or zeroed—slice strictly to the returned count so padded tail slots never corrupt running covariance statistics.

Keep container bit depth separate from host array `dtype`: reading a `PCM16` WAV into `float32`/`float64` leaves its underlying `16`-bit quantization unchanged, so normalization, integer writeback rounding/clipping, channel interleaving, and endianness must all be logged explicitly. Pinned `CMakeLists.txt` specifies `1.2.2` under `LGPL-2.1-or-later` ([`COPYING`](https://github.com/libsndfile/libsndfile/blob/b9103bd48b6c8fb517ae737fe3baee0c718b804c/COPYING) and [`src/sndfile.c`](https://github.com/libsndfile/libsndfile/blob/b9103bd48b6c8fb517ae737fe3baee0c718b804c/src/sndfile.c)), allowing standalone PCM/WAV builds without external Ogg/Vorbis/FLAC/Opus/MPEG libraries.

Our executed native test writes a `16` kHz, 7-frame, 2-channel signed PCM16 WAV (left channel first per row) testing exact 16-bit full-scale endpoints, signs, channel ordering, and short reads:

```text
-32768      0
     0  32767
 16384 -16384
     1     -1
 12345 -23456
     0   1000
 32767 -32768
```

Written via `sf_writef_short` and re-opened: requesting 3 frames per `sf_readf_short` call returns `3, 3, 1, 0` (`7` frames cumulative); requesting 6 scalar items per `sf_read_short` call returns `6, 6, 2, 0` (`14` items cumulative), matching all `14` code words exactly. Normalized `sf_readf_double` divides by `32768` with `0` error, mapping `+32767` to `0.999969482421875` (`< 1.0`).

Independent readback verification avoids calling `libsndfile` again, instead parsing raw RIFF/WAV headers and little-endian `int16` code words via Python's standard-library `wave` module alongside [failure-mode unit tests](../../../../tests/test_codes_industrial_interfaces.py) and [`industrial_interfaces.json`](../../ch10/reports/industrial_interfaces.json).

### I27: `lib_xcore_math` Block Floating-Point and Reference Kernels

Pinned `lib_voice` requires exact dependency version `lib_xcore_math(v3.0.0)` in [`lib_voice/lib_build_info.cmake`](https://github.com/xmos/lib_voice/blob/c9f1a9bf95cd88c7950adf4bf631c217f900ad25/lib_voice/lib_build_info.cmake), pinned here at `16130be45c4002a1f875a4b06ff68d2065cd8c69` to complement `Q15` FIR (I15) with vector, block floating-point (`BFP`), FFT, DCT, and filter primitives.

Block floating-point stores an array of integer mantissas alongside a **single shared exponent**: for instance, mantissas `[-16384, 8192]` paired with exponent `-15` represent `[-0.5, 0.25]`, so identical mantissa bytes under different exponents represent different physical amplitudes! Meanwhile, bit **headroom** counts available left-shift bits before signed overflow, trading overflow safety against right-shift quantization loss across dynamic block rescaling.

Read from `lib_xcore_math/api/` into `src/bfp/`, `src/fft/`, and `src/filter/`, comparing portable C reference paths (`src/arch/ref/`) against XS3 assembly paths in [`lib_build_info.cmake`](https://github.com/xmos/lib_xcore_math/blob/16130be45c4002a1f875a4b06ff68d2065cd8c69/lib_xcore_math/lib_build_info.cmake) (including pre-tabulated vs. generated FFT twiddle tables). Reconstruct real values $m \cdot 2^e$ on 2-point vectors, near-full-scale pulses, non-SIMD remainder lengths, mixed high/low dynamic range vectors, and FFT round-trips before profiling cycle counts and SRAM on physical xCORE silicon—noting [`README.rst`](https://github.com/xmos/lib_xcore_math/blob/16130be45c4002a1f875a4b06ff68d2065cd8c69/README.rst) (`XTC Tools 15.3.1` and host VPU simulation caveats) and [`LICENSE.rst`](https://github.com/xmos/lib_xcore_math/blob/16130be45c4002a1f875a4b06ff68d2065cd8c69/LICENSE.rst) (`XMOS Public Licence v1`).


### I28: FastEnhancer Explicit Streaming Caches and Single-Channel Denoising

FastEnhancer addresses single-channel neural noise suppression ([pinned `README.md`](https://github.com/aask1357/fastenhancer/blob/f85223bd546b27f39dc0744e0310dcd246f750a4/README.md); commit `f85223bd546b27f39dc0744e0310dcd246f750a4` fetched under [`MIT` `LICENSE`](https://github.com/aask1357/fastenhancer/blob/f85223bd546b27f39dc0744e0310dcd246f750a4/LICENSE) into `codes/chapters/ch00/upstream/_downloads/fastenhancer/`, audited 2026-09-26 without copying leaderboard numbers).

**What the model computes.** In [`models/fastenhancer/default/model.py`](https://github.com/aask1357/fastenhancer/blob/f85223bd546b27f39dc0744e0310dcd246f750a4/models/fastenhancer/default/model.py), `RNNFormerBlock` pairs temporal RNN state recurrence with intra-frame frequency attention, and `ONNXModel.forward` applies power-law magnitude compression to `(real, imag)` STFT spectra, predicts a complex mask, multiplies in complex form, uncompresses, and returns enhanced spectra plus updated recurrent caches (note `noncausal` siblings also exist under `models/`, so verify config causality explicitly!).

**Per-chunk streaming caches.** In [`scripts/export_onnx.py`](https://github.com/aask1357/fastenhancer/blob/f85223bd546b27f39dc0744e0310dcd246f750a4/scripts/export_onnx.py) and [`scripts/test_onnx.py`](https://github.com/aask1357/fastenhancer/blob/f85223bd546b27f39dc0744e0310dcd246f750a4/scripts/test_onnx.py), `Model.forward` takes `wav_in`, STFT overlap cache, iSTFT overlap-add cache, and `cache_in_*` recurrent tensors—zero-initialized at utterance start and fed back chunk-by-chunk (zeroing per chunk breaks streaming continuity, whereas omitting reset between utterances leaks prior speaker/noise state).

Default `test_onnx.py` runs at `16` kHz, `n_fft = 512`, `hop_size = 256` (`16` ms hop), zero-padding the tail to flush caches, concatenating chunks, stripping `n_fft - hop_size = 256` leading samples (`16` ms offline OLA alignment, distinct from live hardware capture+compute+queue latency!), and slicing `:length`.

Moreover, [lines 42–51 of `test_onnx.py`](https://github.com/aask1357/fastenhancer/blob/f85223bd546b27f39dc0744e0310dcd246f750a4/scripts/test_onnx.py#L42-L51) include zero-pad flush hops inside the printed RTF denominator (`(idx + hop_size) / sr` with loop bound `length + n_fft - hop_size`): a `1000`-sample input runs `5` hops (`1280`-sample denominator), printing `0.125` on `10` ms mocked time versus true source-duration RTF $0.01 / (1000 / 16000) = 0.16$ (and `16384` vs. `16000` on `1.0` s input; see [Exercise E10-32](../../../../chapters_en/10_engineering-practice.md#e10-32) and I32's mocked AST wrapper audit).

**Minimal verification protocol.** Once version-matched checkpoints are acquired, pin mono WAV input, checkpoint SHA-256, chunk size, ONNX Runtime version, thread count, and EP; verify whole-utterance vs. cache-preserving chunked parity against `--test-streaming` alongside a per-chunk cache-zeroing counterexample, and profile cold start, steady state, tail drain, and session reset separately (remembering FastEnhancer has neither far-end AEC reference nor multi-channel WPE prediction). Checking the [arXiv record](https://arxiv.org/abs/2509.21867) and [official ICASSP 2026 program (`PaperNum=3646`, `SLP-P2.8` on 2026-05-05)](https://www.cmsworkshops.com/ICASSP2026/view_paper.php?PaperNum=3646&bare=1) on 2026-09-29 confirms conference presentation while keeping `16` kHz original runs separate from the `48` kHz C port in [I31](#faster-enhancer-runtime).


<a id="i29stkdelayl"></a>

### I29: STK `DelayL` Fractional Delay, Magnitude Attenuation, and Cross-Block State

From **STK** (*The Synthesis ToolKit in C++*, pinned at [`6aacd357d76250bb7da2b1ddf675651828784bbc`](https://github.com/thestk/stk/tree/6aacd357d76250bb7da2b1ddf675651828784bbc), [`MIT-STK` license](https://github.com/thestk/stk/blob/6aacd357d76250bb7da2b1ddf675651828784bbc/LICENSE) verified 2026-09-26; sparse `include/src/README/LICENSE` checkout excluding `rawwaves` and audio drivers), we study fixed linear-interpolating fractional delay `DelayL`: `src/DelayL.cpp` allocates the delay buffer, `include/DelayL.h::setDelay` sets read-pointer lag and linear interpolation weights relative to the write pointer, `tick` writes the current input sample *before* reading interpolated output (supporting `0` delay cleanly while advancing pointers once per `tick`, whereas `nextOut` merely reads cached output without advancing!), and `include/Filter.h::clear` zeros buffer history (allocate max delay capacity during initialization rather than inside real-time callbacks).

| Inspected Method | State Variable / Step Verified | Independent Mathematical Criterion |
|---|---|---|
| `DelayL.h::setDelay` | Read pointer trails write pointer in sample units | Delays `0`, `0.5`, `1.0` yield current sample, 2-point average, and 1-sample delayed value |
| `DelayL.h::tick` / `nextOut` | Write-before-read ordering, linear weights, ring wraparound | Matches standalone causal 2-tap FIR sample-for-sample |
| `Filter.h::clear` | Zeros delay-line history (not a normal chunk-boundary operation) | `37`-sample chunked streaming matches whole-buffer; per-chunk `clear()` serves as failure control |

**Executed locally.** [`stk_delay_probe.cpp`](../../ch10/examples/stk_delay_probe.cpp) and [`run_stk_delay_probe.py`](../../ch10/examples/run_stk_delay_probe.py) compile unmodified `Stk.cpp` and `DelayL.cpp` in `/tmp` on a `512`-sample, `16` kHz two-tone input (`500` Hz + `6000` Hz, amplitude `0.18` each, zero initial state, distinct from the `2` s listening WAV).

On Darwin arm64 (AppleClang 21.0.0, `C++17 -O2`, default `StkFloat = double`), max absolute error against independent 2-tap FIR is `0` across delays `0`, `0.5`, and `1.0`, `0` across state-preserving `37`-sample chunks, and `0.1519102855` when calling `clear()` before each `37`-sample chunk ([`stk_delay.json`](../../ch10/reports/stk_delay.json)). Re-run offline via `.venv/bin/python -m codes.chapters.ch10.examples.run_stk_delay_probe --report codes/chapters/ch10/reports/stk_delay_current.json` (`StkFrames` multi-channel overloads and modulated delays were not tested).

Always evaluate magnitude droop alongside phase delay: at `6000` Hz (`16` kHz rate), one `0.5`-sample linear interpolation attenuates amplitude to `≈ 0.382683`, and cascading two `0.5`-sample stages drops amplitude to `≈ 0.146447` despite summing to `1.0` sample of phase delay—see [Appendix B Exercise E17-02 and Figure 41](../../../../chapters_en/17_appendix-guide.md#e17-02).

<a id="rnnoise-stream-contract"></a>

### I30: RNNoise Amplitude Scaling, Frame State, and Demo Head/Tail Framing

Combining Bark-scale pitch DSP with GRU band-gain estimation ([Valin, *A Hybrid DSP/Deep Learning Approach to Real-Time Full-Band Speech Enhancement*, MMSP 2018](https://arxiv.org/pdf/1709.08243)), `xiph/rnnoise` is pinned at `70f1d256acd4b34a572f999a05c87bf00b67730d` (`BSD-3-Clause` in `COPYING`; later upstream model architectures differ in parameter count from the 2018 paper).

#### Why `float` Buffers Do Not Imply Unit-Normalized Amplitude

[`src/denoise.h`](https://github.com/xiph/rnnoise/blob/70f1d256acd4b34a572f999a05c87bf00b67730d/src/denoise.h) defines `FRAME_SIZE = 480` (`10` ms hop at `48` kHz, `960`-sample analysis window; query `rnnoise_get_frame_size()` dynamically), while `rnnoise_process_frame` updates analysis/synthesis/pitch/GRU state, returns VAD probability, and writes enhanced audio to `out`.

Crucially, [`examples/rnnoise_demo.c`](https://github.com/xiph/rnnoise/blob/70f1d256acd4b34a572f999a05c87bf00b67730d/examples/rnnoise_demo.c) casts raw headerless signed 16-bit `short` samples directly to `float` **without dividing by `32768`**! Host pipelines delivering unit-normalized `[-1.0f, 1.0f]` PortAudio floats must multiply by `32768.0f` on entry and scale/clamp on exit—having `float*` on both sides never guarantees matching physical units!

#### Head/Tail Sample Framing and Model Object Lifetime

In `rnnoise_demo.c`, frame `0` (`first = 1`) primes state **without writing output**, subsequent full `480`-sample frames write `480` samples each, and any `< 480`-sample remainder at EOF is dropped without zero-pad draining—processing `1440` samples (`3` frames) thus writes `960` samples (a demo-wrapper convention rather than universal delay compensation). Furthermore, [`include/rnnoise.h`](https://github.com/xiph/rnnoise/blob/70f1d256acd4b34a572f999a05c87bf00b67730d/include/rnnoise.h) requires any custom `RNNModel` (and its backing buffer or `FILE*`) to outlive all `DenoiseState` instances referencing it.

Audited statically on 2026-09-28 and compiled natively with identity `create/process/destroy` stubs on 2026-10-01 (I32, verifying amplitude pass-through and head/tail counts without compiling neural weights or running `autogen.sh` network downloads), full evaluation requires pinning model SHA-256 and testing exact frame accounting separately.

<a id="faster-enhancer-runtime"></a>

### I31: `faster-enhancer.c` Custom Quantized C Runtime and Real-Time Pacing

Included strictly as an embedded runtime case study rather than a new model family, Kim's [*faster-enhancer.c* preprint (2026-07-28)](https://arxiv.org/html/2607.25350v1) is an independent unaffiliated C/SIMD port of `48` kHz `FastEnhancer-Medium`, pinned (§3.1) at [`kdrkdrkdr/faster-enhancer.c`](https://github.com/kdrkdrkdr/faster-enhancer.c/tree/7d78dab11bca854fb200426c621a96d7b37a2983):

```text
7d78dab11bca854fb200426c621a96d7b37a2983
```

#### Architectural Changes Relative to General-Purpose Inference Engines

Preserving the network topology while replacing generic graph engines with fixed-shape C/SIMD kernels, dynamic per-frame `int8` activation quantization, `fp16` recurrent state buffers, and init-time ISA dispatch, §2 clamps `int8` values to symmetric `[-127, 127]` so pairwise multiply-adds satisfy $2 \times 127^2 = 32258 < 32767$ (`INT16_MAX`) across arm64 NEON / x86_64 AVX2 without overflow. Public `include/fe.h` exposes init, `320`-sample frame step (`≈ 6.67` ms at `48` kHz, `704`-sample STFT OLA delay `≈ 14.67` ms per §2.1), reset, and free on a single global single-threaded instance with **no scalar C fallback** (unusable as-is on bare MCUs or multi-stream servers without refactoring).

#### Why Frame-Deadline Pacing Differs from Unpaced File Throughput

Sections 3.1 and 3.4 separate unpaced tight-loop file benchmarks from `6.67` ms deadline-paced frame delivery alongside warmup, repeat, and CPU power telemetry—illustrating why DVFS frequency scaling under real-time pacing can trigger deadline misses invisible in unpaced mean RTF.

#### Licensing, Sparse Source Selection, and Reproduction Gaps

Verified on 2026-09-28 (`MIT` `LICENSE`, `NOTICE` retaining FastEnhancer copyright and third-party kernel-design references; 57 plain-text files under `include/`, `src/`, `tests/`, `docs/` excluding `weights/fe.q8` and test WAVs; CMake `3.20` + C11, MSVC rejected; status in [`SOURCE_STATUS.json`](../SOURCE_STATUS.json), zero compilation or inference run here).

Re-auditing [`include/fe.h` lines 45–67](https://github.com/kdrkdrkdr/faster-enhancer.c/blob/7d78dab11bca854fb200426c621a96d7b37a2983/include/fe.h#L45-L67), [`src/fe_pipeline.c` lines 71–108](https://github.com/kdrkdrkdr/faster-enhancer.c/blob/7d78dab11bca854fb200426c621a96d7b37a2983/src/fe_pipeline.c#L71-L108), and [`tests/fe_run_file.c` lines 271–309](https://github.com/kdrkdrkdr/faster-enhancer.c/blob/7d78dab11bca854fb200426c621a96d7b37a2983/tests/fe_run_file.c#L271-L309) on 2026-10-05 confirms: (1) float weights/biases/embeddings alias pointers directly into the caller's raw weight blob, which must outlive `fe_free()`, while calling `fe_run()` before successful `fe_init()` leaves the output buffer unmodified (stale contents, not zeros!) and model hot-swapping requires `fe_free()` first; (2) `fe_init()` runs `128` silent warmup frames before zeroing streaming state (while dynamic activation scales still recompute every frame) and `sizeof` state excludes the live weight blob; and (3) paced-mode `elapsed / audio_s` includes inter-frame sleep waits (`≈ 1.0` even when compute is fast), distinct from per-frame `ft[t]` service-time quantiles (and excludes hardware capture/playback DMA!).

<a id="industrial-upstream-audit"></a>

### I32: Execution Scope of the Pinned Industrial Interface Diagnostics

Corresponding to Chapter 10 audits `SRC10-01..03`, the [historical 2026-09-28 JSON report](../../ch10/reports/industrial_upstream_interfaces.json) binds the [historical `bde483` script bytes](https://github.com/nanless/microphone-array-tutorial/blob/bde483bcc429a553aeaf8830ca3687ab46831db0/codes/chapters/ch10/examples/audit_industrial_upstream_interfaces.py), while the 2026-10-05 current re-verification appears in [§9](#industrial-current-interface-contracts). Verifying full commits, clean worktrees, and file/license SHA-256 digests before and after execution with strict `NaN`/`Infinity`-free JSON:

| Target Project | Action Executed in This Audit | Unexecuted Scope |
|---|---|---|
| `pystoi` | 12 calls to unmodified Python `pystoi.stoi`, logging exceptions, warnings, and return values | Author MATLAB/Octave parity, internal resampling, natural speech corpora |
| `DeepFilterNet` | Static inspection of pinned Rust `ladspa/src/lib.rs` + standalone Python queue/counter trace | `cargo` Rust compilation, neural checkpoints, PipeWire/LADSPA host execution |
| Other industrial repos | Pinned doc/header verification (native C probes archived in their own separate JSON reports) | Zero claims of live hardware acceptance |

Re-run into the current report inside a pre-existing SciPy-equipped isolated Python environment:

```bash
PYTHONDONTWRITEBYTECODE=1 /private/tmp/masp-appb-pra-venv/bin/python \
  codes/chapters/ch10/examples/audit_industrial_upstream_interfaces.py \
  --report codes/chapters/ch10/reports/industrial_upstream_interfaces_current.json
```

While that `/private/tmp/...` path reflects the local scratch venv used during our audit, standard [offline unit tests](../../../../tests/test_codes_industrial_upstream_interfaces.py) validate report invariants without requiring upstream caches, SciPy, or network access. Primary literature passes opened official PortAudio, ALSA, ITU-T P.835 / P.566 (verifying in-force vs. pre-published status directly), ITU-R BS.1770-5, EBU R128, ViSQOL v3, FastEnhancer ICASSP 2026, and `faster-enhancer.c` documents directly (recording honestly that Cees Taal's personal STOI PDF link timed out during that pass, so Python STOI findings rely strictly on pinned `pystoi` source and runtime probes).

<a id="industrial-controlled-contracts"></a>

#### 2026-10-01: Bounded Execution of Five Pinned Source Contracts

The 2026-10-01 [historical 5-item contracts report](../../ch10/reports/industrial_contracts.json) binds [script bytes at commit `bde483`](https://github.com/nanless/microphone-array-tutorial/blob/bde483bcc429a553aeaf8830ca3687ab46831db0/codes/chapters/ch10/examples/audit_industrial_contracts.py) (with 2026-10-05 re-verification in [§9](#industrial-current-interface-contracts)), leaving the 2026-09-28 report untouched while verifying `origin`, `HEAD`, cleanliness, file/license SHA-256 and Git blobs before and after compiling/extracting original routines in `/tmp`:

| Pinned Upstream Target | Executed Code & Stimuli | Independent Verification Result | Explicit Test Doubles & Unexecuted Scope |
|---|---|---|---|
| CMSIS-DSP `Q15` FIR | 2 unmodified C translation units, 4 taps; half-LSB, `-32768` saturation, asymmetric whole vs. `2+1+4` chunks | Half-LSB yields `[0, 0, -1, -1, 0]`, saturates at `32767`; matches reversed-tap integer convolution | Minimal typedef header & scalar `__SSAT`; zero Cortex-M SIMD/MVE or board cycle counts |
| SpeexDSP Noise Helpers | 4 original C helper bodies + original float `arch.h` macros; 2 bins, rising/falling power & minima rollover | Hard flags `[0, 1]`; noise `[2, 2] -> [3, 2] -> [2.5, 1.75]`; dual running-minima states verified | Custom C driver; excludes full `preprocess.c`, FFT, and Cohen soft-MCRA/IMCRA |
| WebRTC Legacy VAD Entry | Original rate/frame validator + `Process` router bodies; 4 rates × `10/20/30` ms + 6 invalid inputs | All 12 valid pairs pass; `44.1` kHz, `15` ms, `NULL` handle/audio, uninit, and `0` length return `-1`; 3 calls per rate | 4 `CalcVad*` stubs count calls only; zero GMM classification or speech ROC benchmarks |
| RNNoise Demo Loop | Complete unmodified `rnnoise_demo.c` `main`; `1440, 1457, 479, 0`-sample raw PCM16 | Writes `960, 960, 0, 0` samples (drops first `480` & trailing `17`); passes raw `1000.0f` instead of `1000/32768` | Identity copy stubs for `create/process/destroy`; zero GRU weights or denoising |
| FastEnhancer CLI Wrapper | Unmodified `test_onnx.py::main` AST; `1000, 16000, 256, 1, 0` samples with in-memory sink | `256`-sample hops, continuous cache threading, strips `256` then slices `:length`; `1000`-sample denominator is `1280` | Mocked loader, ORT session, deterministic `0/0.01` s clock, and RAM sink; `0`-length true RTF logged as `null` |

In the Speex helper check (`beta = 0.25`, `beta_1 = 0.75`, tol `1e-6`), exact equality $0.4 \times 10 = 4$ against minimum `4` leaves flag `0`, whereas minimum `3.99` sets flag `1`; subsequently, power `1 < 2` updates bin 1 downward to `1.75` despite its active flag, while running minima update from `[5, 7]` / tmp `[6, 4]` + current `[3, 8]` $\to$ `[3, 7]` / `[3, 4]` and rollover `[9, 2]` $\to$ `[3, 2]` / `[9, 2]`.

In the FastEnhancer wrapper check, deterministic `0` and `0.01` s mock timestamps pair with a ramp-emitting mock session incrementing recurrent state by `+1` per hop into an in-memory sink; `0`-sample input still runs one zero-padded hop and prints `0.625` upstream, whereas our true source-duration RTF logs `null` (division by zero source duration). Re-run offline via:

```bash
.venv/bin/python -B -m codes.chapters.ch10.examples.audit_industrial_contracts \
  --report codes/chapters/ch10/reports/industrial_contracts_current.json
.venv/bin/python -B -m unittest tests.test_codes_industrial_contracts -v
```

Omitting `--report` prints JSON to `stdout` and cleans `/tmp` afterward; if optional `_downloads/` checkouts are absent on a fresh clone, unit tests skip live upstream execution cleanly while failing hard on any version/digest mismatch.

<a id="distributed-network-deployment"></a>

## 8. Distributed Node Observations, Filter Updates, and Network Timelines

Verified on 2026-10-04 across original papers, pinned upstream code, and official network RFCs/documentation, this section connects Chapter 10's clock interfaces (§10.2) and multi-node engineering (§10.7) with [Extended Topic II (Chapter 13)](../../../../chapters_en/13_distributed-enhancement.md)—extending Chapter 5's linear MWF objective to address what nodes exchange over links, when local vs. broadcast filters update, and how packet timestamps map onto acoustic sample timelines.

### I33: Signal Model and Convergence Guarantees of Original DANSE

Align notation first between Bertrand & Moonen's 2010 two-part IEEE TSP papers—[Part I: Sequential Node Updating](https://homes.esat.kuleuven.be/~abertran/reports/09-65.pdf) (`10.1109/TSP.2010.2052612`) and [Part II: Simultaneous and Asynchronous Node Updating](https://homes.esat.kuleuven.be/~abertran/reports/09-178.pdf) (`10.1109/TSP.2010.2052613`)—and Chapter 13 of this textbook (where $K$ in Part I/II is NOT node count!):

| Mathematical Quantity | DANSE Part I / II Notation | This Book (Chapter 13) |
|---|---|---|
| Number of sensor nodes | $J$ | $U$ |
| Broadcast / compressed dimension per node | $K$ | $Q$ |
| Shared latent target source dimension | $Q$ | $S$ |

Thus upstream's dimension condition $K = Q$ corresponds to $Q = S$ here: on a shared frequency bin and sample timeline, node $u$ observes $\vec{x}_u$ and targets $\vec{d}_u = \mathbf{B}_u \vec{s}$ for a common latent source vector $\vec{s}$ (requiring nonsingular square $\mathbf{B}_u$ when $Q = S$). Centralized LMMSE solves $\mathbf{R}_{xx}\mathbf{W}_u^\star = \mathbf{R}_{x d_u}$; in DANSE, node $u$ stacks its raw local channels $\vec{x}_u$ with $Q$-dimensional compressed broadcasts from peer nodes, solves a local LMMSE system, and updates its future outgoing broadcast compressor from the local sensor block of its own weight matrix—broadcasting linear combinations of *raw local observations only* rather than re-broadcasting fused outputs containing remote signals ([Part I §§II, III-A, IV-A, Theorems III.1 & IV.1](https://homes.esat.kuleuven.be/~abertran/reports/09-65.pdf)).

Sequential convergence theorems assume nonsingular $\mathbf{R}_{xx}$, full-rank shared-subspace transforms $\mathbf{B}_u$, **exact stationary second-order statistics**, and round-robin updates—proving that *adaptive* $Q = S$ compression attains centralized LMMSE, NOT that arbitrary *fixed* beamformers are lossless (§13.8 provides an explicit spatially correlated noise counterexample), nor guaranteeing optimality under $Q < S$ compression loss, finite-snapshot VAD errors, or rank-deficient pseudoinverse regimes ([Part I §IV-C: *DANSE Under Rank Deficiency*](https://homes.esat.kuleuven.be/~abertran/reports/09-65.pdf)).

Under **simultaneous** node updates, peer broadcast compressors change concurrently. In [Part II §IV-B Eqs. (22)–(28), pp. 5297–5298](https://homes.esat.kuleuven.be/~abertran/reports/09-178.pdf), $rS\text{-DANSE}^+$ introduces an auxiliary subspace optimization over $\mathbf{G}$, and Theorem IV.1 requires diminishing step sizes $0 < \alpha_i \le 1$, $\alpha_i \to 0$, $\sum_i \alpha_i = \infty$ alongside exact statistics; simplified §IV-C $rS$-DANSE (omitting $\mathbf{G}$ optimization and using constant $\alpha = 0.5$ or $0.7$) is an empirical relaxation that cannot invoke Theorem IV.1. Meanwhile, Part II §IV-D's parameter-asynchronous DANSE permits irregular filter update instants across nodes while still requiring persistent per-node updates and **sample-level clock synchronization** (Footnote 8).

### I34: Filter Update Ordering and Paper Hyperparameters in WOLA-DANSE Source Code

Section III of the [IWAENC 2010 WOLA-DANSE paper](https://homes.esat.kuleuven.be/~abertran/reports/IWAENC10.pdf) embeds frequency-domain SDW-MWF inside root-Hann, 50%-overlap WOLA filterbanks with exponential covariance averaging and zero link delay, while §IV runs experiments at `32` kHz, `512`-point window, $\mu = 5$, $\alpha = 0.5$, and oracle VAD—hyperparameters that cannot be inferred from demo script defaults.

In pinned author [`WOLA_DANSE1.m` v1.4](https://github.com/AlexanderBertrandLab/Old_Code/blob/a24b73fcd2dc028659535d07bb08068b06108616/WOLA_DANSE1.m), `fs` has no default (`16` kHz in demo usage, giving default `L = 512 * fs / 16000 = 512`), while actual defaults set `mu = 1`, `alpha = 0.7`, simultaneous updates enabled, `2` s covariance / `0.2` s external-weight half-lives, and `3` s equivalent frame accumulation per VAD class before enabling filter updates.

Static control-flow inspection ([current distributed contracts report](../../ch13/reports/upstream_distributed_contracts_current.json), zero MATLAB execution) documents four interface specifics: (1) oracle VAD inspects only sample `1` of each frame, and each frame **broadcasts current `Wext` $\to$ smooths `Wext` toward the *old* target $\to$ computes the new target**, so a new target first influences transmitted broadcasts two broadcast frames later; (2) even in sequential round-robin token mode, internal weights `Wint` update every frame while the token gates only external broadcast target updates; (3) symmetric `hann` omits point-wise COLA normalization and tail draining; and (4) `inv()` applies zero diagonal loading while `eig()` picks the algebraically largest real eigenvalue followed by `abs()`, distinct from noise-whitened GEVD or PSD projection. Converting such a zero-delay script into networked code requires explicit sample sequence numbers, broadcast compressor epochs, jitter buffers, and missing-packet rules.

### I35: What Tree Topologies, Low-Rank GEVD Statistics, and SRO Compensation Each Change

| Method & Primary Reference | Core Algorithmic Modification | Validity Preconditions & Scope Boundaries |
|---|---|---|
| [TI-DANSE, §§III–IV, Theorems 1/2](https://homes.esat.kuleuven.be/~abertran/reports/15-87.pdf) (*IEEE TSIPN* `3(1):130–144`, **2017**; DOI string contains `2016`) | Broadcasts gauge-normalized $\mathbf{P}_u = \mathbf{W}_{uu}\mathbf{G}_u^{-1}$, sums in-network along a tree, and subtracts self-contribution so local dimension becomes $M_u + Q$ | Requires invertible $\mathbf{G}_u$, connected tree, exact same-sample summation, and common target; keep output optimality distinct from parameter coordinate convergence |
| [GEVD-DANSE, §III Eqs. (9)–(17), §IV](https://ftp.esat.kuleuven.be/SISTA/abertran/reports/GEVD_DS_TSP2016.pdf) (2016) | Performs noise-metric GEVD on $(\mathbf{R}_{yy}, \mathbf{R}_{nn})$, retaining chosen rank with gains $1 - 1/\lambda$ before distributed updates | Requires valid nonsingular noise statistics; retaining rank $< S$ targets centralized low-rank GEVD-MWF rather than full-rank LMMSE (distinct from `eig()` in `WOLA_DANSE1.m`) |
| [TI-GEVD-DANSE, §§III–IV](https://ftp.esat.kuleuven.be/sista/pdidier/eusipco2024/eusipco24_pdidier_submission_v1.pdf) (EUSIPCO 2024 author draft) | Combines local GEVD inside tree TI-DANSE with shared coordinate gauge normalization to suppress filter scale drift | Identical gauge transform must hit covariances and filters together (never clip node weights independently!); empirical under-rank plots are not formal proofs |
| [TI-DANSE+ Conf. v1, §§III–V](https://arxiv.org/html/2506.20001v1) (2025) | Updating root keeps per-neighbor branch partial sums separate (local DOFs scale with root degree) plus max-updating-node tree selection | Stationary-statistics convergence demonstration without formal proof; iteration savings do not convert directly to hardware wall time |
| [TI-DANSE+ Extended v2, §§III-C–G, §IV](https://arxiv.org/html/2506.02797v2) (2026-03-10 preprint) | Adds Theorems 1/2 plus communication, tree-selection, and covariance-reset analysis | Mixing stale SCMs across broadcast epochs corrupts coordinates; §III-G low-rank GEVD lacks an equivalent proof, and §IV-D centralized filter application is not live WOLA networking |
| [SRO-Aware DANSE, §IV-A–C Eqs. (11)–(20)](https://arxiv.org/html/2211.02489v2) (*IEEE OJSP* 2023) | Estimates SRO from cross-frame coherence drift, accumulates fractional phase shifts + integer buffer slips, and introduces sample-streaming WOLA convolutive approximation | Assumes stationary nodes/sources, constant SRO, slow broadcast drift, and lossless links; convolutive WOLA approximation still incurs analysis framing latency |

Note bibliographic precision: TI-DANSE appeared in *IEEE TSIPN* Vol. 3, No. 1, pp. 130–144, **2017** (`10.1109/TSIPN.2016.2623095`), and because Eq. (11) in the EUSIPCO 2024 TI-GEVD draft writes a target eigenvalue difference inconsistent with Eq. (12)'s gain, our whitening examples follow the [2016 GEVD paper](https://ftp.esat.kuleuven.be/SISTA/abertran/reports/GEVD_DS_TSP2016.pdf) ($\lambda - 1$ and $1 - 1/\lambda$).

Keep three SRO regimes distinct: Chapter 10's regression uses exact synthetic timestamps; 2023 [SRO-DANSE §IV-B](https://arxiv.org/html/2211.02489v2) tracks blind inter-node coherence drift across stationary geometry (noting that `100` ppm at `16` kHz drifts `1.6` samples/s = `1` sample every `0.625` s, whereas `1` full `512`-sample hop takes `320` s—never treat `320` s as a window over which phase drift is negligible!); and [DWACD (ICASSP 2022 §§II–III)](https://ris.uni-paderborn.de/download/33807/48990/gburrek_icassp22.pdf) accommodates moving talkers across short locally stationary active windows while separating acoustic ToF/STO from crystal SRO ([bounded 3-helper audit](04_source_reproduction.md#distributed-reproduction)).

Finally, updating a peer's broadcast compressor changes the coordinate basis of received signals: [TI-DANSE+ 2026 §III-G1](https://arxiv.org/html/2506.02797v2) recommends accumulating SCMs strictly within a fixed compressor epoch, resetting across updates, and waiting for sufficient speech+noise and noise-only support—whereas under a known *invertible* gauge change within the same subspace, old covariances can be transformed algebraically without discarding history ([Exercise E13-25](../../../../chapters_en/13_distributed-enhancement.md#e13-25)).

### I36: RTP Timestamps, Reference Clocks, and Dante Platform Boundaries

Separate packet NIC arrival time, first-sample ADC timestamp, and scheduled playout time. In [RFC 3550 §5.1 & §6.4.1](https://www.rfc-editor.org/rfc/rfc3550.html#section-5.1), RTP sequence numbers increment per packet while RTP timestamps mark the first audio sample on a stream's media clock starting from a **random offset**—so raw RTP timestamps across two nodes cannot be subtracted directly without RTCP Sender Report `(NTP, RTP)` mappings, RFC 3550 interarrival jitter is packet network jitter rather than crystal ADC ppm, and sequence-number gaps require payload/timestamp inspection to determine exact lost sample counts.

Similarly, [RFC 7273 §§4–6](https://www.rfc-editor.org/rfc/rfc7273.html) separates timestamp reference clocks (`ts-refclk`) from media sample clocks (`mediaclk`): disciplining an OS wall clock via PTP never proves that an asynchronous USB/I2S audio ADC crystal is hardware-locked to that PTP clock!

Audinate's [Dante Controller Clock Synchronization](https://dev.audinate.com/GA/dante-controller/userguide/webhelp/content/clock_synchronization.htm) and [Dante Latency](https://dev.audinate.com/GA/dante-controller/userguide/webhelp/content/latency.htm) pages document leader/follower clocking (PTPv1 by default, plus PTPv2 when AES67/RTP is enabled; hardware PLLs or external word clock on dedicated Dante silicon vs. host computer clocking on Dante Virtual Soundcard) and deterministic playout latency (`max(tx_setting, rx_setting)`, with larger minimums on general-purpose PCs)—platform-specific guarantees that never transfer automatically to commodity Wi-Fi UDP streams (no Dante hardware or PTP network was connected here).

### I37: Compute Payload Bitrates First Before Budgeting Link Bandwidth

With $U$ nodes, $M_u$ local mics per node, and $Q$ **real time-domain PCM** broadcast channels at rate $f_s$ and $b$ bits/sample, raw node payload is $f_s b M_u$ bit/s vs. compressed payload $f_s b Q$ bit/s (`768` vs. `256` kbit/s for `16` kHz, PCM16, $M_u = 3, Q = 1$).

Across $U = 4$ nodes, pairwise unicast transmits $U(U-1) = 12$ streams (`3.072` Mbit/s), true single-send broadcast/multicast emits $U = 4$ streams (`1.024` Mbit/s), and tree up-sum + down-broadcast sends $2(U-1) = 6$ directed edge streams (`1.536` Mbit/s, [Exercise E13-23](../../../../chapters_en/13_distributed-enhancement.md#e13-23), plus $8 Q^2$ bytes per `complex64` $Q \times Q$ matrix update).

Overlapped STFT streams (`512` FFT, `256` hop = `62.5` frames/s at `16` kHz, `257` one-sided `complex64` bins) carry $257 \times 8 \times 62.5 = 128500$ B/s = `1.028` Mbit/s per channel before framing metadata ([Exercise E13-15](../../../../chapters_en/13_distributed-enhancement.md#e13-15)).

Adding explicit L3/L4 packet headers to `10` ms (`160`-sample = `320` B) mono PCM16 packets—`12` B base RTP header ([RFC 3550 §5.1](https://www.rfc-editor.org/rfc/rfc3550.html#section-5.1)), `8` B UDP header ([RFC 768](https://www.rfc-editor.org/rfc/rfc768.html)), and `20` B optionless IPv4 header ([RFC 791 §3.1](https://www.rfc-editor.org/rfc/rfc791.html#section-3.1))—yields `360` B $\times 100\text{ pkt/s} = 288$ kbit/s (`+12.5%` over `256` kbit/s payload, still excluding Ethernet/VLAN/802.11 MAC contention/retransmissions/SRTP/PTP, and always evaluated alongside compression subspace loss, quantization noise, and jitter-buffer packet-drop rates).


<a id="industrial-current-interface-contracts"></a>

## 9. Chapter 10 Current Upstream Interface Re-Verification and Runtime Boundaries

On 2026-10-05, four bounded Chapter 10 interface suites were re-executed on the pinned upstream checkouts, saving standalone current reports—[five C/AST contracts](../../ch10/reports/industrial_contracts_current.json), [`pystoi` & DeepFilterNet queue audit](../../ch10/reports/industrial_upstream_interfaces_current.json), [three-library native C probe](../../ch10/reports/industrial_interfaces_current.json), and [STK `DelayL` native probe](../../ch10/reports/stk_delay_current.json)—while keeping all four historical reports bound to their original dates, tool digests, and inputs (printing to `stdout` by default and rejecting any in-repo destination other than each script's designated `*_current.json` path).

Maintaining identical measurement semantics to the historical three-library C run, all builds execute in external `/tmp` directories with pre- and post-run `origin`, `HEAD`, license, and blob verification while keeping `selection_mismatch` statuses separate from bounded helper checks (recording real scalar CMSIS FIR and Speex helpers alongside explicit WebRTC VAD, RNNoise demo, and FastEnhancer stubs without running neural weights, ARM intrinsics, or live audio hardware).

Native build reports tabulate *pre-checked candidate source files* (`137 / 19 / 362` blob-verified regular files across the three C libraries) separately from *actual compiler inputs* (`104` `.d` dependency files, `149` upstream translation inputs, `2` generated config headers, and `190` system SDK headers)—never equating pre-checked file counts with executed algorithm counts!

Uncompiled symlinks in `libsoxr` and local CRLF/blob-differing Win32 `.def` export files in `libsamplerate` are itemized explicitly with exclusion reasons, failing hard if any `.d` dependency ever touches an unverified file rather than mutating local caches; similarly, `stk` pre-checks `212` members while compiling strictly `2` upstream `.cpp` files plus headers.

`pystoi` preserves all 12 exception, `1e-5` sentinel, and `finite_zero_energy_invalid` cases without treating finite zero-energy outputs as valid intelligibility scores.

DeepFilterNet performs static Rust control-flow inspection and standalone Python counter arithmetic only—confirming that from initial count `0`, strict `t_proc_change > 1000` is first satisfied on host callback `1002` (assuming `rtf < 0.5`, delay, and queue conditions also hold), or on subsequent stable callback `1001` following an underrun reset + trailing increment, keeping host callback quantum distinct from model `hop` without compiling Rust or running a LADSPA host.

`faster-enhancer.c`'s 57-file plain-text sparse selection has its blobs and MIT/NOTICE licenses verified while weights and compilation remain unexecuted.

Finally, cumulative `rtf` and optional `frame_service_rtf` are defined separately: cumulative RTF divides total processing time over a declared window by total newly arrived source audio duration, requiring producers to log window bounds, frame counts, and whether pacing waits enter the numerator—since stateless JSON schema validators check numeric types and bounds without proving cross-report clock, scope, or reset equivalence.
