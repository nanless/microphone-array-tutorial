# Echo Cancellation, Dereverberation, and Speech Separation: Algorithms, Industrial Implementations, and Source Code Research

🌐 [简体中文](../research/02_aec_wpe_separation.md) ｜ **English**

Re-verification date for AEC playout gain notifications, historical echo tails, current source contracts, and LMPAN/DiffVQE manuscripts: 2026-10-04; re-verification date for SGMSE+/StoRM active sampler dispatch and zero-peak boundaries: 2026-10-02; re-verification date for BSS/GSS/CSS pinned source contracts, ArrayDPS tap points and hyperparameter layers, and neural method provenance: 2026-10-01; re-verification date for WPE primary literature, NARA multi-channel history ordering, and NeMo pinned source code: 2026-10-01; re-verification date for AEC upstream APA functions, APM channel boundaries, and neural Kalman taxonomy: 2026-10-01; initial verification date for AEC primary papers, interface diagnostics, and recent candidates: 2026-09-28; original baseline verification date: 2026-09-22; re-verification date for AEC industrial interfaces and paired dataset metadata: 2026-09-23; re-verification date for BSS/GSS source entry points and experiments: 2026-09-24; WPE and cACGMM source entry points expanded and verified offline on 2026-09-29; derivation exercises plus Stream.FM and TF-Locoformer static source audits completed on 2026-09-26; Chapter 15 target speaker extraction plain-text sparse selection verified on 2026-11-17. Corresponding main text: [Chapter 6](../../../../chapters_en/06_aec.md), [Chapter 7](../../../../chapters_en/07_wpe-dereverberation.md), [Chapter 8](../../../../chapters_en/08_speech-separation.md), and [Chapter 15](../../../../chapters_en/15_target-extraction-diarization.md). Organized by signal-processing task, this guide details algorithm mechanics, exact upstream file locations, and reproducible experiments; canonical commit hashes and license boundaries are governed by [SOURCES.lock.json](../SOURCES.lock.json).

An "external implementation" entry indicates where source code executing that algorithm resides—not that this textbook has trained, compiled, or benchmarked every mode of that external system. Results actually executed locally in this repository are itemized explicitly, including the [teaching baseline test suite](../../../../tests/test_codes_aec_wpe_sep_track.py) and the independent `nara_wpe` parity check in W01; entries without attached local execution logs specify experimental reproduction designs. Upstream code, model checkpoints, and speech datasets are each governed by their own independent license terms.

## 1. Define Signal Interfaces and Target Problems First

Acoustic Echo Cancellation (AEC) exploits a known loudspeaker playback reference to estimate and subtract the resulting acoustic echo; Weighted Prediction Error (WPE) uses delayed history of the microphone signals themselves to predict and subtract late reverberation; Blind Source Separation (BSS) relies on statistical independence, spectral low-rank structure, or spatial covariance models to separate multiple unknown sources; and supervised neural separation/extraction further requires clean training targets or target-speaker conditioning clues. Lumping all of these distinct modules under the generic label "noise reduction" obscures their fundamentally different input contracts and acceptance criteria.

| Task | Required Inputs | Mandatory Persistent State or Conditions | Minimal Acceptance Targets |
|---|---|---|---|
| AEC | Synchronized microphone and playback reference | Reference ring buffer, adaptive filter taps, bulk delay, step-size/freeze control | Far-end single-talk cancellation, double-talk near-end preservation, echo-path change |
| WPE | Phase-preserving single- or multi-channel STFT | Delayed regressor history, prediction taps, time-varying power & weighted covariance | Late-tail suppression, early-reflection/direct-path preservation, startup & prolonged silence |
| BSS | Multi-channel mixture and declared source count | Demixing matrices or source SCMs, projection-back scale, cross-frequency permutation alignment | Separation gain, scale restoration, local optima sensitivity, rank-deficient mixtures |
| GSS | Multi-channel recording and speaker activity diarization | Activity masks, spatial cACGMM posteriors, speech/noise SCMs, reference mic | Diarization errors, context padding, cross-frequency cluster consistency |
| Neural Separation / TSE | Mixture audio, plus enrollment clue when conditioned | Checkpoint weights, utterance/running normalization, temporal context, output slots | Out-of-domain mixtures, missing sources, target over-suppression, strict causality |
| CSS | Long unsegmented recording and chunking schedule | Overlapped windows, cross-chunk permutation stitching, fixed output stream count | Speaker turn-taking, silence gaps, duplicate leakage across slots, long-run RAM |

Even array axis ordering can never be guessed from variable names alone: whereas general modules in this repository use `(channels, freq_bins, frames)`, our teaching WPE and mask-based SCM kernels use `(freq_bins, channels, frames)`, and ESPnet's `DNN_WPE` expects `(batch, frames, channels, freq_bins)`. Always annotate axis permutations explicitly and test them on asymmetric two-channel signals rather than relying on coincidental shape compatibility.

<a id="aec"></a>

## 2. AEC: From Adaptive Filtering to Complete Audio Processing Pipelines

### A01 LMS, NLMS, and Leaky Updates

Normalized Least Mean Squares (NLMS) normalizes each gradient step by the instantaneous energy of the reference regressor vector, stabilizing adaptation across varying playback levels. Leaky NLMS additionally shrinks filter weights slightly toward zero on every step—bounding weight drift in poorly excited frequency bands at the cost of a small steady-state bias. Both require that the reference tap genuinely capture the loudspeaker signal generating the acoustic echo. Derivations and runnable baseline code appear in [§6.1.2](../../../../chapters_en/06_aec.md#sec-6-1-2) and `nlms()` in [`aec.py`](../../ch06/core/aec.py).

Recommended code-reading trace: `reference delay line -> echo estimate -> a priori error -> energy normalization denominator -> freeze gate -> tap update`. A minimal verification experiment feeds a known short FIR echo path: verify weight convergence under near-end silence first, then inject near-end speech mid-utterance while freezing updates. For a failure control, shift pure bulk delay beyond the FIR filter length $L$—confirming that running more iterations can never recover taps outside the model support. When sweeping leakage factors, report parameter bias alongside output residual power.

Executed [Exercise E06-11](../../../../chapters_en/06_aec.md#e06-11) cleanly separates *statistical expectation bias* from *finite-sample online trajectory drift*: on a deterministic 4-sample sequence where reference and near-end signals are strictly orthogonal, the exact batch least-squares tap is `0.8`, yet four online NLMS steps initialized right at `0.8` drift to `0.89375` (whereas setting the near-end signal to $0.25\times$ the reference shifts the batch LS solution itself to `1.05`). A single finite-sample trajectory endpoint must never be misdiagnosed as asymptotic estimator bias.

### A02 FDAF, Multi-Tap Frequency-Domain NLMS, MDF/PBFDAF, and Overlap-Save

Single-block Frequency-Domain Adaptive Filtering (FDAF) replaces long time-domain convolution with FFT multiplication, while Multi-Delay Filtering (MDF) / Partitioned-Block FDAF (PBFDAF) splits a long echo tail into $P$ shorter partitions so algorithmic block latency stays small while covering long room tails. Fixing a single DFT bin, whenever the $P$ taps of "frequency-domain multi-tap NLMS" correspond to the current and $P-1$ past reference blocks and share the exact error spectrum and normalization denominator of [Eq. (6-5)](../../../../chapters_en/06_aec.md#sec-6-2), the two updates are component-wise identical rather than distinct algorithms—though identical terminology across papers never guarantees identical step-size control ([Páez Borrallo & García Otero, 1992](https://doi.org/10.1016/0165-1684(92)90077-A) is retained as a historical citation without inferring step-size heuristics from its abstract alone; constrained overlap-save derivations below follow Shynk's full text).

Never treat "multi-tap frequency-domain NLMS" as an unambiguous specification without checking three distinct tap domains: (1) sample-wise time-domain NLMS with $L$ sample delays; (2) partitioned overlap-save PBFDAF holding $P$ reference-block spectra per bin to represent an exact $L$-tap linear FIR; and (3) per-bin cross-frame filtering in a windowed STFT filterbank. Only (2) matches this book's PBFDAF baseline block-for-block when block size, history, valid-region error windowing, power smoothing, step size, time-domain gradient projection, and freeze scheduling all match; by contrast, windowed STFT subband filtering omits cross-band aliasing terms induced by sample-level FIR paths even under perfect-reconstruction analysis/synthesis windows—as derived originally in [Avargel & Cohen (2007), §II Eqs. (7)–(11)](https://webee.technion.ac.il/Sites/People/IsraelCohen/Publications/TASL_May2007.pdf) and illustrated in Chapter 6 via a 2-point DFT unit-sample delay counterexample.

Denoting the $2N$-point DFT matrix by $\mathbf{F}$ and the valid-half selection matrix by $\mathbf{S}$, the exact overlap-save error spectrum is $\mathbf{E} = \mathbf{F}\mathbf{S}^\top(\vec{d} - \mathbf{S}\mathbf{F}^{-1}\mathbf{Q})$—never unwindowed $\mathbf{D} - \mathbf{Q}$. While constrained FFT correlation without per-bin power normalization computes the exact block-LMS gradient on valid FIR taps (still distinct from sample-by-sample in-block weight updates), per-bin preconditioning by $\sum_p |X_p|^2 + \delta$ in Eq. (6-5) alters the time-domain gradient direction and re-couples DFT bins upon time-domain zero-padding projection—nor does that scalar sum invert the full $P \times P$ cross-partition reference correlation matrix of exact block RLS; see [Shynk (1992), §IV, Fig. 5 & Eq. (30)](https://course.ece.cmu.edu/~ece792/handouts/Shynk92.pdf), [Schneider & Kellermann (2016), §§3–4](https://doi.org/10.1186/s13634-015-0302-2), [Soo & Pang (1990)](https://doi.org/10.1109/29.103078), pinned [SpeexDSP `libspeexdsp/mdf.c`](https://gitlab.xiph.org/xiph/speexdsp/-/blob/8e29a256ef0235ebbe7fcb8417b5ac7731eb8307/libspeexdsp/mdf.c), and [Worked Example 6-4](../../../../chapters_en/06_aec.md#sec-6-2).

When reading `speex_echo_state_init_mc()` and `speex_echo_cancellation()`, trace partition count $P$, reference spectrum ring history, overlap-save output slicing, zero-padded error FFT, and time-domain gradient projection—checking block size, FFT scaling, initialization, instantaneous vs. smoothed bin power, regularization floor, DTD/leakage step control, and projection cadence item-by-item before comparing two PBFDAF implementations (and remembering that Worked Example 6-2's complexity speedup compares against sample-domain time NLMS, not a renamed frequency-domain variant).

Minimal experiments should freeze adaptation first to verify bit-accurate sample-domain linear convolution equivalence across partitions, then enable adaptation and inject two distinct faults: (a) omitting valid-half output slicing (which immediately corrupts linear convolution with circular wraparound), and (b) disabling time-domain gradient projection (unconstrained MDF/FDAF, a legitimate algorithmic variant that alters weight trajectories without omitting valid-output slicing; cf. [Avargel & Cohen, 2007 §II](https://webee.technion.ac.il/Sites/People/IsraelCohen/Publications/TASL_May2007.pdf)).

This repository provides the self-contained [partitioned PBFDAF kernel](../../ch06/core/aec_partitioned.py), [two-block hand-check & system-identification demo](../../ch06/aec_partitioned_demo.py), and [unit tests](../../../../tests/test_codes_aec_partitioned.py): at $N = P = 2$, block 1's unconstrained partition-0 candidate exhibits a wrap-around 3rd time tap of $1/3$ that projection zeros out while partition 1 retains its valid 1st tap of $1/3$, alongside narrowband sine vs. broadband white-noise held-out identification controls proving that driving training residual to zero on a pure tone fails to identify the broadband room path (note also that Worked Example 6-2's `224 / 736` complex-op counts tally core FFTs and complex multiplies while excluding $|X|^2$ recomputations—at least `+32` ops/sample—and Python diagnostic IFFTs).

### A03 AUMDF, Partition Constraint Scheduling, and Proportional Adaptation

Comments in SpeexDSP's `mdf.c` specify **Alternately Updated MDF (AUMDF)**, which rotates time-domain gradient-projection FFT/IFFT pairs across partitions to halve transform cost while applying proportional adaptation weights across partitions. Because proportional allocation steers larger steps toward dominant early partitions, fast convergence on sparse impulsive paths never proves equal convergence speed on diffuse reverberant tails.

While reading `proportional adaptation rate` and `MDF / AUMDF` blocks in `mdf.c`, distinguish *"whether a partition undergoes time-domain projection on this frame"* from *"whether that partition receives a frequency-domain gradient update"*. Compare equal-length, equal-energy sparse vs. dense FIR paths and test abrupt shifts in the dominant reflection lag—without mislabeling Speex's partition-level proportional weighting as classic sample-wise PNLMS/IPNLMS.

### A04 RLS, Kalman, FDKF, and Partitioned State-Space AEC

**Keep estimation objectives distinct.** Real time-domain Recursive Least Squares (RLS) in [Eqs. (6-16)–(6-18)](../../../../chapters_en/06_aec.md#sec-6-2) propagates the inverse exponentially weighted reference correlation matrix via the Woodbury / matrix inversion lemma to solve regularized weighted least squares exactly—requiring finite initialization $\mathbf{P}_{-1} = \delta^{-1}\mathbf{I}$ to match the explicit penalty $\lambda^{n+1}\delta\|\vec{w} - \vec{w}_{-1}\|^2$ on early samples, as verified in [`aec_rls.py`](../../ch06/core/aec_rls.py), [`aec_rls_demo.py`](../../ch06/aec_rls_demo.py), and [`test_codes_aec_rls.py`](../../../../tests/test_codes_aec_rls.py).

Standard $L$-tap RLS stores an $L \times L$ inverse covariance state at $O(L^2)$ operations per sample, diverges under prolonged zero reference when $\lambda < 1$ divides $\mathbf{P}$ repeatedly without excitation, and absorbs near-end speech into filter taps during unfrozen double-talk (while square-root, QR, and Householder RLS variants improve conditioning at distinct operation counts; note [MathWorks `dsp.RLSFilter` docs](https://www.mathworks.com/help/dsp/ref/dsp.rlsfilter-system-object.html) document RLS variants while MathWorks' [AEC tutorial](https://www.mathworks.com/help/audio/ug/acoustic-echo-cancellation-aec.html) actually uses FDAF). Note also that our teaching `aec_rls.py` runs an explicit $O(L^3)$ Cholesky positive-definiteness check on every accepted sample for numerical safety—so its Python runtime must never represent optimized $O(L^2)$ C RLS, while [`aec_rls_kalman_comparison.py`](../../ch06/aec_rls_kalman_comparison.py) logs true prior echo prediction error separately from double-talk total residual.

**Fast RLS is distinct from standard frequency-domain NLMS.** Whereas PBFDAF normalizes by scalar bin power, [Schneider & Kellermann's GFDAF analysis (2016, §§3–4)](https://doi.org/10.1186/s13634-015-0302-2) connects generalized frequency-domain algorithms to *approximate block RLS*, and structured fast-RLS families—such as [Cioffi & Kailath's Fast Transversal Filters (1984)](https://doi.org/10.1109/TASSP.1984.1164334) and [Slock & Maouche's FFT/displacement-rank RLS (1994)](https://doi.org/10.1016/0165-1684(94)90018-3)—carry their own rescue/stability preconditions.

**Kalman AEC introduces an explicit state-space prior.** [Eqs. (6-20)–(6-22)](../../../../chapters_en/06_aec.md#sec-6-2) model the time-varying acoustic path as a latent state with process covariance $\mathbf{Q}$ (path drift uncertainty), measurement variance $\Psi$ (near-end speech + background noise power), and state-error covariance $\mathbf{P}$ (distinct from RLS's inverse regressor correlation, coinciding algebraically only under $\mathbf{A}=\mathbf{I}, \Psi=1, \mathbf{P}_n^- = \mathbf{P}_{n-1}/\lambda$, i.e., state-dependent $\mathbf{Q}_n = (\lambda^{-1}-1)\mathbf{P}_{n-1}$). Crucially, a spike in innovation power alone cannot distinguish **double-talk** (which demands *raising* $\Psi$ to slow adaptation) from **abrupt echo-path change** (which demands *raising* $\mathbf{Q}$ to speed re-convergence). Two standalone runnable demos verify exact rational arithmetic: [single-bin complex scalar Kalman](../../ch06/aec_kalman_scalar_demo.py) ([tests](../../../../tests/test_codes_aec_kalman_scalar_demo.py): $K = 3/7$, posterior $22/35$, variance $3/35$; $\Psi: 0.4 \to 4 \implies K = 3/16$; $X = i \implies K = -i/2$) and [2-tap real matrix Joseph-form Kalman](../../ch06/aec_kalman_matrix_demo.py) ([tests](../../../../tests/test_codes_aec_kalman_matrix.py): $\mathbf{P}^- = \operatorname{diag}(1, 4), \vec{x} = [1, 1]^\top, \Psi = 1 \implies \vec{K} = [1/6, 2/3]^\top$ with off-diagonal posterior $-2/3$, whereas forcing diagonal $\mathbf{P}$ distorts step-2 innovation variance from $9/2$ to $19/6$).

**From full-matrix Kalman to FDKF / PBFDKF.** [Enzner & Vary (*Signal Processing* 2006, §§2–4)](https://doi.org/10.1016/j.sigpro.2005.09.013) ([author note](https://homepages.ruhr-uni-bochum.de/gerald.enzner/StateSpaceFDAF.html)) formulate state-space FDAF as a diagonalized frequency-domain Kalman approximation unifying linear subtraction and uncertainty-driven post-filtering, extended to partitioned blocks in [Kuech, Mabande & Enzner (ICASSP 2014)](https://doi.org/10.1109/ICASSP.2014.6853806) and multi-channel/multi-reference covariance block structures (`VD-PBFDKF` vs. cross-speaker `SD-PBFDKF`, both dropping cross-partition blocks) in [Zhu et al. (Interspeech 2021, §§2–3, Figs. 4–5)](https://www.isca-archive.org/interspeech_2021/zhu21d_interspeech.pdf)—always specify whether "diagonal" means diagonal across bins, partitions, speaker channels, or scalar entries.

**Pinned upstream source entry points and verified gaps** (per [`SOURCES.lock.json`](../SOURCES.lock.json)):

| Upstream Source & Entry Point | Implemented State / Mechanics | Verified Gaps & Boundaries |
|---|---|---|
| [`pyroomacoustics` v0.10.0 `adaptive/rls.py`](https://github.com/LCAV/pyroomacoustics/blob/0dd39f2614b7fc44b2cc63dbe7d60f4641068890/pyroomacoustics/adaptive/rls.py) | MIT; general real scalar RLS and `BlockRLS`, fetched locally | `BlockRLS` updates once at block end; lacks bulk-delay search, DTD, and RES |
| [`pyaec` `time_domain_adaptive_filters/rls.py`](https://github.com/ewan-xu/pyaec/blob/5b9c02c57075d790b7df8652884618189d49bbc4/time_domain_adaptive_filters/rls.py) | Apache-2.0; time-domain Woodbury RLS + Kalman + `FDKF`/`PFDKF` demos | Fetched; `rls()`/`kalman()` drop final $N$ samples, while `fdkf()`/`PFDKF()` crash on removed `np.complex` (see A18) |
| [`MetaAF` `metaaf/optimizer_rls.py`](https://github.com/adobe-research/MetaAF/blob/56c4665bdc51c2e0595a7c0cd9b1266408adceff/metaaf/optimizer_rls.py) | JAX per-bin complex matrix RLS optimizer; core fetched | [`zoo/aec/`](https://github.com/adobe-research/MetaAF/tree/56c4665bdc51c2e0595a7c0cd9b1266408adceff/zoo/aec) includes Kalman/RLS baselines under separate non-commercial Adobe Research License (see A13; unexecuted) |
| [`echocatzh/PFDKF` `pfdkf.py`](https://github.com/echocatzh/PFDKF/blob/7c8c86b5691966c330015d8e0960db0733b4844f/pfdkf.py) | MIT; independent third-party partitioned frequency-domain Kalman demo | Fetched & 1 initial block executed; defaults to `res=True` post-filtering (see A18); not Kuech et al.'s official code |
| [`Subband_Kalman_AEC` `main.m`, `saf_kalman.m`](https://github.com/changxuding/Subband_Kalman_AEC/tree/f0c4f7030769d94dea2da3c814837448f171d422) | MIT MATLAB code; standard, square-root, information, and SRIF subband Kalman | Plain-text `.m` files verified; omits `.mat`/WAV assets, requires MATLAB, and includes optional post-filter |
| [`nay0648/bssaec2020` `SimulatedExperiment/`](https://github.com/nay0648/bssaec2020/tree/a3f52249ee61f19e2369823f65c6c31392dcf042/SimulatedExperiment) | Author MATLAB simulation scripts for Aux-ICA / weighted RLS | No explicit redistribution license (`index_only`); calls direct matrix inversion per sample rather than $O(L^2)$ Woodbury or C++ challenge code |

**Do not over-extrapolate industrial citations.** [Wang et al. (AEC Challenge 2021, §§2, 4.1, Tables 1–2)](https://arxiv.org/pdf/2102.08551) report an internal C++/SSE2 cascade (`GCC-PHAT delay -> 20 ms / 10 ms / 320-DFT 5-tap wRLS -> 1-frame-lookahead Deep-FSMN RES`) averaging `0.61` ms/frame (`0.19` ms delay+wRLS, `0.42` ms Deep-FSMN) on an Intel i5-8350U while noting occasional double-talk over-suppression—proving challenge execution without releasing proprietary C++ source or isolating wRLS-only scores. Similarly, [DSP Concepts Audio Weaver `SbKalmanAEC4RefV1` documentation](https://documentation.dspconcepts.com/awe-designer/8.D.2.7/sbkalmanaec4refv1) documents a commercial WOLA subband Kalman module (up to 4 references, exposing `resetP`, `masterReset`, `freeze`, `P`, `PXTerm`, with default residual noise suppression enabled) without publishing redistributable C source, while WebRTC AEC3 and SpeexDSP use linear AUMDF/partitioned adaptive filters rather than RLS/Kalman. Fair head-to-head benchmarking must separate synthetic known-component ERLE from real-recording total power changes, log linear vs. post-filtered taps separately, and never equate numerically equal $\lambda$ and $(Q, \Psi)$ as "matched hyperparameters."

<a id="apa-original-audit"></a>

**APA: Constraining Each Path Update Across $P$ Recent Regressors.** Originally proposed by [Ozeki & Umeda (*Electronics and Communications in Japan, Part I* `67(5):19–27`, 1984)](https://doi.org/10.1002/ecja.4400670503) (publisher metadata verified 2026-10-01 without full-text PDF access, so all derivations follow [Chapter 6 §6.2](../../../../chapters_en/06_aec.md#sec-6-2)), Affine Projection Algorithm (APA) stacks the $P$ most recent causal regressor vectors into a $P \times P$ Gram system—accelerating convergence on colored speech relative to $P=1$ NLMS without peeking at future samples or maintaining RLS's infinite-horizon $L \times L$ covariance inverse.

Pinned [`pyaec` `apa.py`](https://github.com/ewan-xu/pyaec/blob/5b9c02c57075d790b7df8652884618189d49bbc4/time_domain_adaptive_filters/apa.py) (`Apache-2.0` third-party teaching code; note `pyroomacoustics` v0.10.0 contains no APA implementation) hard-codes $0.01\mathbf{I}$ Gram loading with explicit `np.linalg.inv` and zeros buffers on every call without streaming state, freeze, DTD, or RES. Loaded directly in full without AST surgery or patching by [`audit_upstream_apa.py`](../../ch06/examples/audit_upstream_apa.py) ([2026-10-01 historical report](../../ch06/reports/upstream_apa.json); [2026-10-04 current report](../../ch06/reports/upstream_apa_current.json); exact rational adjugate verification in [`test_codes_ch06_upstream_apa.py`](../../../../tests/test_codes_ch06_upstream_apa.py)):

| Deterministic Input (`mu = 1` Throughout) | Actual Upstream Prior Residual Returned | Boundary Verified |
|---|---|---|
| 6-sample zero reference, `d = [1..6]`, $N = P = 2$ | `[1, 2, 3, 4]` | Truncates final $N=2$ samples (`len - N` loop bound) |
| `x = [1, 2, 3, 4, 5], d = 2x`, $N = P = 1$ | $[2, 4/101, 6/(101\times 401), 8/(101\times 401\times 901)]$ | Scalar NLMS limit with hard-coded `0.01` regularization |
| `x = [1, 1, 0, 0, 0], d = [1, 3, 2, 0, 0]`, $N = P = 2$ | $[1, 203/101, 40502/1040401]$ | $2\times 2$ Gram update using shared prior weights; final 2 samples truncated |
| 5-sample complex `x = i, d = 1`, $N = P = 1$ | `[1, 1, 1, 1]` + `ComplexWarning` | Real `np.zeros` buffers silently discard imaginary parts |

```bash
.venv/bin/python -m codes.chapters.ch06.examples.audit_upstream_apa --report codes/chapters/ch06/reports/upstream_apa_current.json
```

### A05 DTD: Binary Freeze Decisions vs. Continuous Step-Size Control

Double-Talk Detection (DTD) protects near-end speech from corrupting adaptive echo-path weights rather than synthesizing an echo replica directly. Whereas Geigel level ratios, cross-correlation, and coherence detectors output binary freeze gates, pinned SpeexDSP [`mdf.c`](https://gitlab.xiph.org/xiph/speexdsp/-/blob/8e29a256ef0235ebbe7fcb8417b5ac7731eb8307/libspeexdsp/mdf.c) explicitly eschews binary DTD in favor of continuous leakage-based learning-rate modulation ([Valin, *IEEE TASLP* 2007](https://people.xiph.org/~jm/papers/valin_taslp2006.pdf)), and WebRTC AEC3's `subtractor.cc` computes separate spectral update gains for its `refined` and `coarse` filters—neither can be reduced to a binary `if (dtd) freeze` flag.

**Executed: 4-state synthetic DTD confusion matrix and abrupt path-change counterexample.** In [Chapter 6 §6.4.1](../../../../chapters_en/06_aec.md#sec-6-4-1), [`double_talk.py`](../../ch06/core/double_talk.py) and [`aec_dtd_demo.py`](../../ch06/aec_dtd_demo.py) evaluate `16` kHz, `160`-sample frames across `20` silence, `40` far-end single-talk (including recovery), `20` near-end single-talk, and `20` double-talk frames: gating NLMS lowers double-talk true-echo residual RMS from `0.1492` (unfrozen) to `0.0812`, yet on an abrupt far-end-only path jump $[0.6, 0.25] \to [0.1, 0.75]$, **all `20/20` frames are misclassified as double-talk**—freezing adaptation right when re-convergence is needed.

**Stopping far-end playback does not immediately silence room echo tails.** When the current render frame goes silent while microphone energy persists, reverberant room decay from *earlier* render frames may still dominate. Indeed, pinned WebRTC AEC3's [`residual_echo_estimator.cc` line 307](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/aec3/residual_echo_estimator.cc#307) and exponential [`reverb_model.cc`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/aec3/reverb_model.cc) propagate historical render power across frames with positive decay—proving statically why current-frame render energy $= 0$ never implies zero acoustic echo tail.

### A06 AEC3 Bulk Delay, Linear Subtraction, Residual Suppression, and Comfort Noise

Read WebRTC's [`modules/audio_processing/aec3/`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/aec3/) in dataflow order: `echo_canceller3.cc -> block_processor.cc -> render_delay_buffer.cc / render_delay_controller.cc -> echo_path_delay_estimator.cc -> subtractor.cc -> residual_echo_estimator.cc -> suppression_gain.cc -> comfort_noise_generator.cc`. Internally, [`echo_path_delay_estimator.cc`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/aec3/echo_path_delay_estimator.cc) downmixes and downsamples capture audio to run time-domain matched filters against the render ring buffer before [`render_delay_controller.cc`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/aec3/render_delay_controller.cc) aligns buffer delay, while `subtractor.cc` runs dual `refined` + `coarse` partitioned filters at that aligned delay (not two delay hypotheses).

Applications enter through [`api/audio/audio_processing.h`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/api/audio/audio_processing.h#509) using `10` ms PCM frames (`ProcessReverseStream()` $\to$ `set_stream_delay_ms()` $\to$ `ProcessStream()`; interleaved `int16` strictly at `8, 16, 32, 48` kHz `NativeRate` with matching stream rates/layouts vs. deinterleaved `float` supporting `44.1` kHz and asymmetric layouts, distinct from internal `64`-sample subband blocks in `aec3_common.h`). `set_stream_delay_ms()` clamps inputs outside `[0, 500]` ms with `kBadStreamParameterWarning` in [`audio_processing_impl.cc`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/audio_processing_impl.cc).

**Playout-volume notifications do NOT rescale reference samples.** Calling `CreatePlayoutVolumeChange()` queues an integer level change in [`audio_processing_impl.cc` line 816](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/audio_processing_impl.cc#816) that sets boolean `gain_change`, forwarded by [`subtractor.cc` line 138](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/aec3/subtractor.cc#138) into [`refined_filter_update_gain.cc` line 52](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/aec3/refined_filter_update_gain.cc#52)—where the `gain_change` handler is an **empty `TODO`**, accompanied only by ERLE state resets in [`aec_state.cc` line 157](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/aec3/aec_state.cc#157).

Meanwhile, setting `Config::EchoCanceller::export_linear_aec_output = true` enables `GetLinearAecOutput()` to export the `16` kHz (`160`-sample / `10` ms) linear subtraction tap (`num_proc_channels()` channels in `audio_processing_impl.cc`), while [`block_processor.cc`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/aec3/block_processor.cc) bypasses capture processing before first render arrival and resets delay/buffer state on render underrun/overrun:

| Observed Symptom | First Telemetry to Inspect | Diagnostic Experiment & Layer to Fix |
|---|---|---|
| Heavy echo right from $t=0$ | Render-before-capture call order, `set_stream_delay_ms` warning/readback, `GetLinearAecOutput()` | Shift broadband reference by known fixed lags; fix reference tap point, HAL timestamps, or tail coverage *before* tuning RES |
| Gradual degradation over minutes | Sliding-window render-mic lag trajectory, hardware timestamps, FIFO level & XRUNs | Constant large offset = bulk delay; steady linear slope = clock SRO requiring shared clock or ASRC (`HasClockdrift()` does NOT resample PCM); step jump = dropped buffer |
| Quiet final output yet loud echo in `GetLinearAecOutput()` | Side-by-side `GetLinearAecOutput()` vs. `ProcessStream()` + near-end speech clipping | Aggressive RES is masking broken linear cancellation; fix reference linearity, delay, or adaptation gating first |

### A07 RES/NLP and Near-End Speech Protection

Residual Echo Suppression (RES) / Nonlinear Processing (NLP) applies spectral gain attenuation based on estimated uncancelled echo (`residual_echo_estimator.cc` + `suppression_gain.cc` in AEC3; `speex_preprocess_ctl(st, SPEEX_PREPROCESS_SET_ECHO_STATE, echo_st)` + `speex_echo_get_residual()` in Speex [`preprocess.c`](https://gitlab.xiph.org/xiph/speexdsp/-/blob/8e29a256ef0235ebbe7fcb8417b5ac7731eb8307/libspeexdsp/preprocess.c) per [Speex Manual §6.2](https://www.speex.org/docs/manual/speex-manual/node7.html)). Always archive linear canceller output separately from post-RES output, testing quiet near-end fricatives overlapping residual echo alongside [Exercises E06-04–E06-06](../../../../chapters_en/06_aec.md#sec-6-9) ([`exercises_enhancement.py`](../cross_chapter/exercises_enhancement.py), [`test_codes_enhancement_round2.py`](../../../../tests/test_codes_enhancement_round2.py), [`test_codes_exercises_enhancement.py`](../../../../tests/test_codes_exercises_enhancement.py): sign-inverted replica quadruples echo energy $+6$ dB; stationary mic noise caps raw I/O power ratio even under $100\%$ echo cancellation; cubic loudspeaker nonlinearity on `500` Hz generates `1500` Hz unmodellable by any linear FIR).

**Interface and state comparison across three implementations** (BSD-3-Clause notices for WebRTC and SpeexDSP in [`THIRD_PARTY_EN.md`](../THIRD_PARTY_EN.md); note [`xiph/rnnoise`](https://github.com/xiph/rnnoise) is single-input NS without playback reference, NOT an AEC):

| Implementation | Input Format & Persistent State | Bulk Delay, Double-Talk, and Residual Handling |
|---|---|---|
| Textbook `NLMSState` | Equal-length 1D float arrays; retains regressor delay line and FIR taps; external `freeze` flag | Zero internal delay search, DTD, RES, or SRC; exact arithmetic baseline under known delay and oracle freeze |
| WebRTC AEC3 via APM | `10` ms render/capture PCM frames; persistent APM, render delay buffer, and dual `refined`/`coarse` states | Caller sets `stream_delay_ms`, internal matched-filter lag tracker; `GetLinearAecOutput()` exports `16` kHz linear tap alongside post-RES final output |
| SpeexDSP AUMDF | `spx_int16_t` capture & playback frames (`10–20` ms frame, `100–500` ms tail per [`speex_echo.h`](https://gitlab.xiph.org/xiph/speexdsp/-/blob/8e29a256ef0235ebbe7fcb8417b5ac7731eb8307/include/speex/speex_echo.h)); `speex_echo_state_init_mc()` for multi-channel | Synchronous `speex_echo_cancellation()` adds zero extra playback buffer, whereas async `playback/capture` helpers copy only mono `frame_size` scalars ([`mdf.c` L687–733](https://gitlab.xiph.org/xiph/speexdsp/-/blob/8e29a256ef0235ebbe7fcb8417b5ac7731eb8307/libspeexdsp/mdf.c#L687-733)); continuous leakage step control; optional RES requires `SPEEX_PREPROCESS_SET_ECHO_STATE` |

Real paired recordings from [Microsoft AEC Challenge pinned `datasets/README.md`](https://github.com/microsoft/AEC-Challenge/blob/6c633d0a9d2a143a0e364899b91b06f127315b18/datasets/README.md) ([Dataset Licenses](https://github.com/microsoft/AEC-Challenge/blob/6c633d0a9d2a143a0e364899b91b06f127315b18/README.md#dataset-licenses) separate from code MIT) pair Windows loopback `*_lpb.wav` with microphone `*_mic.wav` without clean echo/near-end ground truth, governed by four protocol rules: (1) vary one condition per synthetic test (single-talk, double-talk, near-end only, silence, path step, bulk delay, dropped render frame, or two-clock SRO); (2) respect each library's exact frame/rate/gain contract; (3) save both linear and post-RES taps, reporting total I/O power ratios on real pairs distinctly from ground-truth ERLE; and (4) measure double-talk distortion, recovery time, CPU, RAM, and latency on aligned sample intervals without ranking mismatched output taps together.

**Executed: Synchronous SpeexDSP AEC on a real far-end single-talk pair (2026-09-23).** On GUID `-0AcvGNEdEK-DQGxWmtq2Q_farend_singletalk_{lpb,mic}.wav` from pinned commit `6c633d0a9d2a143a0e364899b91b06f127315b18` (`16` kHz mono PCM16, `188320` and `188480` samples, SHA-256 `9b204ad5473726526d14830103e53647897699ef89d49624964b7f2449040426` and `6b4c3e01b969c5cad91f248ff967cfa03df6554f3a060d6e5b06b7d20341bba6`, sliced to common `[0, 188320)` without shifting; diagnostic NCC peak on `[16000, 64000)` is `+498` samples $\approx 31.1$ ms, $r = 0.260$), floating-point `Release` `libspeexdsp.dylib` (`8e29a256ef0235ebbe7fcb8417b5ac7731eb8307`, macOS arm64 / AppleClang 21, dylib SHA-256 `c3e70172a3a9bf60b60bfd0bddfd58550a1899fef09078a9c7d5270d4a12105d`) runs `speex_echo_cancellation()` with `160`-sample frames, `4096`-sample tail (`256` ms), explicit `16` kHz rate set/readback, `[0, 48000)` warmup, and `[48000, 188320)` scoring window ($10\log_{10}(\sum d[n]^2 / \sum e[n]^2)$):

| Reference Fed to Identical Synchronous Speex Core | Scoring-Window Digital Power Change | Interpretation Boundary |
|---|---:|---|
| Matched `lpb` reference | `5.57` dB | Core AUMDF output (`1.62` dB below zero-ref output; 8 contiguous 1 s sub-blocks: `4.55, 4.99, 5.42, 6.41, 5.19, 5.79, 5.61, 6.38` dB); NOT clean-echo ground-truth ERLE |
| All-zero reference | `3.95` dB | Internal DC notch / pre-emphasis / de-emphasis alters power even with zero echo subtraction |
| Intentionally `1` s-delayed (`16000`-zero-padded) reference | `−0.23` dB | Misaligned reference slightly raises output energy on this clip |

Reproduce via [`aec_real_pair_experiment.py`](../../ch06/examples/aec_real_pair_experiment.py) (which rejects `131`-byte Git LFS pointer stubs and restricts optional `--output-wav` strictly outside the repository):

```bash
mkdir -p codes/chapters/ch00/upstream/_downloads/aec-challenge/datasets/real
curl -fL 'https://media.githubusercontent.com/media/microsoft/AEC-Challenge/6c633d0a9d2a143a0e364899b91b06f127315b18/datasets/real/-0AcvGNEdEK-DQGxWmtq2Q_farend_singletalk_lpb.wav' -o codes/chapters/ch00/upstream/_downloads/aec-challenge/datasets/real/-0AcvGNEdEK-DQGxWmtq2Q_farend_singletalk_lpb.wav
curl -fL 'https://media.githubusercontent.com/media/microsoft/AEC-Challenge/6c633d0a9d2a143a0e364899b91b06f127315b18/datasets/real/-0AcvGNEdEK-DQGxWmtq2Q_farend_singletalk_mic.wav' -o codes/chapters/ch00/upstream/_downloads/aec-challenge/datasets/real/-0AcvGNEdEK-DQGxWmtq2Q_farend_singletalk_mic.wav
cmake -S codes/chapters/ch00/upstream/_downloads/speexdsp -B /private/tmp/speexdsp-aec-20260923 -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DUSE_FIXED_POINT=OFF
cmake --build /private/tmp/speexdsp-aec-20260923 --parallel 4
.venv/bin/python codes/chapters/ch06/examples/aec_real_pair_experiment.py --speex-library /private/tmp/speexdsp-aec-20260923/libspeexdsp.dylib
```

**Second real-pair check on official double-talk (2026-09-23).** On GUID `-2jLGNCgf0WDpKMY2iup7g` `doubletalk_{lpb,mic}.wav` (SHA-256 `790203e68cc8cc9b02becc4efefbde0ddf962fff4d223608b00ffc8b2682ed71` and `4e7e35290d763d2d22173976c2cfd46739c9afe76729a0c283c1b69609f22041`; `196186` and `196320` samples truncated to `196160`; full-clip RMS `0.05574` and `0.02760` with `1` full-scale mic sample; `[0, 3)` s warmup, `[3, 8)` s scoring RMS `0.03860`), matched, all-zero, and `1` s-late references yield total I/O power changes of `2.234`, `2.204`, and `2.200` dB (`0.029` dB difference vs. zero-ref; [`aec_doubletalk_experiment.py`](../../ch06/aec_doubletalk_experiment.py)):

```bash
.venv/bin/python -m codes.chapters.ch06.aec_doubletalk_experiment --speex-library /private/tmp/speexdsp-aec-20260923/libspeexdsp.dylib
```

**Uncontrolled teaching PBFDAF comparison on the same two real pairs (2026-09-24).** Running [`aec_pbfdaf_real_pair_compare.py`](../../ch06/examples/aec_pbfdaf_real_pair_compare.py) with fixed synthetic-demo hyperparameters $N = 160, 2N = 320, P = 26, L = 4160$ (matching Speex's $\lceil 4096 / 160 \rceil = 26$ partition capacity), $\mu = 0.3, \delta = 0.05$, and **zero DTD/freeze/delay control** produces negative `float64` power ratios (meaning unfrozen step adaptation *increases* residual energy on unaligned/double-talk real audio, with `1331` and `11` samples exceeding PCM16 full scale on `1` s-late controls):

| Official Real Pair & Scoring Window | Reference Condition | Teaching PBFDAF `float64` Total Power Change |
|---|---|---:|
| Far-end single-talk `[48000, 188320)` | Matched pair | `−6.700` dB |
| Same | All-zero reference | `0.000` dB |
| Same | Intentionally `1` s late | `−18.017` dB |
| Double-talk `[48000, 128000)` | Matched pair | `−6.363` dB |
| Same | All-zero reference | `0.000` dB |
| Same | Intentionally `1` s late | `−1.665` dB |

```bash
.venv/bin/python -m codes.chapters.ch06.examples.aec_pbfdaf_real_pair_compare --pair both
```

**Controlled near-end injection checks: strictly synthetic vs. semi-synthetic (2026-09-23).** Using [Eq. (6-36)](../../../../chapters_en/06_aec.md#sec-6-8) to compare two fresh runs differing on microphone input by a known injected near-end signal $s$ (`output_delta` $\Delta$, projection gain $g_\Delta$, and unaligned relative squared error $E_\Delta$), the synthetic harness loads `aec_far.wav`, `aec_near.wav` (non-zero on `[19200, 28800)`), and `aec_microphone.wav` ($d_0 = d_1 - s$ exact in PCM16; `32`-tap NLMS $\mu = 0.4, \varepsilon = 10^{-8}$):

| Stimulus & Implementation | Playback Reference | $g_{\Delta}$ | $E_{\Delta}$ | Diagnostic Meaning |
|---|---|---:|---:|---|
| Synthetic; NLMS with oracle freeze | Matched | `1.000` | `0` | Identical frozen weights across both arms reproduce injected $s$ sample-for-sample |
| Synthetic; NLMS unfrozen | Matched | `0.936` | `0.322` | Near-end $s$ perturbs adaptive weights between the two arms |
| Synthetic; Speex synchronous core | Matched | `0.878` | `0.208` | Reflects DC/pre-emphasis filtering, continuous adaptation shift, and PCM16 quantization |
| Synthetic; Speex synchronous core | All-zero | `0.879` | `0.204` | Near-identical deviation under zero reference proves pre/de-emphasis contributes |

The semi-synthetic harness superimposes `[48000, 128000)` of `nearend_singletalk_mic.wav` (SHA-256 `192291e973b7dee519102f30ae75e4209712ce55d02f55a9020a7c3b37a1d0e0`, peak sum `13633` PCM16 counts, unclipped) onto the far-end single-talk mic (`[0, 3)` s warmup, `[3, 8)` s injection scoring, `[8, 11)` s post-injection recovery; [`aec_controlled_doubletalk.py`](../../ch06/aec_controlled_doubletalk.py) and [`test_codes_aec_controlled_doubletalk.py`](../../../../tests/test_codes_aec_controlled_doubletalk.py)):

| Semi-Synthetic Stimulus, Speex Synchronous Core | $g_{\Delta}$ | $E_{\Delta}$ | `[8, 11)` s Post-Injection Output Diff RMS |
|---|---:|---:|---:|
| Matched reference | `0.922` | `0.140` | `111.04` PCM16 counts |
| All-zero reference | `0.923` | `0.128` | `5.62` PCM16 counts |

```bash
.venv/bin/python -m codes.chapters.ch06.aec_controlled_doubletalk \
  --speex-library /private/tmp/speexdsp-aec-20260923/libspeexdsp.dylib \
  --include-real-hybrid
```

**WebRTC AEC3 local build & offline execution (2026-09-23).** Syncing full WebRTC commit `0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e` via `depot_tools` (`db1dc923aa3c34fa748015e8c7a435d955b4c105`) inside git-ignored `webrtc_aec3_checkout/` (halting `runhooks` once compiler/symbol tools completed before multi-GB video test-asset downloads) and building `rtc_tools:audioproc_f` ([`rtc_tools/BUILD.gn`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/rtc_tools/BUILD.gn)) under Xcode 27.0 (`27A266a`) with host `libc++` and macOS `15.4` SDK (`use_custom_libcxx=false mac_sdk_path="/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk"`, bypassing bundled `libc++` `INFINITY` and `ld64.lld` `arm64e.x1` `.tbd` incompatibilities; binary SHA-256 `82f8eebaf574eb10fdee695ee7ea05a54418300a045a3ffb0de8f4c78f4592d3`):

```bash
PATH=/path/to/depot_tools:$PATH DEPOT_TOOLS_UPDATE=0 \
  buildtools/mac/gn gen out/aec3 \
  --args='is_debug=false rtc_include_tests=true rtc_enable_protobuf=true target_cpu="arm64" use_custom_libcxx=false mac_sdk_path="/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk"'
PATH=/path/to/depot_tools:$PATH DEPOT_TOOLS_UPDATE=0 \
  autoninja -C out/aec3 -j 4 rtc_tools:audioproc_f
```

Driven by [`aec3_offline_compare.py`](../../ch06/examples/aec3_offline_compare.py) with `--custom_call_order_file` set to `rc` (overriding [`wav_based_simulator.cc`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/modules/audio_processing/test/wav_based_simulator.cc)'s default `cr` order) and `--fixed_interface=true --aec=1 --agc=0 --agc2=0 --ns=0 --hpf=0 --ts=0 --stream_delay=0`:

```bash
.venv/bin/python -m codes.chapters.ch06.examples.aec3_offline_compare \
  --audioproc codes/chapters/ch00/upstream/_downloads/webrtc_aec3_checkout/src/out/aec3/audioproc_f \
  --output-dir /private/tmp/aec3-comparison-new-run
```

| Real Paired Recording & Scoring Interval | AEC3 Linear Output (`--linear_aec_output`) | AEC3 Final Output (`--o`) |
|---|---:|---:|
| Far-end single-talk, `[48000, 188320)` (`[3, 11.77)` s) | `7.685` dB | `15.731` dB |
| Double-talk, `[48000, 128000)` (`[3, 8)` s) | `3.384` dB | `3.851` dB |

(Output WAV SHA-256s: linear `6022c2a5354c9f760a3cac54e8e8ac9c63c1aec7bd16f3f2275f7257f27abc4e` & `9ca8a7e81bc7c07e8b4cf337832de20ff0cea5b94e3e937d7bbcb870baed8fab`; final `cece4d1134285b5fe251a999d5b8915642753b1be19e63a72882d1e2e624944d` & `adcdd3b2d98181d0a9ae4c7d9a85181aa39c9bc2a50d4b51acb9a428df60c5b4`—total power ratios only, NOT ground-truth ERLE.)

**Executed: Same-input synthetic ground-truth PCM interface comparison (2026-09-24).** [`aec_same_input_truth.py`](../../ch06/examples/aec_same_input_truth.py) feeds a `16` kHz, `6` s noiseless stimulus (seed `20260924`, 5-tap path `[0.70, 0, -0.25, 0, 0.12]`, exact PCM near-end injection on `[3, 4)` s; single-talk scored on `[2, 3)` s, injection on `[3, 4)` s; supports `--build-manifest` schema v1 for custom binaries) across all four output taps after aligning each tap's independently impulse-measured bulk output delay (without delay alignment, AEC3 linear $g_\Delta$ drops to `0.199`):

| Output Tap | Impulse-Measured Fixed Output Delay (samples) | Aligned Echo-In / Total-Out Power Ratio (dB) | Injected Increment Gain $g_\Delta$ | Increment Error $E_\Delta$ |
|---|---:|---:|---:|---:|
| Textbook NLMS linear residual | 0 | 79.27 | 1.000 | `≈ 0` |
| Speex synchronous core output | 0 | 37.12 | 0.982 | 0.027 |
| AEC3 exported linear tap | 64 | 60.58 | 0.964 | 0.241 |
| AEC3 final post-processed output | 128 | 30.63 | 0.193 | 0.763 |

(Read-only macOS CoreAudio enumeration and the short probe logged in [`AEC_DEVICE_PROBE_2026-09-23.md`](../../ch06/reviews/AEC_DEVICE_PROBE_2026-09-23.md) confirmed built-in mic ID `71` and speaker ID `76` share `main` / `bltn` clock domain without loopback endpoints or verified acoustic echo coupling, so live hardware AEC remains marked untested.)

### A08 PNLMS/IPNLMS, Subband AEC, and Nonlinear Echo-Path Models

Proportional step allocation (PNLMS/IPNLMS), subband decomposition, and Volterra/Hammerstein/Wiener cascades address distinct bottlenecks—sparse impulse responses, colored regressor conditioning, and loudspeaker/amplifier nonlinearities ([§§6.2–6.3](../../../../chapters_en/06_aec.md#sec-6-2) & [§6.3](../../../../chapters_en/06_aec.md#sec-6-3)).

In this repository, [`aec_ipnlms.py`](../../ch06/core/aec_ipnlms.py) implements real IPNLMS per Eq. (6-11) (matching NLMS step-for-step at $\kappa = -1$ when NLMS regularization equals $L\delta$, and exposing the zero-initialization $\kappa = +1$ zero-share boundary in [`test_codes_aec_ipnlms_subband.py`](../../../../tests/test_codes_aec_ipnlms_subband.py)); [`aec_subband.py`](../../ch06/core/aec_subband.py), [`aec_crossband_demo.py`](../../ch06/aec_crossband_demo.py), and [`test_codes_aec_crossband.py`](../../../../tests/test_codes_aec_crossband.py) contrast exact 2-band Haar cross-band matrices, diagonal in-band NLMS, and full $2\times 2$ cross-band adaptive NLMS per Eq. (6-10) alongside [`aec_ipnlms_subband_demo.py`](../../ch06/aec_ipnlms_subband_demo.py) and [`aec_advanced_exercises.py`](../../ch06/aec_advanced_exercises.py) (`E06-08–E06-10`). Contrast also pinned `pyroomacoustics` [`adaptive/subband_lms.py`](https://github.com/LCAV/pyroomacoustics/blob/0dd39f2614b7fc44b2cc63dbe7d60f4641068890/pyroomacoustics/adaptive/subband_lms.py), critical-sampling cross-filters in [Gilloire & Vetterli (1992)](https://doi.org/10.1109/78.149989), and full-band tap updates in [Lee & Gan's NSAF (2004)](https://doi.org/10.1109/LSP.2004.833445), while keeping Figure 32's [four analytical 2-band Haar cross-term WAVs](05_exercises_and_audio.md#sec-u-b41d17b2cd) and [seven short-path linear-AEC WAVs](05_exercises_and_audio.md#sec-u-042e52c4d2) strictly separate from real-device claims ([`AEC_DEVICE_PROBE_2026-09-23.md`](../../ch06/reviews/AEC_DEVICE_PROBE_2026-09-23.md)).

### A09 Multi-Reference Stereophonic AEC and Non-Identifiability

Multi-reference AEC (`speex_echo_state_init_mc(frame_size, filter_length, nb_mic, nb_speakers)`—never swap `nb_mic` and `nb_speakers`) sums adaptive paths across $K$ loudspeaker channels: whenever stereo playback channels are collinear or strongly correlated, total echo residual may vanish on the current pan position yet leave individual speaker-to-mic paths non-identifiable, causing abrupt echo bursts the instant stereo panning shifts ([§6.6](../../../../chapters_en/06_aec.md#sec-6-6)).

### A10 DTLN-aec: Dual-Input Two-Stage Representation and Streaming State

[`breizhn/DTLN-aec`](https://github.com/breizhn/DTLN-aec/tree/9d24e128b4f409db18227b8babb343016625921f) (`MIT` code; TFLite checkpoints/audio licensed separately; [ICASSP 2021 paper](https://doi.org/10.1109/ICASSP39728.2021.9413510)) feeds `16` kHz mono `mic.wav` and `lpb.wav` (`512`-sample window, `128`-sample hop, `384`-sample zero-pad head/tail sliced back to common length) through an STFT magnitude-masking LSTM stage followed by a learned 1D-conv feature-domain masking stage—preserving LSTM states across hops *within* `process_file()` while resetting at each file boundary (and clamping only positive peaks `if np.max(out) > 1`, as verified via our mock-interpreter AST test in A18 without loading TensorFlow/TFLite weights).

### A11 NKF-AEC: Linear AEC with Neural Kalman Gain Estimation

[`fjiang9/NKF-AEC`](https://github.com/fjiang9/NKF-AEC/tree/8ac58fb8fb9ced48579f9aa310745c54f98d7e1f) (`src/nkf.py`, [paper](https://arxiv.org/abs/2207.11388)) uses a complex RNN to predict per-bin multi-tap Kalman gains inside a strictly *linear* `16` kHz AEC assuming pre-compensated bulk delay; lacking an explicit redistribution license file, it is retained as `index_only` without copying code or checkpoints, and kept distinct from ASRU 2023 `NeuralKalman`.

### A12 NeuralKalman, Deep Adaptive AEC, and DeepVQE

[NeuralKalman (ASRU 2023)](https://doi.org/10.1109/ASRU57964.2023.10389780), [Deep Adaptive AEC (ICASSP 2022)](https://minjekim.com/wp-content/uploads/icassp2022_hzhang.pdf), and [DeepVQE (Interspeech 2023)](https://www.isca-archive.org/interspeech_2023/ristea23_interspeech.html) inject neural networks at fundamentally different points in the AEC loop (with no license-verified official author repositories confirmed as of our audit; third-party reimplementations in [Seidel et al. *IEEE SPM* `41(6):24–38`, Nov. 2024](https://doi.org/10.1109/MSP.2024.3449557) / [author v1 §III & Table I](https://arxiv.org/html/2501.16367v1#S3) clarify the architectural taxonomy without standing in as original author releases):

| Method | Exact Modification Relative to Linear Baseline | Output Tap & Structural Assumptions | Minimal Ablation & Source Status |
|---|---|---|---|
| `NKF` | Neural network predicts multi-tap Kalman gain; raw render spectrum & linear state recursion retained | Linear subtraction residual; OLA, multi-tap, $A = 1$, pre-aligned reference | Sweep path jump & double-talk gain response; A11 author repo `index_only` (unlicensed) |
| `NeuralKalman` | Learns reference nonlinearity, state update/transition, while retaining analytic gain formula | OLS echo-cancelled output (1-tap in survey comparison); domain-dependent reference transform | Ablate learned reference nonlinearity vs. learned state transition; survey code $\ne$ ASRU author release |
| `Deep Adaptive AEC` | Learns reference feature transform & update gain around multi-tap linear filter state | Survey evaluates NR-free variant (`OLA`, $A = 1$) | Swap learned gain vs. learned reference transform independently; no verified author release |
| `DTLN-aec` | End-to-end two-stage mask estimation (STFT + learned conv basis) without explicit linear FIR taps | Final neural output; `16` kHz mono, stateful across `128`-sample hops | A10 MIT code fetched; mock AST probe $\ne$ TFLite inference |
| `Haubner CTF-AEC` | DNN predicts frequency-selective step-size & error-normalization masks for linear CTF filter | Linear CTF residual; `center_stft=True` batch-reset training script | A14 `BSD-4-Clause` author code fetched without checkpoints |
| `DeepVQE` | Joint neural attention alignment + echo, noise, and reverberation suppression | Joint post-processed output (no separable linear echo tap) | Paper architecture index only; avoid partial NS-only third-party clones |

### A13 Meta-AF: Meta-Learned Adaptive Optimizers and Directory-Scoped Licensing

[`adobe-research/MetaAF`](https://github.com/adobe-research/MetaAF/tree/56c4665bdc51c2e0595a7c0cd9b1266408adceff) ([paper](https://arxiv.org/abs/2204.11942)) meta-learns GRU/FGRU complex frequency-domain filter optimizers (`metaaf/filter.py`, `core.py`, `optimizer_gru.py`, `optimizer_fgru.py`, `meta.py` under `UIUC/NCSA` open-source license, paired with task recipes in `zoo/aec/aec.py`, `optimizer_kf.py`, `aec_baselines.py`, `aec_eval.py`, `zoo/__config__.py`, `zoo/metrics.py` under [`zoo/LICENSE`](https://github.com/adobe-research/MetaAF/blob/56c4665bdc51c2e0595a7c0cd9b1266408adceff/zoo/LICENSE) §§1.1 & 1.3 permitting non-commercial academic/teaching research redistribution with original notices—fetched into git-ignored `_downloads/` without JAX/Haiku/SpeexDSP-Python execution or checkpoints).

### A14 DNN Step-Size Control for Linear CTF-AEC

[Haubner et al.'s official repository](https://github.com/ThomasHaubner/e2e_dnn_ad_control_for_lin_aec/tree/7a003133d742698de7acba9510d9586d7d57a584) (*IEEE/ACM TASLP* Vol. 32, pp. 227–238, **2024**, [`10.1109/TASLP.2023.3325923`](https://doi.org/10.1109/TASLP.2023.3325923), `BSD-4-Clause`, zero checkpoints included, expecting user `train_data.h5`/`test_data.h5`) traces from `main_train.py` and [`libPython/class_frontend.py`](https://github.com/ThomasHaubner/e2e_dnn_ad_control_for_lin_aec/blob/7a003133d742698de7acba9510d9586d7d57a584/libPython/class_frontend.py) into `update_filter()` in [`class_aec_ctf.py`](https://github.com/ThomasHaubner/e2e_dnn_ad_control_for_lin_aec/blob/7a003133d742698de7acba9510d9586d7d57a584/libPython/class_aec_ctf.py)—retaining an explicit linear Convolutive Transfer Function (CTF) echo filter modulated by DNN step-size and error-normalization masks, while setting `center_stft=True` and re-initializing filter states per sequence in `forward_batch()`.

### A15 Integrated Multi-Channel AEC/NR and Cascade Ordering

[`Arnout-Roebben/Integrated_AEC_NR`](https://github.com/Arnout-Roebben/Integrated_AEC_NR/tree/23c6b567c7863a8ee9bafd38bad0d3ff2f25e185) (`MIT` MATLAB R2024a code, [*IEEE TASLPRO* 2025](https://doi.org/10.1109/TASLPRO.2025.3648802), [`ReadMe.md`](https://github.com/Arnout-Roebben/Integrated_AEC_NR/blob/23c6b567c7863a8ee9bafd38bad0d3ff2f25e185/ReadMe.md)) compares `MWF`, `MWFext`, `AEC -> NR`, `NR -> AEC`, and `NRext -> AEC -> PF` inside [`Util/Process/process.m`](https://github.com/Arnout-Roebben/Integrated_AEC_NR/blob/23c6b567c7863a8ee9bafd38bad0d3ff2f25e185/Util/Process/process.m)—noting critically that `process.m` receives pre-separated oracle signal components $s, n, e_s, e_n, l_s, l_n$ and forms **oracle VAD masks directly from clean speech `s_f` and clean echo `es_f`** (and excluding `Audio/sig.mat` speech data from MIT code licensing).

### A16 Controllers, Residual Suppressors, and Nonlinearity Evidence in Original Papers

Audited directly from primary PDFs on 2026-09-28 without cross-paper ranking:
- **PercepNet AEC** ([Valin et al., 2021 §§3–5](https://arxiv.org/pdf/2102.05245)): pairs Speex-derived `150` ms `16` kHz AUMDF (`8` kHz / `400` ms coarse delay estimator; peak-energy partition always projected, remaining partitions alternated) with a 32-band ERB gain + comb-filter post-processor (`2` lookahead frames = `40` ms algorithmic latency) whose neural inputs concatenate mic, linear residual, and **unaligned** far-end $f(n)$ rather than delayed $f(n-D)$.
- **DVSS Learned Step-Size NLMS** ([Ivry, Cohen & Berdugo, §§3, 4.3, Table 4](https://israelcohen.com/wp-content/uploads/2022/04/Deep_Adaptation_Control_for_Acoustic_Echo_Cancellation.pdf)): predicts continuous NLMS step sizes every `8` ms (`16` ms STFT, `96` ms past context, `150` ms FIR, `80/10/10` h at `16` kHz, path jump at `4.5–5.5` s), reporting `3.4` s / `95%` convergence (`<-10` dB misalignment sustained) vs. `7.9` s / `58%` for DTD-free constant-$\mu = 3\times 10^{-5}$ NLMS.
- **Nonlinear Volterra / Cascade & Amplifier Clipping**: [Stenger & Rabenstein (NSIP 1999, §§2–3.1)](https://www.eurasip.org/Proceedings/Ext/NSIP99/Nsip99/papers/146.pdf) separate memoryless clipping from symmetric 2nd-order Volterra kernels (cf. `svf()` factored 2nd-order cascade in [`pyaec` `volterra.py`](https://github.com/ewan-xu/pyaec/blob/5b9c02c57075d790b7df8652884618189d49bbc4/nonlinear_adaptive_filters/volterra.py), which likewise truncates `M` tail samples); [Nollett & Jones (NSIP 1997, §6 vs. §5)](https://www.iwaenc.org/proceedings/1997/nsip97/pdf/author/ns970521.pdf) measure `19.5` dB ERLE (`30+200`-tap LNL cascade) vs. `11.1` dB (`230`-tap NLMS, `+8.4` dB gain) on hard-clipped synthetic vowels through an inexpensive PC speaker; and [Shah et al. (2015)](https://doi.org/10.1109/TASLP.2015.2425955) ([institutional citation](https://scholarsmine.mst.edu/ele_comeng_facwork/6620/)) indexes speaker voltage/current sensing without quoting unverified abstract dB figures.

### A17 Industrial Source Code: Which Architectural Layer Owns Default Constants?

- **WebRTC AEC3 internal blocks vs. APM frames**: [`echo_canceller3_config.h`](https://webrtc.googlesource.com/src/+/0467d2b91cc20b9b001c2bbb73d43ea6b2491f3e/api/audio/echo_canceller3_config.h) defaults both `refined` and `coarse` filters to `13` subband blocks (`12` initial blocks over `2.5` s) = $13 \times 64 = 832$ samples (`52` ms linear tail at `16` kHz—far shorter than APM's `500` ms `set_stream_delay_ms` clamp), plus `5`-block default delay, `4×` decimation, `5` matched filters, and `32`-sample headroom.
- **SpeexDSP partition ceiling & DC/pre-emphasis**: [`mdf.c`](https://gitlab.xiph.org/xiph/speexdsp/-/blob/8e29a256ef0235ebbe7fcb8417b5ac7731eb8307/libspeexdsp/mdf.c) rounds `4096 / 160` up to `26` partitions (`4160` taps), rotates projection FFTs rather than skipping gradient updates, and runs DC notch + pre-/de-emphasis even on zero reference.
- **SOF wrapper vs. external AEC library**: `src/audio/google/` in SOF (`BSD-3-Clause`; `Kconfig` defaults: `48` kHz, $\le 2$ mic / $\le 2$ ref, `200` KB arena, `44` ms nominal delay, headroom `4`) ships `google_rtc_audio_processing_mock.c`—which literally *adds* mic + reference to test DMA plumbing—while real builds link an external `third_party/lib/google_rtc_audio_processing` binary not included in SOF source.
- **XMOS `lib_voice` fixed-point BFP AEC**: [`lib_voice`](https://github.com/xmos/lib_voice/tree/c9f1a9bf95cd88c7950adf4bf631c217f900ad25) (`v1.1.0`, XMOS Public Licence v1, `lib_xcore_math` `v3.0.0`) specifies `240`-sample frame advance, `512`-point window (`257` real FFT bins), and `10` main / `5` shadow partitions (with `16` unused taps per partition) in `lib_voice/api/aec/aec_defines.h`.

### A18 Actual Diagnostics and Preserved Failures on Pinned Upstream AEC Interfaces

Executed under Python `3.13.12` and NumPy `2.5.3` on 2026-09-28 ([`audit_aec_upstream_interfaces.py`](../../ch06/examples/audit_aec_upstream_interfaces.py), [historical report](../../ch06/reports/aec_upstream_interfaces.json), [offline unit tests](../../../../tests/test_codes_aec_upstream_interfaces.py)) and re-verified on 2026-10-04 into the [current report](../../ch06/reports/aec_upstream_interfaces_current.json) ([source-contract tests](../../../../tests/test_codes_ch06_source_contracts.py), preserving historical bytes bound to commit `621d727a63475e3ee8ad2a29d518889c42e92571`):

1. **`pyaec` original functions**: On 6-sample zero reference, `d = [1..6]`, $N = 2$, both `rls()` and `kalman()` return `[1, 2, 3, 4]`—truncating the final $N=2$ samples via `nIters = min(...) - N` (and initializing `P = delta * I` rather than $\delta^{-1}\mathbf{I}$), while `fdkf()` and `PFDKF()` raise `AttributeError` on removed `np.complex`.
2. **`echocatzh` original `PFDKF`**: With $N = 1, M = 2$, zero weights, and `x = d = [1, 1]` calling `filt()` *without* `update()`, `res=False` returns error `[1, 1]` and echo `[0, 0]`, whereas default `res=True` returns error `≈ [0.5, 0.5]` and second output `≈ [0.5, 0.5]` despite **all-zero linear weights**—proving the second return value is `d - final_output` rather than linear echo.
3. **`DTLN-aec` AST `process_file()` with mock interpreter**: On `768`-sample zero mic, `640`-sample zero ref, stage-1 unit mask, and stage-2 constant `±0.5` blocks (`4×` overlap $\to \pm 2.0$), `if np.max(out) > 1` scales `+2.0` down to `+0.99` while leaving **`-2.0` unscaled**.

```bash
.venv/bin/python codes/chapters/ch06/examples/audit_aec_upstream_interfaces.py
.venv/bin/python -m codes.chapters.ch06.examples.audit_aec_upstream_interfaces --report codes/chapters/ch06/reports/aec_upstream_interfaces_current.json
```

### A19 2024–2026 Candidates, Evaluation Standards, and Inclusion Boundaries

Audited on 2026-09-28, 2026-10-01, and 2026-10-04:
- **Adopted with verified source**: **Haubner 2024 CTF-AEC** (A14) and **Seidel et al. Nov. 2024 *IEEE SPM* Neural Kalman Survey** ([`10.1109/MSP.2024.3449557`](https://doi.org/10.1109/MSP.2024.3449557) / [arXiv v1](https://arxiv.org/html/2501.16367v1)).
- **Align-ULCNet** ([EUSIPCO 2025, pp. 406–410, §II](https://eusipco2025.org/wp-content/uploads/pdfs/0000406.pdf)): dual-encodes and soft-aligns Kalman residual + playback reference before masking Kalman residual (no license-verified AEC author repo; generic single-channel ULCNet repo does not substitute).
- **E2E-AEC** ([2026-01 preprint v1 §2](https://arxiv.org/html/2601.16774v1), confirmed in [official ICASSP 2026 program `PaperNum=5855`, `AASP-L9.2` on 2026-05-07](https://www.cmsworkshops.com/ICASSP2026/view_paper.php?PaperNum=5855&bare=1)): `20` ms / `10` ms unidirectional GRU trained stage-wise using offline GCC-PHAT pseudo-labels and optional VAD output gating (no official checkpoints/configs released).
- **MSA-EchoLite** ([2026-08-04 preprint §II](https://arxiv.org/html/2608.03650v1)): feeds PBFDAF echo estimate + mic into Bark-compressed frequency/time blocks (unverified publication/code).
- **Device-Calibrated Stereo Path Priors** ([Wang et al., Interspeech 2025, pp. 773–777, `10.21437/Interspeech.2025-1572`](https://www.isca-archive.org/interspeech_2025/wang25o_interspeech.html), [§§2–4 PDF](https://www.isca-archive.org/interspeech_2025/wang25o_interspeech.pdf)): initializes `VD-MCSSFDAF` Eq. (17) priors and Eq. (19) process noise from pre-calibrated speaker-to-mic paths on fixed-geometry devices (connects directly to [Exercise E06-13](../../../../chapters_en/06_aec.md#e06-13); no official code).
- **LMPAN** ([Interspeech 2026, pp. 4945–4949, `10.21437/Interspeech.2026-191`](https://www.isca-archive.org/interspeech_2026/liu26b_interspeech.html); [author v1 §§2–3](https://arxiv.org/html/2607.02062v1)): pairwise soft-aligns `(mic, ref, NLMS_residual)` into GTCRN with training-only WavLM loss (`32` ms / `16` ms framing vs. `"100 lags = 1 s"` text discrepancy; [`alibaba/unified-audio`](https://github.com/alibaba/unified-audio) still lists AEC as `Developing`).
- **DiffVQE** ([Interspeech 2026, pp. 4971–4976, `10.21437/Interspeech.2026-2337`](https://www.isca-archive.org/interspeech_2026/lugo26_interspeech.html); [author v2 §§2–4](https://arxiv.org/pdf/2605.08189v2)): 1-step score correction atop a discriminative predictor, yet non-causal due to `512/128` bidirectional GRUs and offline GCC-PHAT ([demo page](https://ifnspaml.github.io/DiffVQE-Demo/) / [repo](https://github.com/ifnspaml/DiffVQE-Demo) ship audio samples only).
- **Evaluation Tools & ITU-T Standards**: [`ifnspaml/EC-Evaluation-Toolbox`](https://github.com/ifnspaml/EC-Evaluation-Toolbox/tree/aec8873325ed8e4d93eadc53e5c0af84be4db4fd) lacks a license (`index_only`, missing `speechlightning`); [ITU-T P.862 (PESQ)](https://www.itu.int/rec/T-REC-P.862/en) is marked **Superseded/Deleted (2024-01-05)** whereas [ITU-T P.863 (POLQA)](https://www.itu.int/rec/T-REC-P.863/en) (`2018` + `2020` amendment) remains in force under commercial licensing; and [Microsoft AEC Challenge](https://github.com/microsoft/AEC-Challenge) remains at its ICASSP 2023 4th edition (preserving our 5 locally modified audio worktree files untouched).

<a id="wpe"></a>

## 3. WPE and Joint Convolutive Filtering

### W01 Offline WPE and MIMO-WPE

Variance-normalized delayed linear prediction uses delayed multi-channel STFT history to predict and subtract current-frame late reverberation while time-varying inverse-power weights prevent loud speech onsets from dominating the regression. In multi-channel MIMO-WPE, the per-bin regressor dimension equals $\text{channels} \times \text{taps}$, becoming ill-conditioned under insufficient frames or collinear channels. Our teaching baseline lives in [`dereverberation.py`](../../ch07/core/dereverberation.py), paired with upstream [`nara_wpe` `wpe.py`](https://github.com/fgnt/nara_wpe/blob/a166779cca2088817e330481bd20af1a2c598555/nara_wpe/wpe.py) (`build_y_tilde()`, `wpe_v6()`, etc.; [Yoshioka & Nakatani, *IEEE TASL* 2012](https://doi.org/10.1109/TASL.2012.2210879)).

Trace lag ordering, prediction delay $\Delta$, `"valid"` vs. `"full"` covariance support, complex conjugation, and power floors—verifying single-tap closed-form normal equations first, then duplicating a mono channel into two identical channels, and injecting ultra-short clips, silent bins, and singular SCMs.

Executed locally on macOS arm64 (Python `3.13.12`, NumPy `2.5.3`, `complex128`) on 2026-09-22 and re-verified on 2026-09-23, [`compare_wpe_reference.py`](../../ch07/examples/compare_wpe_reference.py) compares unmodified `nara_wpe` `0.0.11` `wpe_v6()` against `offline_wpe()` (`seed=20260922`, shape `(4, 2, 192)`, `2` taps, `delay=3`, `statistics_mode="valid"`, zero power smoothing, zero diagonal loading; 2-mic, 1-mic, channel-swapped, and scaled inputs across `1` and `3` iterations): excluding the first `4` startup frames (bypassed in our baseline vs. zero-padded in upstream), **max absolute complex difference is `3.17e-15`** and **max relative $L_2$ difference is `5.76e-16`**.

Meanwhile, Exercise E07-04 proves offline non-causality on two 4-frame sequences sharing an identical prefix yet differing on the final frame—shifting frame-1 output from `0` to `-1/9` because batch covariance pools future frames—and [Exercise E07-06](../../../../chapters_en/07_wpe-dereverberation.md#e07-06) derives from spherical Gaussian NLL ($M\log\lambda$) why shared variance equals the cross-channel mean squared magnitude (`[1+j, 2] -> 3`, doubling channel 1 amplitude $\to 9$; cf. diagonal vs. scaled-identity models in §§V-A/B and tail-energy caveats in §VI of [Yoshioka & Nakatani 2012 full text](https://www.audiolabs-erlangen.de/resources/aps-w23/papers/sap_Yoshioka2012.pdf)).

**Pinned `0.0.11` statistical conventions.** `build_y_tilde()` stacks taps outer and channels inner; `statistics_mode="valid"` accumulates from `delay + taps - 1` while predicting all frames; default `wpe` aliases `wpe_v7`, `_stable_solve()` falls back to `lstsq` without relative trace loading, and `get_power_inverse()` floors power at `1e-10 * max(power)` across **the entire passed array**—meaning batching all bins in `wpe_v6`/`v7` vs. looping bin-by-bin in `wpe_v8` yields different weights on high-dynamic-range spectra (verified concretely in W11).

### W02 Frame-Recursive Online WPE

In `nara_wpe/wpe.py`, `OnlineWPE` maintains recursive inverse covariance, filter taps, power, and history ring buffers (`step_frame()` takes `(freq_bins, channels)` and predicts using *prior* taps before updating state; see [Caroselli et al., Interspeech 2017 §3](https://www.isca-archive.org/interspeech_2017/caroselli17_interspeech.pdf) and [Drude et al., ITG 2018 §4.3 Eqs. (11)–(17)](https://groups.uni-paderborn.de/nt/pubs/2018/ITG_2018_Drude_Paper.pdf)).

**Critical `0.0.11` lag-index shift across three same-named `delay` APIs.** On 1-bin, 1-mic, 1-tap sequence `[10, 20, 30, 40]` at 0-based frame $t = 3$ (`X_3 = 40`) with `delay = 1` and unit tap `1`:
- Offline `build_y_tilde()` taps $X_{t-1} = 30 \implies \text{residual } 10$;
- Stateless `online_wpe_step()` taps $X_{t-2} = 20 \implies \text{residual } 20$;
- Stateful `OnlineWPE.step_frame()` reads its ring buffer *before* inserting $X_t$, tapping $X_{t-3} = 10 \implies \text{residual } 30$.

**Multi-tap, multi-channel flattening order also differs across all three.** Feeding channel 0 = `[10, 20, 30, 40, 50, 60]` and channel 1 = `[100, 200, 300, 400, 500, 600]` with `2` taps and `delay = 1` at 0-based frame `5` (`X_5 = [60, 600]`), [`audit_upstream_wpe_contracts.py`](../../ch07/examples/audit_upstream_wpe_contracts.py) and [`upstream_wpe_contracts.json`](../../ch07/reports/upstream_wpe_contracts.json) verify:

| `nara_wpe` `0.0.11` Entry Point | Flattened Regressor Vector at Frame 5 | Axis & Temporal Ordering |
|---|---|---|
| `build_y_tilde()` | `[50, 500, 40, 400]` | Tap-major (channels inner); newest to oldest |
| `online_wpe_step()` | `[40, 30, 400, 300]` | Channel-major (taps inner); newest to oldest |
| `OnlineWPE._get_prediction()` (2D frame input) | `[20, 30, 200, 300]` | Channel-major (taps inner); **oldest to newest** |

Setting output-0 filter weights to `[1, 2, 3, 4]` and output-1 weights to `0` at current sample `60` yields exact hand-computed predictions:

$$
\begin{aligned}
p_{\mathrm{stateless}}&=40+2\times30+3\times400+4\times300=2500,\\
p_{\mathrm{state}}&=20+2\times30+3\times200+4\times300=1880.
\end{aligned}
$$

Hence stateless residual is $60 - 2500 = -2440$ whereas stateful `OnlineWPE` residual is $60 - 1880 = -1820$—matching unmodified upstream output identically.

Furthermore, [`compare_online_wpe_reference.py`](../../ch07/examples/compare_online_wpe_reference.py) (`seed=20260922`, `(48, 2, 1)`, `2` taps, `delay=2`, $\alpha = 0.95$) confirms `2.39e-15` max diff against our standalone MIT-attributed NumPy adapter, `0` diff across outer chunking vs. `2.56` jump when re-instantiating `OnlineWPE` mid-stream at frame `24`, while [`wpe_silence_boundary.py`](../../ch07/examples/wpe_silence_boundary.py) and [`chapter07_online_wpe_silence.json`](../../ch07/reports/chapter07_online_wpe_silence.json) prove that **prolonged zero input drives `OnlineWPE`'s inverse covariance `inv_cov /= alpha` to `Inf`/`NaN` overflow** despite power-inverse flooring.

### W03 Block-Online WPE and API Constraints

In `nara_wpe/wpe.py`, `OnlineWPE.step_block()` expects pre-stacked history windows with `_get_prediction()` asserting `block_shift == 1`, whereas [`nara_wpe/tf_wpe.py`](https://github.com/fgnt/nara_wpe/blob/a166779cca2088817e330481bd20af1a2c598555/nara_wpe/tf_wpe.py) implements distinct `block_wpe_step()` (default `2` s block, `forgetting_factor=0.7` weighting `0.3 * old + 0.7 * new` with per-block zero-padded head and intra-block future leakage) and `recursive_wpe()` (`only_use_final_filters=True` re-filters the entire utterance non-causally; static audit only, TF1 not executed). Executed [`wpe_temporal_contract.py`](../../ch07/wpe_temporal_contract.py) perturbs frame `30` by $+5 + 2j$ on `(48, 2, 1)`: causal `OnlineWPE.step_frame` exhibits exact `0` diff on frames `< 30`, $\sqrt{29}$ at frame `30`, and `≈ 1.261` thereafter (bit-identical across `48×1`, `7+13+5+23`, and `30+18` outer loops), whereas batch `wpe_v6` leaks `0.210` backward onto frames `< 30`.

### W04 DNN-WPE: Mask/Power Networks Coupled to Analytical Linear Solvers

[ESPnet `DNN_WPE`](https://github.com/espnet/espnet/blob/be79590bb2ff26ffb01bc825c5f68cb9418b7f0d/espnet2/enh/layers/dnn_wpe.py) (`dnn_wpe.py` + `mask_estimator.py` + `wpe.py`, mapping `(B, T, C, F) -> (B, F, C, T)`) multiplies predicted masks onto observed power in iteration 1 (`eps=1e-6` floor; [Kinoshita et al., Interspeech 2017 §5.1.2](https://www.isca-archive.org/interspeech_2017/kinoshita17_interspeech.pdf) uses `2` s blocks, $\alpha = 0.7$, `512/128` STFT at `16` kHz, `delay=3`, `37/10` taps, `3` iterations, and a unidirectional LSTM that still ingests `±5` lookahead/context frames). Static AST audit of pinned ESPnet uncovers three silent quirks: (1) constructor args `diagonal_loading`, `diag_eps`, and `use_torch_solver` are stored on `self` **yet never read or passed in `forward()`**, leaving `get_filter_matrix_conj()` hard-coded to `1e-10 * I` + `.inverse()`; (2) `MaskEstimator` line `82` calls out-of-place `masked_fill()` and **discards its return value** without zeroing padded mask frames prior to covariance estimation; and (3) `normalization=True` divides by temporal mask sum without an `eps` guard.

<a id="nemo-wpe"></a>

**NVIDIA NeMo: Same-named WPE carries incompatible mask, iteration, and regularization contracts.** Auditing pinned `v2.4.0` commit `2381f42f6979449b5b99538f8f80135831009b51` (`Apache-2.0` [LICENSE](https://github.com/NVIDIA-NeMo/Speech/blob/2381f42f6979449b5b99538f8f80135831009b51/LICENSE), 11 plain-text files fetched; [25.07 API docs](https://docs.nvidia.com/nemo-framework/user-guide/25.07/nemotoolkit/audio/api.html), [`masking.py` L1028–1063](https://github.com/NVIDIA-NeMo/Speech/blob/2381f42f6979449b5b99538f8f80135831009b51/nemo/collections/audio/modules/masking.py#L1028-L1063), [`multichannel.py` L754–1034](https://github.com/NVIDIA-NeMo/Speech/blob/2381f42f6979449b5b99538f8f80135831009b51/nemo/collections/audio/parts/submodules/multichannel.py#L754-L1034), [`enhancement.py` L99–134](https://github.com/NVIDIA-NeMo/Speech/blob/2381f42f6979449b5b99538f8f80135831009b51/nemo/collections/audio/models/enhancement.py#L99-L134); 12 AST assertions in [`upstream_wpe_contracts.json`](../../ch07/reports/upstream_wpe_contracts.json)):

| Contract Dimension | Pinned NeMo `WPEFilter` / `MaskBasedDereverbWPE` Behavior | Cross-Framework Porting Consequence |
|---|---|---|
| Mask-to-power formula | Squares $(m \cdot \lvert X\rvert)^2$ (`m` clamped to $[10^{-10}, 1]$) before averaging across channels | Squares the mask ($m^2 \lvert X\rvert^2$), unlike ESPnet's linear $m \lvert X\rvert^2$ |
| Multi-iteration target | Calls `self.filter(input=output, power=power)`—re-filtering the *previous residual* | Diverges from standard WPE (which re-estimates power from residual while always predicting *original* mixture $X$) |
| Startup & length masking | Left-pads zeros (`full` startup) and zeros weights past `input_length` *before* covariance sum | Differs from both NARA `"valid"` and ESPnet post-solve zeroing |
| Power & matrix loading | Inverse weight `1 / (mean_power + eps)` (`eps=1e-8`); matrix adds `diag_reg * trace(Q) + eps` (`diag_reg=1e-6`) | `trace(Q)` is NOT divided by `C * K` matrix dimension |
| Linear solver | Direct `torch.linalg.solve` (`complex128`); `diag_reg=None/0` skips matrix loading | No `lstsq` fallback on singular matrices |

Our standalone rational NumPy check confirms both differences: magnitude `2` with mask `0.5` yields power `1` under $(m|X|)^2$ vs. `2` under $m|X|^2$; on `[1, 2, 3, 4]` (`1` tap, `delay=1`, zero-padded head, unit initial power $\implies w_1 = 20/14 = 10/7$, residual `[1, 4/7, 1/7, -2/7]`), iteration 2 predicting original $X$ yields $w_2 = 146/101$, output `[1, 56/101, 11/101, -34/101]`, whereas re-filtering the residual yields $w_2 = 28/103$, output `[1, 216/721, -9/721, -234/721]`.

### W05 WPD: Simultaneous Spatial Filtering and Delayed Reverberation Cancellation

Weighted Power minimization Distortionless response (WPD; [Nakatani & Kinoshita, *IEEE SPL* 2019](https://doi.org/10.1109/LSP.2019.2911179) / [arXiv](https://arxiv.org/abs/1812.08400), ML formulation in [EUSIPCO 2019](https://arxiv.org/abs/1908.02710)) stacks current-frame channels and delayed history into a single distortionless optimization: [Exercise E07-05](../../../../chapters_en/07_wpe-dereverberation.md#e07-05) solves a 4D block system yielding `[0.4, 0.6, -0.2, 0]` (objective `0.6` vs. `2/3` when history taps are forced to zero), and [Exercise E07-07](../../../../chapters_en/07_wpe-dereverberation.md#e07-07) ([`enhancement_structure_exercises.py`](../cross_chapter/enhancement_structure_exercises.py)) verifies Schur-complement equivalence ($\mathbf{G} = \mathbf{C}^{-1}\mathbf{B}^\mathsf{H}$, $\mathbf{S} = \mathbf{A} - \mathbf{B}\mathbf{C}^{-1}\mathbf{B}^\mathsf{H}$) between WPD and factored `WPE -> weighted MPDR` under matched weights and loading ([Boeddeker et al., ICASSP 2020 §4 & Appendix](https://arxiv.org/pdf/1910.13707); note in [ESPnet `beamformer.py`](https://github.com/espnet/espnet/blob/be79590bb2ff26ffb01bc825c5f68cb9418b7f0d/espnet2/enh/layers/beamformer.py), `wpd` calls `get_WPD_filter_with_rtf()` whereas `wpd_souden` calls `get_WPD_filter_v2()`, and `beamformer.py` vs. `beamformer_th.py` differ).

### W06 AR-FastMNMF and Joint Dereverberation/Separation

AR-FastMNMF ([ICASSP 2021](https://ieeexplore.ieee.org/document/9414857); [`sekiguchi92/SoundSourceSeparation`](https://github.com/sekiguchi92/SoundSourceSeparation/tree/897fe87fea3d85a243d8a3fd36c2232bb0548ad3) `src/` and `FastBSSD.py`) optimizes autoregressive reverberation taps jointly with diagonalizable multi-source MNMF—governed by an academic-institution-only research license that precludes commercial product inclusion.

### W07 MINT, Cepstral Liftering, and Late Reverberant Spectral Variance

- **Known-RIR MINT** ([Miyoshi & Kaneda, 1988](https://www.kanedayyy.jp/asp/pub/sig_ac/MINT1988paper.pdf)): exact BézoutFIR inversion requires known RIRs with zero common polynomial zeros across channels and degrades sharply under near-common zeros or sensor noise ([Exercises E07-15–E07-17](../../../../chapters_en/07_wpe-dereverberation.md#e07-15) verify QR/SVD conditioning, regularized MINT, and quantized PCM scoring).
- **Homomorphic Cepstral Liftering** ([Oppenheim, Schafer & Stockham, *Proc. IEEE* `56(8):1264–1291`, 1968, §VI *Echo Removal*](https://dsp-group.mit.edu/wp-content/uploads/2024/11/nonnlinearfiltering_1968.pdf), [`10.1109/PROC.1968.6570`](https://doi.org/10.1109/PROC.1968.6570)): real power cepstrum $\operatorname{IFFT}(\log|X|^2) = 2\operatorname{IFFT}(\log|X|)$ discards phase completely ([Exercise E07-14](../../../../chapters_en/07_wpe-dereverberation.md#e07-14)) and concentrates only simple periodic echoes rather than dense room RIRs.
- **Late Reverberant Spectral Variance** ([Habets, Gannot & Cohen, *IEEE SPL* `16(9):770–773`, 2009, §§II–V](https://israelcohen.com/wp-content/uploads/2018/05/SPL_Sep2009.pdf)): models diffuse late decay parameterized by $T_{60}(f)$ and DRR ([Exercise E07-18](../../../../chapters_en/07_wpe-dereverberation.md#e07-18): at $T_{60} = 0.6$ s and $H = 8$ ms hop, power decays by $0.831764\times$/frame = `-60` dB over `75` frames; note `SET_DEREVERB` in pinned SpeexDSP `preprocess.c` has its energy update commented out and is non-functional).

### W08 Industrial Evidence, BTK20 C++ Source, and GPU-GSS WPE Tap Points

- **Google Home 2017** ([Li et al., Interspeech 2017 §§2–3, Fig. 1](https://www.cs.cmu.edu/~chanwook/MyPapers/b_li_interspeech_2017.pdf); [Caroselli et al. 2017 §3](https://www.isca-archive.org/interspeech_2017/caroselli17_interspeech.pdf)) documents streaming 2-mic WPE ahead of fCLP + acoustic modeling without releasing production code.
- **NTT Evaluation Package** ([official page](https://www.rd.ntt/cs/team_project/media/signal/wpe/), [evaluation-only license](https://www.rd.ntt/cs/team_project/media/signal/wpe/licence.html), [8-mic `512/128` config](https://www.rd.ntt/cs/team_project/media/signal/wpe/config.html)) ships closed MATLAB p-code (`index_only`, distinct from MIT `nara_wpe`).
- **BTK20 C++** ([`dereverberation.cc`](https://github.com/kkumatani/distant_speech_recognition/blob/feff19ec8bcb770f6530fe280dc3ccafc2f5984a/btk20_src/dereverberation/dereverberation.cc)): inclusive `[lowerN, upperN]` bounds (`upperN - lowerN + 1` taps, channel-grouped, zero-padded past `lowerN`), independent per-output residual power (floor `1e-6`), `loadDb` + `diagonal_bias` GSL Cholesky solve over **the entire buffered utterance** before resetting stream readback (batch despite `FeatureStream` naming).
- **GPU-GSS WPE** ([`gss/core/enhancer.py`](https://github.com/desh2608/gss/blob/10fad18cae85e2e4342c77421abc70c9c5da23ed/gss/core/enhancer.py) & [`desh2608/wpe`](https://github.com/desh2608/wpe/tree/bd2857b5b8de36df4f436a93574c088bea142042)): defaults to `15` s context, `10` taps, `delay=2`, `3` iterations, `1024/256` STFT (`bf_drop_context` trims beamformer masks *after* WPE has already used full context, and OOM frequency-chunking alters global `1e-10 * max` power floors).

### W09 Two Fetched Neural WPE Extensions: MetaAF-WPE and TSO-VACE-WPE

- **MetaAF-WPE** ([`zoo/wpe/wpe.py`](https://github.com/adobe-research/MetaAF/blob/56c4665bdc51c2e0595a7c0cd9b1266408adceff/zoo/wpe/wpe.py); [paper §VII-C & Fig. 13](https://arxiv.org/pdf/2204.11942)): replaces closed-form normal equations with GRU/FGRU tap updates (`w * d_input` unconjugated convention, root-Hann `512/256`, `5` taps, `delay=2`; note higher SRR in Fig. 13 can accompany *lower* STOI, and `nara_wpe_eval.py::fit_nara()` passing 3D windows into `OnlineWPE.step_frame()` actually matches single-frame calls bit-for-bit in W11).
- **TSO-VACE-WPE** ([Yang & Chang, *IEEE/ACM TASLP* 2022](https://doi.org/10.1109/TASLP.2022.3205752) / [arXiv §II-A, Fig. 1](https://arxiv.org/pdf/2112.13569); [`dreadbird06/tso_vace_wpe`](https://github.com/dreadbird06/tso_vace_wpe/tree/10ee77dd020d58af508feb77251f9353208cb33a), `MIT` code fetched excluding `models/`/`data/`): synthesizes a virtual 2nd channel from mono input (`1024/256` STFT, `delay=3`, `30/15` real/virtual taps, PSD net observing real mic only, per-bin `1e-10 * max` floor unprotected on all-zero bins).

### W10 Distributed Arrays and 2024–2026 Dereverberation Candidates

[Lohmann et al. (ICASSP 2023, §§4–5)](https://arxiv.org/pdf/2301.07649) decompose inter-node TDOA into integer frames, integer samples, and fractional samples (`9` mics, `16` kHz, `1024/256`, noiseless):

| Candidate & Reference | Modified Step & Inclusion Decision |
|---|---|
| [Lohmann et al. 2024 Mic Subset Selection, §§3–4](https://arxiv.org/pdf/2401.08486) | Adds group-lasso sparsity on non-reference mic taps before re-running WPE; `index_only` (no verified author release) |
| [Lohmann et al. 2024 Reference Selection, §§4–5](https://arxiv.org/pdf/2411.03168) | Scores candidate reference channels via scale-invariant output sparsity (requires candidate WPE runs); unexecuted |
| [Stream.FM (*IEEE TASLPRO* 2026)](https://doi.org/10.1109/TASLPRO.2026.3696215) | Replaces linear WPE with causal multi-stage flow matching; audited statically in N14 |

<a id="two-stage-dereverb"></a>

**Explicitly tunable early-reflection training targets.** [Lemercier et al. (*EURASIP JASM* 2023, `10.1186/s13636-023-00285-8`)](https://doi.org/10.1186/s13636-023-00285-8) ([arXiv §§2.4, 3, 4](https://arxiv.org/pdf/2204.02978); [UHH page](https://www.inf.uni-hamburg.de/en/inst/ab/sp/publications/tasl2022-dereverberation.html) linking unlicensed inference-only [`sp-uhh/2sderev`](https://github.com/sp-uhh/2sderev)) places loss directly on DNN-WPE output followed by a Wiener post-filter—retaining early `40` ms RIR reflections for hearing-aid targets vs. `16` ms for cochlear-implant targets while masking startup frames from training loss.

<a id="pnp-wpe-candidate"></a>

**PnPWPE.** [2023 preprint §§II–III, Eqs. (22)–(38)](https://arxiv.org/pdf/2312.02773) couples noisy WPE with Regularization-by-Denoising (RED) priors via variable splitting (`index_only`, no verified author code/weights).

<a id="known-arrival-control"></a>

**Known-arrival-delay source-index control.** Following [Lohmann et al. ICASSP 2023 §4](https://arxiv.org/html/2301.07649v1), [Exercise E07-19](../../../../chapters_en/07_wpe-dereverberation.md#e07-19) and [Figure 70](../../../../figures_en/fig70_prediction_delays.png) ([`delay_audio/MANIFEST.json`](../../ch07/delay_audio/MANIFEST.json), [Exercises Guide §57](05_exercises_and_audio.md#known-arrival-audio)) demonstrate on 6 synthetic PCM WAVs how uniform prediction delay on an early-arriving auxiliary mic copies and cancels the direct target completely, whereas arrival-compensated lag indexing preserves it.

<a id="rir-encoder-candidate"></a>

**Contrastive RIR embeddings for generative restoration (`rir-encoder`).** [Khanagha & Gerkmann (2026 preprint §3)](https://arxiv.org/html/2606.09557v1) modulate SGMSE+ via FiLM using contrastive RIR embeddings; 8 MIT-licensed plain-text files (`20717` B) from pinned [`sp-uhh/rir-encoder`](https://github.com/sp-uhh/rir-encoder/tree/8ee0ba7e083c8e38170dc13632cc9937a702d0f0) (`8ee0ba7e083c8e38170dc13632cc9937a702d0f0`) were fetched on 2026-10-04, revealing three static boundaries: `ReverbDataset._stft` sets `center=False` whereas extraction `emb_stft` defaults to `center=True`; `AdaptiveAvgPool1d(1)` averages across full padded time axes despite passing `lengths` to Conformer; and per-worker RIR label dicts + import-time side effects in `train.py` preclude blind import.

<a id="usdnet-plus-candidate"></a>

**USDnet++ weak supervision from WPE/WPD.** [Pang et al. (Interspeech 2026, pp. 4490–4494, §§3–4, Eqs. (3)–(15))](https://www.isca-archive.org/interspeech_2026/pang26_interspeech.pdf) uses WPE/WPD outputs as auxiliary unsupervised targets alongside FCP (with optional future taps and two-pass dual-network inference; `index_only`, no verified author code).

### W11 Pinned WPE Source Diagnostics and Literature Audit Scope

Re-executed on 2026-10-04 into three dedicated `*_current.json` reports—[`upstream_wpe_contracts_current.json`](../../ch07/reports/upstream_wpe_contracts_current.json), [`wpe_upstream_interfaces_current.json`](../../ch07/reports/wpe_upstream_interfaces_current.json), and [`chapter07_online_wpe_silence_current.json`](../../ch07/reports/chapter07_online_wpe_silence_current.json)—while preserving byte-identical 2026-09-28 ([`wpe_upstream_interfaces.json`](../../ch07/reports/wpe_upstream_interfaces.json), [`test_codes_wpe_upstream_interfaces.py`](../../../../tests/test_codes_wpe_upstream_interfaces.py)) and 2026-10-01 ([`upstream_wpe_contracts.json`](../../ch07/reports/upstream_wpe_contracts.json), [`test_codes_ch07_upstream_wpe_contracts.py`](../../../../tests/test_codes_ch07_upstream_wpe_contracts.py), bound to commit `34d494dd8074169a92f79bf90b9ff2bbcf25220f`) historical reports:
- **Prolonged-silence overflow**: `OnlineWPE` first turns non-finite at frame **`1024`** ($\alpha = 0.5$) and frame **`13838`** ($\alpha = 0.95$), failing to recover on subsequent non-zero speech.
- **Executed `wpe_v6/v7` vs. `wpe_v8` dynamic-range counterexample**: On 2-bin, 1-mic, 6-frame real `complex128` input with bin 0 = `[1, 2, 1, 4, 2, 1]` and bin 1 = `1e-8 * [1, 3, 2, 1, 5, 2]` (`1` tap, `delay=1`, `1` iter, `"valid"`), global-floor `v6/v7` clamps all six bin-1 inverse powers to `6.25e8` yielding tap **`13/20`**, whereas per-bin `v8` yields exact tap **`2940/5693`**. Moreover, helper `get_power()` on `arange(1, 13).reshape(2, 1, 6)` with `psd_context=1` averages across *frequency* (`[25, 34, 45, 58, 73, 90]` on both rows) and raises `ValueError` on `psd_context=(0, 2)`.

```bash
.venv/bin/python -m codes.chapters.ch07.examples.audit_upstream_wpe_contracts --report codes/chapters/ch07/reports/upstream_wpe_contracts_current.json
.venv/bin/python -m codes.chapters.ch07.examples.audit_wpe_upstream_interfaces --report codes/chapters/ch07/reports/wpe_upstream_interfaces_current.json
.venv/bin/python -m codes.chapters.ch07.examples.wpe_silence_boundary --upstream --report codes/chapters/ch07/reports/chapter07_online_wpe_silence_current.json
```

<a id="bss"></a>

## 4. Blind Source Separation, Spatial Mixture Models, and GSS

### B01 FDICA and the Frequency Permutation Problem

Frequency-Domain ICA ([`ssspy` `bss/fdica.py`](https://github.com/tky823/ssspy/blob/38b9389e8b1914422561f1936d9b28d042d62d2c/ssspy/bss/fdica.py) + `algorithm/` permutation aligners) estimates independent demixing matrices per DFT bin, requiring inter-frequency permutation alignment distinct from utterance-level PIT evaluation (note our teaching `pit_permutation()` enumerates $K!$ up to $K \le 8$, whereas larger $K$ uses $O(K^3)$ Hungarian assignment on additive pairwise losses).

### B02 AuxIVA, IP/ISS Updates, and Projection-Back Scale Restoration

Auxiliary-function IVA ([Ono, WASPAA 2011, pp. 189–192](https://doi.org/10.1109/ASPAA.2011.6082320); [`ssspy` `bss/iva.py`](https://github.com/tky823/ssspy/blob/38b9389e8b1914422561f1936d9b28d042d62d2c/ssspy/bss/iva.py) implementing `IP1, IP2, ISS1, ISS2, IPA`) couples all frequency bins via joint source spherical/Laplace norms ([Exercises E08-04–E08-07](../../../../chapters_en/08_speech-separation.md#sec-8-8)).

**Executed: Blind `AuxLaplaceIVA` estimation on pinned `ssspy` (`38b9389e8b1914422561f1936d9b28d042d62d2c`, Apache-2.0).** [`reproduce_auxiva_reference.py`](../../ch08/examples/reproduce_auxiva_reference.py) (`seed=2468`, two 8-segment modulated Gaussian sources at `8` kHz, `4` s, full-rank mixing `[[1, 0.9], [0.8, 1]]`, `256/64` periodic Hann STFT round-trip diff `2.7e-15`, `30` IP iterations + mic-0 projection-back) raises mean SI-SDR from `-0.02` dB to **`37.25` dB** (`SI-SDRi = 37.27` dB, output/ref RMS ratios `0.999, 0.999`, internal contrast `568.99 -> -8.69`), whereas feeding deterministic harmonic textbook WAVs yields `0.14` dB (`RMS` ratios `0.536, 1.674`—illustrating breakdown under strongly correlated harmonic envelopes), and rank-1 mixing raises `LinAlgError: Singular matrix`.

**BSS Eval v3/v4 vs. zero-mean SI-SDR boundaries.** Whereas [Chapter 8 SI-SDR](../../../../chapters_en/08_speech-separation.md#sec-8-5) permits strictly a single scalar projection without time shifts, [Vincent, Gribonval & Févotte (2006, §§II–III)](https://doi.org/10.1109/TSA.2005.858005) decompose `e_spat + e_interf + e_artif` via `512`-tap FIR filters—differing sharply between [`mir_eval` `0.8.2` (`MIT`)](https://github.com/mir-evaluation/mir_eval/blob/8db0b3812e2032544c1fc00d02d4256cab043f3d/mir_eval/separation.py) ([LICENSE](https://github.com/mir-evaluation/mir_eval/blob/8db0b3812e2032544c1fc00d02d4256cab043f3d/LICENSE.txt)) and [`museval` `0.4.1` BSS Eval v4 (`MIT`)](https://github.com/sigsep/sigsep-mus-eval/blob/b2a2246f343eeb3c17dd69882db3eb73b6714617/museval/metrics.py) ([v4 docs](https://github.com/sigsep/sigsep-mus-eval#bsseval-v4), [LICENSE](https://github.com/sigsep/sigsep-mus-eval/blob/b2a2246f343eeb3c17dd69882db3eb73b6714617/LICENSE)), where `bsseval_sources_version=False` places spatial distortion `e_spat` in the **denominator** rather than numerator.

**Executed: Three `projection_back` contracts on identical inputs + PRA `init_eig=True` crash.** On correlated outputs `[1, 0, 1]` and `[0, 1, 1]` with ref `[1, 2, 3]`, `ssspy` joint least squares and matrix-inverse row 0 both recover exact `[1, 2]`, whereas `pyroomacoustics` per-source scalar fitting returns **`[2, 2.5]`** (and returns **conjugated** `2 - j` on `[1, j]` vs. `(2 + j)[1, j]`, [Exercise E08-11](../../../../chapters_en/08_speech-separation.md#e08-11), [B12 report](../../ch08/reports/separation_upstream_interfaces.json)). Moreover, calling PRA `0.10.0` `auxiva(..., init_eig=True, n_iter=0)` on 2-mic, 1-source input crashes immediately with `ValueError` assigning `(2, 1)` eigenvectors into `(1, 2)` source rows.

### B03 OverIVA, FIVE, and Overdetermined Target Extraction

When $M > N_{\text{src}}$, OverIVA ([Scheibler & Ono, 2019 §§2–3](https://arxiv.org/html/1905.07880v2)) pairs non-Gaussian/time-varying target rows with an orthogonal stationary Gaussian background subspace—avoiding energy-PCA's fatal pitfall of discarding weak targets below strong background noise (`diag(9, 4, 1)` counterexample and exact `W=[1, 0.5], C=[[2, 1], [1, 2]] -> J=0.8, U=[0.8, -1], det=-1.4` verification in [Exercise E08-23](../../../../chapters_en/08_speech-separation.md#e08-23) and [`overiva_teaching.py`](../../ch08/core/overiva_teaching.py)). Static audit of the 4 MIT-licensed plain-text files fetched on 2026-10-04 from [`onolab-tmu/overiva` (`1cb3189112889ebfe2ccbbde0f55b1db6020fc48`, license SHA `a872a2a232076ca5034f53a33709efecfe61c7481dbda84b595549f7015f03bc`)](https://github.com/onolab-tmu/overiva/tree/1cb3189112889ebfe2ccbbde0f55b1db6020fc48) notes `(frames, bins, mics)` input, column-oriented `(bins, mics, targets)` `W` unscaled by projection-back, unguarded `kwargs.pop('proj_back')` in `auxiva_pca.py`, and separate GPL-3.0 licensing on [`fakufaku/piva`](https://github.com/fakufaku/piva/tree/7fa273e9aa597aba57067aef2e7b6c999dadef26) (which also implements single-target `FIVE`).

### B04 ILRMA Spectral Low-Rank NMF and Spatial Updates

Independent Low-Rank Matrix Analysis ([Kitamura et al., *IEEE/ACM TASLP* 2016](https://doi.org/10.1109/TASLP.2016.2577880); [`pyroomacoustics` `bss/ilrma.py`](https://github.com/LCAV/pyroomacoustics/blob/v0.10.0/pyroomacoustics/bss/ilrma.py); [`ssspy` `bss/ilrma.py`](https://github.com/tky823/ssspy/blob/38b9389e8b1914422561f1936d9b28d042d62d2c/ssspy/bss/ilrma.py) where `T @ V` uses `(sources, bins, bases) @ (sources, bases, frames)` despite swapped docstring labels; [Exercise E08-08](../../../../chapters_en/08_speech-separation.md#e08-08)) exposes a concrete **axis-scaling bug in pinned `pyroomacoustics` `0.10.0`**: on a `(9, 3, 2)` test spectrum (`PCG64(8)`, `np.random.seed(7)`, `2` bases, `1` iter, `proj_back=False`, scale factors `≈ 0.721573, 0.351573`), upstream scales **input mic columns `W[:, :, s]`** instead of source rows and returns the *pre-scaled* signal spectrum—leaving returned `Y` and `W @ X` inconsistent by **`3.7390` max abs diff** ([B12](#separation-upstream-audit)).

### B05 Full-Rank Spatial Covariance MNMF

Multichannel NMF ([Sawada et al., *IEEE TASLP* 2013 §§IV–V](https://www.kecl.ntt.co.jp/icl/signal/sawada/mypaper/IEEEtaslp2013sawada.pdf), [`10.1109/TASL.2013.2239990`](https://doi.org/10.1109/TASL.2013.2239990); [`ssspy` `bss/mnmf.py::GaussMNMF`](https://github.com/tky823/ssspy/blob/38b9389e8b1914422561f1936d9b28d042d62d2c/ssspy/bss/mnmf.py)) models full-rank diffuse/underdetermined source SCMs reconstructed via multichannel Wiener filtering—requiring multi-seed dispersion reporting because mixture-sum consistency $\sum_n \hat{\vec{y}}_n = \vec{x}$ holds trivially for any Wiener partition even when sources fail to separate.

### B06 FastMNMF, FastMNMF2, and Joint Diagonalization

FastMNMF ([Sekiguchi et al., EUSIPCO 2019 §§II–III](https://www.eurasip.org/Proceedings/Eusipco/eusipco2019/Proceedings/papers/1570533283.pdf) / [arXiv](https://arxiv.org/abs/1903.03237)) and directivity-constrained FastMNMF2 ([*IEEE/ACM TASLP* 2020](https://doi.org/10.1109/TASLP.2020.3019181); MIT implementations in [`pyroomacoustics` `fastmnmf.py`](https://github.com/LCAV/pyroomacoustics/blob/v0.10.0/pyroomacoustics/bss/fastmnmf.py) & [`fastmnmf2.py`](https://github.com/LCAV/pyroomacoustics/blob/v0.10.0/pyroomacoustics/bss/fastmnmf2.py), distinct from academic-institution-only [`SoundSourceSeparation` LICENSE](https://github.com/sekiguchi92/SoundSourceSeparation/blob/897fe87fea3d85a243d8a3fd36c2232bb0548ad3/LICENSE)) jointly diagonalize source SCMs (`G_NFM` per-bin vs. frequency-shared `g_NM`): our executed `0`-iteration check ([B12](#separation-upstream-audit)) proves both **mutate caller `W0 = 2I` in-place to `I`** without copying (`Q_FMM = W0`) and return `(mics, frames, bins, sources)` on `mic_index="all"`, while static reading notes [`fastmnmf2.py` L143–147](https://github.com/LCAV/pyroomacoustics/blob/0dd39f2614b7fc44b2cc63dbe7d60f4641068890/pyroomacoustics/bss/fastmnmf2.py#L143-L147) caches `tmp1/tmp2` before refreshing `lambda/Y` under `accelerate=False` (and indexes block-diagonal [Distributed FastMNMF 2026 preprint §§3.1–3.2](https://arxiv.org/html/2605.19388v1) as synchronous centralized multi-array modeling without released code).

### B07 TRINICON Broadband Time-Domain Block Adaptation

Implemented in [`pyroomacoustics` `bss/trinicon.py`](https://github.com/LCAV/pyroomacoustics/blob/v0.10.0/pyroomacoustics/bss/trinicon.py) following [Aichner et al. (*Signal Processing* 2006, §§1–2 & Table 1)](https://www.buchner-net.com/lnt2006_19.pdf) ([DOI](https://doi.org/10.1016/j.sigpro.2005.06.022); [Buchner ICASSP 2004](https://www.buchner-net.com/publications.html)), our executed `j_max=0` frozen-filter check (`filter_length=4, block_length=8, n_blocks=2` $\implies$ hop `4`, ramp `0.001 * n`) uncovers an **axis-modulo padding bug**: while `16` samples return `16` exact delayed samples, **`17` samples return `19` samples with output index `16` zeroed (`0.0` instead of exact `0.014`)** because upstream computes remainder modulo `n_mics` (`axis=0`) rather than `n_samples` (`axis=1`) ([B12](#separation-upstream-audit)).

### B08 cACGMM: Directional Unit-Vector Clustering, Activity Gating, and Underflow

Complex Angular Central Gaussian Mixture Models ([Ito et al., EUSIPCO 2016 Eq. (14)](https://www.eurasip.org/Proceedings/Eusipco/Eusipco2016/papers/1570256519.pdf), [`10.1109/EUSIPCO.2016.7760429`](https://doi.org/10.1109/EUSIPCO.2016.7760429); [`pb_bss` `cacgmm.py`](https://github.com/fgnt/pb_bss/blob/10acc347fc9ea21e3d312806a0bd751d0d0af183/pb_bss/distribution/cacgmm.py); [Chapter 8 §8.4's two-frame worked example](../../../../chapters_en/08_speech-separation.md#sec-8-4) and [Exercise E08-09](../../../../chapters_en/08_speech-separation.md#e08-09), [`enhancement_structure_exercises.py`](../cross_chapter/enhancement_structure_exercises.py), [`test_codes_enhancement_structure.py`](../../../../tests/test_codes_enhancement_structure.py)) cluster unit-normalized directional vectors while constructing beamformer SCMs from *unnormalized* STFT observations. Crucially, executing AST-extracted `gss/cacgmm/utils.py` posterior math reveals two engineering behaviors in GPU-GSS: (1) clamping posteriors to `[1e-10, 1 - 1e-10]` without re-normalizing gives `[ε, ε, 1 - ε]` (sum `1 + ε`) when only background is active; and (2) subtracting max log-density across **all** components *before* activity masking causes **total underflow to `[0, 0, 0]`** when an *inactive* component has log-density `1000` vs. active components at `0` (where exact active-support normalization gives `[0, 3/8, 5/8]`).

### B09 GSS: Diarization Activity, WPE, Guided cACGMM, and MVDR Beamforming

Guided Source Separation ([Boeddeker et al., CHiME-5 2018 §3.2, Fig. 3, Eqs. (5)–(10)](https://www.isca-archive.org/chime_2018/boeddecker18_chime.pdf); [Raj et al. GPU-GSS, Interspeech 2023](https://www.isca-archive.org/interspeech_2023/raj23_interspeech.pdf); [`desh2608/gss`](https://github.com/desh2608/gss/tree/10fad18cae85e2e4342c77421abc70c9c5da23ed)) adds an always-active noise class so all-speaker-silent frames retain positive denominator mass, verified across:
- **Executed E-step mislabeling check** ([`gss_activity_error_demo.py`](../../ch08/gss_activity_error_demo.py), [`guided_activity_posterior()`](../../ch08/core/separation.py)): priors `[0.4, 0.4, 0.2]`, relative densities `[8, 1, 1]` $\implies$ true active posterior `16/17` vs. `0` when missed, and `2/3` false target posterior on silence when falsely marked active.
- **Executed teaching `cACGMM -> SCM -> MVDR` 5-WAV suite** ([Chapter 8 §8.4.5](../../../../chapters_en/08_speech-separation.md#sec-8-4-5), [`gss_teaching.py`](../../ch08/core/gss_teaching.py), [`gss_audio/MANIFEST.json`](../../ch08/gss_audio/MANIFEST.json), [`STATE.npz`](../../ch08/gss_audio/STATE.npz); read-only replay `.venv/bin/python -m codes.chapters.ch08.examples.gss_teaching_demo --check`): raises `[3200, 23200)` SI-SDR from `1.280328` dB (PCM `1.280339` dB) to **`13.761977` dB** (PCM **`13.762022` dB**), falling back cleanly to ref mic 0 when target activity is 100% missed.
- **Pinned GPU-GSS `bss_iterations_post=1` & power-SCM tap points** ([`audit_upstream_separation_contracts.py`](../../ch08/examples/audit_upstream_separation_contracts.py)): when default `bss_iterations_post=1` is non-zero, final `cur.predict(Obs.T)` **omits the activity mask**, while `bf_block()` rightly passes unnormalized post-WPE `Obs` into `get_power_spectral_density_matrix()` (verified on `[[1, 2, 0], [0, 0, 3]]` with weights `[1, 1, 0]` $\to \operatorname{diag}(2.5, 0)$ vs. `1e-12` scaled weights hitting `1e-10` floor $\to \operatorname{diag}(0.05, 0)$).
- **GPU-GSS vs. NeMo [`MaskEstimatorGSS`](https://github.com/NVIDIA/NeMo/blob/2381f42f6979449b5b99538f8f80135831009b51/nemo/collections/audio/modules/masking.py) comparison** ([`upstream_separation_contracts.json`](../../ch08/reports/upstream_separation_contracts.json)):

| Porting Dimension | Pinned `desh2608/gss` (GPU-GSS) | Pinned NeMo `MaskEstimatorGSS` |
|---|---|---|
| Pipeline scope | Utterance + context WPE, cACGMM, and MVDR | Batch `(B, M, F, T)` layer returning `(B, K, F, T)` masks only |
| Always-active noise class | Appended automatically inside GSS entry | Caller must append noise row explicitly in `activity` |
| Directional normalization | Zero-norm divisor guard; raw WPE spectrum fed to MVDR | $x / (\|x\| + \epsilon)$ (sub-unit norm on tiny inputs) |
| Posterior gating math | Clamped to `[1e-10, 1 - 1e-10]`; unmasked final `predict` when `post > 0` | Shifts by max log-pdf across *all* components before masking + denominator `eps` |

### B10 GPU-GSS Batching Parameters and Memory Boundaries

GPU-GSS ([`README.md`](https://github.com/desh2608/gss/blob/10fad18cae85e2e4342c77421abc70c9c5da23ed/README.md); [Raj et al. 2023 §4 & Table 4](https://www.isca-archive.org/interspeech_2023/raj23_interspeech.pdf)) batches frequency bins and speaker segments (`max-batch-duration`, `max-segment-length`, `context-duration=15` s, `20` BSS iters across CHiME-6's **six 4-mic Kinect arrays = 24 channels**, never misquoted as a "6-mic array"), measuring cumulative job-hour speedups rather than single-utterance streaming latency.

<a id="mffca"></a>

### B11 mfFCA: Multi-Frame Full-Rank Spatial Covariance Modeling

Multi-frame Full-rank Spatial Covariance Analysis ([Sawada et al., *IEEE/ACM TASLP* `31:3589–3602`, 2023, §§III–IV](https://www.kecl.ntt.co.jp/icl/signal/sawada/mypaper/IEEEtaslp2023sawada.pdf), [`10.1109/TASLP.2023.3313446`](https://doi.org/10.1109/TASLP.2023.3313446); 9 evaluation-licensed files fetched from [`nttcslab-sp/mfFCA` `1d6b422fc56f5f9ef612fd87ae1402e1311c5c3c`](https://github.com/nttcslab-sp/mfFCA/tree/1d6b422fc56f5f9ef612fd87ae1402e1311c5c3c) under non-redistributable internal-evaluation `LICENCE.txt`) stacks shifted STFT frames in `fca.py::mfFCA` (`(frames, bins, mics)`)—noting that positive shifts in `_shifted_matrix()` read **future frames** (non-causal), `optimizationEM()` floors spatial matrices/powers at `param_floor=1e-4`, and `main_synthetic.ipynb` / CuPy `main_speech.py` remain unexecuted locally.

<a id="separation-upstream-audit"></a>

### B12 Pinned Separation Interface Diagnostics: Execution Scope and Offline Verification

Executed via [`audit_separation_upstream_interfaces.py`](../../ch08/examples/audit_separation_upstream_interfaces.py) ([historical report](../../ch08/reports/separation_upstream_interfaces.json)) and expanded on 2026-10-01 / 2026-10-02 via [`audit_upstream_separation_contracts.py`](../../ch08/examples/audit_upstream_separation_contracts.py) ([report](../../ch08/reports/upstream_separation_contracts.json), [`test_codes_ch08_upstream_separation_contracts.py`](../../../../tests/test_codes_ch08_upstream_separation_contracts.py)—recording honestly how `ssspy` `CACGMM` import failed on missing `packaging` in `/private/tmp/masp-ch04-pra-venv` on 2026-10-01 before passing 0-iter posterior + 1-step shape update to `0` error inside `packaging`-equipped `/private/tmp/masp-appb-pra-venv` on 2026-10-02):

| Audited Upstream Target | Observed Pinned Behavior | Invalid Over-Extrapolation to Avoid |
|---|---|---|
| PRA `AuxIVA` `init_eig=True` | Raises `ValueError` on 2-mic, 1-source input | Assuming default identity init fails |
| `projection_back` APIs | `ssspy` joint LS `[1, 2]` $\ne$ PRA scalar `[2, 2.5]` (conjugated) | Assuming same-named functions are interchangeable |
| PRA `ILRMA` return contract | Scales mic axis & returns pre-scaled `Y` (`3.7390` diff vs. `W @ X`) | Blaming mathematical ILRMA theory |
| PRA `TRINICON` tail padding | `17`-sample input zeros valid index `16` (`0.0` vs. `0.014`) | Claiming adaptive separation quality was benchmarked |
| PRA `FastMNMF` / `FastMNMF2` | Mutates caller `W0` in-place even at `n_iter=0` | Mistaking `0`-iter Wiener sum identity for source separation |
| GPU-GSS posterior helper | `[1e-10, 1-1e-10]` clamping & inactive-max log-pdf underflow | Claiming full CuPy CHiME-6 WER failure |

```bash
PYTHONDONTWRITEBYTECODE=1 /private/tmp/masp-ch04-pra-venv/bin/python -B \
  -m codes.chapters.ch08.examples.audit_upstream_separation_contracts
.venv/bin/python -B -m unittest discover -s tests \
  -p 'test_codes_ch08_upstream_separation_contracts.py' -v
```

<a id="neural"></a>

## 5. Neural Separation, Target Speaker Extraction, and Continuous Output

**2026-10-04 current source-contract runs archived separately.** Re-executing all four Chapter 8 audit scripts produced [`separation_upstream_interfaces_current.json`](../../ch08/reports/separation_upstream_interfaces_current.json), [`upstream_separation_contracts_current.json`](../../ch08/reports/upstream_separation_contracts_current.json), [`streamfm_source_audit_current.json`](../../ch08/reports/streamfm_source_audit_current.json), and [`tflocoformer_source_audit_current.json`](../../ch08/reports/tflocoformer_source_audit_current.json)—leaving all four historical JSON reports byte-identical and bound to historical tool commit `7b80fab1d1012f104cdf441423218bf0898ce918`.

Before running, each current tool validates destination report paths, verifies official `origin`, isolated Git root, pinned `HEAD`, every loaded source/license Git blob and SHA-256, dependency hashes, and working-tree cleanliness before and after execution, keeps `source_selection_mismatch` statuses separate from used-blob verification, and compiles PRA/`ssspy` `.py` files directly from raw source bytes (bypassing stale `.pyc` caches and rejecting pre-imported unbound modules; verified via [negative-control tests](../../../../tests/test_codes_ch08_source_contracts.py)):

```bash
/private/tmp/masp-appb-pra-venv/bin/python -m codes.chapters.ch08.examples.audit_separation_upstream_interfaces --report codes/chapters/ch08/reports/separation_upstream_interfaces_current.json
/private/tmp/masp-appb-pra-venv/bin/python -m codes.chapters.ch08.examples.audit_upstream_separation_contracts --report codes/chapters/ch08/reports/upstream_separation_contracts_current.json
.venv/bin/python -m codes.chapters.ch08.examples.streamfm_source_audit --report codes/chapters/ch08/reports/streamfm_source_audit_current.json
.venv/bin/python -m codes.chapters.ch08.examples.tflocoformer_source_audit --output codes/chapters/ch08/reports/tflocoformer_source_audit_current.json
```

### N01 Conv-TasNet: Learned Short-Window Encoder, TCN, and Causal Normalization

Conv-TasNet ([Luo & Mesgarani, *IEEE/ACM TASLP* 2019](https://doi.org/10.1109/TASLP.2019.2915167) / [arXiv §II-C Eqs. (12)–(14)](https://arxiv.org/pdf/1809.07454); [`asteroid` `models/conv_tasnet.py`](https://github.com/asteroid-team/asteroid/blob/fce87469132760fbab41c20616ea0f0e079aad38/asteroid/models/conv_tasnet.py) + `masknn/norms.py`) highlights an important abbreviation collision: whereas the paper writes `cLN` for *cumulative* causal normalization across past frames, Asteroid's `masknn/norms.py` binds `cLN = ChanLN` (instantaneous per-frame channel-wise norm) and `cgLN = CumLN` (cumulative running norm), while default non-causal configs use utterance-global `gLN`.

### N02 DPRNN: Intra-Chunk, Inter-Chunk, and Overlap-Add Dual-Path Modeling

Dual-Path RNN ([Luo et al., ICASSP 2020](https://doi.org/10.1109/ICASSP40776.2020.9054266) / [arXiv §2.1](https://arxiv.org/pdf/1910.06379); [`asteroid` `models/dprnn_tasnet.py`](https://github.com/asteroid-team/asteroid/blob/fce87469132760fbab41c20616ea0f0e079aad38/asteroid/models/dprnn_tasnet.py) + `masknn/recurrent.py`) segments long feature sequences into overlapped 3D chunks—noting statically that Asteroid's intra-chunk RNN is *always* bidirectional (`bidirectional=False` toggles only inter-chunk RNN) and default `gLN` still reads full-utterance statistics, nor does internal OLA chunking solve cross-utterance CSS speaker tracking.

### N03 SepFormer: Dual-Path Self-Attention and Full-Utterance Context

SepFormer ([Subakan et al., ICASSP 2021](https://doi.org/10.1109/ICASSP39728.2021.9413901) / [arXiv §2](https://arxiv.org/pdf/2010.13154); [`speechbrain` `lobes/models/dual_path.py`](https://github.com/speechbrain/speechbrain/blob/89ead74d163463d30c62329a09cfdb4c54f5abc1/speechbrain/lobes/models/dual_path.py) + `recipes/WSJ0Mix/separation/`) replaces RNNs with intra-/inter-chunk Transformers—requiring explicit audits of attention lookahead masks, global normalization, quadratic memory scaling on long recordings, and cross-window permutation stitching.

### N04 TF-GridNet: Spectral, Temporal, and Full-Band Cross-Frame Modeling

Pinned [ESPnet `tfgridnet_separator.py`](https://github.com/espnet/espnet/blob/be79590bb2ff26ffb01bc825c5f68cb9418b7f0d/espnet2/enh/separator/tfgridnet_separator.py) ([Wang et al., ICASSP 2023](https://doi.org/10.1109/ICASSP49357.2023.10094992)) explicitly marks `TFGridNet` as offline, taking `(batch, samples, mics)` with fixed `n_imics`—and divides raw input directly by `torch.std(..., dim=(1, 2))` **without an `eps` zero-variance guard**, creating `NaN`/`Inf` risk on digital silence or constant inputs alongside full-utterance non-causality.

### N05 S4M: Readable Backbone Files Do Not Equal a Turnkey Training Recipe

[`JusperLee/S4M`](https://github.com/JusperLee/S4M) (`MIT` code; [Interspeech 2023 paper](https://doi.org/10.21437/Interspeech.2023-696)) ships architectural files (`S4M.py`, `s4.py`, `s4d.py`, `lssl.py`, `sashimi.py`) **without** dataset preparation scripts, training loops, or checkpoints—and structured state-space layers inside a non-causal encoder/decoder do not automatically yield causal streaming.

### N06 SPMamba and Bidirectional Selective State-Space Modeling

[`JusperLee/SPMamba`](https://github.com/JusperLee/SPMamba/tree/f939f60a10db8a66aa69ec09685684307af47412) (`Apache-2.0` code: `audio_train.py`, `audio_test.py`, `configs/`, `look2hear/`; [arXiv preprint](https://arxiv.org/abs/2404.02063)) swaps **bidirectional** Mamba blocks into a TF-GridNet backbone—linear sequence complexity does NOT imply causal zero-lookahead streaming, and custom `mamba-ssm` CUDA kernels require GPU availability checks.

### N07 Mamba-TasNet and Dual-Path Mamba

[`xi-j/Mamba-TasNet`](https://github.com/xi-j/Mamba-TasNet/tree/a35c692f27213781a11b1606c375cda1e1f0fb62) (`GPL-3.0` license—distinct redistribution obligations from MIT/Apache; `train_wsj0mix.py`, `hparams/WSJ0Mix/`, `modules/`, `inference.ipynb`) links two distinct preprints in its README—[*Speech Slytherin* / Mamba-TasNet](https://arxiv.org/abs/2407.09732) and [*Dual-Path Mamba*](https://arxiv.org/abs/2403.18257)—which must neither be conflated with each other nor share checkpoints with `SPMamba`.

### N08 SpeakerBeam: Enrollment Conditioning and Evaluation-Only Licensing

[`BUTSpeechFIT/speakerbeam`](https://github.com/BUTSpeechFIT/speakerbeam/tree/91af02cc617afa35fedfbdbf32533012cd0a8672) ([*IEEE JSTSP* 2019](https://doi.org/10.1109/JSTSP.2019.2922820); 8 files fetched into git-ignored evaluation cache: `src/models/td_speakerbeam.py`, `adapt_layers.py`, `base_models_informed.py`, `system.py`, `src/datasets/librimix_informed.py`, `egs/libri2mix/`—note `LICENSE.txt` is a **restrictive internal-evaluation license** rather than open source, and its appendix still lists legacy filename `convtasnet_informed.py` renamed to `td_speakerbeam.py` in the tree) conditions `TDConvNetInformed` on an auxiliary enrollment utterance—requiring disjoint enrollment vs. test utterances plus wrong-speaker, cross-mic, and absent-target failure checks.

### N09 CSS: Cross-Chunk Permutation Stitching and Long-Form State

Continuous Speech Separation ([LibriCSS paper](https://arxiv.org/abs/2001.11482); [`microsoft/NOTSOFAR1-Challenge`](https://github.com/microsoft/NOTSOFAR1-Challenge) `MIT` code / `CC-BY-4.0` data, where [`run_inference.py`](https://github.com/microsoft/NOTSOFAR1-Challenge/blob/6f58e08b008f7530ba4141f0aeb02447c70b6fd7/run_inference.py) triggers full dev-set/model downloads even in debug single-session mode, and legacy `css/css_with_conformer/separate.py::Separator` defines `device` while its own `__main__` passes `device_id`) stitches overlapped separation windows across arbitrarily long meetings:
- As proved in Exercise E08-06, two consecutive chunks each scoring `20` dB under local PIT drop to **`−3.69` dB** when stitched with swapped slots.
- Per [MeetEval (CHiME-7 2023, §3.2 Eq. (3) & §3.4 Eq. (5))](https://www.isca-archive.org/chime_2023/neumann23_chime.pdf), across four 1-word turns `A:a, B:c, A:b, B:d`, hypothesis streams `H0 = "a d", H1 = "c b"` score **`ORC-WER = 0/4`** yet **`cpWER = 2/4`** under optimal global speaker bijection.

### N10 SGMSE+, StoRM, and Iterative Score-Based Diffusion

Score-based generative enhancement ([`sp-uhh/sgmse`](https://github.com/sp-uhh/sgmse) commit `1961cf...` & [`sp-uhh/storm`](https://github.com/sp-uhh/storm) commit `257e96...`; commits and `MIT` code licenses in [`THIRD_PARTY_EN.md`](../THIRD_PARTY_EN.md) & [`SOURCES.lock.json`](../SOURCES.lock.json); [Appendix B bridge](04_source_reproduction.md#appendix-b-reproduction)) requires multi-step predictor-corrector/ODE sampling—where static audit of `sgmse/model.py` lines `426–467` vs. `720–777` proves `ScoreModel.enhance()` reads `self.sde.sampler_type` on OUVE/SBVE paths (ignoring its same-named public argument), whereas `StochasticRegenerationModel.enhance()` dispatches via its `sampler_type` parameter, and **both divide by whole-utterance peak magnitude without a zero-peak guard** (`NaN` risk on digital silence).

### N11 ArrayDPS and Likelihood-Guided Diffusion Separation

[`ArrayDPS/ArrayDPS` (`750ac2b7c75458f4ca5bad203dafda528f575e55`, `MIT` code)](https://github.com/ArrayDPS/ArrayDPS/tree/750ac2b7c75458f4ca5bad203dafda528f575e55) ([author v3 §§2–3 & Appendices B/C](https://arxiv.org/html/2505.05657v3)) couples an `8` kHz single-speaker diffusion prior with per-step FCP relative-RIR guidance initialized by `2048/256` `100`-iter Gaussian IVA (`separate.py` line `414` $\to$ `src/sampler_spatial_v1_reverb_iva_8kHz.py::Sampler`, NOT generic `src/sampler.py`):

| Configuration Layer | FCP FFT / Hop | FCP Relative Power Floor | Past / Future Tap Frames |
|---|---|---|---|
| Paper Appendix C.3 | `512 / 64` | `0.001` | Multi-frame per paper |
| `Sampler.__init__` defaults | `512 / 128` | `0.01` | `0 / 0` |
| `separate.py` CLI defaults | `512 / 128` | `0.01` | `20 / 0` |
| `README.md` separation recipe | `512 / 128` | `0.001` | `13 / 0` |

Crucially, `separate.py` line `161` computes oracle SDR (`max_snr`) against ground-truth `sources[:, 0, :]` and uses it on lines `180–186` to gate adaptive retry budgets (`2/3/4` at thresholds `18/13/10` dB by default, or `5/5/5` at `35/14/10` dB in README)—an evaluation-driver stopping rule kept distinct from blind ML sample selection in §4.2 ([`audit_upstream_separation_contracts.py`](../../ch08/examples/audit_upstream_separation_contracts.py), [`upstream_separation_contracts.json`](../../ch08/reports/upstream_separation_contracts.json)).

### N12 AudioSep: Natural-Language Queries and Mono Downmixing

[`audio-agi/audiosep` (`944583f18b84589dc965de3ad77525c945334252`, `MIT` code)](https://github.com/audio-agi/audiosep/tree/944583f18b84589dc965de3ad77525c945334252) ([*IEEE/ACM TASLP* 2024](https://doi.org/10.1109/TASLP.2024.3520017)) conditions separation on CLAP text embeddings via `pipeline.py`—automatically resampling and **downmixing inputs to `32` kHz mono** (discarding multi-channel spatial phase cues and targeting sound-class descriptions rather than voice-biometric speaker identity).

### N13 Online SpatialNet: Causal Layers vs. Stateless Top-Level `forward()`

In [`Audio-WestlakeU/NBSS` (`cc42fc8ad2e6642c09b8f4169a85b4766dc22b7e`, `MIT` code)](https://github.com/Audio-WestlakeU/NBSS/tree/cc42fc8ad2e6642c09b8f4169a85b4766dc22b7e), [`models/arch/OnlineSpatialNet.py`](https://github.com/Audio-WestlakeU/NBSS/blob/cc42fc8ad2e6642c09b8f4169a85b4766dc22b7e/models/arch/OnlineSpatialNet.py) ([*IEEE SPL* `31:2295–2299`, 2024](https://doi.org/10.1109/LSP.2024.3418714), comparing masked attention, Retention, and `mamba-ssm`/`causal-conv1d` CUDA temporal blocks) includes a bottom-of-file causal prefix test, yet top-level **`OnlineSpatialNet.forward()` passes `state=None` to every sublayer without exposing recurrent/conv caches as inputs or outputs**—demonstrating why *"causal whole-utterance graph"*, *"linear $O(T)$ complexity"*, and *"state-preserving chunk-by-chunk streaming API"* are three distinct properties.

<a id="streamfm"></a>

### N14 Stream.FM: Multi-Stage Flow Matching and Per-Frame Solver States

Targeting causal single-channel speech restoration/dereverberation with few NFE steps ([Stream.FM author v3, 2026-04-21, §§II–III & Supp. S.VIII](https://arxiv.org/html/2512.19442v3); published in [*IEEE TASLPRO* `34:3087–3101`, 2026](https://doi.org/10.1109/TASLPRO.2026.3696215); pinned `AGPL-3.0` commit [`ab2700c1154acc5c2ce67a5344182028336413f5`](https://github.com/sp-uhh/streamfm/tree/ab2700c1154acc5c2ce67a5344182028336413f5) in git-ignored `_downloads/streamfm/`), Stream.FM separates **physical audio frame time $t$** from **ODE flow time $\tau \in [0, 1]$** (training with coupled Gaussian draw $Z$ across $X_0 = Y + \sigma_y Z$ and $X_1 = X + \sigma_x Z$ so velocity target is $X_1 - X_0$, while multi-stage Runge–Kutta in `solve_ode()` invokes the network once per *stage* rather than once per step `N`):

| Pinned Source File | Audited Construct | Streaming Integration Consequence |
|---|---|---|
| `config/streamfm_derev.yaml` | `FlowModel`, `16` kHz mono, `512/256` root-Hann STFT, `256` bins (Nyquist dropped & zero-restored) | `32` ms window / `16` ms hop framing baseline before compute/DMA |
| `sgmse/model.py::FlowModel` | Coupled noise $Z$, velocity target, `enhance_from_features()` | Top-level `enhance_from_features()` runs batch `solve_ode()` offline |
| `sgmse/util/solvers.py` | `get_butcher_tableau()`, `_solver_step_rk()`, `solve_ode()` | Total NFE = `steps × stages_per_step` |
| `sgmse/backbones/streaming_unet.py` | `CausalNCSNpp` & sublayer `init_state()` / `forward_step()` | Per-stage state buffers must never overwrite each other across $\tau$ stages |
| `fit_rk_scheme.py` & `LRK5_streamfm_derev.yaml` | Learned Runge–Kutta tableau optimization | Standalone solver-fitting stage (unexecuted) |

**Executed static AST audit (`streamfm_source_audit.py` $\to$ [`streamfm_source_audit.json`](../../ch08/reports/streamfm_source_audit.json)) uncovers two direct streaming-method mismatches in untouched upstream code:**
1. **`CausalResnetBlockBigGANpp` state tuple length mismatch**: `init_state()` (line `672`) returns a **5-tuple**, whereas `forward_step()` (line `683`) unpacks `state` into **6 variables**.
2. **`CausalConv2d.forward_step()` missing attributes**: lines `899–900` read `self.depthwise_separable` and `self.pointwise_conv`, neither of which is assigned in `CausalConv2d.__init__` (they belong to sibling class `CausalDecoupledConv2d`), while conv ring buffers slice `x[0, :, :, 0]` assuming `batch == 1, time == 1` (`re-run: .venv/bin/python -m codes.chapters.ch08.examples.streamfm_source_audit`).

Next-step experimental verification (unexecuted locally; zero AGPL code copied into our repo) must resolve those state contracts first, verify batch-causal vs. frame-by-frame equivalence on identical random $Z$, check impulse/silence/tail draining and future-perturbation invariance per RK stage, and measure full-stage device latency rather than single-forward GPU throughput.

<a id="tflocoformer"></a>

### N15 TF-Locoformer: Local Conv-SwiGLU, Global Attention, and Complex Spectral Interfaces

Replacing RNNs in TF-GridNet with local `SwiGLUConvDeconv1d` modules paired with rotary/NoPE self-attention ([Saijo et al., IWAENC 2024 §2, Fig. 1, Eqs. (1)–(8)](https://arxiv.org/pdf/2408.03440); pinned `Apache-2.0` MERL commit [`7a615460d347ff7334a13dbb831d16280da72cdc`](https://github.com/merlresearch/tf-locoformer/tree/7a615460d347ff7334a13dbb831d16280da72cdc) audited 2026-09-26), TF-Locoformer maps complex spectra directly to $N$ separated complex spectra:

| Inspected Entry Point | Actual Tensor I/O Contract | Key Verification Items |
|---|---|---|
| `standalone/tflocoformer_separator.py::TFLocoformerSeparator` | Complex `[B, T, F] -> [B, N, T, F]` | Pure spectral module (no internal STFT/iSTFT); `B` is batch, NOT mics |
| `TFLocoformerBlock` / `LocoformerBlock` | Alternating frequency-path & time-path modeling | Macaron FFN residual scaling (`1.0` in code vs. `1/2` in paper Eqs. (2) & (4)), `RMSGroupNorm` |
| `espnet2/enh/separator/tflocoformer_separator.py` | Complex spectrum $\to$ list of complex spectra | Docstring says waveform, actual `forward` returns complex spectra |
| `egs2/whamr/enh1/conf/tuning/train_enh_tflocoformer.yaml` | `8` kHz, `256/64` STFT (`32/8` ms), `6` layers, dim `128`, `4` heads/groups, FFN `192`, kernel `8` | Differs from paper Table 1 dimensions; WHAMR recipe needs 2-src dereverb modifications |
| `tests/test_tflocoformer.py` & `..._load_pretrained_weights.py` | 3D `[B, T, F]` with explicit `norm_type="rmsgroupnorm"` | Leaves default `norm_type` and 4D input branch untested |

Tracing dimensions on $B=2, T=50, F=129$, dim `128`, $N=2$: 3D complex input stacks real/imag to `[2, 2, 50, 129]` $\to$ `3×3` conv `[2, 128, 50, 129]` $\to$ utterance-wide `GroupNorm(1, 128)` $\to$ `tf_order="ft"` frequency attention on `[100, 129, 128]` ($B H T F^2$ score elements) & unmasked non-causal time attention on `[258, 50, 128]` ($B H F T^2$ score elements, non-streaming even under `pos_enc="nope"`) $\to$ deconv `[2, 2, 2, 50, 129]` $\to$ complex `[2, 2, 50, 129]`.

**Executed static AST audit ([`tflocoformer_source_audit.py`](../../ch08/examples/tflocoformer_source_audit.py) $\to$ [`tflocoformer_source_audit.json`](../../ch08/reports/tflocoformer_source_audit.json), [`test_codes_tflocoformer_source_audit.py`](../../../../tests/test_codes_tflocoformer_source_audit.py)) catches two concrete bugs in both standalone and ESPnet wrappers:**
1. **Typo in default `norm_type` argument**: `__init__` (line `65` / `73`) sets default `norm_type="rmsgrouporm"` (missing `'n'`), whereas `LocoformerBlock` (line `296` / `320`) asserts `norm_type in ("layernorm", "rmsgroupnorm")`.
2. **Broken 4D input branch**: asserting `input.shape[1] == 1` followed by `transpose(1, 2)` rejects `[B, T, 1, F]` (`T != 1`) while turning `[B, 1, T, F]` into `[B, 2*T, 1, F]`, failing the `2`-channel encoder conv whenever $T \ne 1$. (Note also `separate.py` runs `mix /= mix.std(axis=-1)` without zero-variance protection on `cuda:0`; re-run read-only via `.venv/bin/python -m codes.chapters.ch08.examples.tflocoformer_source_audit`.)

### N16 Candidate Screening and Literature Audit Boundaries

Audited on 2026-09-28: **mfFCA 2023** is adopted in [B11](#mffca); **SAM Audio 2025** ([paper](https://arxiv.org/pdf/2512.18099), [`facebookresearch/sam-audio`](https://github.com/facebookresearch/sam-audio) under custom `SAM License` with text/visual/span prompts and multi-repo dependencies) and **FLASeformer 2025** ([Interspeech 2025](https://www.isca-archive.org/interspeech_2025/wang25j_interspeech.pdf), linear attention without verified pinned author repo) are retained as literature pointers without local execution (noting also that Ono 2011 and Kitamura 2016 DOIs were verified while formulas were cross-checked directly against pinned `ssspy`/`pyroomacoustics` code and independent algebra).

### N17 DiCoW v1: Diarization-Conditioned Target-Speaker ASR

[DiCoW (*Computer Speech & Language*, Vol. 95, 2026)](https://doi.org/10.1016/j.csl.2025.101841) performs **target-speaker ASR** (outputting *text transcripts* rather than separated waveforms), conditioned on 4-class frame-level diarization masks (`silence, target, non-target, overlap`) injected via `TargetSpeakerAmplifier` in [`src/models/whisper_ctc.py`](https://github.com/BUTSpeechFIT/TS-ASR-Whisper/blob/0ea6679d44405f5ff39188030123524686c198e9/src/models/whisper_ctc.py) of [`TS-ASR-Whisper` v1 `0ea6679`](https://github.com/BUTSpeechFIT/TS-ASR-Whisper/tree/0ea6679d44405f5ff39188030123524686c198e9) ([`src/main.py`](https://github.com/BUTSpeechFIT/TS-ASR-Whisper/blob/0ea6679d44405f5ff39188030123524686c198e9/src/main.py), [`src/train.py`](https://github.com/BUTSpeechFIT/TS-ASR-Whisper/blob/0ea6679d44405f5ff39188030123524686c198e9/src/train.py)) and [`DiCoW` v1 `e9326bd`](https://github.com/BUTSpeechFIT/DiCoW/tree/e9326bd536bf632e823357438b210102903ba620)—where [`example.py`](https://github.com/BUTSpeechFIT/DiCoW/blob/e9326bd536bf632e823357438b210102903ba620/example.py) ([`modeling_dicow.py`](https://github.com/BUTSpeechFIT/DiCoW/blob/e9326bd536bf632e823357438b210102903ba620/modeling_dicow.py)) loads `BUT-FIT/DiCoW_v1` whereas [`app.py`](https://github.com/BUTSpeechFIT/DiCoW/blob/e9326bd536bf632e823357438b210102903ba620/app.py) in that exact same commit loads `BUT-FIT/DiCoW_v2`. Keep code (`Apache-2.0`), [DiCoW v1](https://huggingface.co/BUT-FIT/DiCoW_v1)/[v2](https://huggingface.co/BUT-FIT/DiCoW_v2) model cards (`CC BY 4.0`), gated Hugging Face access on [`pyannote/speaker-diarization-3.1`](https://huggingface.co/pyannote/speaker-diarization-3.1) & [`segmentation-3.0`](https://huggingface.co/pyannote/segmentation-3.0) (`MIT` + user gate), and later non-commercial `CC BY-NC 4.0` [`DiariZen` weights](https://github.com/BUTSpeechFIT/DiariZen/blob/main/MODEL_LICENSE) ([DiCoW license note](https://github.com/BUTSpeechFIT/DiCoW#license)) strictly distinct.

### N18 FlowSep: Text-Queried Latent Flow Matching for Sound Events

[FlowSep (ICASSP 2025)](https://doi.org/10.1109/ICASSP49660.2025.10890129) ([`Audio-AGI/FlowSep` `d8164db58bd461ef5bb6df8ffd372b536ee6afb4`](https://github.com/Audio-AGI/FlowSep/tree/d8164db58bd461ef5bb6df8ffd372b536ee6afb4): [`lass_inference.py`](https://github.com/Audio-AGI/FlowSep/blob/d8164db58bd461ef5bb6df8ffd372b536ee6afb4/lass_inference.py), [`lass_config/2channel_flow.yaml`](https://github.com/Audio-AGI/FlowSep/blob/d8164db58bd461ef5bb6df8ffd372b536ee6afb4/lass_config/2channel_flow.yaml), external [Zenodo checkpoint](https://zenodo.org/records/13869712) and [AudioCaps](https://audiocaps.github.io/) VAE/BigVGAN dependencies) targets text-queried general audio separation rather than speaker identity; lacking any `LICENSE`/`COPYING` file in the pinned tree as of 2026-09-28, it is kept strictly `index_only` (`fetch_enabled=false`) with zero files copied locally.


<a id="meco-candidate"></a>

### N19 MeCo: One-Step Mean-Flow Post-Separation Correction

[MeCo (2026 preprint v1 §§2–4)](https://arxiv.org/html/2606.09677v1) applies a `1-NFE` mean-velocity displacement step conditioned on preliminary separated estimates plus the multi-channel mixture. On 2026-10-04, 16 MIT-licensed core/config/entry plain-text files (`95305` B, root MIT SHA `acf88dd4e67d3f5d50a9ec79e6f035339ae1bc072afffed56a73387d1fb1679b`) were fetched from [`rlaehghks5/MECO` `375ac4dee2a8e193aaa48553289499a6f43d93e6`](https://github.com/rlaehghks5/MECO/tree/375ac4dee2a8e193aaa48553289499a6f43d93e6)—deliberately excluding GoogleResearch `Apache-2.0` helpers and two `NVIDIA Source Code License-NC` CUDA kernels:
- [`flowmse/odes.py`](https://github.com/rlaehghks5/MECO/blob/375ac4dee2a8e193aaa48553289499a6f43d93e6/flowmse/odes.py) and [`odesolvers.py`](https://github.com/rlaehghks5/MECO/blob/375ac4dee2a8e193aaa48553289499a6f43d93e6/flowmse/sampling/odesolvers.py) initialize from `estimate + 0.487 * z` (complex Gaussian noise) and step dimensionless flow time $1.0 \to 0.03$ ($\Delta\tau = 0.97$, NOT `30` ms), while [`model_MeCo.py`](https://github.com/rlaehghks5/MECO/blob/375ac4dee2a8e193aaa48553289499a6f43d93e6/flowmse/model_MeCo.py) uses central finite-difference JVP, gradient clipping, and interval-weighted `xr` loss:

| Audited Condition | Static Fact in Pinned Source ([`data_module_multichannel.py`](https://github.com/rlaehghks5/MECO/blob/375ac4dee2a8e193aaa48553289499a6f43d93e6/flowmse/data_module_multichannel.py), [`evaluate.py`](https://github.com/rlaehghks5/MECO/blob/375ac4dee2a8e193aaa48553289499a6f43d93e6/evaluate.py)) | Invalid Assumption to Avoid |
|---|---|---|
| Input channel count | Multi-channel backbone hard-codes `6` complex / `12` real channels (`1` state + `1` est + `4` mix) | Assuming `forward` accepts mono or arbitrary mic counts out-of-the-box |
| STFT & training crop | `16` kHz, `510` FFT, `64` hop, `center=True`, `256` frames = `16320` samples (`1.02` s) | Equating `4` ms hop to system latency or paper's `4` s training crop |
| Evaluation framing | Utterance-peak normalization & padding to multiples of `64` frames | Treating `1 NFE` as causal frame-streaming |
| RNG determinism | `evaluate.py` declares `--seed 1234` **without ever consuming it** | Assuming `--seed` makes `0.487 * z` deterministic |
| Silence protection | Training divides by `max(|est|)` with zero floor; `evaluate.py` adds `1e-12` | Treating training and eval silence guards as identical |
| Upstream pairing | Expects pre-computed `(clean, estimate, mixture)` triplets on disk | Expecting MeCo to run upstream BSS or speaker assignment |
| Proxy metrics | `evaluate.py` computes `SI-SIR`/`SAR` using `est - clean` as noise proxy | Mistaking `est - clean` proxy residual for true multi-source BSS Eval decomposition |

<a id="sis-candidate"></a>

### N20 SIS: Speaker Identity Supervision During Training Only

[MERL Technical Report TR2026-137 §§2–3](https://www.merl.com/publications/docs/TR2026-137.pdf) uses an auxiliary utterance from the same speaker solely as a **training-time** contrastive supervision signal alongside PIT (requiring zero enrollment clue at inference time, unlike SpeakerBeam TSE; cf. mixture-consistency vs. mask trade-offs in [Exercises E08-20 & E08-22](../../../../chapters_en/08_speech-separation.md#e08-22)); inspecting [`merlresearch/sis_sep`](https://github.com/merlresearch/sis_sep) on 2026-10-04 found only a placeholder announcement `README.md` with zero code or license.

<a id="target-extraction-diarization"></a>

### N21 Extended Topic IV: Conditional Target Extraction, Multi-Channel Diarization, and Far-Field Multimodal Front Ends

Corresponding to [Chapter 15 (Extended Topic IV)](../../../../chapters_en/15_target-extraction-diarization.md) and [Exercises & Audio §62](05_exercises_and_audio.md#tse-exercises-audio), whereas blind separation emits unattributed output slots, Target Speaker Extraction (TSE) and Target-Speaker VAD (TS-VAD) anchor speaker identity via external enrollment/spatial/visual clues—demanding explicit verification of silent-target SI-SDR singularity, Fresnel near-field spherical wavefront curvature identifiability, overlap diarization leakage, and asynchronous audio-visual shutter alignment across five self-contained NumPy kernels:
- [`conditional_tse.py`](../../ch15/core/conditional_tse.py): Feature-wise Linear Modulation (FiLM), zero-mean SI-SDR singularity on silent targets, dual-regime target-present/absent energy-floor loss, and single-speaker zero-overlap fricative-preservation bypass gating.
- [`distance_bubble.py`](../../ch15/core/distance_bubble.py): Fresnel spherical-wave steering vectors, broadside co-directional near-field target vs. far-field interferer LCMV "distance bubble" beamforming, and white-noise-gain (WNG) collapse at small apertures.
- [`diarization_tsvad.py`](../../ch15/core/diarization_tsvad.py): Theoretical missed-speech DER lower bound of single-label frame clustering on overlapped meetings, joint speaker-embedding cosine + spatial prior TS-VAD posterior updates, and exclusive soft-gated GSS beamforming.
- [`av_alignment.py`](../../ch15/core/av_alignment.py): Asynchronous video frame-rate and camera shutter-exposure integration alignment, plus confidence-gated acoustic fallback under lip occlusion.
- [`eend_enrollment_robustness.py`](../../ch15/core/eend_enrollment_robustness.py): Fixed-slot-cap DER lower bounds under multi-talker overlap, EEND-EDA attractor existence gating, and crosstalk-contaminated enrollment purification via weighted exclusive pooling.

Upstream reference [`wesep` (`c2437a64d314baa3b7a49421d0030c608adf5a32`, `Apache-2.0`)](https://github.com/wenet-e2e/wesep/tree/c2437a64d314baa3b7a49421d0030c608adf5a32) in [`SOURCES.lock.json`](../SOURCES.lock.json) has 6 plain-text files fetched (`runtime/separate/separate_engine.h`, `examples/spatial/mc_librimix/local/build_spatial_cues.py`, `examples/audio/librimix/local/build_audio_cues.py`, etc.; `source_verified`, `execution="not_run"`, zero checkpoints/LibriMix data downloaded or C++/PyTorch inference compiled).

## 6. Common Logging Checklist for Reproduction Experiments

| Logging Category | Mandatory Fields to Record | Ambiguity If Omitted |
|---|---|---|
| Input | File SHA-256, sample rate, mic/source counts, reference mic, target signal definition | Whether two runs solve the same task |
| Framing & Time | STFT window/hop, center padding, reference bulk delay, lookahead, chunk stitching | Algorithmic distortion vs. bulk shift/truncation |
| State | Initialization, cross-chunk persistence, freeze gates, resets, fault recovery | Suitability for continuous streaming |
| Numerics | Precision (`fp32`/`fp64`/`complex128`), floors, diagonal loading, stopping rules, RNG seed | Ill-conditioned matrix breakdown vs. seed variance |
| Model | Upstream commit SHA, config YAML, checkpoint SHA-256, utterance/running normalization | Exact network reproducibility |
| Metrics | Output tap point, active scoring intervals, delay alignment, PIT/ORC/cpWER rule, input baseline | Cross-table metric comparability |
| Compute & RAM | Target CPU/GPU/DSP, thread count, batch size, timer boundaries, peak resident memory | Real-time hardware budget compliance |
| Licensing | Separate licenses for code, checkpoints, speech corpora, and third-party subcomponents | Redistribution and commercial deployment rights |

Scale experimental complexity in four deliberate stages: `analytical toy arrays -> synthetic regression -> pinned public recordings -> long-form target-hardware capture`. Synthetic checks prove exact adherence to mathematical equations, public corpora enable cross-lab reproducibility, and physical hardware loops validate reference taps, crystal drift, OS drivers, loudspeaker nonlinearities, and human interaction—none of the three tiers substitutes for the others.
