# AGENTS_EN.md — Long-Term Engineering & Authoring Guidelines

🌐 [简体中文](AGENTS.md) ｜ **English**

> Read this document in full before modifying anything in this repository. It governs both technical/pedagogical quality and how investigations, edits, verifications, and handoffs are conducted.
> Whenever rules appear to conflict, resolve them strictly according to the priority hierarchy in Section 2.1 rather than applying lower-priority clauses mechanically.

## 1. Project Mission and Repository File Map

### 1.1 Target Audience and Pedagogical Promise

This repository is an open textbook and runnable engineering laboratory on **Microphone Array Signal Processing: From Fundamentals to Industrial Practice**, serving beginners, academic researchers, and audio systems engineers (available in both Chinese and English editions). Starting from mathematical and signal-analysis foundations (Chapter 0) and why multi-microphone arrays unlock spatial degrees of freedom (Chapter 1), the book progresses through sound source localization (`DOA`), beamforming, acoustic echo cancellation (`AEC`), blind dereverberation (`WPE`), speech separation, and multi-target source tracking, then bridges into real-time hardware deployment, objective/subjective evaluation, Pareto system selection, and four frontier chapters on acoustic camera imaging, bandwidth-constrained distributed wireless arrays, spatial/wearable arrays, and conditioned target extraction with multichannel speaker diarization.

Every core chapter delivers four things simultaneously:

1. **Exact definitions paired with physical intuition**;
2. **Zero-gap mathematical derivations** that never skip sign, conjugation, dimension, or normalization steps;
3. **Hand-verifiable numerical examples** backed by pure-NumPy code;
4. **Explicit operating assumptions, failure boundaries, and engineering selection criteria**.

Readers do not need prior background in array signal processing. Self-contained mathematical derivations and 16 foundational exercises are consolidated in Chapter 0 (`00_math-foundations.md`), while the master notation table, bilingual glossary, and compact formula lookup live in Appendix A (`16_appendix-symbols-math.md`); subsequent chapters include just the minimal prerequisite setup needed on the spot and link directly to the corresponding foundational sections.

### 1.2 Authoritative Source Files vs. Generated Artifacts

Before editing any file, determine whether it is an **authoritative source file** or a **generated artifact**. Whenever a generated artifact exhibits an issue, fix its authoritative source file or build script and regenerate—never hand-patch compiled outputs, PNGs, WAVs, or JSON manifests.

| Category | Path | Role | Modification Rule |
|---|---|---|---|
| Tutorial chapter sources | `chapters/*.md` (Chinese), `chapters_en/*.md` (English) | Overview, 16 main/frontier chapters (Ch 0–15), and 2 appendices | Edit directly; verify cross-chapter links and heading parity after editing |
| Pedagogical Python sources | `core/`, `examples/`, and chapter runners under `codes/chapters/{ch00,ch01..ch15,appendix_a,appendix_b}/` | Each shared numerical kernel lives once in the chapter where it is first derived; Ch 0 math kernels/exercises and cross-chapter runners live in `ch00/` | Edit the single authoritative implementation; synchronize unit tests, `COVERAGE*.md`, and referencing chapters (legacy `codes.array_tutorial` and `codes.examples` paths are retired) |
| Ch 0 math kernels & figures | `codes/chapters/ch00/core/chapter00_math.py`, `codes/chapters/ch00/chapter00_exercises.py`; `scripts/make_math_figures.py` → Figs 84–87 & `ch00/reports/figure84`–`87*.json` | 16 foundational math exercises (`E00-01`–`E00-16`) and 4 numbered figures | Binds 4 authoritative sources with fixed PRNG seeds; covers spherical geometry, $\mathrm{SO}(3)$, phasors/COLA, subspaces/SVD/Woodbury, Rayleigh/GEV/Khatri–Rao, Wirtinger/NNLS/MM, complex Gaussians/CRLB/Wiener/Kalman, Hungarian assignment, and Bessel/Legendre/spherical harmonics; generates zero WAVs |
| Shared file I/O primitives | `codes/chapters/ch00/io_contracts.py` | Regular parent-chain & exact-membership validation, strict finite JSON parsing, metadata typing, and atomic temp-file report writing | Shared across audio generators and upstream fetchers; each caller owns its schema, scoring, and provenance rules; binds SHA-256 when participating in generation (preflight checks do not claim immunity to concurrent TOCTOU races or power-loss crashes) |
| Language switch line | `scripts/language_bars.py` → one `🌐` line at the top of each of 68 Chinese/English Markdown pair members (`chapters/`↔`chapters_en/`, `research/`↔`research_en/`, root and subdirectory `*.md`↔`*_EN.md`) | GitHub cannot run scripts inside Markdown, so every document opens with one link to its counterpart; the line is repository navigation only and `build_markdown_helpers.shield_math` strips it before the website and PDF are rendered (they have their own language pills) | After adding or renaming a paired document run `python scripts/language_bars.py`; `quality_check_en.py` re-verifies read-only through the same source; never hand-edit the line |
| Main 109-WAV audio layout | `codes/chapters/ch00/audio/MANIFEST.json`, per-chapter `audio/*.wav`; `scripts/code_layout.py` group-to-chapter map | 27 groups, 109 mathematically synthesized WAVs under a single manifest partitioned across chapters | Regenerate from updated sources; manifest binds 20 source files and logs owning chapter, source SHA-256, and WAV SHA-256; never hand-edit individual WAVs |
| Known noise-weight audio | `codes/chapters/appendix_a/core/weighted_audio.py`, `examples/generate_weighted_audio.py` → `appendix_a/weighted_audio/` | 5 `16 kHz`, `32,000`-sample mono/stereo synthetic WAVs + standalone manifest | Shared `700 Hz` target & `3500`/`4000 Hz` interferers, known fixed weights summing to 1; separates analytic, float, and readback PCM integer metrics over `28,800`-sample steady window; 4 bound sources, gain 1, strict 6 regular members, `--check` read-only replay (not blind noise estimation or stochastic independence proof) |
| Equal-DRR spectral audio | `codes/chapters/appendix_b/core/response_audio.py`, `examples/generate_response_audio.py` → `appendix_b/response_audio/` | 5 `16 kHz`, `32,002`-sample mono synthetic WAVs + standalone manifest | Shared `2000`/`4000 Hz` source across two short FIRs with identical reflection energy; gain 1, full 2-sample tail; separates analytic, float, and integer PCM over `28,800`-sample window; 4 bound sources, 6 regular members, `--check` read-only replay |
| Appendix A upstream solver contract | `codes/chapters/appendix_a/examples/audit_upstream_solver_contracts.py` → `appendix_a/reports/upstream_solver_contracts_current.json` | Compiles & executes pinned `pb_bss` solver module from verified source bytes alongside standalone NumPy controls | Retains 4 legacy cases, 2 NumPy cases, and 4 shape/complex-batch controls; separates lock/status/sparse policy from touched blobs; writes in-repo only to `_current.json` (preserving legacy report bytes); defaults to `stdout` without touching `.pyc` caches |
| Appendix B TAC source contract | `codes/chapters/appendix_b/examples/audit_tac_contracts.py` → `appendix_b/reports/tac_contracts_current.json` | 10 AST structural invariants across 4 pinned author TAC files with zero PyTorch execution | Separates lock/status/sparse policy from verified blobs across 4 bound runner sources; defaults to `stdout`, writes in-repo only to `_current.json` without touching legacy `tac_contracts.json` |
| Fig 60 numerical report | `scripts/make_figures.py::fig_weighted_noise` → `appendix_a/reports/figure60_weighted_noise.json` | Three fixed weight sets, analytic MSE/NMSE, and three readback PCM integer `E/D` ratios | Regenerates with figure over identical window/reference without fitting gain/delay; spectral-line orthogonality holds strictly on the declared integer-period window |
| Fig 61 numerical report | `scripts/make_figures.py::fig_equal_drr_response` → `appendix_b/reports/figure61_equal_drr_response.json` | Short RIR taps, two-tone analytic power gains, and two readback PCM integer `E/D` ratios | Regenerates with figure; names RIR DRR separately from source-weighted reflection power |
| STFT convolution audio | `codes/chapters/ch02/core/stft_convolution.py`, `examples/generate_stft_convolution.py` → `ch02/stft_audio/` | 3 shared-source, full-tail, shared-gain synthetic WAVs + standalone manifest | Separates float vs. PCM; compares round-trip against dry source and framewise MTF against full linear convolution (`--check` read-only, separate from main 109 WAVs) |
| Digital sweep deconvolution audio | `codes/chapters/ch02/core/deconvolution.py`, `examples/generate_sweep_audio.py` → `ch02/sweep_audio/` | 4 `16 kHz` mono synthetic WAVs (`32,000`-pt excitation, three `32,160`-pt responses) + manifest | Shared gain 1, full 160-pt tail; retains full `65,536`-pt IR, `161`-pt support, and out-of-support leakage in float & PCM; 5 bound sources, 5 regular members, `--check` read-only |
| Ch 2 upstream model contract | `codes/chapters/ch02/examples/audit_upstream_models.py` → explicit `reports/upstream_model_contracts.json` | Bounded calls on pinned `PRA`/`Acoular`/`doatools`/`pyfar` + doc/license check on `speed-of-sound-in-air` | Defaults to `stdout`; separates used blob identity from worktree sparse status without overwriting legacy `upstream_models.json` |
| Coherent reflection audio | `codes/chapters/ch04/core/reflection_audio.py`, `examples/generate_reflection_audio.py` → `ch04/reflection_audio/` | 4 `16 kHz`, `32,012`-sample mono/stereo synthetic WAVs + manifest | Same-source `4 kHz` direct (`8`/`8` samples) + specular reflection (`12`/`11` samples), full 12-pt tail, gain 1; `27,200`-pt steady window separates analytic (`0.1I` virtual noise), float, and noise-free PCM `E/D`; 5 sources, 5 members, `--check` read-only |
| Known unitary focusing audio | `codes/chapters/ch04/core/focus_audio.py`, `examples/generate_focus_audio.py` → `ch04/focus_audio/` | 4 `16 kHz`, `32,024`-sample mono/4-ch synthetic WAVs + manifest | Known `1`/`3 kHz` unitary permutation over `[2400, 29600)` window; 6 bound sources, `/32768` PCM decode (`--check` read-only; not blind CSSM) |
| Distributed known-stats audio | `codes/chapters/ch13/core/distributed.py`, `core/distributed_audio.py`, `examples/generate_distributed_audio.py` → `ch13/distributed_audio/` | 17 `16 kHz`, `32,000`-sample mono/4-ch synthetic WAVs + manifest | 9 bound sources, gain 1; separates `28,800`-pt steady window and `800`-pt packet-drop window across analytic, float, and PCM `E/D`; caches effective receiver weights post-broadcast; 18 regular members, `--check` read-only |
| Ch 4 upstream DOA & SAID contracts | `ch04/examples/audit_upstream_doa.py`, `audit_said_compression.py` → explicit `upstream_doa_contracts.json`, `said_compression_contracts.json` | Bounded `PRA`/`doatools` DOA probes and single-file SAID compression audit | Defaults to `stdout`; preserves upstream bugs & `not_run` statuses without overwriting legacy reports; 4 legacy `reproduce_*` scripts write only outside repo |
| Distributed figures & contracts | `scripts/make_figures.py::fig_distributed_*` → `ch13/reports/figure65`–`67*.json`; `audit_upstream_distributed_contracts.py` → `upstream_distributed_contracts_current.json` | Fixed steering, post-broadcast updates, SRO/dropouts, static MATLAB & 3 `paderwasn` AST helpers | Logs solver calls separately from effective weights; MATLAB unexecuted, DWACD statically checked; preserves MIT/GPL conflict & `not_run` boundaries |
| Wearable spatial array audio & figs | `codes/chapters/ch14/core/wearable_audio.py`, `examples/generate_wearable_audio.py` → `ch14/wearable_audio/`; `scripts/make_wearable_figures.py` → Figs 76–79 & `ch14/reports/figure76`–`79*.json` | 6 `16 kHz`, `32,000`-sample mono/stereo synthetic WAVs + manifest | Gain 1, `[1600, 30400)` (`28,800`-pt) analytic/float/integer PCM; 6 bound sources, 7 regular members, `--check` read-only (not measured ear-canal/wind-tunnel recordings) |
| Conditioned TSE & distance-bubble audio & figs | `codes/chapters/ch15/core/tse_audio.py`, `examples/generate_tse_audio.py` → `ch15/tse_audio/`; `scripts/make_tse_figures.py` → Figs 80–83 & `ch15/reports/figure80`–`83*.json` | 6 `16 kHz`, `32,000`-sample mono/3-ch synthetic WAVs + manifest | Gain 1, `[1600, 30400)` (`28,800`-pt) analytic/float/integer PCM `E/D`; 6 bound sources, 7 regular members, `--check` read-only |
| Acoustic imaging snapshot audio | `codes/chapters/ch12/core/imaging_audio.py`, `examples/generate_imaging_audio.py` → `ch12/imaging_audio/` | 5 `24 kHz`, `48,004`-sample mono/stereo synthetic WAVs + manifest | Gain 1, 4-sample tail; 20 snapshots × `[252, 2160)` (`38,160` pts) separating analytic, float & PCM CSMs; 6 regular members, `--check` read-only |
| Imaging figures & Acoular contract | `scripts/make_figures.py::fig_imaging_*` → `ch12/reports/figure62`–`64*.json`, `figure75_imaging_objectives.json`; `audit_upstream_imaging_contracts.py` → `upstream_imaging_contracts_current.json` | 2-grid, 441-grid spherical PSF/DAMAS/CLEAN-SC, objective controls & pinned Acoular probes | Preserves Acoular `r_diag=False` CLEAN-SC over-allocation, CMF half-triangle scaling & endpoint conventions without modifying upstream or legacy reports |
| Author DAMAS Python contract | `ch12/examples/audit_damas_author_contracts.py` → explicit `damas_author_contracts_current.json` | 9 bounded Python controls on author GPLv3 `damas.py` (8 matches, 1 zero-CSM `NameError`) | Defaults to `stdout`; verifies upstream file & 3 direct helpers before/after without running MATLAB/MEX or overwriting legacy reports |
| Derivative-constrained audio | `codes/chapters/ch05/core/derivative_audio.py`, `examples/generate_derivative_audio.py` → `ch05/derivative_audio/` | 4 `16 kHz`, `32,002`-sample mono/3-ch synthetic WAVs + manifest | 1-sample causal reference, 8 bound sources, fixed seed, gain 1; separates analytic, float & PCM over `[2400, 29600)` (`--check` read-only) |
| Binaural beamforming audio | `codes/chapters/ch05/examples/binaural_beamforming_demo.py` → `ch05/binaural_beamforming_audio/` | 4 `16 kHz`, `32,000`-sample/ear stereo synthetic PCM16 WAVs + manifest | Gain 1, `27,200`-pt steady window, 6 bound sources, 5 regular members (`--write` to generate, `--check` read-only) |
| Per-bin phase reference audio & Fig 68 | `codes/chapters/ch05/core/phase_audio.py`, `examples/generate_phase_audio.py` → `ch05/phase_audio/`; `scripts/make_beamforming_figures.py` → Fig 68 & JSON | 3 `16 kHz`, `32,000`-sample mono two-tone phase controls + manifest | Gain 1, `[2400, 29600)` (`27,200` pts), 5 bound sources, 4 regular members (`--check` read-only) |
| Affine Projection (`APA`) audio | `codes/chapters/ch06/core/apa_audio.py`, `examples/generate_apa_audio.py` → `ch06/apa_audio/` | 6 `16 kHz`, `32,013`-sample mono synthetic WAVs + manifest | Colored AR(1) reference / 16-tap FIR, frozen after `[0, 24000)` training, scored on `8,000`-pt holdout (13-pt tail excluded from score); 8 sources, gain 1 (`--check` read-only) |
| Playback gain & tail-echo audio & Fig 69 | `ch06/core/reference_audio.py`, `examples/generate_reference_audio.py` → `ch06/reference_audio/`; `scripts/make_reference_figures.py` → Fig 69 & JSON | 6 `16 kHz`, `33,600`-sample mono synthetic WAVs + manifest | Frozen 2-tap FIR, gain step $1\to 2$ & 1,600-sample post-stop tail; separates step/steady/tail windows across 8 bound sources, 7 regular members (`--check` read-only) |
| Ch 6 upstream AEC contracts | `ch06/examples/audit_aec_upstream_interfaces.py`, `audit_upstream_apa.py` → `reports/*_current.json` | Pinned `pyaec` / `echocatzh` probes & DTLN AST mock check | Preserves tail truncation, `np.complex` failure & complex warnings separately from historical `621d727` reports |
| Known arrival & history audio & Fig 70 | `ch07/core/prediction_delays.py`, `core/delay_audio.py`, `examples/generate_delay_audio.py` → `ch07/delay_audio/`; `scripts/make_delay_figures.py` → Fig 70 & JSON | 6 `16 kHz`, `33,024`-sample mono/stereo synthetic WAVs + manifest | `32,256`-pt source, 384-pt propagation vs. 768-pt history padding, gain 1, `29,184`-pt scoring window; 5 sources, 7 members (`--check` read-only) |
| Ch 7 upstream WPE contracts | Two WPE audits + `wpe_silence_boundary.py` → three `*_current.json` reports | Full `nara_wpe` NumPy, bounded AST checks & long-silence blowup/recovery | Separates sparse-selection status from verified blobs while keeping historical `34d494` reports untouched |
| Known-path MINT inverse audio | `codes/chapters/ch07/core/mint_teaching.py`, `examples/mint_teaching_demo.py` → `ch07/mint_audio/` | 6 `16 kHz`, `32,512`-sample mono/stereo synthetic WAVs + manifest | Two sparse FIR scenes (well-conditioned vs. near-common-zero) + post-path noise, gain 1, `27,200`-pt window + 512-pt tail; 5 sources (`--check` read-only) |
| Multi-frequency geometry audio | `codes/chapters/ch03/core/geometry_audio.py`, `examples/generate_geometry_audio.py` → `ch03/geometry_audio/` | 3 `32 kHz` dual-tone mono/6-ch synthetic WAVs + manifest | Gain 1, steady-window per-bin LS complex ratio across 6 bound sources (`--check` read-only) |
| Known-direction baseline audio | `codes/chapters/ch03/core/baseline_calibration.py`, `core/baseline_audio.py`, `examples/generate_baseline_audio.py` → `ch03/baseline_audio/` | 6 `16 kHz`, `32,000`-sample mono/stereo `500 Hz` synthetic WAVs + manifest | `27,200`-pt steady window (850 cycles), `<1 ms` principal delay bound, 6 sources, 7 regular members (`--check` read-only) |
| Ch 3 upstream coarray contract | `codes/chapters/ch03/examples/audit_upstream_coarray.py` → explicit `reports/upstream_coarray_contracts.json` | 7 bounded `doatools` coarray (`DA`/`SS`) & perturbation methods | Defaults to `stdout`; preserves upstream list/non-square input bugs without overwriting legacy reports |
| Appendix B room simulation assets | `codes/chapters/appendix_b/examples/room_srp_exercise.py` → `appendix_b/room_audio/` | 6 positions × (dry, direct-only, full room) = 18 white-noise WAVs, manifest, `ROOM_RESULTS.png` & `RESULTS.json` | Generated in isolated `pyroomacoustics==0.10.0` env into fresh dir; `check_room_assets.py` verifies 21 regular members, 7 sources & 81-tap/40-pt fractional-delay metadata read-only |
| Binaural ITD/ILD cue audio | `codes/chapters/ch01/core/binaural_cues.py`, `examples/generate_binaural_cues.py` → `ch01/binaural_audio/` | 5 stereo synthetic WAVs + manifest (left channel first, gain 1) | Binds 4 sources; separates float vs. PCM scoring (`--check` read-only; synthetic panning, not measured HRTFs) |
| Directional spectral cue audio | `codes/chapters/ch01/core/spectral_cues.py`, `examples/generate_spectral_cues.py` → `ch01/spectral_audio/` | Two `32,000`-pt mono sources + two `32,001`-pt stereo outputs (`16 kHz`) + manifest | 2-tap FIR across `1`/`7 kHz` source tilts, gain 1, 1-sample tail, `28,800`-pt window; 4 sources, 5 members (`--check` read-only) |
| Ch 1 `libmysofa` interpolation contract | `codes/chapters/ch01/examples/audit_libmysofa_interpolation.py` → explicit `reports/libmysofa_interpolation.json` | Full upstream `libmysofa` interpolation function across 5 structural controls | Retains midpoint attenuation under common delay vs. per-direction control; defaults to `stdout` |
| Oracle T-F mask audio | `codes/chapters/ch08/core/mask_representation.py`, `examples/mask_representation_demo.py` → `ch08/mask_audio/` | 6 `16 kHz`, `32,000`-sample mono synthetic WAVs + manifest | Full-record FFT oracle bins + shared envelope, gain 1, `27,200`-pt window, 5 bound sources (`--check` read-only) |
| Prescribed-slot CSS audio & Fig 71 | `codes/chapters/ch08/core/css_audio.py`, `examples/generate_css_audio.py` → `ch08/css_audio/`; `scripts/make_css_figures.py` → Fig 71 & JSON | 4 `16 kHz`, `32,000`-sample stereo synthetic WAVs + manifest | Prescribed slot order with $-2\times$ block-2 jump fit over `6,400`-pt overlap; scores `27,200`-pt main & `10,400`-pt trailing windows; 6 sources, 5 members (`--check` read-only) |
| Controlled GSS teaching audio | `codes/chapters/ch08/core/gss_teaching.py`, `examples/gss_teaching_demo.py` → `ch08/gss_audio/` | 5 synthetic WAVs, `STATE.npz` & manifest (7 regular files total) | Rejects symlinks/extra files; `--check` replays in memory without disk writes |
| Dual-scenario selection audio | `codes/chapters/ch11/core/selection_audio.py`, `examples/generate_selection_audio.py` → `ch11/scenario_audio/` | 2 scenarios × (target, mixture, 3-tap, 9-tap) = 8 `16 kHz`, `32,008`-pt WAVs + manifest | Shared gain `0.8`, 8-pt tail, group-delay-aligned `[1600, 30400)` window; 4 sources, 9 regular members (`--check` read-only) |
| Noise-mismatch audio | `codes/chapters/ch10/core/noise_mismatch.py`, `examples/generate_noise_mismatch.py` → `ch10/noise_audio/` | 6 `16 kHz`, `32,000`-sample mono synthetic WAVs + manifest | Fixed-preamble, target-polluted & oracle-variance controls at gain 1 across two `6,400`-pt windows (`--check` read-only) |
| PCM-to-Kalman tracking audio | `codes/chapters/ch09/core/tracking_audio.py`, `examples/chapter09_tracking_audio.py` → `ch09/tracking_audio/` | 2 continuous-motion synthetic WAVs + per-frame GCC/gating/KF manifest | Separates float vs. readback PCM and `state_time` vs. availability time; 9 bound sources, 3 regular members (`--check` read-only) |
| Track lifecycle Fig 72 & JSON | `codes/chapters/ch09/core/lifecycle.py`, `examples/tracking_lifecycle_demo.py`, `scripts/make_tracking_figures.py` → Fig 72 & `codes/chapters/ch09/reports/figure72_tracking_lifecycle.json` | 197 frames over readback stereo PCM under two age-clock protocols | Binds 12 sources & 3 inputs; integer half-sample ticks, monotonic non-reused track IDs, zero extra WAVs |
| Moving-source free-field audio | `codes/chapters/ch09/core/moving_source.py`, `examples/moving_source_audio.py` → `ch09/moving_audio/` | 3 free-field moving-source WAVs + ground-truth trajectory manifest | 5 bound sources, 4 regular members (`--check` read-only) |
| DEMAND real environmental audio | `codes/chapters/ch02/core/real_recordings.py`, `examples/prepare_real_recordings.py` → `ch02/real_audio/` | 16-ch `NRIVER` 10 s excerpt + 3 mono derivatives (4 WAVs, 8 regular members total) | Pinned Zenodo SHA-256 & CC BY-SA 3.0 attribution/license preserved (`--check` offline read-only) |
| External upstream lockfiles | `codes/chapters/ch00/SOURCES.lock.json`, `ARCHIVE_SOURCES.lock.json`, `THIRD_PARTY*.md` | Official URLs, 40-char commit hashes / archive SHA-256s, licenses & boundaries | Record only verified primary sources; never infer redistribution rights from public visibility |
| Source verification status | `codes/chapters/ch00/SOURCE_STATUS.json`, `ARCHIVE_SOURCE_STATUS.json` | Offline verification of local `_downloads/` checkouts against current lockfiles | Generated strictly via `--verify --report`; never hand-edit failures to passing |
| Historical lockfile snapshots | `codes/chapters/ch00/source_snapshots/`, `codes/chapters/ch00/core/source_history.py` | Verbatim lockfile/status snapshots bound by historical contract reports | Preserves exact historical bytes and verifies referenced project records read-only |
| Fig 45 numerical report | `scripts/make_figures.py::fig_stft_convolution` → `codes/chapters/ch02/reports/figure45_stft_convolution.json` | 6-point full convolution, 4-point per-frame MTF & 8-point zero-padded controls with COLA denominator | Regenerated with figure; 8-pt control retains 4-pt synthesis window (not a full CTF filterbank) |
| Fig 13 Monte Carlo report | `scripts/make_figures.py` → `codes/chapters/ch04/reports/figure13_gcc_reverb.json` | Per-trial events, peak contrast & counts across 9 conditions × 150 trials | Regenerated simultaneously with figure; preserves shared PRNG seeds, failures & Wilson intervals |
| Upstream research manuals | `codes/chapters/ch00/research/*.md` (Chinese), `research_en/*.md` (English) | Per-algorithm upstream entrypoints, axis/sign discrepancies, industrial configs & contract audits | Cross-linked with chapters and `COVERAGE*.md`; separate proposed checks from executed runs |
| Repository `README` sources | `README.md`, `README_EN.md`, `scripts/README*.md`, `codes/chapters/README*.md`, `codes/chapters/ch00/{research*,upstream}/README*.md`, `codes/chapters/ch02/real_audio/README.md`, `scripts/vendor/mathjax-3.2.2/README*.md` | Maintained documentation entrypoints across Chinese and English editions | Keep shared counts/facts synchronized; never touch `ch02/real_audio/` member count (its English translation lives inside `codes/chapters/README_EN.md`) |
| Course Home & GitHub cover SVGs | `scripts/course_home/{index.html,home.css,home.js,export_portable.py}` → `site/index.html` & `assets/readme/course-{intro.html,cover.svg,cover-still.svg}` | Single-source course landing page, portable HTML & animated/static GitHub cover SVGs | Included in site input digest; `--check` verifies portable exports read-only (cover SVGs are not counted among the 87 numbered figures) |
| Homepage LOCATA room case | `scripts/course_home/prepare_locata.py` → `assets/readme/cases/locata/{before,after,mvdr,focus}.wav`, `levels.svg`, `MANIFEST.json`, `NOTICE.txt` | Pinned LOCATA dev Task 2 DICIT recording: SRP scan, DSB, loaded MVDR + postfilter & B-bypass | Binds 6 upstream files + MVDR kernel SHA-256 at shared export gain (`--check` read-only) |
| Interactive 4-Mic Browser Lab | `assets/readme/lab/{README.md,README_EN.md,array-lab.html,export-lab.mjs}` → `single.wav`, `left.wav`, `right.wav` & `delay-and-sum.svg` | Bilingual offline zero-permission 4-mic browser lab | Single embedded `<script id="lab-model">` kernel verified byte-for-byte via `node assets/readme/lab/export-lab.mjs --check` |
| Fig 20 WPE report | `scripts/make_figures.py::fig_wpe` → `codes/chapters/ch07/reports/figure20_wpe.json` | Parameters, frame selection, energy denominators & spectral-band errors on shared signal | Regenerated with figure; never conflate plotted spectrogram range with metric integration band |
| Fig 59 selection scenarios report | `scripts/make_figures.py::fig_selection_scenarios` → `codes/chapters/ch11/reports/figure59_selection_scenarios.json` | Dual-scenario readback PCM errors, analytic controls, scene weight $q$ & minimax bounds | Regenerated from standalone scenario manifest; preserves integer NMSE denominators without averaging dB |
| Ch 11 MeetEval contracts | `ch11/examples/audit_meeting_scoring_interfaces.py`, `audit_meeting_kernel_contracts.py` → `reports/*_current.json` | 32 Python module imports + 2 native C++ kernels across 13 cases | Defaults to `stdout`; writes in-repo only to `*_current.json` while preserving legacy reports |
| Fig 74 selection evidence report | `scripts/make_selection_figures.py` → `codes/chapters/ch11/reports/figure74_selection_evidence.json` | DS & 3 loaded MVDR weights on shared true covariance under `5°` mismatch & hard gates | Reuses Ch 5 beamformer kernel across 8 bound sources; unmeasured device latency remains undecided |
| Figs 57/58 engineering reports | `scripts/make_figures.py::fig_engineering_limits`, `fig_noise_mismatch` | Fig 57 stores exact soft-update fractions, finite EDC tail & blocking schedule; Fig 58 reads Ch 10 noise manifest | Regenerated from authoritative sources; separates float decomposition from readback PCM |
| Fig 56 tracking Fisher report | `scripts/make_figures.py::fig_tracking_information` → `codes/chapters/ch09/reports/figure56_tracking_information.json` | Instantaneous bearing range-scale nullspace & correlated-error conditioning | Regenerated with figure; pure geometric model omitting acoustic propagation delay |
| Ch 9 upstream tracking contracts | `ch09/examples/audit_upstream_tracking_contracts.py`, `audit_tracking_upstream_interfaces.py` → `ch09/reports/*_current.json` | Native ODAS C, SAF C `step` counting stubs & FilterPy/Stone Soup probes | Separates sparse-selection status from verified blobs without modifying legacy reports |
| Known dead-mic audio & Fig 73 | `ch10/core/channel_selection.py`, `core/channel_audio.py`, `examples/generate_channel_audio.py`; `scripts/make_channel_figures.py` → 6 WAVs & Fig 73 JSON | Known dead mic 0 excised across observations, steering & covariance axes | 8 audio sources, 9 figure sources, 7 regular members, `27,200`-pt window (`--check` read-only) |
| Ch 10 industrial contracts | 4 industrial/SRC/delay probes → `ch10/reports/*_current.json` | Bounded C/Python upstream calls vs. worktree sparse status | Defaults to `stdout`; preserves all 4 legacy historical reports intact |
| Figure plotting scripts | 12 `scripts/make_*_figures.py` scripts (`make_figures.py`, `make_aec_figures.py`, `make_beamforming_figures.py`, `make_reference_figures.py`, `make_delay_figures.py`, `make_css_figures.py`, `make_tracking_figures.py`, `make_channel_figures.py`, `make_selection_figures.py`, `make_wearable_figures.py`, `make_tse_figures.py`, `make_math_figures.py`) | Generate all 87 numbered pedagogical figures (`figures/` and `figures_en/`) | Fix plotting scripts and re-render whenever a figure needs updating; never edit PNGs manually |
| Site & PDF build scripts | `scripts/build_site.py`, `build_pdf.py`, `build_markdown_helpers.py`, `inline_layout.js`, `pdf_math.js`, `heading_aliases.py`, `legacy_sequential_anchors.json`, `chapter_identity.py`, `chapter_numbering.json` | Compile static course editions (`site/`, `site_en/`), combined HTML, and tagged A4 PDFs (`dist/*.pdf`) while preserving deep-link aliases | Shared script edits invalidate build digests for both HTML and PDF outputs |
| Vendored PDF MathJax/SRE | `scripts/vendor/mathjax-3.2.2/` | Pinned MathJax 3.2.2 SVG, SRE 4.0.6 English rules, Apache-2.0 licenses & legacy CHTML WOFFs | Preserve verbatim upstream bytes, SHA-256 checks & licenses |
| Generated figures | `figures/fig*.png`, `figures_en/fig*.png` | 87 numbered pedagogical figures | Generated strictly by plotting scripts |
| Compiled static websites | `site/**`, `site_en/**` | 19 current chapter pages + 8 legacy compatibility pages + 6 research manual pages (33 HTMLs) + staged audio copies | Generated by `build_site.py`; never edit directly |
| Dated QA audit snapshots | `reviews/2026-10-05-*.json`, `reviews/2026-10-05-human-listening-register.csv` | Dated technical audio audit, unexecuted human-listening plan/register, optional test audit & PDF formula/font controls | Dated historical evidence; never claim they pre-certify future edits |
| Compiled monograph outputs | `dist/combined*.html`, `dist/microphone-array-tutorial*.pdf` | Combined HTML intermediates and tagged A4 monograph PDFs | Generated by `build_pdf.py`; never edit directly |

`chapters/` (and `chapters_en/`) contains **19 Markdown files** (Overview, 16 chapters `00`–`15` including 4 frontier chapters, and 2 appendices `16`–`17`), accompanied by **6 research manual Markdown files**, **87 numbered figures**, and 1 room-simulation plot in Appendix B. Keep file counts, figure counts, chapter titles, and build entrypoints synchronized across `AGENTS*.md`, `README*.md`, `00_overview.md`, and `scripts/README*.md`.

Chapter numbering & legacy PDF topic identities: Chapter 0 (`00_math-foundations.md`) uses `E00-*` and tags `(0-k)`; Chapters 12–15 (Frontier Topics I–IV) use `E12`–`E15` and `(12-k)`–`(15-k)`; Appendices A and B (Chapters 16–17) use `E16`/`E17` and `(16-k)`/`(17-k)`. To preserve `ch-0` (Overview) and previously published PDF named destinations, `scripts/chapter_identity.py` and `chapter_numbering.json` map stable internal PDF topic IDs (`ch-18` Ch 0 Math, `ch-14` Ch 12 Imaging, `ch-15` Ch 13 Distributed, `ch-16` Ch 14 Wearable, `ch-17` Ch 15 TSE, `ch-12` Ch 16 Appendix A, `ch-13` Ch 17 Appendix B) alongside 8 legacy HTML routes and fragment aliases.

### 1.3 Primary Build & Verification Commands

Run from the repository root:

```bash
.venv/bin/python codes/chapters/ch00/examples/generate_audio_samples.py  # Main 109 WAVs + manifest
.venv/bin/python -m codes.chapters.ch01.examples.generate_binaural_cues  # 5 binaural WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch01.examples.generate_spectral_cues  # 4 spectral-cue WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch02.examples.generate_stft_convolution  # 3 STFT convolution WAVs
.venv/bin/python -m codes.chapters.ch02.examples.generate_sweep_audio  # 4 sweep deconvolution WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch03.examples.generate_geometry_audio  # 3 multi-freq geometry WAVs
.venv/bin/python -m codes.chapters.ch03.examples.generate_baseline_audio  # 6 baseline calibration WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch04.examples.generate_reflection_audio  # 4 coherent-reflection WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch04.examples.generate_focus_audio  # 4 unitary-focusing WAVs
.venv/bin/python -m codes.chapters.ch05.examples.generate_derivative_audio  # 4 derivative-constraint WAVs
.venv/bin/python -m codes.chapters.ch05.examples.generate_phase_audio  # 3 phase-reference WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch05.examples.binaural_beamforming_demo --write  # 4 stereo WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch06.examples.generate_apa_audio  # 6 colored-ref APA WAVs
.venv/bin/python -m codes.chapters.ch06.examples.generate_reference_audio  # 6 playback-gain/tail WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch07.examples.mint_teaching_demo  # 6 MINT inverse WAVs
.venv/bin/python -m codes.chapters.ch07.examples.generate_delay_audio  # 6 arrival/history WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch08.examples.mask_representation_demo  # 6 oracle-mask WAVs
.venv/bin/python -m codes.chapters.ch08.examples.generate_css_audio  # 4 CSS polarity/gain WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch08.examples.gss_teaching_demo  # 5 GSS WAVs + STATE.npz
.venv/bin/python -m codes.chapters.ch09.examples.chapter09_tracking_audio  # 2 PCM-to-KF tracking WAVs
.venv/bin/python -m codes.chapters.ch09.examples.moving_source_audio  # 3 moving-source WAVs
.venv/bin/python codes/chapters/ch02/examples/prepare_real_recordings.py --check  # Offline check of 4 DEMAND WAVs
.venv/bin/python -m codes.chapters.ch10.examples.generate_noise_mismatch  # 6 noise-mismatch WAVs
.venv/bin/python -m codes.chapters.ch10.examples.generate_channel_audio  # 6 dead-mic WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch11.examples.generate_selection_audio  # 8 dual-scenario WAVs (--check read-only)
.venv/bin/python -m codes.chapters.appendix_a.examples.generate_weighted_audio  # 5 weighted-noise WAVs (--check read-only)
.venv/bin/python -m codes.chapters.appendix_b.examples.generate_response_audio  # 5 equal-DRR WAVs (--check read-only)
.venv/bin/python -m codes.chapters.appendix_b.examples.check_room_assets  # Read-only check of 21 room members
.venv/bin/python -m codes.chapters.ch12.examples.generate_imaging_audio  # 5 imaging snapshot WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch13.examples.generate_distributed_audio  # 17 distributed WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch14.examples.generate_wearable_audio  # 6 wearable WAVs (--check read-only)
.venv/bin/python -m codes.chapters.ch15.examples.generate_tse_audio  # 6 TSE distance-bubble WAVs (--check read-only)
.venv/bin/python scripts/make_figures.py      # Figs 1–24, 28–31, 35–67, 75
.venv/bin/python scripts/make_aec_figures.py  # Figs 25–27, 32–34
.venv/bin/python scripts/make_beamforming_figures.py  # Fig 68
.venv/bin/python scripts/make_reference_figures.py  # Fig 69
.venv/bin/python scripts/make_delay_figures.py  # Fig 70
.venv/bin/python scripts/make_css_figures.py  # Fig 71
.venv/bin/python scripts/make_tracking_figures.py  # Fig 72
.venv/bin/python scripts/make_channel_figures.py  # Fig 73
.venv/bin/python scripts/make_selection_figures.py  # Fig 74
.venv/bin/python scripts/make_wearable_figures.py  # Figs 76–79
.venv/bin/python scripts/make_tse_figures.py  # Figs 80–83
.venv/bin/python scripts/make_math_figures.py  # Figs 84–87
.venv/bin/python -m codes.chapters.ch00.chapter00_exercises  # Ch 0 16 math foundation exercises
.venv/bin/python scripts/build_site.py        # Build static course website
.venv/bin/python scripts/build_pdf.py         # Build combined HTML & tagged A4 PDF
.venv/bin/python -m unittest discover -s tests -p 'test_codes*.py' -v  # Pedagogical code regression suite
```

Prose-only edits do not require re-rendering unchanged figures, whereas edits to plotting scripts require regenerating and visually inspecting affected PNGs, and edits to chapters, navigation, or builders require rebuilding affected HTML and PDF outputs. Offline PDF math rendering uses vendored MathJax 3.2.2 SVG + SRE 4.0.6 English rules without network access (whereas the static website loads MathJax online). Whenever headless Chrome, optional libraries, network access, or permissions are unavailable on a host, report the unexecuted checks and reasons honestly rather than claiming full verification.


## 2. Rule Priority, Authorization, and Scope Control

### 2.1 Rule Priority Hierarchy

Three **non-negotiable hard constraints** must be satisfied simultaneously on every task:

- **Obey system/platform safety and data-integrity boundaries**: Never perform unauthorized or destructive operations, and never corrupt user work or Git history.
- **Obey the user's task scope and authorization boundaries**: Authorization defines what actions may be taken, never permission to fabricate facts, equations, numbers, or citations.
- **Guarantee factual and mathematical accuracy**: Never write known errors, unverified assertions, or invented figures merely to fill a table or match requested phrasing.

Once all three hard constraints hold, trade off quality objectives in the following strict order:

1. **Primary-source provenance and reproducibility**: State exact assumptions and citations for key claims; specify configurations, conventions, and PRNG seeds for simulations.
2. **Reader comprehension and cohesive human textbook narrative**: Break down every concept, comparison, derivation, hand calculation, figure, and table step by step—weaving definitions, intuition, algebra, worked examples, and failure boundaries into natural academic prose rather than bolting cold boilerplate disclaimers onto the ends of paragraphs.
3. **Cross-book consistency and artifact synchronization**: Keep terminology, notation, figure numbers, navigation, and compiled HTML/PDF outputs aligned.
4. **Typesetting and local stylistic polish**: Ensure crisp figures and balanced paragraph spacing without altering technical meaning or substituting prose polish for fact-checking.

Lower-priority rules never override higher-priority rules: never invent a benchmark number just to make a sentence "concrete," and never replace a precise domain term with an imprecise colloquialism just for stylistic variety.

### 2.2 Authorization Boundaries

- Reading files and running read-only checks do not mutate repository state. When the user requests edits, running existing local builds and test suites is part of standard verification; pure review tasks must never overwrite generated artifacts under the guise of "checking."
- Keep edits strictly within the user's requested scope. If out-of-scope issues are noticed, log their file paths, line numbers, and impact without rewriting unrelated chapters unprompted.
- Never open external issues/PRs, send external messages, deploy websites, or `git commit` / `git push` unless explicitly requested by the user.
- Inspect `git status` / `git diff` before editing. Pre-existing user or collaborator changes belong to them—coordinate if scopes overlap, and never overwrite or revert their work.
- When user authorization for edits, builds, commits, or pushes has already been granted for the task, do not repeatedly re-prompt per tool or retry; verify sandbox/approval capabilities across active subagents separately rather than assuming older restricted subagents inherited updated permissions.

### 2.3 Handling Uncertainty

Internal review ledgers and scratch notes may mark items `to-be-verified` (`TBD`) to prevent guessing, whereas published Markdown, HTML, and PDF deliverables must **never** contain placeholder tokens (`TODO`, `TBD`, `link pending`, `number pending`).

Whenever primary evidence is incomplete, state bounded, well-sourced conclusions clarifying what is established versus what remains open—or omit non-essential assertions altogether. Never delete defect records, weaken quality gates, or invent numbers to simulate completion.

### 2.4 Controlling Change Scope

Before editing, state the objective, target source files, and verification plan. Keep each commit focused on a coherent theme without bundling unrelated cleanups.

Whenever a term or equation changes, trace every direct consumer: definition site, cross-references, master notation table (Appendix A), figure captions, Overview, and exercise solutions. When adding, deleting, or moving numbered exercises, figures, tables, or equations, verify every cross-reference by **semantic content**—never assume uniform offset shifting or check link existence alone.

Published `sec-x-y` and `sec-N` URL fragments are public deep-link contracts. Before renumbering headings, record the topic owned by each existing anchor so old IDs are never reassigned to unrelated topics, and register historical aliases in both the website and combined PDF while updating visible prose cross-references to current section numbers. Ensure every HTML/PDF `id` remains unique.

## 3. Tiered Workflows and Multi-Agent Collaboration

### 3.1 Classify by Risk and Blast Radius First

Line count alone does not measure blast radius—a single sign change in a steering-vector equation affects the entire book and cannot be treated as a minor tweak.

| Tier | When Applicable | Typical Examples |
|---|---|---|
| **Minor Fix** | Local phrasing, typos, formatting, or broken local links that do not alter technical claims, equations, figures, or section structure | Fixing a typo, standardizing punctuation, repairing a relative path |
| **Chapter-Level** | Modifies technical content, derivations, worked examples, figures, terminology, or pedagogical order within one or a few related chapters | Fixing an HPBW derivation, adding an anchor hand calculation, re-rendering figures |
| **Book-Wide** | Alters overall chapter structure, unifies notation/terminology globally, sweeps prose style, updates build pipelines, or audits all chapters | Reordering chapters, book-wide formula audit, redesigning site navigation |

When in doubt, escalate to the higher tier. Any change touching core equations, physical units, citations, experimental numbers, or generator scripts requires at least **Chapter-Level** handling.

### 3.2 Minor Fix Workflow

1. Read the target passage and surrounding context to confirm no symbol collisions, duplicate definitions, or downstream cross-references are affected;
2. Apply the minimal edit without rewriting adjacent paragraphs;
3. Inspect `git diff`, links, and Markdown syntax;
4. Rebuild and spot-check affected HTML/PDF pages if published outputs change;
5. Report exact modified locations and executed checks. (Escalate immediately to Chapter-Level if any technical conclusion must change.)

### 3.3 Chapter-Level Workflow

1. **Audit**: Read the target chapter and every passage referencing it. Where multi-agent review is available, run at least one independent read-only review pass focusing on math derivations, pedagogical flow, or code reproducibility.
2. **Plan**: Deduplicate findings, bound scope, and verify primary equations and external numbers before editing.
3. **Edit**: Partition write ownership strictly across mutually exclusive files (at most one agent writing any file at a time; shared indexes, Appendix A notation, `README*`, and build scripts are coordinated centrally by the lead agent).
4. **Verify**: Re-read context, check cross-chapter links, re-run affected plotting scripts and tests, rebuild HTML/PDF, and inspect rendered figures, equations, tables, navigation, and bookmarks. (In single-agent environments, run two distinct passes—technical accuracy first, cross-chapter consistency & rendering second—without mislabeling self-review as independent parallel review.)

### 3.4 Book-Wide Workflow

When the user explicitly specifies a sequential per-chapter multi-subagent protocol (e.g., three read-only reviewers per chapter covering math/pedagogy, primary sources/code, and experiments/figures, completing full chapter verification before advancing), follow that sequential gate strictly without batching chapters prematurely.

Otherwise, execute book-wide updates in two phases:
1. **Read-Only Audit Phase**: Partition non-overlapping domains (curriculum/navigation; foundations/geometry/Appendix A; DOA/beamforming/AEC/WPE/separation/tracking; engineering/selection/exercises; plotting/site/PDF build toolchain). Subagents submit structured defect records without modifying files; the lead agent deduplicates and prioritizes a single consolidated ledger.
2. **Mutually Exclusive Write & Unified Build Phase**: Each writer re-reads assigned files prior to editing to ensure zero concurrent drift, never overwriting collaborator changes. Once all source edits settle, run all figure generators, test suites, website builds, and PDF builds from a single coordinator, verifying every resolved defect item by item. Short chapter execution notes belong in `codes/chapters/README*.md` rather than proliferating one-off Markdown files.

### 3.5 Defect Record Format, Severity, and Lifecycle States

```text
ID | File & Line/Figure | Category | Severity | Evidence | Proposed Fix | Owner | Status | Verification
```

- **P0 Blocker**: Build failure, missing critical artifacts, or data-safety risk (blocks release unconditionally).
- **P1 Critical**: Incorrect core equation, numerical result, physical unit, citation, or technical conclusion that would mislead readers (must be fixed prior to release).
- **P2 Major**: Skipped derivation step, missing failure boundary, weak example, or broken deep-link navigation impairing comprehension or reproduction.
- **P3 Minor**: Local formatting, phrasing, or minor stylistic inconsistency that does not alter technical meaning.

Valid states: `Identified`, `In Scope`, `Modified` (source edited, awaiting artifact check), `Verified` (both source and compiled outputs verified), `Deferred` (with explicit reason and prerequisites; never allowed for P0/P1 in a release), and `Out of Scope`.

### 3.6 Ten Review Dimensions

1. **Completeness** (definitions, comparisons, selection criteria, cross-links);
2. **Logical order** (no use-before-definition, circularity, or silent symbol redefinition);
3. **Expository depth** (zero skipped algebra steps while keeping tangential trivia bounded);
4. **Equations & numbers** (derivations, units, approximations, metric denominators, code parity);
5. **Examples & boundaries** (reproducible hand calculations, explicit operating assumptions & failure modes);
6. **Figures & tables** (axis labels, units, grayscale-safe legends, effective font sizes, model fidelity);
7. **Structure & navigation** (heading hierarchy, TOC, prev/next links, sidebar `nav`, PDF bookmarks);
8. **Rendering & reproducibility** (MathJax rendering, table scrolling, code blocks, deterministic scripts);
9. **Accessibility** (keyboard navigation, visible focus, `aria-current`, table `<th scope="col">`, tagged PDF `/StructTreeRoot`, formula `ActualText`);
10. **Prose style** (absence of AI boilerplate, corporate jargon, anthropomorphic slang, or unbacked absolutes per Section 6).

## 4. Evidence Tiers, Citations, Currency, and Metric Conventions

### 4.1 Four Evidence Tiers

| Tier | Source Type | What It Can Support | Usage Rule |
|---|---|---|---|
| **A: Primary Evidence** | Official standards, specifications, original papers, official corpus manuals, verified upstream source commits | Equation provenance, standard clauses, algorithm definitions, benchmark numbers, software behavior | Preferred; cite exact section, equation, table, version, or 40-char commit hash |
| **B: Authoritative Exposition** | Classic monographs, peer-reviewed surveys, university lecture notes | Consensus definitions, historical context, physical intuition, taxonomy | Valid for exposition, never a substitute for primary benchmark tables (cite edition/pages) |
| **C: Project / Vendor Docs** | Repository `README`s, vendor datasheets, application notes, author technical posts | Implementation-specific parameters, APIs, vendor claims, reproduction steps | Supports only claims about that specific project/product with version & verification date |
| **D: Discovery Leads** | Search snippets, aggregators, forum posts, unverified secondary summaries | Locating keywords and tracing back to Tier A–C sources | **Never** cite as evidence for technical claims or list in references |

### 4.2 Mandatory Citation Scope

- Cite original papers, standards, or classic textbooks whenever introducing a named algorithm, non-obvious equation, standard limit, corpus size, hardware spec, or software default.
- Pin paper benchmark citations to exact table/figure/section numbers rather than linking bare landing pages.
- Label in-repo derivations and experiments explicitly as textbook derivations or simulations pointing to the exact script/module.
- Prefer DOIs, standards bodies, archival proceedings, and official repositories with verification dates.

### 4.3 Temporal Currency & Inclusion Gate for Recent Literature

- Label preprints explicitly as **preprints** rather than peer-reviewed publications; bind accepted-in-press papers to venue and year without guessing volume/page numbers.
- Bind "recent" or "current" claims to explicit verification dates or publication years, and bind challenge references (`CHiME`, `DNS`, `AEC`, `LOCATA`, `DCASE`, `NOTSOFAR`) to exact edition, track, split, and leaderboard status.
- **Four-part gate before adding a recent paper**: the prose must explain (1) what concrete failure of prior baselines it solves, (2) what exact mathematical/architectural step changed, (3) what assumptions it relies on, and (4) how a reader can run a minimal verification or apply it in system selection. If two or more are missing, keep it out of main chapter headings.

### 4.4 Eight-Point Minimum Specification for Any Performance Number

1. **Dataset** (name, version, split, real recording vs. synthetic simulation);
2. **Acoustic conditions** ($f_s$, physical device count vs. per-device mic count vs. selected channels, close-talk/optical reference role, source count, SNR, SER, $T_{60}$, noise field);
3. **Algorithm configuration** (window, hop, filter taps, guard delay $\Delta$, iterations, causal vs. lookahead, oracle priors);
4. **Metric definition** (exact variant, unit, direction—never conflate `SI-SDR` with `SI-SDRi`, or absolute WER points with relative WER reduction: a drop from `20%` to `15%` is **5 percentage points** or **25% relative reduction**, never "5% improvement");
5. **Baseline parity** (identical data, preprocessing, $f_s$, and scoring script);
6. **Statistical aggregation** (mean, median, single deterministic case, or confidence interval with trial count $N$, failure count, denominator, and shared-seed disclosure);
7. **Compute & latency context** (CPU/MCU/GPU model, threads, batch size, precision, timing scope; compute real-time factor `RTF` strictly against **newly advanced hop duration $H/f_s$**, never window length $N/f_s$ or multiplied by channel count $M$);
8. **Exact source pointer** (paper table/figure, standard clause, or in-repo script path).

Additional measurement rules:
- Multi-scenario selection fixes scenario weights, tokenization, and failure handling upfront; candidates whose confidence intervals straddle a hard threshold or lack hardware measurements are classified as **insufficient evidence** rather than passing.
- Compare exact integer/rational counts before float casting where tie-breaking matters, and preserve log-probabilities or explicit underflow flags rather than reporting tiny positive densities as mathematical zero.
- State whether bulk propagation delay, algorithmic group delay, and scalar gain are aligned or intentionally retained, and distinguish linear-subtractor output from post-NLP output in AEC comparisons.
- Separate unlabeled spatial set error (`OSPA`/`GOSPA`), track ID switches (`IDSW`/`MOTA`), and coasting uncertainty during dropouts in multi-target tracking.


## 5. Equations, Worked Examples, Figures, Terminology, and Beginner Clarity

### 5.1 Eight-Step Equation Verification Protocol

1. **State physical assumptions upfront** (far-field vs. near-field, narrowband vs. broadband, stationarity, statistical independence, linearity, causality).
2. **Pin sign & coordinate conventions** (DFT sign, delay sign, azimuth/elevation origin, steering conjugation, inner product; highlight any departure from earlier chapters).
3. **Define every symbol on first appearance** (meaning, scalar/vector/matrix type, dimensions, units; distinguish transpose $^\top$, conjugate $^*$, and Hermitian transpose $^H$).
4. **Write out non-trivial intermediate steps** (never hide sign flips, dimension reductions, normalizations, or Woodbury/Wirtinger steps behind "it is easy to show").
5. **Audit mathematical & dimensional consistency** (matrix shapes, summation bounds, non-negativity of probabilities/powers, Hermitian PSD property of covariances, matching units on additive regularization $\varepsilon\mathbf{I}$ vs. relative scale $\alpha$).
6. **Check limiting & degenerate cases** (zero noise, $M=1$, $M=2$, anechoic $T_{60}\to 0$, exact steering match, $\alpha\to 0$ or $\alpha\to\infty$).
7. **Recompute independently** on a small hand-calculable numerical example sharing identical conventions across prose, plotting scripts, and `core/` code.
8. **Tagging & cross-referencing**: Assign equation tags such as `\tag{5-1}` strictly sequentially within each chapter to core definitions, closed-form solvers, objective functions, state recursions, and metrics; verify that every prose reference ("Eq. (5-1)") matches the target equation's actual assumptions and variables.

Additional equation rules:
- Never place multi-line display equations inside `>` blockquotes where `>` could leak into TeX; keep display math as top-level blocks.
- State explicit bin ranges whenever writing $f_k = k f_s / N$ (distinguishing non-negative `rfft` bins, negative upper-half full-DFT bins, and even-$N$ Nyquist bins under fractional delays).
- Distinguish zero population mean $\mathbb{E}[x]=0$ from zero finite-sample mean $\bar{x}=0$, and use $\|\vec{z}\|_2^2 = \vec{z}^H\vec{z}$ over $\mathbb{C}^M$.
- Verify decibel conventions explicitly: $10\log_{10}(P_2/P_1)$ for power/energy/variance/PSD ratios; $20\log_{10}(|A_2|/|A_1|)$ only when comparing amplitudes across identical impedance/proportionality constants; specify whether multi-block `ERLE` averages dB values or pools linear energy first; never strip reference suffixes from `dB SPL`, `dBFS`, or `dBm`.

### 5.2 Four-Tier Worked Example Taxonomy

| Tier | Pedagogical Purpose | Minimum Requirements |
|---|---|---|
| **Anchor Hand Calculation** | Walks the reader through a core algorithm for the first time | Full inputs, units, formula substitution, key intermediate values, rounding rule, final answer, "what this shows," and "when this breaks" (at least one per core algorithm) |
| **Auxiliary Example** | Clarifies a local step or parameter trade-off | Focuses on the step at hand without omitting intermediate quantities that determine sign, order of magnitude, or conclusion |
| **Engineering Case** | Illustrates parameter sizing, compute/memory cost, or field failure modes | Hardware, dataset, acoustic setup, software version, metric, and generalization limits |
| **Counterexample / Boundary Case** | Refutes an over-broad intuition or exposes algorithmic breakdown | Pinpoints which assumption changed, which step broke, and what cannot be inferred from the failure |

### 5.3 Objective Standards for Figures and Tables

All figures are generated deterministically by `scripts/make_*_figures.py` with fixed PRNG seeds:
- Label every axis with quantity name and unit (or explicit dimensionless definition).
- Pair colors with distinct line styles, markers, or direct callouts so grayscale printing remains legible.
- Verify effective rendered font sizes at A4 body width (~`165 mm`, $\ge 7\,\mathrm{pt}$ for conclusion-bearing labels) and web body width (~`860 px`, $\ge 11\,\mathrm{CSS\ px}$)—split crowded figures rather than shrinking text.
- State array geometry, $f_s$, frequency, SNR, $T_{60}$, and normalization in captions, and keep `(a)(b)(c)` panel scales consistent when comparing the same physical quantity.
- Ensure system block diagrams match the exact signal/reference/mask/covariance dataflow of the generating code (never omit playback loopback in AEC or collapse activity/mask/SCM/steering into one vague arrow), and ensure RIR direct-vs-reverberant splits use identical filter lengths and high-pass settings.
- Tables must specify column names, units, consistent precision, distinct markers for `missing` / `N/A` / `not_run`, and fit A4 width without horizontal clipping.

### 5.4 Terminology, Abbreviations, and Notation

- Expand every technical acronym on its first occurrence in each standalone chapter (`Minimum Variance Distortionless Response (MVDR)`), using consistent terminology thereafter.
- Disambiguate overloaded terms explicitly (e.g., per-bin PHAT magnitude normalization vs. spatial covariance whitening; OS device `suspend` vs. VAD tail `hangover`).
- Keep notation consistent across adjacent chapters and synchronized with Appendix A (`16_appendix-symbols-math.md`).

### 5.5 Beginner-Friendly "Zero-Gap" Exposition Standards

- Explain the physical or geometric meaning in plain language **before** introducing any matrix factorization, Wirtinger derivative, Sherman–Morrison–Woodbury update, subspace projection, or constrained optimization, then write out every algebraic step explicitly.
- Motivate new concepts from observable acoustic phenomena, minimal numerical examples, or bounded analogies before presenting general equations.
- Walk the reader through every equation, `(a)(b)(c)` figure panel, and multi-row comparison table immediately after presenting it.
- Structure heading hierarchies (`h2`/`h3`/`h4`) so every major algorithm, experiment, or troubleshooting guide is directly reachable from the sidebar TOC and PDF bookmarks.

### 5.6 Paragraph and Block-Level Structure

- Separate Markdown paragraphs, lists, display equations, tables, code blocks, and blockquotes with explicit blank lines; indent continuation paragraphs/math inside ordered list items by 4 spaces and verify compiled `<ol>/<li>` structure in HTML.
- Keep definitions, assumptions, display equations, and immediate symbol glosses adjacent; break long multi-clause enumerations into bulleted lists or tables.

## 6. Academic Prose Style Guidelines

### 6.1 Three-Step Self-Edit Pass

1. **Deletion check**: Cut throat-clearing warmups ("It is worth noting that...") whenever deleting them preserves 100% of the factual and causal content.
2. **Referent check**: Replace vague pronouns ("this," "it," "the above approach") whenever the antecedent is ambiguous with the exact signal, matrix, module, or metric.
3. **Concretization check**: Replace corporate/colloquial metaphors ("fallback safety net when inversion blows up") with observable actions and conditions ("switch to delay-and-sum weights whenever the condition number of $\hat{\mathbf{R}}_{xx}$ exceeds $10^4$").

### 6.2 Anti-Patterns to Eliminate

- Formulaic AI transitions, meta-commentary labels ("In plain English:", "One-sentence takeaway:"), internal sprint jargon ("deliverable," "closed-loop alignment"), anthropomorphic slang ("feeding masks into MVDR," "the room sabotages the array"), unbacked marketing superlatives ("universally optimal," "completely eliminates reverberation"), and stacked unit-test-style negative disclaimers at the end of every paragraph ("this is not real speech, not blind estimation, not an industrial ranking").
- Avoid mistranslations or formatting artifacts: distinguish *spatial aliasing* from spectral aliasing, *speech absence* from missing files, *finite support* (non-zero tap range) from software support, leading-zero decimals (`0.5`, never `.5`), and clean Markdown links without citation-tool title residue.

### 6.3 Acceptable Idioms and Technical Usage

- Sequential transitions ("First," "Next," "Finally") are appropriate when describing genuine ordered steps of an algorithm; use numbered lists whenever steps are complex.
- Standard technical acronyms (`DOA`, `STFT`, `MVDR`, `WPE`) are retained after expanding them on first occurrence in each chapter.
- Code identifiers, CLI commands, API symbols, and established international algorithm names remain in their exact canonical form.
- Domain terms with precise meanings—"signal chain / pipeline" for module ordering, "theoretical bound" for a formula/standard ceiling, and "baseline" for a reproducible reference experiment—are encouraged.
- Contrastive constructions ("not A, but B"), bounded physical analogies, and occasional rhetorical questions are acceptable when informative and non-repetitive.

### 6.4 Positive Writing Examples

- **Clean transition**: "The previous section estimated inter-microphone time delays. This section solves for the 3D source coordinates from multiple delay pairs."
- **Bounded quantitative claim**: "For an equal-weight delay-and-sum beamformer, once the target wavefront is time-aligned and sensor noises are mutually uncorrelated with equal variance, doubling the number of microphones increases white-noise gain by approximately `3 dB`."
- **Direct causal link**: "Sub-sample time-delay resolution finer than $1/f_s$ requires peak interpolation across adjacent lag bins."
- **Figure walkthrough**: "The red double-headed arrow spans the two `−3 dB` crossings, marking the half-power beamwidth (`HPBW`)."
- **Analogy with explicit boundary**: "Diagonal loading can be intuited as penalizing excessive weight norm. That analogy captures numerical conditioning; algebraically, the algorithm adds $\varepsilon\mathbf{I}$ directly to $\hat{\mathbf{R}}_{xx}$ prior to inversion."

### 6.5 Three-Tier Information Routing & Global Mathematical Notation Canon

1. **Three-Tier Information Routing** (keep machine-audit ledgers out of human textbook prose):
   - **Main Textbook Chapters (`chapters/*.md`, `chapters_en/*.md`)**: Physical intuition, modeling assumptions, zero-gap derivations, hand calculations, figure walkthroughs, physical failure boundaries, and engineering selection criteria.
   - **Upstream Research Manuals (`codes/chapters/ch00/research/*.md`, `research_en/*.md`)**: Exact third-party file/line pointers, tensor axis discrepancies, upstream bug analyses, and reproduction commands.
   - **Machine Manifests, JSON Reports & Unit Tests (`MANIFEST.json`, `reports/*.json`, docstrings, `tests/`)**: 40-char Git SHAs, `source_selection_mismatch` flags, exact directory member counts, and IEEE 754 `float64` extreme-exponent guards ($10^{200}$).
2. **Global Mathematical Notation Canon**:
   - **Continuous vs. discrete brackets**: Continuous-time signals use parentheses $x(t)$; discrete-time sequences and integer-lag correlations use square brackets $x[n]$, $R[l]$ (never reuse `n` for both noise signal and sample index—write additive noise as $v_m[n]$ or $\vec{v}$).
   - **Imaginary unit, transposes & expectation**: Upright imaginary unit `\mathrm{j}` (never italic `j` or `i`); real transpose `^\top` (never `^T` or `^\mathsf{T}`); Hermitian transpose `^H`; expectation `\mathbb{E}[\cdot]`; base-10 log `\log_{10}`, natural log `\ln`; upright text subscripts `\mathrm{...}` (never deprecated `\rm`).
   - **Vectors, matrices & covariance subscripts**: Column vectors use `\vec a, \vec w, \vec x`; matrices use upright bold `\mathbf R, \mathbf A, \mathbf I`; spatial covariance matrices use double subscripts `\mathbf R_{xx}, \mathbf R_{ss}, \mathbf R_{nn}` (finite-snapshot estimate `\hat{\mathbf R}_{xx}`), strictly distinct from scalar cross-correlation $R_{x_1 x_2}(\tau)$ / $R[l]$.
   - **Frequency & equation tags**: Continuous frequency $f$ ($\mathrm{Hz}$), analog angular frequency $\omega = 2\pi f$ ($\mathrm{rad/s}$), digital angular frequency $\Omega_{\mathrm{d}} = 2\pi f/f_s$ or $2\pi f/f_s$; strictly single-dollar (`$`) and double-dollar (`$$`) TeX delimiters (never LaTeX backslash-paren/bracket delimiters, never backtick-wrapped dollar math, and never Chinese punctuation inside math mode); monotonic chapter equation tags `\tag{c-k}` starting at `1`.

## 7. Verification Matrix and Release Gates

| Change Type | Mandatory Verification | Pass Criteria |
|---|---|---|
| Prose, structure, terminology | Re-read modified & adjacent sections; search for stale terms/placeholders; verify internal links | Zero broken sentences, duplicate definitions, or stale terms; all targets exist |
| Equations, derivations, hand calcs | Recompute independently; audit units, conjugation & signs; check a limiting case | Every intermediate step reproducible; dimensions match; limiting case agrees with physics |
| Numbers, tables, experimental claims | Recompute via script or standalone math; verify rounding across prose, tables & figures | Reproducible within stated precision under identical conventions |
| Citations, standards, upstream versions | Verify against primary paper/standard/official repo at pinned version & date | Reachable primary source with explicit version & operating conditions |
| Figure plotting scripts | Re-render target PNGs; inspect full-size & scaled axes, units, legends, font sizes & colors | Deterministic execution; matches prose model; zero clipping or illegible text |
| Markdown → static website | Run `build_site.py`; check logs, page count, internal links, anchors, images & MathJax | Clean build; zero missing pages, dead links, missing images, or unrendered TeX |
| Combined HTML & monograph PDF | Run `build_pdf.py`; inspect cover, TOC, edited pages, wide tables, formulas, figures & bookmarks | Clean build; zero overflow, clipping, blank formulas, or mis-targeted bookmarks |
| Navigation, TOC, metadata | Cross-check titles, order, prev/next links & bookmarks across Markdown, website & PDF | 1-to-1 title/order parity and reachable destinations across all three surfaces |
| Semantic cross-references | Parse equation tags, section numbers, figure numbers & section links | Unique numbering; figure numbers match filenames; section links land on exact subsection |
| Web & PDF accessibility | Check `lang`, `alt` text, heading hierarchy, `<th scope="col">`, keyboard focus, `/StructTreeRoot`, `Formula` `ActualText` & bookmark targets | Supported checks pass; uncertified screen-reader/PDF-UA items stated honestly |
| Build scripts or dependencies | Run affected builders plus regression tests (`unittest discover` & `quality_check.py`) | Zero exit code, complete outputs, zero unintended dependency/file mutations |

Current full-build baselines enforced by explicit constants in `scripts/quality_check.py`: **87 numbered PNGs** + **1 room plot**, **19 current tutorial pages** + **8 legacy compatibility pages** + **6 research manual pages** (**33 HTML pages**), and **19 chapter-level + 172 section-level + 827 subsection-level = 1,018 PDF bookmarks**. Quality gates must compare against independent explicit baselines rather than self-comparing dynamically generated counts, verify PNG `SourceScriptDigest` provenance, check sparse/blank PDF pages, and enforce independent expected values in numerical unit tests (checking Hermitian/PSD properties, scaling invariance, complex-input rejection on real APIs, overflow guards, and degenerate zero-energy/all-dropout edge cases).


## 8. Engineering Discipline

### 8.1 Workspace Hygiene and Pre-Existing Changes

- Read target files in full and inspect workspace diffs before making edits.
- Respect user-owned modifications—never overwrite, revert, or opportunistically reformat unrelated files.
- Never run destructive history/data commands unless explicitly requested and scoped by the user.
- Partition large multi-file efforts across mutually exclusive file sets, reconciling shared notation and cross-references at the end.

### 8.2 Authoritative Sources and Deterministic Regeneration

- Always modify source files first and regenerate downstream artifacts (`site/`, `dist/`, `figures/`, `MANIFEST.json`).
- Numbered figures are generated exclusively through the 12 registered plotting scripts in `scripts/` (`make_figures.py`, `make_aec_figures.py`, `make_beamforming_figures.py`, `make_reference_figures.py`, `make_delay_figures.py`, `make_css_figures.py`, `make_tracking_figures.py`, `make_channel_figures.py`, `make_selection_figures.py`, `make_wearable_figures.py`, `make_tse_figures.py`, `make_math_figures.py`) with fixed PRNG seeds and verified `SourceScriptDigest` PNG metadata. Full releases must also regenerate all 87 numbered figures from current sources; any diffs must trace to defect records or intentional edits.

### 8.3 Virtual Environments, Dependencies, and Reproducibility

- Use the repository's primary `.venv` (`requirements.txt`: `numpy`, `matplotlib`, `Markdown`, `pypdf`) without silently upgrading packages or adding heavyweight frameworks.
- Install optional external dependencies (`scipy`, `pyroomacoustics==0.10.0`, C FFTW) strictly inside isolated `/tmp` virtual environments, recording exact package versions and commit hashes while keeping `.venv` clean.
- Pin formal PDF cover dates via `--build-date YYYY-MM-DD` or `SOURCE_DATE_EPOCH`.

### 8.4 Network Access, Credentials, and External Actions

- Verify equations, standards, datasheets, and challenge rules against primary sources rather than search snippets.
- Never expose tokens, keys, or credentials in docs, logs, or responses, and never create external issues/PRs/uploads without explicit user instruction.

### 8.5 Pedagogical Code, Third-Party Upstream Sources, Dual Licensing, and Supply-Chain Rules

- **Dual-License Model (`LICENSE` + `LICENSE-DOCS`)**: Original software, numerical cores, scripts, and tests in this repository are licensed under the **MIT License** (`LICENSE`), while original textbook prose, upstream research manuals, pedagogical figures, and compiled HTML/PDF editions are licensed under **CC BY-NC-SA 4.0** (`LICENSE-DOCS`). Third-party components (`NumpyOnlineWPE011`, vendored MathJax 3.2.2 / SRE 4.0.6), public acoustic datasets (`DEMAND` under CC BY-SA 3.0, `LOCATA` excerpts under CC BY 4.0), and local Git-ignored `_downloads/` checkouts retain their respective upstream licenses and are never relicensed.
- Classify every method in `codes/chapters/ch00/COVERAGE*.md` into exactly one status (`Runnable Local Baseline`, `External Reference Implementation`, `Principle Index`, or `Explicitly Excluded`), keeping pedagogical implementations dependency-light and tested against independent hand/analytic expectations.
- Pin every Git upstream in `SOURCES.lock.json` by official HTTPS URL, full 40-char commit hash, verification date, code license, and sparse-selection rules (`fetch_enabled` records local acquisition policy rather than a legal license grant; `source_selection_mismatch` on legacy checkouts is preserved honestly without overwriting user caches).
- Track non-Git source tarballs/ZIPs separately in `ARCHIVE_SOURCES.lock.json` (`ARCHIVE_SOURCE_STATUS.json`), verifying SHA-256 digests and path/size safety prior to extraction without folding them into the Git project count.
- Scrub outer Git environment variables inside acquisition tools, separate source-code licenses from checkpoint/corpus licenses, verify time-indexed tap/buffer/padding/silence behavior when auditing upstream APIs, guard array-slice view aliasing and C pointer strides/types, and verify DC/Nyquist scaling in zero-padded real-signal FFT interpolation.

### 8.6 Chapter Exercises and Audio Assets

- Assign stable IDs (`E00-01`–`E17-13`) to all 352 exercises, keeping chapter prose, worked solutions, CLI return values, unit tests, `COVERAGE*.md`, `codes/chapters/README*.md`, and `05_exercises_and_audio.md` synchronized 1-to-1 (adding exercises never inflates the 388-row algorithm count; importing a module must never trigger experiments, disk writes, or network calls).
- Specify time units, block sizes, update instants, and exact freeze semantics (weights vs. covariance vs. reference history vs. state prediction) on stateful exercises; enforce explicit ordered token lists in WER scorers.
- Keep one shared export gain per comparison group (never per-file peak normalization), separate analytic expectations, float components, and readback `PCM16` integer metrics (`E`, `D`, `E/D`), preserve CC BY-SA 3.0 / CC BY 4.0 attributions alongside staged real-audio copies, default HTML `<audio>` players to non-autoplay with volume-down guidance, and ensure every `--check` mode is strictly read-only.

## 9. Completion Handoff Format

Structure completion reports concisely in order (omitting inapplicable items):

1. **Executive Summary** (2–3 sentences stating what was completed and whether any items remain open);
2. **Key Changes by Theme** (what changed and why);
3. **Clickable File List** (authoritative source files and generated artifacts);
4. **Verification Results** (exact commands/checks executed and pass/fail/not-run outcomes—never write a bare "tested");
5. **Unverified Items & Risks** (missing host tools/permissions and exact follow-up commands);
6. **Primary Sources Consulted** (when external references were used).

## 10. Book-Wide Research, Probabilistic Tracking, and Formula Rendering Boundaries

### 10.1 Research & Page-by-Page QA Ledgers

- Keep items as `Modified` until final compiled HTML/PDF artifacts are inspected, promoting to `Verified` only with concrete page numbers or report paths.
- Synchronize counts separately across algorithms (`388`), exercises (`352`), figures (`87` + `1` room plot), published WAVs (`283`), indexed Git projects (`129`), fetched checkouts (`112`), and archive sources (`2`).
- Verify PDF bookmarks via named destinations (`id` → page & viewport coordinates) rather than text matching that could land on TOC lines.
- Distinguish original paper-author repositories from maintainer re-implementations, preserve unmodified upstream failures separately from isolated compatibility shims, and save post-E-step posteriors separately from post-M-step parameters in iterative EM/BSS routines.

### 10.2 Tracking Timestamps, Circular Statistics, and Upstream State

- Log `state_time`, `measurement_time`, and frame-end data availability time separately so downstream consumers never double-extrapolate already predicted states.
- Wrap circular innovations to $[-180^\circ, 180^\circ)$, flag undefined circular means when resultant length vanishes, preserve exact zero-prior support in particle likelihood shifts, include between-mode mean spread in mixture/IMM variances, and distinguish PHD intensity mass, cardinality distributions, and Bernoulli existence probabilities.

### 10.3 TeX Source Hygiene and Visual Formula Auditing

- Guard string replacements against host-language backslash escapes (`\a`, `\t`, `\r`, `\b`) corrupting TeX commands such as `\approx`, `\tau`, `\rho`, or `\beta`; reject raw C0/DEL control bytes outside `\n`, `\r`, and `\t`.
- When capturing scrolling browser screenshots for visual QA, await layout, webfont readiness, and paint completion after every scroll/jump rather than only on initial page load, visually verifying rendered glyphs rather than trusting non-zero container widths alone.
