# Algorithm Implementation, Industrial Deployment, and Source Code Research Manual

🌐 [简体中文](../research/README.md) ｜ **English**

**[← Back to Book Overview](../../../../chapters_en/00_overview.md)** · **[Chapter-by-Chapter Teaching Code Index](../../README_EN.md)** · **[Algorithm Coverage Matrix](../COVERAGE_EN.md)** · **[129-Item Upstream Version Lockfile](../SOURCES.lock.json)**

This manual picks up immediately after the mathematical derivations in the main chapters and systematically dissects the engineering interface questions that inevitably arise when moving "from blackboard equations to real open-source libraries and embedded firmware": which tensor dimension corresponds to sensors versus frequency bins, how cross-frame states are initialized and frozen, where echo and voice-activity control signals enter the pipeline, which file contains the upstream source entry point, and how to construct minimal controlled experiments to verify implementation correctness.

All material is compiled from original papers, maintainer documentation, and official repositories; verification dates for individual entries are recorded in each thematic guide. See [SOURCES.lock.json](../SOURCES.lock.json) for the pinned version manifest and [COVERAGE_EN.md](../COVERAGE_EN.md) for the per-algorithm coverage matrix. Algorithms sharing the same name across different projects may adopt different normalizations, spectral axis orders, and output reference channels; one must never compare output arrays from two functions directly without first checking their conventions.

This manual comprises this landing page together with the 5 thematic guides listed below, forming 6 research web pages in total; combined with the 19 current tutorial pages (the overview page plus 18 content pages spanning Chapters 0–15 and the 2 appendices), the site contains 25 current content pages. The site additionally retains 8 legacy-route compatibility pages on the same topics, generating 33 HTML files in total. Compatibility pages do not count as new chapters or independent research content.

The lockfile registers 129 Git upstream projects, of which 112 have local fetching enabled with existing local workspaces. As of 2026-11-17, 89 items pass offline verification under the full sparse-selection policy (`source_verified`), while 22 previously fetched workspaces still follow an older exclusion rule and are reported as selection mismatches (`source_selection_mismatch`) and retained as-is; the real recording files in AEC Challenge contain local modifications and remain recorded as a verification failure (1 item). Another 17 items retain source-index metadata only (`index_only`). A selection mismatch neither means that source code was not fetched nor indicates an algorithm execution failure. As a compiled dependency of MeetEval, `kaldialign` has only had its source fetched and has neither been built nor executed. The four-file selection from the TAC authors has been fetched into an independent git-ignored directory; its README-declared non-commercial license, static structural audit, and unexecuted original network boundary are documented in [Appendix B Source Contracts](04_source_reproduction.md#tac-contracts).

The Spatial & Tracking guide begins with binaural cues, SAF HRIR/HRTF conventions, and `libmysofa` data interfaces, then examines SBL, RobustSBL, online subband GSC, SMP-PHAT, and multi-frame association; the AEC, WPE & Separation guide examines online WPE interface comparisons, DNN adaptive control, and long-recording separation; the Industrial Deployment guide examines resampling, loudness, STOI/ESTOI, ViSQOL, frame-wise file I/O, block floating-point arithmetic, and STK fractional-delay state management. The AEC, WPE & Separation guide additionally audits the time-frequency architecture, tensor shapes, and configuration boundaries of the pinned TF-Locoformer source code.

This inventory still includes unexecuted and principle-only entries; exact conditions and method-level outcomes are governed by the individual thematic guides.

Executed interface experiments across three industrial libraries and synthetic-input SBL experiments are accompanied by original driver scripts and [execution records](04_source_reproduction.md).
The reports distinguish nominal checks, state-reset counterexamples, and runs reaching iteration limits; the actual compiled execution of SMP-PHAT separately documents a numerical failure caused by negative-lag index conversion.
These small-scale interface comparisons do not constitute full-dataset paper reproductions or product-level device acceptance.

The AEC, WPE & Separation guide also documents a synchronous AEC experiment running pinned SpeexDSP on a real playback-loopback/microphone recording pair alongside zero-reference and misaligned-reference controls; this is recorded strictly as an observation of that interface on that specific clip, without extrapolating to ground-truth ERLE, double-talk protection, or hardware acceptance.

Readers can also recompute the Chapter 5 same-input analytical beamformer comparisons, Chapter 6 synthetic ground-truth AEC interface comparisons, Chapter 7 online WPE temporal boundaries, Chapter 8 pinned AuxIVA blind estimation and GSS activity mislabeling experiments, and the Appendix B six-position room simulation.

The room exercise is accompanied by [18 white-noise synthetic WAV files and one result figure](../../appendix_b/room_audio/MANIFEST.json), managed separately from the 109 samples in the main manifest. Every item specifies fixed inputs, runnable entry points, and failure or interpretation boundaries; analytical and synthetic figures never stand in as performance numbers for natural speech, physical hardware, or complete meeting separation pipelines.

| Reading Material | Corresponding Chapters | Core Questions Addressed |
|---|---|---|
| [Spatial Processing and Source Tracking](01_spatial_and_tracking.md) | Chapters 1–5, 9, 12, 14 | Broadband localization, robust beamforming, spherical and wearable spatial arrays, industrial acoustic imaging, data association, and multi-target state estimation |
| [Echo Cancellation, Dereverberation, and Speech Separation](02_aec_wpe_separation.md) | Chapters 6–8, 15 | Adaptive step-size control, online states, blind separation and target speaker extraction, multi-channel diarization, neural architectures, and inference dependencies |
| [Industrial Deployment and Evaluation](03_industrial_deployment.md) | Chapters 10–11, 13, Appendix B | Device audio callbacks, clock drift and resampling, embedded fixed-point processing, inference runtimes, and objective scoring |
| [Reproducing Upstream Code After Fetching](04_source_reproduction.md) | Whole Book | Commit verification, input adaptation, minimal controlled experiments, independent numerical checks, and result logging |
| [Exercises and Audio Experiments](05_exercises_and_audio.md) | Chapters 0–11, Extended Topics I–IV, Appendices A/B | Recomputation, licensing, and boundaries for 352 code exercises, 109 main synthetic WAVs, 30 standalone synthetic sets (170 WAVs), and 4 real-recording excerpts/derivatives; see the [Overview Complete Audio Table](../../../../chapters_en/00_overview.md#audio-assets) for per-set entry points |

"Experiments" in this research manual fall into two strictly separated categories. Executed experiments explicitly record inputs, environments, and measured outputs; suggested reader experiments state concrete steps and pass/fail criteria without pre-filling performance numbers. Having upstream source code locally does not imply that installation, training, inference, or product acceptance has been completed on this machine.

## Starting from a Concrete Problem

- Localization peak points in the opposite direction: check the sign convention of time delay, array sensor coordinates, FFT phase sign, and steering-vector conjugation before debating algorithmic resolution.
- Works offline but fails on hardware: check audio clock drift, reference tap point, block size, cross-frame history state, and input gain first; these conditions simultaneously affect multiple downstream modules.
- High separation score but degraded ASR: verify output permutation, reference channel selection, scale normalization, time alignment, and scoring regions before diagnosing speech distortion.
- Linear solver returns finite arrays that violate constraints or disagree with hand calculations: read [Appendix A Solver Contracts](04_source_reproduction.md#appendix-solver-contracts) first, checking data types, numerical rank, residuals, and physical objectives separately.
- Considering replacing a module with a new paper method: align the paper's input assumptions, training supervision, causality, and output target before deciding whether to replace an existing block.
- Looking to reproduce neural research directions in Appendix B: pick one output target from [Eight Task Categories and Version Boundaries](04_source_reproduction.md#appendix-b-reproduction), then verify source commits, checkpoints, datasets, and timing scopes; static structural evidence for TAC does not substitute for running the original network.
- Wanting to recompute a worked answer first or audition a time delay or clipping artifact: select a stable exercise ID or audio group from [Exercises and Audio Experiments](05_exercises_and_audio.md) and inspect waveforms and parameters alongside Figure 29. Synthetic examples are neither natural speech nor substitutes for formal listening tests.

## Common Order for Reading Upstream Source Code

1. Read the entry-point function and configuration schema; record input/output shapes, physical units, sampling rate, reference channel, and state lifetime.
2. Locate the four execution paths: state initialization, nominal update, missing-observation/silence freeze, and state reset.
3. Map the code operations back to the textbook objective function or recursion; optimizers, diagonal loading, power floors, and projections are integral parts of the implementation.
4. Read maintainer unit tests, distinguishing shape smoke tests, numerical accuracy tests, audio perceptual tests, and hardware-in-the-loop tests.
5. Run a side-by-side check on a small-scale textbook input or independent analytical example. Running without crashing and producing numerically trustworthy outputs are two distinct conclusions.

When interpreting results, distinguish the following four evidence tiers according to the question the evidence actually answers:

| Evidence Tier | Supported Conclusion | Existing Records or Verification Entry Points |
|---|---|---|
| Teaching Baseline | Formulas and numerical values are reproducible under specified models, inputs, and boundaries | [Coverage Matrix](../COVERAGE_EN.md), [Exercises and Audio Experiments](05_exercises_and_audio.md) |
| Source Verification | Pinned commits, local files, and selected sparse-checkout scopes match the lockfile | [Git Source Status](../SOURCE_STATUS.json), [Archive Source Status](../ARCHIVE_SOURCE_STATUS.json) |
| Method Execution | Interface behavior or numerical outputs under a specified environment, input, and configuration | [Actual Execution Records and Unexecuted Scopes](04_source_reproduction.md) |
| Device Acceptance | Real-time deadlines, resource budgets, acoustic performance, and fault recovery on target hardware meet specifications | [Device Inspection Steps and Measurement Boundaries](03_industrial_deployment.md); governed strictly by explicitly logged hardware measurements |

Local upstream source status is generated by the fetch tooling; see [Upstream Source Fetching Guide](../upstream/README_EN.md).
Licenses and dependencies belong to each respective upstream project; this manual never merges source files from disparate projects into a single relicensed software package.

Chapter 4 additionally provides the [pinned `doatools` ESPRIT interface comparison](../../ch04/reports/doatools_esprit_reference.json) and [per-trial statistics for Figure 13](../../ch04/reports/figure13_gcc_reverb.json). Bias under the upstream default weighting, unweighted execution, and the independent reference implementation are recorded separately; see Section 14 of the Spatial & Tracking guide for complete numbers, applicability scopes, and source entry points.

Chapter 12 (Extended Topic I) provides the [Acoustic Imaging chapter](../../../../chapters_en/12_acoustic-imaging.md), [layer-by-layer source contracts](01_spatial_and_tracking.md#imaging-contract-audit), and [18 exercises plus five-WAV experiments](05_exercises_and_audio.md#sec-47-1). Outputs represent reference sound-pressure contributions and integrated region quantities; upstream algorithms, teaching examples, and industrial operating conditions are verified separately.

Extended Topic II (Distributed Arrays) cross-links [node objectives and compression conditions](../../../../chapters_en/13_distributed-enhancement.md), [original papers and industrial audio networks](03_industrial_deployment.md#distributed-network-deployment), [pinned source contracts](04_source_reproduction.md#distributed-reproduction), and [20 exercises plus 17 WAVs](05_exercises_and_audio.md#distributed-exercises-audio). Original statistical algorithms, known-ground-truth synthetic audio, source fetching, and full-device execution are documented separately.

Extended Topic III (Spatial and Wearable Arrays) provides the [spherical harmonic capture, binaural matching, and cross-modal fusion chapter](../../../../chapters_en/14_spatial-wearable-arrays.md) alongside [15 exercises and six-WAV experiments](05_exercises_and_audio.md#wearable-exercises-audio). Extended Topic IV (Target Speaker Extraction and Multi-Channel Diarization) provides the [clue conditioning, near-field distance-bubble, and TS-VAD/GSS chapter](../../../../chapters_en/15_target-extraction-diarization.md) alongside [14 exercises and six-WAV experiments](05_exercises_and_audio.md#tse-exercises-audio).

The [WASN synchronization and LibriCSS public evaluation sources](04_source_reproduction.md#overview-source-entrypoints) recommended in the Book Overview have now had bounded sparse selections fetched at pinned commits. Source identity verification, static inspection of original interfaces, and unexecuted dependency/scoring/hardware boundaries are preserved separately; source file counts never stand in for industrial acceptance.


## Mapping Between Documentation, Commands, and Results

- Learning entry points and chapter numbering: [Repository Home](../../../../README_EN.md), [Book Overview](../../../../chapters_en/00_overview.md), and [Chapter 0 Mathematical Foundations](../../../../chapters_en/00_math-foundations.md). Chapters 0–15 form the main text (with Chapters 12–15 serving as Extended Topics I–IV), and Chapters 16–17 are Appendices A and B.
- 352 stable code exercises and 30 standalone synthetic audio sets (32 manifests in total, including the main 109-WAV manifest and 1 real-recording set): see the [Teaching Code Directory](../../README_EN.md#audio-inventory) and the Exercises & Audio guide in this manual; exercise counts, WAV counts, and algorithm counts are tallied separately.
- From source assets to HTML/PDF builds: see the [Build & Publishing Guide](../../../../scripts/README_EN.md#audio-regeneration), which verifies authentic generator sources and manifests before publishing site copies.
- Fetching pinned upstream sources: see the [Upstream Fetching Guide](../upstream/README_EN.md); fetching neither installs nor executes upstream code. Current sparse-selection status and bounded original-method execution results are recorded separately.

Reports cited across these guides serve as evidence for their respective pinned versions and fixed inputs; they do not automatically imply that an entire upstream system is runnable in the current environment. Before re-running any audit, first verify the actual identities of the reporting script, upstream source files, and runtime dependencies; write new outputs only to designated current-report paths or ordinary paths outside the repository, and never overwrite digests or retroactively mark historical reports as passing. Where formal human listening tests or assistive-technology evaluations lack actual logged records, they remain explicitly marked as unexecuted.
