# On-Demand Acquisition of Third-Party Reference Implementations

🌐 [简体中文](README.md) ｜ **English**

**[简体中文](README.md)**

`fetch_upstreams.py` reads `SOURCES.lock.json` in the parent directory (`codes/chapters/ch00/`), downloads only entries marked eligible for automated acquisition, and checks out each working tree at its pinned commit hash. By default, repositories are placed under `_downloads/` in this directory, which is ignored by Git.

`_downloads/` is a local, on-demand cache of third-party upstream source trees; upstream `README` files and other Markdown documents inside `_downloads/` are not part of this textbook or its research manuals, nor are they distributed in this repository's Git commits (see root `LICENSE` and `LICENSE-DOCS` for the dual-license terms governing original repository code and documentation). Both acquisition tools accept `--destination` to target a new directory outside the repository for fresh downloads or offline verification; pass the exact same `--destination` path when verifying. This flag does not migrate existing caches. Existing working trees may contain local edits or research data and are never overwritten or deleted merely to tidy a directory.

The acquisition scripts fetch only the registered primary source repositories—they do not install Python/C++ dependencies, execute build scripts, or recursively assemble multi-repository build environments. After fetching, always inspect the `LICENSE`, `NOTICE`, and dependency documentation at the pinned commit. Acquiring upstream source files does not mean a project has been compiled, executed, or validated against published paper benchmarks.

```bash
.venv/bin/python codes/chapters/ch00/upstream/fetch_upstreams.py --list
.venv/bin/python codes/chapters/ch00/upstream/fetch_upstreams.py --project nara_wpe
.venv/bin/python codes/chapters/ch00/upstream/fetch_upstreams.py --project wasn-platform
.venv/bin/python codes/chapters/ch00/upstream/fetch_upstreams.py --project libricss
.venv/bin/python codes/chapters/ch00/upstream/fetch_upstreams.py --all --report tmp/source-acquisition.json
.venv/bin/python codes/chapters/ch00/upstream/fetch_upstreams.py --verify --report tmp/source-verification.json
```

`--all` processes only entries with `fetch_enabled: true`. When an existing checkout already matches the official remote URL, pinned commit hash, clean worktree state, required entrypoints, and full sparse-selection policy recorded in the lockfile, it is reused directly; otherwise, the directory is preserved untouched and reported as an error. If an individual project fails, the tool continues processing remaining projects and exits with a non-zero status code while recording every failure in the JSON report. `--verify` performs strictly offline local checks without network access; `execution: not_run` in the resulting report explicitly indicates that the upstream code itself was not executed.

Verification reports record the requested source scope alongside both `observed_sparse_patterns` (the rules actually configured in the checkout) and `expected_sparse_patterns` (the current policy). A sparse checkout must match the exact ordered sequence of inclusion paths plus all 32 extension-exclusion rules; missing a rule, adding extra exclusions or inclusions, inserting comments, or altering rule order yields `source_selection_mismatch` without automatically rewriting sparse rules or deleting cached files. Twenty-two older checkouts created on 2026-10-02 still carry the earlier 20-rule exclusion list and are therefore no longer counted as passing under the current full-selection policy, while historical run reports preserve the exact conditions in effect when they were generated. All Git subprocesses scrub outer repository-discovery and configuration environment variables so commands cannot accidentally redirect into this textbook repository; the tool also ignores user-global and system-level Git configs. If a network proxy is required, configure standard environment-level HTTP/HTTPS proxy variables rather than relying on global Git `insteadOf` URL rewrites.

Fresh checkouts enable sparse-checkout mode to omit common model-weight, audio, shared-library, NumPy binary, and archive file extensions; the complete list lives in `OMITTED_EXTENSIONS` inside `fetch_upstreams.py`. Submodules are not fetched recursively, and `GIT_LFS_SKIP_SMUDGE=1` prevents automatic Git LFS binary smudging. Inline model coefficients or Git-tracked binary assets outside those extensions may still exist inside source files, so sparse checkout must never be misrepresented as a guarantee that all third-party non-code assets were excluded. Legacy non-sparse checkouts are left intact and labeled `existing_checkout_preserved`; such a checkout can pass verification under that preservation policy only when its lockfile entry does not mandate an explicit `source_paths` subset, and doing so does not prove that models or datasets are absent from that worktree.

Both acquisition scripts share the standard-library primitives in [`io_contracts.py`](../io_contracts.py) for strict JSON parsing and pre-write report path validation. Lockfile parsing rejects duplicate keys, non-finite numbers (`NaN`/`Infinity`), and floating-point overflow. Whenever `--report` is supplied, the output path is preflight-checked before any download, cache creation, or source extraction begins: paths containing lexical `..`, symlinks in the parent chain, directories, special files, hard-linked targets, or paths colliding with acquisition scripts, shared I/O modules, `_downloads/`, archive caches, or either lockfile are rejected. Standard macOS system symlinks (`/tmp`, `/var`, `/etc` resolving under `/private`) remain supported. Reports are written only when `--report` is explicitly passed—first serializing finite JSON, then writing a temporary file in the same directory and atomically replacing the regular report file; if preflight validation or replacement fails, any pre-existing report remains intact. These checks enforce local non-concurrent safety contracts and do not claim immunity against concurrent TOCTOU filesystem races. Each newly generated report binds the SHA-256 digests of the active acquisition script and shared I/O module under `report_sources`, kept distinct from the lockfile digest and `execution: not_run`; historical reports are never retroactively rewritten to mimic current tool runs.

Isolated local checkouts are well suited for studying large frameworks or copyleft implementations. Redistribution, commercial integration, pretrained model weights, and acoustic corpora remain governed by each upstream project's own license terms. For example, the WebRTC checkout contains only the primary `src` repository rather than a full `depot_tools` build workspace, and the same applies to other multi-repo or submodule-heavy projects. Entries without an established license or restricted to evaluation/academic research retain their explicit index-only or restricted status and must never be lumped together as general open-source software.

On 2026-10-04, bounded text selections of 11 files from `CN-UPB/WASN` and 27 files from `LibriCSS` (source, documentation, and license files) were fetched. The former documents synchronization entrypoints and known-SRO resampling, while the latter documents meeting data preparation and continuous evaluation scoring; neither constitutes an installable end-to-end pipeline. Exact commit hashes, license notes, static interface observations, and non-execution boundaries are documented in [the two Overview source entrypoints](../research_en/04_source_reproduction.md#overview-source-entrypoints).

Before updating any pinned commit in a lockfile, inspect ownership of any local modifications, verify in a fresh destination directory, and never `git add` `_downloads/` into this repository. For algorithm-by-algorithm source reading workflows, see the [Source Reproduction Manual](../research_en/04_source_reproduction.md).


## Official Source Archives with Pinned Digests

Some official upstream releases are distributed as versioned source tarballs or ZIP archives rather than Git repositories. `fetch_archives.py` reads the dedicated [`ARCHIVE_SOURCES.lock.json`](../ARCHIVE_SOURCES.lock.json) lockfile, verifies the SHA-256 digests of both the archive and its official descriptor before extraction, and extracts only the registered subset of regular source files. Both `tar.xz` and `zip` formats are supported: compressed byte digests are verified first, followed by checks on member paths, symlinks, duplicate entries, size ceilings, and file selections. For instance, HARKTOOL5 3.5.0 is extracted into `_downloads/harktool5-3.5.0/` while the raw archive and `.dsc` descriptor are cached in `_downloads/.archive-cache/`; neither location is tracked by Git.

```bash
.venv/bin/python codes/chapters/ch00/upstream/fetch_archives.py --list
.venv/bin/python codes/chapters/ch00/upstream/fetch_archives.py --project harktool5-3.5.0
.venv/bin/python codes/chapters/ch00/upstream/fetch_archives.py --verify --report codes/chapters/ch00/ARCHIVE_SOURCE_STATUS.json
```

Digest verification happens strictly prior to archive parsing or disk extraction. The tool enforces compressed and uncompressed byte limits, rejects path traversal (`..` or absolute paths), symlinks, special device nodes, and duplicate archive members, and only verifies (never overwrites or repairs) pre-existing extracted directories. Offline `--verify` compares every local file byte-for-byte against the verified archive's registered selection, checks for missing or unexpected extra files/directories, and records per-file SHA-256 digests alongside `execution: not_run`. Plain HTTP upstream URLs rely on the pinned SHA-256 digests in the lockfile for content integrity; cryptographic PGP signature verification of release descriptors is not claimed here. Digest provenance and verification dates are recorded in `ARCHIVE_SOURCES.lock.json`.

Archive-based upstreams are counted separately and are never folded into the Git repository total. For HARK's custom license restrictions and the distinction between Infineon configuration headers and proprietary algorithm cores, see [Third-Party Attributions](../THIRD_PARTY_EN.md); neither acquisition script installs, compiles, or executes upstream binaries.


The Vo RFS multi-target tracking source archive is fetched and verified via the same archive workflow:

```bash
.venv/bin/python codes/chapters/ch00/upstream/fetch_archives.py --project vo-rfs-tracking-updated
.venv/bin/python codes/chapters/ch00/upstream/fetch_archives.py --verify --report codes/chapters/ch00/ARCHIVE_SOURCE_STATUS.json
```

This author-provided ZIP archive is fetched strictly for academic research reading as permitted by its notice, retaining 62 selected `.m`/`.c` source files and original headers. Its descriptor is a pinned snapshot of the official author download page rather than a Debian `.dsc`, official checksum file, or cryptographic signature; the SHA-256 digests in the lockfile were computed locally upon acquisition. The full raw ZIP stays solely inside the Git-ignored `.archive-cache/`, while the extracted selection excludes precompiled MEX binaries and `.mat` data; neither MATLAB nor MEX binaries were executed. Keep these two archive entries distinct from the Git repository count.


## Chapter 2 High-Precision Environmental Calibration Reference

The author repository for the humid-air speed-of-sound model is pinned at commit `5c7e6652cf2fd7b4c52e83d8fa3782b3c56e61de`. `fetch_upstreams.py` permits literal ASCII spaces inside sparse-checkout paths while continuing to reject wildcards, newlines, absolute paths, and `..` traversal. All 43 LabVIEW `.vi` files under `source labview files/`, plus `LICENSE` and `ReadMe.txt`, were fetched according to the registered selection, omitting precompiled `.exe` binaries, the separately licensed paper PDF, and any LabVIEW runtime. `.vi` files are graphical LabVIEW source files; fetching them on disk does not imply that block diagrams were inspected in LabVIEW or that the model was executed. Physical assumptions and runtime prerequisites are documented in [Chapter 2 Environmental Parameter Study](../research_en/01_spatial_and_tracking.md#sec-1-2).

```bash
.venv/bin/python codes/chapters/ch00/upstream/fetch_upstreams.py --project speed-of-sound-in-air
```

## Chapter 9 Bounded Text Selections for LOCATA

```bash
.venv/bin/python codes/chapters/ch00/upstream/fetch_upstreams.py --project locata-io
.venv/bin/python codes/chapters/ch00/upstream/fetch_upstreams.py --project locata-eval
```

These two pinned commits fetch 15 and 16 explicitly licensed text files, respectively, strictly for upstream source inspection, excluding unattributed helper scripts, audio recordings, and `.mat` evaluation blobs. Individual ODC-BY 1.0 file headers and the 2-clause BSD license for `munkres` are preserved verbatim; see [`SOURCES.lock.json`](../SOURCES.lock.json) for the full file-by-file manifest. Because omitted third-party helpers are required by the original top-level entrypoints, these selections do not form a standalone runnable evaluator; static coordinate, windowing, and scoring findings are detailed in [Spatial & Tracking Study §61.4](../research_en/01_spatial_and_tracking.md#tracking-model-lifecycle-sphere).


## Chapter 12 (Topic I) Original Author DAMAS / CMF Source Code

```bash
.venv/bin/python -B codes/chapters/ch00/upstream/fetch_upstreams.py --project damas-author
```

Twenty-one GPLv3 text files (74,510 bytes total) from pinned author commit `61987952e2237e6b088a169ee891dd96576f2565` are fetched into an isolated Git-ignored worktree, preserving the full GPLv3 `LICENSE` and source notices while excluding three `.mat` files and precompiled MEX binaries. Our bounded contract audit executes only small Python controls and does not run MATLAB, MEX binaries, or the large-scale paper experiments. Differences among matrix products, full-CSM vs. scan-domain NNLS objectives, and original Python edge cases are analyzed in [Source Reproduction Manual §13](../research_en/04_source_reproduction.md#imaging-author-source).

## Chapter 5 openMHA Binaural Fixed-Beamformer Selection

```bash
.venv/bin/python -B codes/chapters/ch00/upstream/fetch_upstreams.py --project openmha
```

Official tag `v4.18.1` resolves to commit `0b9f087e514f9cc0a901fcde9db17de35f9ec438`. The local worktree fetches 10 pinned text files totaling 5,904,110 bytes: root `COPYING` and `README.md`, the 4-microphone binaural example configuration with its precomputed filter-coefficient table, and the C++/header/Makefile sources for `acSteer` and `steerbf`. Full AGPLv3 license text and file-level notices remain inside the isolated Git-ignored worktree. The coefficient table serves as precomputed input to that configuration; sample WAVs, raw HRIR corpora, and the openMHA runtime are omitted, so the checkout is not an executable openMHA installation.

Acquisition verified the pinned commit, clean worktree state, exact 10-file sparse selection, and raw file bytes; plugins were neither compiled nor executed, and no human listening tests were conducted. When reading this example, start with channel routing in `beamformer.cfg`, then trace filter-weight loading order and output conjugation across `acSteer` and `steerbf`.

During the first full-table verification after adding `openmha`, two entries failed: the sole change in the `sof` worktree was an untracked macOS Finder `.DS_Store` file, which was moved intact outside the repository without touching tracked files or sparse rules. Re-running full verification yielded 76 passing (subsequently 89 passing after the 13 pure-text additions on 2026-11-17), 22 legacy sparse-rule mismatches (`source_selection_mismatch`), 17 `index_only` entries, and 1 local-modification failure in the `aec-challenge` audio cache, still exiting with status code `1` as designed. All 115 pre-existing lockfile records and status entries remained consistent with historical snapshots.

## How Acquisition Entrypoints Map to Verification Reports

The primary Git lockfile registers 129 upstream repositories, while archive sources are tracked separately in `ARCHIVE_SOURCES.lock.json`; the active sparse-selection policy enforces 32 excluded extensions. Total indexed projects, fetched local checkouts, checkouts passing current sparse-policy verification, and methods actually executed are read separately from the lockfiles, `SOURCE_STATUS.json` / `ARCHIVE_SOURCE_STATUS.json`, and per-chapter `reports/*.json` contract files—download counts are never conflated with algorithmic coverage.

```bash
.venv/bin/python codes/chapters/ch00/upstream/fetch_archives.py --verify --report tmp/archive-verification.json
```

The command above verifies the default local Git-ignored archive cache and writes a fresh finite-JSON report to `tmp/`; custom cache locations require passing the matching `--destination` flag and fail honestly if files are missing or altered. Routine spot checks should write to `tmp/` so historical method reports remain untouched. Whenever the repository-wide `SOURCE_STATUS.json` or `ARCHIVE_SOURCE_STATUS.json` needs updating, regenerate it via actual `--verify --report` execution rather than hand-editing status fields.

Download failure, `source_selection_mismatch` on the whole worktree, exact blob verification of specific source files used by a probe, runtime method failure, and `not_run` are distinct outcomes. Bounded function-level probes may record verified blob hashes for the exact upstream files and direct dependencies they touch, without ever upgrading a worktree-wide `source_selection_mismatch` to `passed`. Local edits and `.pyc` bytecode inside shared caches are left undisturbed; isolated virtualenvs or compatibility shims live outside the repository and are logged separately so patched behavior is never misattributed to unmodified upstream code.

Following the current chapter sequence, acoustic-imaging upstream contract probes live under Chapter 12 (`codes.chapters.ch12.examples`) and distributed-array probes live under Chapter 13 (`codes.chapters.ch13.examples`), while Appendices A and B (Chapters 16 and 17) retain module prefixes `appendix_a` and `appendix_b`. Exact CLI commands and optional environment prerequisites are listed in the [Chapter Code Index](../../README_EN.md) and the [Build Toolchain Guide](../../../../scripts/README_EN.md#optional-verification). Run all commands in this directory and across every `README` from the repository root.
