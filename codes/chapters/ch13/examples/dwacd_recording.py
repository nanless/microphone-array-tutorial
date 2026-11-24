r"""Run the pinned original DWACD and resampler on the chapter 13 recordings.

Use the repository's upstream source fetcher for the fixed paderwasn checkout.
In a separate directory, check out the paderbox revision used by paderwasn::

    git clone https://github.com/fgnt/paderbox.git /tmp/paderbox-src
    git -C /tmp/paderbox-src checkout 7fed5b44be2effcedb7a26778ada6c5668b1d6bd

The reference run used this isolated environment::

    uv venv --python 3.13 /tmp/dwacd-env
    uv pip install --python /tmp/dwacd-env/bin/python \
        numpy==2.5.3 scipy==1.18.1 setuptools==81.0.0 wheel==0.48.0 Cython==3.3.1
    uv pip install --python /tmp/dwacd-env/bin/python --no-build-isolation \
        /tmp/paderbox-src

Then run from the tutorial repository root::

    /tmp/dwacd-env/bin/python -B -m codes.chapters.ch13.examples.dwacd_recording \
        --paderbox-source /tmp/paderbox-src --libriwasn LibriWASN_200_0L.zip

This older paderbox imports the former scipy.signal.hamming name while loading
an MFCC module. The runner supplies its alias during import in newer SciPy.
DWACD and the resampler run after the alias is removed, using the original code.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import scipy
from scipy import signal

from codes.chapters.ch04.core import upstream_contracts as source
from codes.chapters.ch13.examples.clock_recovery import (
    BLOCK, FS, _libriwasn_inputs, _relative_lag, _synthetic_signals,
)


PADERBOX_URL = "https://github.com/fgnt/paderbox.git"
PADERBOX_REVISION = "7fed5b44be2effcedb7a26778ada6c5668b1d6bd"
CALIBRATION = 16 * BLOCK
PADERWASN_MODULES = (
    "paderwasn/__init__.py",
    "paderwasn/synchronization/__init__.py",
    "paderwasn/synchronization/sro_estimation.py",
    "paderwasn/synchronization/sync.py",
    "paderwasn/synchronization/time_shift_estimation.py",
    "paderwasn/synchronization/utils.py",
)


def _evaluate(reference, remote, estimator, compensate, true_ppm=None):
    active = np.ones(reference.size, dtype=bool)
    full = estimator(remote, reference, active, active)
    prefix = slice(0, CALIBRATION)
    calibration = estimator(remote[prefix], reference[prefix], active[prefix], active[prefix])
    first_estimate = estimator.settling_time - 1
    calibrated_ppm = float(np.median(calibration[first_estimate:]))
    corrected = compensate(remote, calibrated_ppm)
    # Use the same full blocks before and after, away from the resampler's tail.
    comparisons = []
    for start in range(CALIBRATION, reference.size - BLOCK, BLOCK):
        block = slice(start, start + BLOCK)
        before = _relative_lag(reference[block], remote[block])
        after = _relative_lag(reference[block], corrected[block])
        if before is not None and after is not None:
            comparisons.append((start + BLOCK // 2, before, after))
    centers, before, after = np.asarray(comparisons).T
    return {
        "true_ppm": true_ppm,
        "prefix_fitted_ppm": calibrated_ppm,
        "prefix_ppm_after_settling": calibration[first_estimate:].tolist(),
        "full_ppm_after_settling": full[first_estimate:].tolist(),
        "first_estimate_reference_end_sample": first_estimate * estimator.seg_shift
        + estimator.temp_dist + estimator.seg_len,
        "heldout_center_and_lags": comparisons,
        "heldout_lag_slope_ppm_before": float(np.polyfit(centers, before, 1)[0] * 1e6),
        "heldout_lag_slope_ppm_after": float(np.polyfit(centers, after, 1)[0] * 1e6),
    }


def _try_activity_window(estimator, reference, remote):
    # Use just one window pair so a later active pair cannot fill the output.
    length = estimator.temp_dist + estimator.seg_len
    ref_activity = np.ones(length, dtype=bool)
    remote_activity = np.arange(length) >= estimator.temp_dist
    rate = estimator(remote[:length], reference[:length], remote_activity, ref_activity)
    remote_delayed_count = int(np.sum(remote_activity[:estimator.seg_len]))
    return {
        "settling_pairs": estimator.settling_time,
        "strict_activity_threshold_samples": estimator.src_activity_th,
        "remote_delayed_waveform_activity_samples": remote_delayed_count,
        "remote_original_delayed_activity_slice_samples": length,
        "remote_original_delayed_activity_true_samples": int(remote_activity.sum()),
        "remote_delayed_waveform_passes_threshold": remote_delayed_count > estimator.src_activity_th,
        "original_estimate_ppm": rate.tolist(),
    }


def _check_paderbox(checkout):
    if (source.git(checkout, "remote", "get-url", "origin") != PADERBOX_URL
            or source.git(checkout, "rev-parse", "HEAD") != PADERBOX_REVISION
            or source.git(checkout, "status", "--porcelain", "--untracked-files=all")):
        raise ValueError("expected a clean original paderbox checkout at the pinned revision")
    license_file = source.ordinary_file(checkout / "LICENSE")
    original = source.git(checkout, "rev-parse", f"{PADERBOX_REVISION}:LICENSE")
    if source.git(checkout, "hash-object", "--no-filters", "--", str(license_file)) != original:
        raise ValueError("paderbox license differs from the pinned revision")
    return source.sha(license_file)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--libriwasn", type=Path, required=True)
    parser.add_argument("--paderbox-source", type=Path, required=True, help="clean checkout of the installed paderbox revision")
    parser.add_argument("--paderwasn", type=Path, default=source.CACHE / "paderwasn")
    args = parser.parse_args()
    paderbox = args.paderbox_source.resolve()
    paderbox_license = _check_paderbox(paderbox)
    identity = source.verify_project("paderwasn", args.paderwasn, PADERWASN_MODULES)
    samples, recording = _libriwasn_inputs(args.libriwasn)
    synthetic = list(_synthetic_signals())
    import_alias = not hasattr(signal, "hamming")
    sys.path.insert(0, identity["checkout"])
    try:
        paderbox_record = source.package_preflight("paderbox", paderbox)
        was = source.package_preflight("paderwasn", Path(identity["checkout"]))
        with source.source_imports("paderwasn", was), source.source_imports("paderbox", paderbox_record):
            if import_alias:
                signal.hamming = signal.windows.hamming
            try:
                from paderwasn.synchronization.sro_estimation import DynamicWACD
                from paderwasn.synchronization.sync import compensate_sro
            finally:
                if import_alias:
                    del signal.hamming

            synthetic_results = {f"{ppm:+d}_ppm": _evaluate(ref, sig, DynamicWACD(), compensate_sro, ppm)
                                 for ppm, ref, sig, _ in synthetic}
            pairs = {device: _evaluate(samples[recording["reference"]], sig, DynamicWACD(), compensate_sro)
                     for device, sig in samples.items() if device != recording["reference"]}
            _, ref, sig, _ = synthetic[0]
            activity_test = _try_activity_window(DynamicWACD(settling_time=1), ref, sig)
    finally:
        sys.path.pop(0)
    source.check_unchanged(identity)
    source.package_unchanged(was)
    source.package_unchanged(paderbox_record)
    if _check_paderbox(paderbox) != paderbox_license:
        raise ValueError("paderbox license changed during execution")
    result = {
        "paderwasn": {**identity, "executed_python": was["executed_python_sources"]},
        "paderbox": {"origin": PADERBOX_URL, "revision": PADERBOX_REVISION,
                     "license_sha256": paderbox_license,
                     "executed_python": paderbox_record["executed_python_sources"]},
        "python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
        "scipy_hamming_alias_during_import_only": import_alias,
        "sample_rate_hz": FS,
        "calibration_samples": [0, CALIBRATION],
        "libriwasn": {**recording, "activity_annotation": {
                          "member": "vad/speaker.csv", "speaker": 1, "interval_samples": [762880, 978688],
                          "sha256": "aa82d478ab9aa37027c5e2b0d67ff950994d9f6f3affff08b7e30c447a06284d"},
                      "activity_masks": "all true within the stated single-speaker interval",
                      "pairs": pairs},
        "synthetic": synthetic_results,
        "activity_window_control": activity_test,
    }
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
