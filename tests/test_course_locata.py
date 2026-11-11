"""Check the distributed LOCATA comparison without the development archive."""

from hashlib import sha256
import json
import unittest
import wave
from xml.etree import ElementTree

import numpy as np

from scripts.course_home import prepare_locata as case


class CourseLocataTests(unittest.TestCase):
    def test_real_pcm_measurements_and_provenance(self):
        report = json.loads((case.OUT / "MANIFEST.json").read_text())
        self.assertEqual(report["mvdr"]["calibration_frames"], 41)
        self.assertGreater(report["mvdr"]["min_model_white_noise_gain_db"], 4)
        self.assertLess(report["mvdr"]["max_model_response_error"], 1e-12)
        self.assertLess(report["mvdr"]["strong_demo"]["model_target_blocking_error"], 1e-12)
        self.assertEqual(
            [(name, item["peak_deg"], item["frames"])
             for name, item in report["localization_windows"].items()],
            [("A_only", 35, 40), ("B_only", 16, 54)],
        )
        self.assertEqual(set(report["source_sha256"]), {
            "scripts/course_home/prepare_locata.py",
            "codes/chapters/ch02/core/conventions.py",
            "codes/chapters/ch02/core/spectral.py",
            "codes/chapters/ch03/core/geometry.py",
            "codes/chapters/ch04/core/covariance.py",
            "codes/chapters/ch04/core/doa.py",
            "codes/chapters/ch05/core/beamforming.py",
        })
        for name, digest in report["source_sha256"].items():
            self.assertEqual(sha256((case.ROOT / name).read_bytes()).hexdigest(), digest)
        self.assertEqual(set(report["files_sha256"]), {*case.FILES, "levels.svg"})
        for name, digest in report["files_sha256"].items():
            self.assertEqual(sha256((case.OUT / name).read_bytes()).hexdigest(), digest)
        chart = ElementTree.parse(case.OUT / "levels.svg").getroot()
        self.assertEqual(len(chart.findall("{http://www.w3.org/2000/svg}polyline")), 3)

        streams = []
        for name in case.FILES:
            with wave.open(str(case.OUT / name)) as stream:
                self.assertEqual((stream.getnchannels(), stream.getsampwidth(),
                                  stream.getframerate(), stream.getnframes()),
                                 (1, 2, 48000, 302528))
                streams.append(np.frombuffer(stream.readframes(stream.getnframes()), dtype="<i2"))
        samples = np.asarray(streams)
        self.assertLess(np.max(abs(samples.astype(np.int64))), 32767)
        for name, window in report["comparison"]["windows"].items():
            start, stop = window["sample_range"]
            fit_start, fit_stop = report["mvdr"]["calibration_samples"]
            self.assertTrue(stop <= fit_start or start >= fit_stop)
            segment = samples[:, start:stop].astype(np.int64)
            taper = np.hanning(stop - start)
            # Independent integration via a two-sided spectral projection and
            # the resulting time-domain squared norm.
            f = abs(np.fft.fftfreq(stop - start, 1 / 48000))
            spectrum = np.fft.fft(segment * taper / 32768, axis=-1)
            spectrum[:, (f < 300) | (f > 3400)] = 0
            projected = np.fft.ifft(spectrum, axis=-1).real
            band = np.sum(projected**2, axis=-1) / np.sum(taper**2)
            for index, filename in enumerate(case.FILES):
                self.assertEqual(int(segment[index] @ segment[index]),
                                 window["pcm_sum_squares"][filename])
                self.assertAlmostEqual(float(10 * np.log10(band[index] / band[0])),
                                       window["band_db_vs_reference"][filename], places=11)
            if name == "B_only":
                self.assertLess(10 * np.log10(band[2] / band[1]), -7)
                self.assertLess(10 * np.log10(band[3] / band[1]), -12)
            if name == "A_only":
                self.assertGreater(10 * np.log10(band[2] / band[1]), -1)
                self.assertGreater(10 * np.log10(band[3] / band[1]), -2)


if __name__ == "__main__":
    unittest.main()
