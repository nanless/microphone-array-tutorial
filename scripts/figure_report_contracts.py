"""Fixed numerical-report identities shared by drawing and publication checks.

The registry is deliberately explicit. Runtime report fields or a directory
listing cannot redefine which figure/report pair is expected.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import struct

from codes.chapters.ch00.io_contracts import validate_parent_chain


@dataclass(frozen=True)
class FigureReportContract:
    number: int
    report_relative: str
    png_filename: str
    producer_relative: str


FIGURE_REPORTS = (
    FigureReportContract(13, 'codes/chapters/ch04/reports/figure13_gcc_reverb.json', 'fig13_gcc_reverb.png', 'scripts/make_figures.py'),
    FigureReportContract(20, 'codes/chapters/ch07/reports/figure20_wpe.json', 'fig20_wpe.png', 'scripts/make_figures.py'),
    FigureReportContract(21, 'codes/chapters/ch09/reports/figure21_tracking.json', 'fig21_tracking.png', 'scripts/make_figures.py'),
    FigureReportContract(45, 'codes/chapters/ch02/reports/figure45_stft_convolution.json', 'fig45_stft_convolution.png', 'scripts/make_figures.py'),
    FigureReportContract(56, 'codes/chapters/ch09/reports/figure56_tracking_information.json', 'fig56_tracking_information.png', 'scripts/make_figures.py'),
    FigureReportContract(57, 'codes/chapters/ch10/reports/figure57_engineering_limits.json', 'fig57_engineering_limits.png', 'scripts/make_figures.py'),
    FigureReportContract(59, 'codes/chapters/ch11/reports/figure59_selection_scenarios.json', 'fig59_selection_scenarios.png', 'scripts/make_figures.py'),
    FigureReportContract(60, 'codes/chapters/appendix_a/reports/figure60_weighted_noise.json', 'fig60_weighted_noise.png', 'scripts/make_figures.py'),
    FigureReportContract(61, 'codes/chapters/appendix_b/reports/figure61_equal_drr_response.json', 'fig61_equal_drr_response.png', 'scripts/make_figures.py'),
    FigureReportContract(62, 'codes/chapters/ch12/reports/figure62_imaging_psf.json', 'fig62_imaging_psf.png', 'scripts/make_figures.py'),
    FigureReportContract(63, 'codes/chapters/ch12/reports/figure63_imaging_model_checks.json', 'fig63_imaging_model_checks.png', 'scripts/make_figures.py'),
    FigureReportContract(64, 'codes/chapters/ch12/reports/figure64_imaging_calibration.json', 'fig64_imaging_calibration.png', 'scripts/make_figures.py'),
    FigureReportContract(65, 'codes/chapters/ch13/reports/figure65_distributed_compression.json', 'fig65_distributed_compression.png', 'scripts/make_figures.py'),
    FigureReportContract(66, 'codes/chapters/ch13/reports/figure66_distributed_updates.json', 'fig66_distributed_updates.png', 'scripts/make_figures.py'),
    FigureReportContract(67, 'codes/chapters/ch13/reports/figure67_distributed_transport.json', 'fig67_distributed_transport.png', 'scripts/make_figures.py'),
    FigureReportContract(68, 'codes/chapters/ch05/reports/figure68_phase_reference.json', 'fig68_phase_reference.png', 'scripts/make_beamforming_figures.py'),
    FigureReportContract(69, 'codes/chapters/ch06/reports/figure69_reference_timing.json', 'fig69_reference_timing.png', 'scripts/make_reference_figures.py'),
    FigureReportContract(70, 'codes/chapters/ch07/reports/figure70_prediction_delays.json', 'fig70_prediction_delays.png', 'scripts/make_delay_figures.py'),
    FigureReportContract(71, 'codes/chapters/ch08/reports/figure71_css_polarity_gain.json', 'fig71_css_polarity_gain.png', 'scripts/make_css_figures.py'),
    FigureReportContract(72, 'codes/chapters/ch09/reports/figure72_tracking_lifecycle.json', 'fig72_tracking_lifecycle.png', 'scripts/make_tracking_figures.py'),
    FigureReportContract(73, 'codes/chapters/ch10/reports/figure73_channel_failure.json', 'fig73_channel_failure.png', 'scripts/make_channel_figures.py'),
    FigureReportContract(74, 'codes/chapters/ch11/reports/figure74_selection_evidence.json', 'fig74_selection_evidence.png', 'scripts/make_selection_figures.py'),
    FigureReportContract(75, 'codes/chapters/ch12/reports/figure75_imaging_objectives.json', 'fig75_imaging_objectives.png', 'scripts/make_figures.py'),
    FigureReportContract(76, 'codes/chapters/ch14/reports/figure76_foa_intensity_diffuseness.json', 'fig76_foa_intensity_diffuseness.png', 'scripts/make_wearable_figures.py'),
    FigureReportContract(77, 'codes/chapters/ch14/reports/figure77_bsm_ls_vs_magls.json', 'fig77_bsm_ls_vs_magls.png', 'scripts/make_wearable_figures.py'),
    FigureReportContract(78, 'codes/chapters/ch14/reports/figure78_own_voice_fusion.json', 'fig78_own_voice_fusion.png', 'scripts/make_wearable_figures.py'),
    FigureReportContract(79, 'codes/chapters/ch14/reports/figure79_hear_through_comb.json', 'fig79_hear_through_comb.png', 'scripts/make_wearable_figures.py'),
    FigureReportContract(80, 'codes/chapters/ch15/reports/figure80_tse_conditioning_modes.json', 'fig80_tse_conditioning_modes.png', 'scripts/make_tse_figures.py'),
    FigureReportContract(81, 'codes/chapters/ch15/reports/figure81_nearfield_distance_bubble.json', 'fig81_nearfield_distance_bubble.png', 'scripts/make_tse_figures.py'),
    FigureReportContract(82, 'codes/chapters/ch15/reports/figure82_diarization_tsvad_gss.json', 'fig82_diarization_tsvad_gss.png', 'scripts/make_tse_figures.py'),
    FigureReportContract(83, 'codes/chapters/ch15/reports/figure83_av_asynchrony_occlusion.json', 'fig83_av_asynchrony_occlusion.png', 'scripts/make_tse_figures.py'),
    FigureReportContract(84, 'codes/chapters/ch00/reports/figure84_complex_phasor_fourier.json', 'fig84_complex_phasor_fourier.png', 'scripts/make_math_figures.py'),
    FigureReportContract(85, 'codes/chapters/ch00/reports/figure85_subspace_projection_svd.json', 'fig85_subspace_projection_svd.png', 'scripts/make_math_figures.py'),
    FigureReportContract(86, 'codes/chapters/ch00/reports/figure86_hermitian_rayleigh_gev.json', 'fig86_hermitian_rayleigh_gev.png', 'scripts/make_math_figures.py'),
    FigureReportContract(87, 'codes/chapters/ch00/reports/figure87_wirtinger_probability_bayesian.json', 'fig87_wirtinger_probability_bayesian.png', 'scripts/make_math_figures.py'),
)
REPORT_BY_NUMBER = {entry.number: entry for entry in FIGURE_REPORTS}


def _text_keyword_count(data: bytes, keyword: bytes) -> int:
    """Count a keyword in all three PNG text chunk forms, without decoding it."""
    if not data.startswith(b'\x89PNG\r\n\x1a\n'):
        raise ValueError('not a PNG file')
    offset, count = 8, 0
    while offset < len(data):
        if offset + 12 > len(data):
            raise ValueError('truncated PNG chunk')
        length = struct.unpack_from('>I', data, offset)[0]
        end = offset + 12 + length
        if end > len(data):
            raise ValueError('truncated PNG payload')
        kind = data[offset + 4:offset + 8]
        if kind in (b'tEXt', b'zTXt', b'iTXt'):
            payload = data[offset + 8:offset + 8 + length]
            if payload.split(b'\0', 1)[0] == keyword:
                count += 1
        offset = end
        if kind == b'IEND':
            if offset != len(data):
                raise ValueError('unexpected bytes after PNG IEND')
            return count
    raise ValueError('PNG IEND missing')


def report_bindings_issues(root: Path, figures_dir: Path) -> list[str]:
    """Read-only report-to-PNG byte bindings; numerical checks are separate."""
    from PIL import Image
    root, figures_dir = Path(root), Path(figures_dir)
    issues = []
    for entry in FIGURE_REPORTS:
        label = f'figure{entry.number}'
        try:
            report = validate_parent_chain(root / entry.report_relative)
            picture = validate_parent_chain(figures_dir / entry.png_filename)
            if not report.is_file() or not picture.is_file():
                raise ValueError('report or PNG is missing')
            expected = hashlib.sha256(report.read_bytes()).hexdigest()
            if _text_keyword_count(picture.read_bytes(), b'NumericalReportDigest') != 1:
                raise ValueError('NumericalReportDigest must occur exactly once')
            with Image.open(picture) as image:
                image.verify()
            with Image.open(picture) as image:
                image.load()
                value = image.info.get('NumericalReportDigest')
                if not isinstance(value, str) or value != expected:
                    raise ValueError('NumericalReportDigest does not match actual report bytes')
        except (OSError, ValueError, TypeError, SyntaxError, Image.DecompressionBombError) as exc:
            issues.append(f'{label}: {exc}')
    return issues
