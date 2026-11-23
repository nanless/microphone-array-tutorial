"""Read-only numeric checks for the 35 published figure reports.

The shared figure_report_contracts registry owns filenames. These checks use
small independent formulas, actual PCM integer sums, and saved-event summary
arithmetic. They do not regenerate complete waveforms or replay upstream code.
Stored random observations are checked as observations, not freshly simulated.
"""
from __future__ import annotations

import cmath
import hashlib
import itertools
import math
from pathlib import Path
import struct
import wave

import numpy as np

from codes.chapters.ch00.io_contracts import strict_json_loads, validate_parent_chain
from scripts import quality_check as quality

ROOT = Path(__file__).resolve().parent.parent


def _close(actual, expected, label="value", *, rtol=3e-11, atol=3e-13):
    """Compare typed finite scalars/arrays; discrete values remain exact."""
    if isinstance(expected, dict):
        if type(actual) is not dict or actual.keys() != expected.keys():
            raise ValueError(label + ": fields differ")
        for key, value in expected.items():
            _close(actual[key], value, label + "/" + key, rtol=rtol, atol=atol)
    elif isinstance(expected, (list, tuple)):
        if type(actual) is not list or len(actual) != len(expected):
            raise ValueError(label + ": array shape differs")
        for i, (left, right) in enumerate(zip(actual, expected)):
            _close(left, right, label + f"/{i}", rtol=rtol, atol=atol)
    elif type(expected) is float:
        if (type(actual) is not float or not math.isfinite(actual)
                or not math.isclose(actual, expected, rel_tol=rtol, abs_tol=atol)):
            raise ValueError(label + ": independent numeric value differs")
    elif type(actual) is not type(expected) or actual != expected:
        raise ValueError(label + ": discrete value/type differs")


def _finite_tree(value, label="report"):
    if isinstance(value, dict):
        if 'real' in value and 'imag' in value:
            _z(value)
        for key, item in value.items():
            if type(key) is not str:
                raise ValueError(label + ": non-string key")
            _finite_tree(item, label + "/" + key)
    elif isinstance(value, list):
        for i, item in enumerate(value):
            _finite_tree(item, label + f"/{i}")
    elif type(value) is float and not math.isfinite(value):
        raise ValueError(label + ": nonfinite number")
    elif value is not None and type(value) not in (str, bool, int, float):
        raise ValueError(label + ": invalid JSON type")


def _check_numeric_schema(number, report):
    actual = {}
    def visit(value, path):
        if isinstance(value, dict):
            for key, item in value.items():
                visit(item, path + '/' + key)
        elif isinstance(value, list):
            for item in value:
                visit(item, path + '/*')
        elif type(value) in (bool, int, float):
            actual.setdefault(path, set()).add({bool:'b', int:'i', float:'f'}[type(value)])
    visit(report, '')
    expected = _NUMERIC_SCHEMAS[number]
    if actual != expected:
        missing = sorted(expected.keys()-actual.keys())
        extra = sorted(actual.keys()-expected.keys())
        wrong = sorted(p for p in actual.keys() & expected.keys() if actual[p] != expected[p])
        raise ValueError(f'figure{number} numeric schema differs: missing={missing[:3]}, extra={extra[:3]}, types={wrong[:3]}')


def _read(path):
    return strict_json_loads(validate_parent_chain(path).read_bytes())


def _pcm(folder, name, *, channels=1, samples=32000, rate=16000):
    path = validate_parent_chain(folder / name)
    with wave.open(str(path), "rb") as stream:
        if (stream.getframerate(), stream.getnchannels(), stream.getsampwidth(),
                stream.getnframes(), stream.getcomptype()) != (rate, channels, 2, samples, "NONE"):
            raise ValueError(name + ": PCM format differs")
        data = stream.readframes(samples)
    if len(data) != samples * channels * 2:
        raise ValueError(name + ": truncated PCM")
    values = struct.unpack("<" + str(samples * channels) + "h", data)
    return [values[c::channels] for c in range(channels)]


def _projections(values, frequencies, *, start, stop, rate=16000):
    """Actual PCM coefficients; no gain/delay compensation is applied."""
    return [2 / ((stop-start)*32768) * sum(
        values[n] * cmath.exp(-2j*math.pi*f*n/rate) for n in range(start, stop))
        for f in frequencies]


def _report13(report):
    config = report["configuration"]
    _close(config, {
        "sample_rate_hz":16000,"signal_samples":16000,"frame_samples":1024,
        "direct_delay_samples":9,"true_gcc_lag_samples":-9,"fft_samples":4096,
        "search_lag_samples":[-64,64],"source_lowpass_hz":6000,"phat_additive_epsilon":1e-10,
        "row_drr_db":[6.0,0.0,-6.0],"column_t60_s":[0.1,0.4,0.8],
        "success_tolerance_samples":1,"wilson_z":1.96,"percentile_method":"numpy default linear",
        "additive_noise":False,"source_peak_normalized":True,"frame_window":"rectangular",
        "frame_start_min":0,"frame_start_max_inclusive":14974,"tail_duration_factor":1.1,
        "amplitude_decay_coefficient":6.91,"contrast_definition":"max(g)/(median(abs(g))+1e-12)"}, "configuration")
    s = report["statistics"]
    if set(s) != {"rate","low","high","median","q1","q3","successes","failures",
                  "exceptions","trials","trial_correct","trial_contrasts","seeds"}:
        raise ValueError("figure13 statistics fields differ")
    _close(s["trials"], 150)
    _close(s["seeds"], [2000+17*i for i in range(150)])
    for key in ("trial_correct", "trial_contrasts"):
        if len(s[key]) != 3 or any(len(row) != 3 for row in s[key]):
            raise ValueError("figure13 condition shape differs")
    for i,j in itertools.product(range(3),repeat=2):
        flags, contrasts = s["trial_correct"][i][j], s["trial_contrasts"][i][j]
        if (len(flags) != 150 or any(type(x) is not bool for x in flags)
                or len(contrasts) != 150 or any(type(x) is not float or not math.isfinite(x) or x < 0 for x in contrasts)):
            raise ValueError("figure13 stored observations differ")
        successes = sum(flags); p = successes/150; z = 1.96
        centre = (p+z*z/300)/(1+z*z/150)
        half = z*math.sqrt(p*(1-p)/150+z*z/(4*150**2))/(1+z*z/150)
        ordered = sorted(contrasts)
        def quantile(fraction):
            position = 149*fraction; k = int(position)
            return ordered[k]+(position-k)*(ordered[min(k+1,149)]-ordered[k])
        expected = {"successes":successes,"failures":150-successes,"exceptions":0,
                    "rate":p,"low":max(0.0,centre-half),"high":min(1.0,centre+half),
                    "median":quantile(.5),"q1":quantile(.25),"q3":quantile(.75)}
        for key,value in expected.items():
            if len(s[key]) != 3 or len(s[key][i]) != 3:
                raise ValueError("figure13 summary shape differs")
            _close(s[key][i][j],value,"statistics/"+key)


def _report20(report):
    _close(report["parameters"], {
        "sample_rate_hz":16000,"duration_seconds":1.4,"seed":21001,"rir_length_samples":12800,
        "direct_delay_samples":160,"decay_parameter_seconds":.6,"drr_db":6.0,
        "fft_size":512,"hop_samples":128,"window":"symmetric Hann (numpy.hanning)",
        "center":False,"padding":False,"taps":10,"delay_frames":6,"iterations":3,
        "power_smoothing":"5 frames centered, truncate and renormalize at boundaries",
        "relative_power_floor":1e-5,"relative_loading":1e-6},"parameters")
    _close(report["measurement"], {
        "frequency_bins":257,"frequency_range_hz":[0,8000],"display_frequency_range_hz":[0,2500],
        "frame_count":172,"last_frame_start_ms":1368.0,"active_frames":96,"quiet_frames":61,
        "other_frames":15,"active_threshold_fraction_of_reference_peak":.1,
        "quiet_threshold_fraction_of_reference_peak":.01,
        "nmse_alignment":"one complex gain over all active frame-frequency coefficients"},"measurement")
    if set(report["results"]) != {"input","output"}:
        raise ValueError("figure20 result roles differ")
    for role,row in report["results"].items():
        if set(row) != {"quiet_energy","total_energy","quiet_energy_ratio_db",
                        "active_scale_aligned_nmse_db","active_alignment"}:
            raise ValueError("figure20 requires actual active-alignment sufficient statistics")
        quiet,total = row["quiet_energy"],row["total_energy"]
        if type(quiet) is not float or type(total) is not float or not 0 < quiet <= total:
            raise ValueError("figure20 energy domain differs")
        _close(row["quiet_energy_ratio_db"],10*math.log10(quiet/total),role+" quiet dB")
        stats=row["active_alignment"]
        if set(stats) != {"reference_energy","estimate_energy","reference_conjugate_dot_estimate"}:
            raise ValueError("figure20 alignment statistics fields differ")
        er,ex=stats["reference_energy"],stats["estimate_energy"]
        cross=stats["reference_conjugate_dot_estimate"]
        if (type(er) is not float or type(ex) is not float or er<=0 or ex<=0
                or type(cross) is not list or len(cross)!=2 or any(type(x) is not float for x in cross)):
            raise ValueError("figure20 alignment statistics domain differs")
        residual=1-(math.hypot(*cross)/math.sqrt(er)/math.sqrt(ex))**2
        if not 0 < residual <= 1:
            raise ValueError("figure20 projection residual outside domain")
        _close(row["active_scale_aligned_nmse_db"],10*math.log10(residual),role+" aligned NMSE")


def _report21(report):
    _close(report["seed"],22001)
    _close(report["parameters"],{"frames":120,"particles":200,"observation_std_deg":7,
        "observation_replacement_probability":.2,"likelihood_clutter_probability":.1,
        "angle_bounds_deg":[0,120],"process_std_deg_per_frame":1,
        "velocity_process_std_deg_per_frame":.25,"resample_ess_below":100.0,
        "missing_frame_interval_half_open":[48,58]},"parameters")
    f=report["frames"]
    if set(f)!={"index","truth_deg","observation_deg","estimate_deg","missing","replacement_draw",
                "ess_before_resampling","resampled","reflection_fraction"} or any(len(v)!=120 for v in f.values()):
        raise ValueError("figure21 frame fields/length differ")
    _close(f["index"],list(range(120)))
    _close(f["truth_deg"],[60+40*math.sin(2*math.pi*i/120) for i in range(120)])
    _close(f["missing"],[48<=i<58 for i in range(120)])
    for i in range(120):
        for key in ("replacement_draw","resampled"):
            if type(f[key][i]) is not bool:raise ValueError("figure21 Boolean flag differs")
        if (f["observation_deg"][i] is None) != f["missing"][i]:
            raise ValueError("figure21 missing observation differs")
        if (type(f["estimate_deg"][i]) is not float or not 0<=f["estimate_deg"][i]<=120
                or type(f["reflection_fraction"][i]) is not float or not 0<=f["reflection_fraction"][i]<=1
                or type(f["ess_before_resampling"][i]) is not float or not 0<f["ess_before_resampling"][i]<=200+1e-10):
            raise ValueError("figure21 saved observation/state domain differs")
        _close(f["resampled"][i],f["ess_before_resampling"][i]<100)
    _close(report["counts"],{"observed":110,"missing":10,
        "replacement_draws":sum(f["replacement_draw"]),
        "observed_replacements":sum(r and not m for r,m in zip(f["replacement_draw"],f["missing"])),
        "resampling":sum(f["resampled"]),
        "reflection_at_least_5_percent":sum(x>=.05 for x in f["reflection_fraction"])},"counts")
    def rmse(values, missing):
        errors=[(v-t)**2 for v,t,m in zip(values,f["truth_deg"],f["missing"]) if m==missing]
        return math.sqrt(math.fsum(errors)/len(errors))
    _close(report["scores"],{"raw_observed_rmse_deg":rmse(f["observation_deg"],False),
        "pf_observed_rmse_deg":rmse(f["estimate_deg"],False),
        "pf_missing_rmse_deg":rmse(f["estimate_deg"],True),
        "minimum_ess":min(f["ess_before_resampling"])},"scores")
    phd=report["phd_illustration"]
    _close(phd["not_a_filter_run"],True)
    grid=phd["grid_deg"];intensity=phd["intensity_targets_per_degree"]
    _close(grid,[-90+180*i/399 for i in range(400)],'fixed PHD illustration grid')
    if len(grid)!=len(intensity) or len(grid)<2 or any(b<=a for a,b in zip(grid,grid[1:])):
        raise ValueError("figure21 PHD grid differs")
    if any(type(v) is not float or v<0 for v in intensity):raise ValueError("figure21 intensity domain differs")
    area=math.fsum((b-a)*(x+y)/2 for a,b,x,y in zip(grid,grid[1:],intensity,intensity[1:]))
    _close(phd["integral_targets"],area,"PHD integral")
    _close(phd["integral_targets"],2.,'two-unit illustrative intensity')
    gaussian = [sum(math.exp(-.5*((angle-center)/std)**2)/(std*math.sqrt(2*math.pi))
                    for center,std in ((-25.,4.),(40.,5.))) for angle in grid]
    gaussian_area = math.fsum((b-a)*(x+y)/2
        for a,b,x,y in zip(grid,grid[1:],gaussian,gaussian[1:]))
    # The saved background lies in [0,.0005); its integral is at most .09.
    # Bound every normalized point without replaying the random observations.
    maximum_background_area = .0005*(grid[-1]-grid[0])
    normalizer_lower = gaussian_area
    normalizer_upper = gaussian_area+maximum_background_area
    for value,base in zip(intensity,gaussian):
        lower = 2*base/(gaussian_area+maximum_background_area)
        upper = 2*(base+.0005)/gaussian_area
        if not lower-3e-13 <= value <= upper+3e-13:
            raise ValueError('figure21 intensity differs from fixed Gaussian/background bounds')
        if value == 0:
            if base > 0:
                raise ValueError('figure21 positive Gaussian requires positive intensity')
            continue
        sample_lower = 2*base/value
        sample_upper = 2*(base+.0005)/value
        if not math.isfinite(sample_lower) or not math.isfinite(sample_upper):
            raise ValueError('figure21 background normalizer ratio is not finite')
        normalizer_lower = max(normalizer_lower,sample_lower)
        normalizer_upper = min(normalizer_upper,sample_upper)
    # One T must generate every point: background_i=T*intensity_i/2-gaussian_i.
    # Area(intensity)=2 then makes Area(background)=T-Area(gaussian).
    # Allow only 64 ulps in this O(2) normalizer for Gaussian/quadrature/ratio
    # evaluation; this is not an absolute tolerance on tiny intensities.
    normalizer_slack = 64*math.ulp(gaussian_area+maximum_background_area)
    if normalizer_lower > normalizer_upper+normalizer_slack:
        raise ValueError('figure21 Gaussian/background has no common normalizer')


def _math_report(number,report):
    expected={
      84:{"quadrature_phasor_magnitude":math.sqrt(2),"quadrature_phasor_phase_deg":-45.,
          "recovered_group_delay_us":125.,"parseval_energy":6.,"circular_aliasing_l2_error":.5,
          "zero_padded_max_abs_error":0.,"hann_50pct_max_abs_error":0.,"squared_hann_50pct_ripple":.5},
      85:{"norm_sq_b":2.,"norm_sq_parallel":4/3,"norm_sq_orthogonal":2/3,
          "condition_number_2":1000.,"unregularized_error_norm":10.,
          "tikhonov_eps0p1_error_norm":math.hypot(.1/100.1,.001/.1001),
          "woodbury_interferer_attenuation_ratio":1/7,"schur_determinant":7.},
      86:{"r_nn_eigenvalues":[1.,3.],"max_gev_sinr":4/3,"unwhitened_matched_sinr":1.,
          "whitening_sinr_gain_db":10*math.log10(4/3),"sparse_unique_virtual_lags":7,
          "khatri_rao_vec_max_abs_error":0.},
      87:{"mvdr_output_power":7/19,"lcmv_output_power":3/8,
          "improper_pseudo_variance_abs":math.sqrt(13)/5,"dft_dc_pseudo_variance":1.,
          "dft_nyquist_pseudo_variance":1.,"wiener_mmse":6/11,
          "kalman_step1_pred_var0":5.,"kalman_step1_post_var0":15/8}}
    _close(report["anchors"],expected[number],"anchors",rtol=3e-12,atol=1e-13)
    if number==87:_close(report["seed"],20260089)


def _report76(report):
    _close(report["anchor_points"],{
        "in_phase_0deg":{"intensity_norm":math.sqrt(2),"energy_density":1.5,"diffuseness":1-math.sqrt(2)/1.5},
        "quadrature_90deg":{"intensity_norm":1/math.sqrt(2),"energy_density":1.,"diffuseness":1-1/math.sqrt(2)},
        "out_of_phase_180deg":{"intensity_norm":0.,"energy_density":.5,"diffuseness":1.},
        "open_sphere_first_null_kr":math.pi},"anchor_points")


def _report80(report):
    _close(report["film_input"],[[1.,-.5,2.],[.5,1.5,-1.]])
    _close(report["film_output"],[[1.5,-1.5,3.5],[.5,1.,-.25]])
    leak=[-45+.5*i for i in range(91)]
    _close(report["leak_ratio_db_grid"],leak)
    _close(report["absent_loss_curves_db"],{
        str(tau):[10*math.log10((10**(db/10)+tau)/tau) for db in leak]
        for tau in (1e-4,1e-3,1e-2)})
    alphas=[.2+.05*i for i in range(41)]
    _close(report["scale_alphas"],alphas)
    _close(report["present_neg_sisdr_db"],[-20.]*41)
    _close(report["present_waveform_nmse_db"],[10*math.log10((a-1)**2+.01*a*a) for a in alphas])
    _close(report["bypass_examples"],{"single_fricative_ratio_db":-20.,"single_fricative_bypass":True,
        "dual_overlap_ratio_db":10*math.log10(.25),"dual_overlap_bypass":False})


def _bubble(spacing, regularization=0.):
    """Independent symmetric three-sensor scalar solution, not the LCMV kernel."""
    a=cmath.exp(-2j*math.pi*3430/343*(math.hypot(.4,spacing)-.4))
    if regularization==0:
        side=1/(2*(a.conjugate()-1))
        return [side,-2*side,side]
    cross=2*a+1
    determinant=(3+regularization)**2-abs(cross)**2
    u=(3+regularization)/determinant; v=-cross/determinant
    return [a*u+v,u+v,a*u+v]


def _report81(report):
    radii=[.15+.025*i for i in range(195)]
    spacings=[.01+.005*i for i in range(29)]
    _close(report["radii_m"],radii)
    _close(report["spacings_m"],spacings)
    _close(report["fresnel_m_vs_spacing"],[80*d*d for d in spacings])
    def gain(w,d,r):
        a=cmath.exp(-2j*math.pi*10*(math.hypot(r,d)-r))
        return abs((w[0].conjugate()+w[2].conjugate())*a+w[1].conjugate())**2
    _close(report["radial_gain_db"],{
        "dsb_d10cm":[10*math.log10(max(gain([1/3]*3,.1,r),1e-6)) for r in radii],
        "bubble_d10cm":[10*math.log10(max(gain(_bubble(.1),.1,r),1e-6)) for r in radii],
        "bubble_d2cm_reg5e4":[10*math.log10(max(gain(_bubble(.02,5e-4),.02,r),1e-6)) for r in radii]},
        "radial_gain_db",atol=2e-9)
    wng=[-10*math.log10(sum(abs(w)**2 for w in _bubble(d))) for d in spacings]
    # The original Gram solve loses ~1e-9 dB at the smallest spacing.
    _close(report["wng_db_vs_spacing"],wng,"WNG",rtol=3e-11,atol=2e-8)
    folder=ROOT/'codes/chapters/ch15/tse_audio'
    ref=_pcm(folder,'tse_near_target_ref.wav')[0][1600:30400]
    denom=sum(x*x for x in ref)
    expected={}
    for role in ('farfield_dsb','bubble_d10cm','bubble_d2cm_illcond'):
        actual=_pcm(folder,'tse_'+role+'.wav')[0][1600:30400]
        expected[role]=sum((x-y)**2 for x,y in zip(actual,ref))/denom
    _close(report["pcm_integer_nmse"],expected,"actual PCM NMSE",rtol=1e-13,atol=1e-16)


def _report82(report):
    reference=[[1,1,1,1,0,0],[0,0,1,1,1,0]]
    _close(report["ref_activity"],reference)
    _close(report["single_label_bound"],{"overlap_frames":2,"min_miss_frames":2,
        "ref_speaker_frames":7,"min_miss_s":1.,"ref_speaker_time_s":3.5,"der_lower_bound":2/7})
    spatial=[[.95,.92,.82,.80,-.40,0.],[-.45,-.40,.78,.81,.94,-.05]]
    embeddings=[[.88,.86,.48,.47,.12,.08],[.10,.14,.46,.49,.89,.06]]
    def posterior(k):
        return [[1/(1+math.exp(-(k*(spatial[c][t]-.25)+4*(embeddings[c][t]-.45))))
                 for t in range(6)] for c in range(2)]
    joint=posterior(3)
    _close(report["posteriors_spectral_only"],posterior(0))
    _close(report["posteriors_joint"],joint)
    _close(report["gss_target_exclusive_mask"],[a*(1-b) for a,b in zip(*joint)])
    _close(report["gss_interference_mask"],[1-a for a in joint[0]])
    hypotheses={"single_bound":[[1,1,1,1,0,0],[0,0,0,0,1,0]],
        "imperfect_2c":[[1,1,1,0,1,0],[0,0,0,1,0,1]],
        "over_100":[[0,0,0,0,1,1],[0,0,1,0,1,1],[1,1,0,0,1,1]],
        "tsvad_joint":[[int(x>=.5) for x in row] for row in joint]}
    if set(report["der_cases"])!=set(hypotheses):raise ValueError("figure82 DER cases differ")
    for name,hyp in hypotheses.items():
        size=max(2,len(hyp)); ref=reference+[[0]*6 for _ in range(size-2)]
        hyp=hyp+[[0]*6 for _ in range(size-len(hyp))]
        scores={perm:sum(ref[c][t]*hyp[perm[c]][t] for c in range(size) for t in range(6))
                for perm in itertools.permutations(range(size))}
        row=report["der_cases"][name];perm=row["optimal_permutation"]
        if (type(perm)is not list or any(type(x)is not int for x in perm)
                or tuple(perm)not in scores or scores[tuple(perm)]!=max(scores.values())):
            raise ValueError("figure82 mapping is not an optimal speaker assignment")
        matched=max(scores.values()); miss=fa=conf=0
        for t in range(6):
            nr=sum(r[t] for r in ref);nh=sum(h[t] for h in hyp)
            common=sum(ref[c][t]*hyp[perm[c]][t] for c in range(size))
            miss+=max(nr-nh,0);fa+=max(nh-nr,0);conf+=min(nr,nh)-common
        der=(miss+fa+conf)/7
        _close(row,{"miss_frames":miss,"false_alarm_frames":fa,"speaker_confusion_frames":conf,
            "correct_matched_frames":matched,"ref_speaker_frames":7,"miss_s":miss*.5,
            "false_alarm_s":fa*.5,"speaker_confusion_s":conf*.5,"ref_speaker_time_s":3.5,
            "recording_duration_s":3.,"der":der,"exceeds_unit_der":der>1,
            "optimal_permutation":perm},"der_cases/"+name)


def _online_map(offset,frames):
    # Integer microsecond ticks keep exact equality boundaries independent of
    # the producer's searchsorted/float implementation.
    decisions=[16000*i+8000 for i in range(frames)]
    arrivals=[40000*m+20000+offset for m in range(40)]
    indices=[max((m for m,t in enumerate(arrivals) if t<=decision),default=-1) for decision in decisions]
    return indices,[i>=0 for i in indices],[arrivals[i]/1e6 if i>=0 else 0. for i in indices]


def _report83(report):
    _close(report["alignment_mode"],"online_last_arrived")
    _close(report["audio_decision_times_s"],[.016*i+.008 for i in range(15)])
    for prefix,offset in (("sync",0),("lag120ms",120000)):
        indices,valid,arrivals=_online_map(offset,15)
        _close(report[prefix+"_video_indices"],indices)
        _close(report[prefix+"_valid_video_mask"],valid)
        _close(report[prefix+"_video_arrival_times_s"],arrivals)
    lags=[5.*i for i in range(33)];loss=[]
    for lag in lags:
        indices,valid,_=_online_map(int(lag*1000),30)
        # Visual activity is sampled at the selected exposure centre, with
        # the consonant beginning at audio frame10's centre (.168s).
        loss.append(100*sum(not valid[i] or (indices[i]*.04+.02)<.168 for i in range(10,18))/8)
    _close(report["lags_ms"],lags)
    _close(report["consonant_loss_percent"],loss)
    alpha=[.92,.90,.88,.85,.55,.35,.12,.05,.05,.08,.65,.88]
    av=[.96,.95,.94,.93,.72,.58,.15,.04,.05,.10,.78,.94]
    _close(report["occlusion_alpha_v"],alpha)
    _close(report["occlusion_mask_av"],av)
    _close(report["occlusion_mask_fused"],[a*v+(1-a)*.82 for a,v in zip(alpha,av)])
    _close(report["architecture_chunk_latency_s"],{"streaming_tse":.032,"chunked_joint":1.2,"offline_gss":4.5})
    _close(report["architecture_latency_scope"],"Assigned teaching budgets; not measured model latency or universal minimum context.")


def _report77(report):
    _close(report["anchor_3dir"],{"ls_magnitudes":[math.sqrt(.5),0.,math.sqrt(.5)],
                                     "magls_magnitudes":[1.,1.,1.]})
    case=report["circular_4khz_72dir"]
    _close(case["parameters"],{"frequency_hz":4000.,"radius_m":.025,"sound_speed_m_s":343.,
        "direction_count":72,"regularization":.05,"iterations":25})
    angles=np.arange(72)*math.pi/36
    mic=np.arange(4)*math.pi/2
    kr=2*math.pi*4000*.025/343
    transfer=np.exp(1j*kr*np.cos(mic[:,None]-angles[None,:]))
    target=(1+.55*np.sin(angles))*np.exp(2j*kr*np.sin(angles))
    magnitude=np.abs(target)
    normal=transfer@transfer.conj().T+.05*np.eye(4)
    ls=np.linalg.solve(normal,transfer@target.conj())
    phasor=np.exp(1j*np.angle(target));previous=math.inf;history=[];weight=np.zeros(4,complex)
    for _ in range(25):
        candidate=np.linalg.solve(normal,transfer@(magnitude*phasor.conj()))
        response=candidate.conj()@transfer
        objective=float(np.sum((np.abs(response)-magnitude)**2)+.05*np.vdot(candidate,candidate).real)
        if objective<=previous:
            weight=candidate; previous=objective
            nonzero=np.abs(response)>0
            phasor=phasor.copy();phasor[nonzero]=response[nonzero]/np.abs(response[nonzero])
        history.append(previous)
    _close(case["objective_history"],history,"MagLS sum objective",rtol=2e-10,atol=2e-12)
    _close(case["monotonically_non_increasing"],True)
    for prefix,w in (("ls",ls),("magls",weight)):
        _close(case[prefix+"_weights"],{"real":w.real.tolist(),"imag":w.imag.tolist()},prefix+" weights",rtol=2e-10,atol=2e-12)
        mags=np.abs(w.conj()@transfer)
        _close(case[prefix+"_magnitudes"],mags.tolist(),prefix+" response",rtol=2e-10,atol=2e-12)
        _close(case[prefix+"_weighted_mag_mse"],float(np.mean((mags-magnitude)**2)),prefix+" mean residual",rtol=2e-10,atol=2e-12)
    _close(case["target_magnitudes"],magnitude.tolist())
    mixing={}
    for beta in (0.,.2,.5,1.):
        a=(1+beta)/2;b=(1-beta)/2
        var=a*a+b*b+.4*a*b;cross=2*a*b+.2*(a*a+b*b)
        mixing['beta_'+str(beta)]={"ic_out":cross/var,"snr_gain_db":-10*math.log10(var)}
    _close(report["reference_mixing_points"],mixing)


def _report78(report):
    _close(report["anchor_250hz_sinr_db"],10*math.log10(1/4+4/.04))
    _close(report["anchor_4000hz_sinr_db"],10*math.log10(1/.05+.01/.04))
    names=('wearable_air_wind_2ch','wearable_inear_occluded','wearable_inear_equalized','wearable_fused_output')
    folder=ROOT/'codes/chapters/ch14/wearable_audio'
    reference=_pcm(folder,'wearable_clean_ref.wav')[0][1600:30400]
    expected={}
    for name in names:
        values=_pcm(folder,name+'.wav',channels=2 if name==names[0] else 1)[0][1600:30400]
        expected[name]=sum((x-y)**2 for x,y in zip(values,reference))/(28800*32768**2)
    _close(report["pcm_reference_mse_ch0"],expected,rtol=1e-13,atol=1e-16)
    low=.02*(1/401)**2+.0002*(200/401)**2
    high=.00025*(80/81)**2+.0002*(10/81)**2
    _close(report["analytic_reference_mse_ch0"],dict(zip(names,[.02025,.00945,.02005,low+high])))


def _report79(report):
    _close(report["comb_notches_2ms_up_to_4500hz"],[250.+500*i for i in range(9)])
    folder=ROOT/'codes/chapters/ch14/wearable_audio'
    for field,name in (("clean_ref","wearable_clean_ref.wav"),("hearthrough_comb","wearable_hearthrough_comb.wav")):
        actual=_pcm(folder,name)[0]
        phasors=_projections(actual,[250,4000],start=1600,stop=30400)
        _close(report['pcm_tone_magnitudes_'+field],[abs(z) for z in phasors],field+" actual PCM phasors",atol=2e-12)


def _numeric_view(value):
    """Drop prose while retaining every numeric/Boolean parameter and its shape."""
    if isinstance(value,dict):
        return {k:v for k,x in value.items() if (v:=_numeric_view(x)) is not None}
    if isinstance(value,list):
        items=[_numeric_view(x) for x in value]
        return items if any(x is not None for x in items) else None
    return value if type(value) in (int,float,bool) else None


def _integer_power(parts,start,stop,reference=None,*,phase=False):
    windows=[v[start:stop] for v in parts];n=stop-start;channels=len(parts)
    energies=[sum(x*x for x in row) for row in windows];denom=n*32768**2
    expected={"scoring_interval_samples":[start,stop],"samples_per_channel":n,"channels":channels,
        "integer_squared_sum_E_per_channel":energies,"integer_denominator_D_per_channel":denom,
        "mean_square_per_channel":[x/denom for x in energies],
        "integer_squared_sum_E_all_channels":sum(energies),
        "integer_denominator_D_all_channels":denom*channels,
        "mean_square_all_channels":sum(energies)/(denom*channels)}
    if reference is not None:
        truth=[v[start:stop] for v in reference]
        refs=[sum(x*x for x in row) for row in truth]
        errors=[sum((x-y)**2 for x,y in zip(row,ref)) for row,ref in zip(windows,truth)]
        if phase:
            expected.update(integer_reference_squared_sum=refs[0],integer_error_squared_sum_per_channel=errors,
                total_reference_mse_per_channel=[x/denom for x in errors],
                normalized_squared_reference_error_per_channel=[x/r for x,r in zip(errors,refs)])
        else:
            cross=[sum(x*y for x,y in zip(row,ref)) for row,ref in zip(windows,truth)]
            expected.update(integer_reference_squared_sum_per_channel=refs,
                integer_reference_error_squared_sum_per_channel=errors,
                integer_output_reference_cross_sum_per_channel=cross,
                reference_error_mean_square_per_channel=[x/denom for x in errors],
                reference_NMSE_per_channel=[x/r for x,r in zip(errors,refs)],
                projection_gain_per_channel=[x/r for x,r in zip(cross,refs)],pcm_decode_divisor=32768)
    return expected


def _report68(report):
    _close(_numeric_view(report['parameters']),{
        'sample_rate_hz':16000,'samples_per_channel':32000,'frequencies_hz':[500.,1500.],
        'source_amplitude_per_frequency':.1,'source_duration_s':2.,'fade_samples':320,'common_export_gain':1.,
        'base_weights_real_imag':[[.5,0.],[.5,0.]],
        'high_weight_multiplier_real_imag':{'reference':[1.,0.],'flip':[-1.,0.],'quadrature':[0.,1.]},
        'propagation_delay_samples':0,'algorithm_delay_samples':0,'scoring_interval_samples':[2400,29600],
        'samples_in_scoring_window':27200,'scoring_cycles_per_frequency':[850,2550]})
    folder=ROOT/'codes/chapters/ch05/phase_audio'
    pcm={role:_pcm(folder,'phase_'+role+'.wav') for role in ('reference','flip','quadrature')}
    reference=pcm['reference'][0];renergy=sum(x*x for x in reference[2400:29600])/32768**2
    def pairs(values):return [[z.real,z.imag] for z in values]
    rphasor=_projections(reference,[500,1500],start=2400,stop=29600)
    for role,high,error in (('reference',1+0j,0.),('flip',-1+0j,.02),('quadrature',-1j,.01)):
        row=report['samples'][role];response=[1+0j,high]
        _close(row['analytic'],{'mean_square':.01,'reference_mean_square':.01,'total_reference_mse':error,
            'normalized_squared_reference_error':error/.01,'response_real_imag':pairs(response),
            'common_pure_delay_possible':role=='reference',
            'delay_control':'preserving 500Hz requires 500*tau integer; then 1500*tau is also integer'})
        _close(row['file'],'phase_'+role+'.wav')
        _close(row['pcm_integer_measurements'],_integer_power(pcm[role],2400,29600,pcm['reference'],phase=True),role+' PCM integers')
        period=[.1*(math.cos(2*math.pi*i/32)+(high*cmath.exp(6j*math.pi*i/32)).real) for i in range(32)]
        base={'scoring_interval_samples':[2400,29600],'samples_per_channel':27200,'channels':1,
            'squared_sum':272.,'reference_squared_sum':272.,'error_squared_sum':27200*error,
            'mean_square':.01,'total_reference_mse':error,'normalized_squared_reference_error':error/.01,
            'phasor_real_imag':pairs([.1+0j,.1*high]),'reference_phasor_real_imag':pairs([.1+0j,.1+0j]),
            'response_real_imag':pairs(response),'real_projection_columns':4,
            'projection_scope':'known integer-cycle cos/sin diagnostic; no correction applied','peak':max(map(abs,period))}
        _close(row['float_measurements'],base,role+' float finite-tone controls')
        values=pcm[role][0];energy=sum(x*x for x in values[2400:29600])/32768**2
        err=sum((x-y)**2 for x,y in zip(values[2400:29600],reference[2400:29600]))/32768**2
        phasor=_projections(values,[500,1500],start=2400,stop=29600)
        base.update(squared_sum=energy,reference_squared_sum=renergy,error_squared_sum=err,
            mean_square=energy/27200,total_reference_mse=err/27200,
            normalized_squared_reference_error=err/renergy,phasor_real_imag=pairs(phasor),
            reference_phasor_real_imag=pairs(rphasor),response_real_imag=pairs([x/y for x,y in zip(phasor,rphasor)]),
            peak=max(map(abs,values))/32768)
        _close(row['pcm_measurements'],base,role+' actual PCM projections')


def _report69(report):
    windows={'step':(16000,17600),'stable':(19200,28800),'tail':(32000,33600)}
    _close(_numeric_view(report['parameters']),{'sample_rate_hz':16000,'source_samples':32000,
        'samples_per_channel':33600,'convolution_tail_samples':1600,'frequencies_hz':[500.,730.],
        'source_amplitudes':[.08,.06],'fade_samples':320,'playback_gain_intervals':[[0,16000,1.],[16000,33600,2.]],
        'true_path_delays_samples':[0,1600],'true_path_coefficients':[.7,.35],'filter_length':1601,
        'nlms_step_size':.5,'nlms_epsilon':1e-8,'common_export_gain':1.,
        'scoring_intervals_samples':{k:list(v) for k,v in windows.items()},
        'scoring_cycles_per_frequency':{'step':[50,73],'stable':[300,438]},'ncc_frame_samples':160,
        'ncc_activity_rms':.001,'ncc_coherence_threshold':.85})
    folder=ROOT/'codes/chapters/ch06/reference_audio'
    names=('early','late','echo','early_residual','wrong_gain_residual','late_residual')
    values={name:_pcm(folder,'reference_'+name+'.wav',samples=33600) for name in names}
    # Evaluate only the three published scoring windows of the sparse formula;
    # this is not complete generation/quantization replay.
    def source(n):
        if n<0 or n>=32000:return 0.
        env=1.
        if n<320:env=math.sin(math.pi*n/(2*319))**2
        elif n>=31680:env=math.sin(math.pi*(31999-n)/(2*319))**2
        return env*(.08*math.cos(2*math.pi*50*(n%1600)/1600)+.06*math.sin(2*math.pi*73*(n%1600)/1600))
    def amplitude(name,n):
        early=source(n);prev=source(n-1600);g=2. if n>=16000 else 1.
        gp=2. if n-1600>=16000 else 1.
        echo=.7*g*early+.35*gp*prev;prediction=.7*early+.35*prev
        return {'early':early,'late':g*early,'echo':echo,'early_residual':echo-prediction,
            'wrong_gain_residual':echo-g*prediction,'late_residual':0.}[name]
    coefficients={'early':(1.,1.),'late':(2.,2.),'echo':(1.75,2.1),
        'early_residual':(.7,1.05),'wrong_gain_residual':(-.35,0.),'late_residual':(0.,0.)}
    for role in names:
        row=report['samples'][role];_close(row['file'],'reference_'+role+'.wav')
        analytic=row['analytic'];_close(analytic['source_constant_envelope_mean_square'],.005)
        for key,multiplier in zip(('step','stable'),coefficients[role]):
            _close(analytic[key],{'signed_source_multiplier':multiplier,'mean_square':.005*multiplier**2})
        _close(analytic['tail']['ideal_zero'],role in ('early','late','wrong_gain_residual','late_residual'))
        for key,(start,stop) in windows.items():
            parts=values[role];win=parts[0][start:stop];n=stop-start
            _close(row['pcm_integer_measurements'][key],_integer_power(parts,start,stop),role+'/'+key+' integer')
            base={'scoring_interval_samples':[start,stop],'samples_per_channel':n,'channels':1,
                'mean_square':sum(x*x for x in win)/(n*32768**2),'peak':max(map(abs,win))/32768}
            _close(row['pcm_measurements'][key],base,role+'/'+key+' actual PCM')
            model=[amplitude(role,i) for i in range(start,stop)]
            base.update(mean_square=math.fsum(x*x for x in model)/n,peak=max(map(abs,model)))
            _close(row['float_measurements'][key],base,role+'/'+key+' sparse-model score')
    echo=values['echo'][0][32000:];late=values['late'][0][32000:]
    if any(late) or not any(echo):raise ValueError('figure69 complete pure-echo tail differs')
    rms=[]
    for i in range(10):
        block=echo[160*i:160*(i+1)];n=160
        variance=(n*sum(x*x for x in block)-sum(block)**2)/(n*n*32768**2)
        rms.append(math.sqrt(variance))
    _close(report['plotted_tail_centered_rms'],rms,'actual PCM centred tail RMS',rtol=1e-13,atol=1e-16)
    _close(report['pcm_tail_activity'],{'input':'actual decoded late-reference and true-echo PCM',
        'frame_samples':160,'activity_rms':.001,'coherence_threshold':.85,'tail_interval_samples':[32000,33600],
        'tail_frames':10,'tail_states':[2]*10,'tail_labels':['near_only']*10,'tail_absolute_ncc':[0.]*10,
        'tail_ncc_defined':[False]*10,'tail_near_end_truth_present':[False]*10,
        'available_after_samples':list(range(32160,33601,160)),
        'interpretation':'zero NCC is an undefined zero-variance placeholder; near_only is activity label, not component truth'})


def _report70(report):
    _close(_numeric_view(report['parameters']),{'sample_rate_hz':16000,'hop_samples':128,
        'source_samples':32256,'samples_per_channel':33024,'source_period_samples':1536,
        'source_period_frames':12,'source_periods':21,'pulse_start_frames_in_period':[0,4],
        'pulse_samples':128,'frequency_hz':500.,'amplitude':.1,
        'observation_delays_samples':{'reference':384,'auxiliary':0},
        'common_auxiliary_history_delay_samples':384,'aligned_auxiliary_history_delay_samples':768,
        'reference_protection_delay_samples':384,'source_time_gaps_samples':{'common':0,'aligned':384},
        'physical_propagation_tail_samples':384,'maximum_history_padding_samples':768,
        'complete_retained_support_samples':[0,33024],'scoring_interval_samples':[1536,30720],
        'scoring_samples_per_channel':29184,'scoring_source_periods':19,'common_export_gain':1.})
    _close(report['scoring_interval_samples'],[1536,30720])
    _close(report['source_time_control'],{'receiver_indices':[6,7,8,9],'reference':[0,1,0,0],
        'common_history':[0,1,0,0],'aligned_history':[1,0,0,0]})
    _close(report['power_decay_control'],{'T60_seconds':.6,'hop_seconds':.008,'hops':list(range(76)),
        'relative_power':[10**(-.08*i) for i in range(76)],
        'wrong_amplitude_factor_used_as_power':[10**(-.04*i) for i in range(76)]})
    folder=ROOT/'codes/chapters/ch07/delay_audio'
    reference=_pcm(folder,'delay_reference.wav',samples=33024)[0][1536:30720]
    denominator=sum(x*x for x in reference)
    _close(denominator,9715561194)
    for role in ('common_residual','aligned_residual'):
        values=_pcm(folder,'delay_'+role+'.wav',samples=33024)[0][1536:30720]
        error=sum((x-y)**2 for x,y in zip(values,reference))
        _close(report['pcm_integer_errors'][role],{'E':error,'D':denominator,'nmse':error/denominator})


def _report71(report):
    _close(_numeric_view(report['parameters']),{'sample_rate_hz':16000,'samples_per_channel':32000,
        'frequencies_hz':[500.,1500.],'amplitudes':[.1,.1],'fade_samples_each_end':320,
        'first_block_samples':[0,19200],'second_block_samples':[12800,32000],
        'overlap_samples':[12800,19200],'overlap_samples_per_channel':6400,
        'overlap_cycles_per_channel':[200,600],'second_block_multiplier':-2.,
        'scoring_intervals_samples':{'primary':[2400,29600],'post_overlap':[19200,29600]},
        'scoring_samples_per_channel':{'primary':27200,'post_overlap':10400},
        'post_overlap_cycles_per_channel':[325,975],'common_export_gain':1.})
    _close(report['waveform_display_interval_samples'],[14800,15120])
    overlap=report['overlap_fit']
    _close(overlap['float_application_gain'],[-.5,-.5])
    # The two observed unexported overlap blocks have a32-sample common period.
    # Quantize only this stated diagnostic period, not the complete asset.
    gains=[]
    for harmonic in (1,3):
        sine=[math.sin(2*math.pi*harmonic*i/32) for i in range(32)]
        previous=[round(.1*x*32768) for x in sine];current=[round(-.2*x*32768) for x in sine]
        gains.append(sum(x*y for x,y in zip(previous,current))/sum(x*x for x in current))
    _close(overlap['decoded_pcm_overlap_application_gain'],gains)
    _close(overlap['matching'],{'status':'matched','reason':'matched','current_indices_for_previous':[0,1],
        'absolute_centered_correlation':[[1.,0.],[0.,1.]],'candidate_mean_scores_identity_swap':[1.,0.],
        'margin':1.,'thresholds':{'minimum_centered_rms':1e-8,'minimum_correlation':.5,'minimum_margin':.05}})
    folder=ROOT/'codes/chapters/ch08/css_audio'
    reference=_pcm(folder,'css_reference.wav',channels=2)
    if set(report['pcm_integer_measurements'])!={'reference','naive','polarity','corrected'}:
        raise ValueError('figure71 output roles differ')
    for role,measurements in report['pcm_integer_measurements'].items():
        values=_pcm(folder,'css_'+role+'.wav',channels=2)
        expected={key:_integer_power(values,*window,reference) for key,window in
            {'primary':(2400,29600),'post_overlap':(19200,29600)}.items()}
        _close(measurements,expected,role+' actual PCM',rtol=1e-13,atol=1e-16)


def _report45(report):
    from fractions import Fraction as F
    x=list(map(F,[1,2,3,4]));h=[F(1),F(0),F(1,2)];window=[F(0),F(1,2),F(1),F(1,2)]
    denominator=[F(0)]*6;outputs={name:[F(0)]*6 for name in ('identity','fft4','fft8_keep4')}
    for start in (-2,0,2,4):
        frame=[(x[start+i] if 0<=start+i<4 else F(0))*window[i] for i in range(4)]
        for i in range(4):
            index=start+i
            if not 0<=index<6:continue
            denominator[index]+=window[i]**2
            outputs['identity'][index]+=frame[i]*window[i]
            outputs['fft4'][index]+=sum(frame[(i-k)%4]*h[k] for k in range(3))*window[i]
            outputs['fft8_keep4'][index]+=sum(frame[i-k]*h[k] for k in range(3) if 0<=i-k<4)*window[i]
    linear=[sum(x[i-k]*h[k] for k in range(3) if 0<=i-k<4) for i in range(6)]
    final={name:[float(value/d) for value,d in zip(values,denominator)] for name,values in outputs.items()}
    for key,value in {'input':[1.,2.,3.,4.],'impulse_response':[1.,0.,.5],'periodic_hann':[0.,.5,1.,.5],
        'hop_samples':2,'frame_starts_samples':[-2,0,2,4],'retained_output_samples':6,
        'synthesis_denominator':list(map(float,denominator)),'linear_full':list(map(float,linear)),
        'wola_outputs':final,'errors_against_linear':{name:[a-float(b) for a,b in zip(final[name],linear)]
            for name in ('fft4','fft8_keep4')}}.items():
        _close(report[key],value,'figure45/'+key,atol=1e-13)


def _source_roles(number):
    common={'scripts/make_figures.py'}
    if number in (62,63,64,75):return set(quality.IMAGING_FIGURE_SOURCES)
    if number in (65,67):return set(quality.DISTRIBUTED_AUDIO_SOURCES)|common
    if number==66:return common|{'codes/chapters/ch13/core/distributed.py','codes/chapters/ch02/core/conventions.py',
        'codes/chapters/ch04/core/covariance.py','codes/chapters/ch00/io_contracts.py'}
    if number==68:return common|{'scripts/make_beamforming_figures.py','codes/chapters/ch05/core/phase_audio.py',
        'codes/chapters/ch05/examples/generate_phase_audio.py','codes/chapters/ch02/core/conventions.py',
        'codes/chapters/ch00/core/audio_samples.py','codes/chapters/ch00/io_contracts.py'}
    if number==69:return common|{'scripts/make_reference_figures.py','codes/chapters/ch06/core/reference_audio.py',
        'codes/chapters/ch06/examples/generate_reference_audio.py','codes/chapters/ch06/core/aec.py',
        'codes/chapters/ch06/core/aec_numeric.py','codes/chapters/ch06/core/double_talk.py',
        'codes/chapters/ch02/core/conventions.py','codes/chapters/ch00/core/audio_samples.py','codes/chapters/ch00/io_contracts.py'}
    if number==70:return common|{'scripts/make_delay_figures.py','codes/chapters/ch07/core/prediction_delays.py',
        'codes/chapters/ch07/core/delay_audio.py','codes/chapters/ch07/examples/generate_delay_audio.py',
        'codes/chapters/ch00/core/audio_samples.py','codes/chapters/ch00/io_contracts.py'}
    if number==71:return common|{'scripts/make_css_figures.py','codes/chapters/ch08/core/css.py',
        'codes/chapters/ch08/core/css_audio.py','codes/chapters/ch08/examples/generate_css_audio.py',
        'codes/chapters/ch02/core/conventions.py','codes/chapters/ch00/core/audio_samples.py','codes/chapters/ch00/io_contracts.py'}
    if number==72:return {'scripts/make_tracking_figures.py','codes/chapters/ch09/core/lifecycle.py',
        'codes/chapters/ch09/examples/tracking_lifecycle_demo.py','codes/chapters/ch09/examples/chapter09_tracking_audio.py',
        'codes/chapters/ch09/core/tracking_audio.py','codes/chapters/ch09/core/moving_source.py',
        'codes/chapters/ch09/core/tracking.py','codes/chapters/ch04/core/doa.py','codes/chapters/ch04/core/covariance.py',
        'codes/chapters/ch02/core/conventions.py','codes/chapters/ch00/core/audio_samples.py','codes/chapters/ch00/io_contracts.py'}
    if number==73:return set(quality.CHANNEL_GENERATING_SOURCES)|{'scripts/make_channel_figures.py'}
    if number==74:return {'scripts/make_selection_figures.py','codes/chapters/ch11/core/selection_physics.py',
        'codes/chapters/ch11/core/selection.py','codes/chapters/ch03/core/geometry.py',
        'codes/chapters/ch05/core/beamforming.py','codes/chapters/ch04/core/covariance.py',
        'codes/chapters/ch02/core/conventions.py','codes/chapters/ch00/io_contracts.py'}
    if number in range(76,80):return set(quality.WEARABLE_AUDIO_SOURCES)|common|{
        'codes/chapters/ch14/core/spatial_parametric.py','codes/chapters/ch14/core/binaural_matching.py',
        'scripts/make_wearable_figures.py'}
    if number in range(80,84):return set(quality.TSE_AUDIO_SOURCES)|common|{
        'codes/chapters/ch15/core/conditional_tse.py','codes/chapters/ch15/core/diarization_tsvad.py',
        'codes/chapters/ch15/core/av_alignment.py','scripts/make_tse_figures.py'}
    if number in range(84,88):return common|{'codes/chapters/ch00/core/chapter00_math.py',
        'codes/chapters/ch00/io_contracts.py','scripts/make_math_figures.py'}
    return None


def _identity(contract,report):
    _finite_tree(report)
    number=contract.number
    exercise_fields = {
        45:('exercise','E02-14'),68:('exercise_id','E05-15'),69:('exercise_id','E06-28'),
        70:('exercise_ids',['E07-18','E07-19']),71:('exercise_ids',['E08-24']),
        72:('exercise_id','E09-19'),73:('exercise_id','E10-26'),74:('exercise_id','E11-20'),
        76:('exercise_ids',['E14-02','E14-04','E14-05']),77:('exercise_ids',['E14-06','E14-07','E14-08']),
        78:('exercise_ids',['E14-09','E14-10','E14-13']),79:('exercise_ids',['E14-11','E14-12','E14-13']),
        80:('exercise_ids',['E15-01','E15-02']),81:('exercise_ids',['E15-03','E15-04','E15-05','E15-12']),
        82:('exercise_ids',['E15-06','E15-07','E15-08','E15-09','E15-10']),83:('exercise_ids',['E15-11']),
        84:('exercise_ids',['E00-03','E00-04','E00-05']),85:('exercise_ids',['E00-06','E00-07','E00-08']),
        86:('exercise_ids',['E00-09','E00-10']),87:('exercise_ids',['E00-11','E00-13','E00-14'])}
    if number in exercise_fields:
        field,expected = exercise_fields[number]
        _close(report[field],expected,field)
    if 'schema_version'in report:_close(report['schema_version'],1,'schema_version')
    if 'figure'in report:_close(report['figure'],number,'figure')
    producer_sha=hashlib.sha256(validate_parent_chain(ROOT/contract.producer_relative).read_bytes()).hexdigest()
    if 'source_script' in report:
        _close(report['source_script'],contract.producer_relative,'source_script')
    if 'script_sha256' in report:
        _close(report['script_sha256'],producer_sha,'script_sha256')
    if number in (13,20):
        _close(report['generator'],contract.producer_relative+'::'+
            ('gcc_reverb_statistics' if number==13 else 'fig_wpe'),'generator')
    roles=_source_roles(number)
    if roles is not None:
        expected={p:hashlib.sha256(validate_parent_chain(ROOT/p).read_bytes()).hexdigest() for p in roles}
        _close(report['source_sha256'],expected,'source_sha256')
        if report['source_sha256'].get(contract.producer_relative)!=producer_sha:
            raise ValueError('independent producer role missing')
    else:
        field={13:'generator_sha256',20:'source_sha256',45:'source_script_sha256'}.get(number,'script_sha256')
        _close(report[field],producer_sha,field)
        if number==20:
            expected={p:hashlib.sha256(validate_parent_chain(ROOT/p).read_bytes()).hexdigest() for p in
                ('scripts/make_figures.py','codes/chapters/ch07/core/dereverberation.py','codes/chapters/ch02/core/conventions.py')}
            _close(report['generator_inputs'],expected,'generator_inputs')
    folders={59:'ch11/scenario_audio',60:'appendix_a/weighted_audio',61:'appendix_b/response_audio',
        65:'ch13/distributed_audio',67:'ch13/distributed_audio',68:'ch05/phase_audio',69:'ch06/reference_audio',
        70:'ch07/delay_audio',71:'ch08/css_audio',**{i:'ch14/wearable_audio' for i in range(76,80)},
        **{i:'ch15/tse_audio' for i in range(80,84)}}
    if number in folders:
        folder=ROOT/'codes/chapters'/folders[number]
        _close(report['audio_manifest_sha256'],hashlib.sha256(validate_parent_chain(folder/'MANIFEST.json').read_bytes()).hexdigest(),'audio manifest bytes')
        if 'wav_sha256'in report:
            names={68:quality.PHASE_AUDIO_WAVS,69:quality.REFERENCE_AUDIO_WAVS,70:quality.DELAY_AUDIO_WAVS,
                71:quality.CSS_AUDIO_WAVS,78:set(quality.WEARABLE_AUDIO_CHANNELS),
                79:{'wearable_clean_ref.wav','wearable_hearthrough_comb.wav'},81:set(quality.TSE_AUDIO_CHANNELS)}[number]
            _close(report['wav_sha256'],{p:hashlib.sha256(validate_parent_chain(folder/p).read_bytes()).hexdigest() for p in names},'actual WAV bytes')


def _real_storage_shape(value):
    """Validate a rectangular JSON numeric array before NumPy can coerce it."""
    if type(value) is list:
        shapes = [_real_storage_shape(item) for item in value]
        if shapes and any(shape != shapes[0] for shape in shapes):
            raise ValueError('complex storage must be rectangular')
        return (len(value),) + (shapes[0] if shapes else ())
    if type(value) not in (int,float) or (type(value) is float and not math.isfinite(value)):
        raise ValueError('complex storage requires finite non-Boolean numbers')
    return ()


def _z(value, *, shape=None):
    if type(value)is not dict or set(value)!={'real','imag'}:raise ValueError('complex fields differ')
    real_shape = _real_storage_shape(value['real'])
    if real_shape != _real_storage_shape(value['imag']) or (shape is not None and real_shape != shape):
        raise ValueError('complex real/imag storage shapes differ from the required role')
    try:
        real = np.asarray(value['real'],dtype=float)
        imag = np.asarray(value['imag'],dtype=float)
    except (ValueError,TypeError,OverflowError) as exc:
        raise ValueError('complex storage exceeds the supported numeric representation') from exc
    if not np.isfinite(real).all() or not np.isfinite(imag).all():
        raise ValueError('complex storage must be finite')
    return real+1j*imag


def _distributed_complex_shapes(report):
    """Check the fixed two-node model's dimensions, including startup receivers."""
    def group(values,shapes,label):
        if type(values) is not list or len(values) != len(shapes):
            raise ValueError(label+': complex member count differs')
        for value,shape in zip(values,shapes):
            _z(value,shape=shape)
    for name,run in report['results']['trajectories'].items():
        if name not in ('round_robin','simultaneous','simultaneous_half'):
            raise ValueError('unknown distributed update schedule')
        group(run['compressions'],[(2,)]*2,'final compressions')
        group(run['cached_outputs'],[(4,)]*2,'final outputs')
        group(run['receiver_coefficients_raw_coordinates'],[(3,)]*2,'final receiver coefficients')
        for epoch,row in enumerate(run['history']):
            for field in ('compressions_before','compressions_after'):
                group(row[field],[(2,)]*2,field)
            for field in ('cached_outputs','outputs_at_last_solve'):
                group(row[field],[(4,)]*2,field)
            # In round-robin epoch 1, node 2 has not yet solved with a peer.
            dimensions = [3,2] if name=='round_robin' and epoch==0 else [3,3]
            group(row['receiver_coefficients_raw_coordinates'],[(n,) for n in dimensions],'receiver coefficients')
            group(row['current_receiver_projections'],[(n,4) for n in dimensions],'receiver projections')
            solutions = row['solutions']
            if type(solutions) is not list or len(solutions) != (1 if name=='round_robin' else 2):
                raise ValueError('distributed solution count differs')
            for solution in solutions:
                group(solution['incoming_compressions'],[(2,)]*2,'incoming compressions')
                for field,shape in {'proposal_local_weights':(2,),'new_compression':(2,),
                        'weights':(4,),'effective_weights_after_broadcast':(4,),
                        'compressed_weights':(3,),'raw_compressed_weights':(3,),
                        'normalized_projection':(3,4),'raw_projection':(3,4)}.items():
                    _z(solution[field],shape=shape)


def _complex_close(actual,expected,label):
    expected=np.asarray(expected)
    _close(actual,{'real':expected.real.tolist(),'imag':expected.imag.tolist()},label)


def _imaging_supplement(number,report):
    data=report['results']
    if number in (62,63):
        a=np.array([[1,1],[1,cmath.exp(-2j*math.pi/3)]],complex)
        for name,b,qinv in (('independent',[17/16,.5],[1.,.25]),('coherent',[21/16,.75],[1.2,.45])):
            case=data['two_cell']['cases'][name];history=[[0.,0.]]
            for _ in range(100):
                first=max(b[0]-.25*history[-1][1],0.)
                history.append([first,max(b[1]-.25*first,0.)])
            _close(case['gs'],{'q':history[-1],'history':history,
                'scan_residual_history':[[x+.25*y-b[0],.25*x+y-b[1]] for x,y in history],
                'iterations':100,'passes_per_iteration':1,'sweep':'forward','relaxation':1.},'two-cell full GS')
            measured=a@np.diag([1.,.25])@a.conj().T if name=='independent' else np.outer(a@np.array([1.,.5]),(a@np.array([1.,.5])).conj())
            fitted=a@np.diag(qinv)@a.conj().T;delta=fitted-measured
            numerator=float(np.sum(abs(delta)**2));denominator=float(np.sum(abs(measured)**2))
            upper=float(abs(delta[0,0])**2+abs(delta[1,1])**2+abs(delta[0,1])**2)
            # Tiny arithmetic residuals are bounded, never treated as positive
            # mathematical underflow for these order-one fixed matrices.
            _close(case['full_csm_residual'],{'frobenius_squared':numerator,
                'reference_frobenius_squared':denominator,'relative_frobenius':math.sqrt(numerator/denominator),
                'upper_unweighted_squared':upper,'squared_underflow':{'frobenius_squared':False,
                'reference_frobenius_squared':False,'upper_unweighted_squared':False}},'full CSM residual')
    if number==63:
        _close(data['rank_one_clean_sc']['relative_peak_stop_factor'],float(64*np.finfo(float).eps))
    if number!=62:return
    spatial=data['spherical_scan'];angles=2*np.pi*np.arange(8)/8
    microphones=np.column_stack([.2*np.cos(angles),.2*np.sin(angles),np.zeros(8)])
    axis=np.linspace(-.3,.3,21);xx,yy=np.meshgrid(axis,axis,indexing='xy')
    grid=np.column_stack([xx.ravel(),yy.ravel(),np.full(441,.6)])
    distance=np.linalg.norm(grid[None,:,:]-microphones[:,None,:],axis=2)
    reference=np.linalg.norm(grid,axis=1)
    transfer=reference/distance*np.exp(-2j*np.pi*4000/343*(distance-reference))
    weights=transfer/np.sum(abs(transfer)**2,axis=0)
    covariance=transfer[:,[215,225]]@np.diag([1.,.25])@transfer[:,[215,225]].conj().T
    scan=np.real(np.sum(weights.conj()*(covariance@weights),axis=0))
    directions=grid/reference[:,None];plane=np.exp(2j*np.pi*4000/343*(microphones@directions.T))/8
    plane_scan=np.real(np.sum(plane.conj()*(covariance@plane),axis=0))
    for source,index in enumerate((215,225)):
        local=np.flatnonzero(np.linalg.norm(grid[:,:2]-grid[index,:2],axis=1)<=.09+1e-12)
        peak=int(local[np.argmax(scan[local])]);plane_peak=int(local[np.argmax(plane_scan[local])])
        _close(spatial['local_peaks'][source],{'source_index':source,'true_grid_index':index,
            'peak_grid_index':peak,'true_position_m':grid[index].tolist(),'peak_position_m':grid[peak].tolist(),
            'peak_power':float(scan[peak]),'position_error_m':float(np.linalg.norm(grid[peak]-grid[index])),
            'plane_peak_position_m':grid[plane_peak].tolist(),'plane_peak_power':float(plane_scan[plane_peak])},'complete local peak')
    changed=spatial['reference_change'];new_distance=np.linalg.norm(grid[[215,225]]-[0,0,.2],axis=1)
    changed_q=np.array([1.,.25])*(reference[[215,225]]/new_distance)**2
    _close(changed['reference_position_m'],[0.,0.,.2]);_close(changed['q'],changed_q.tolist())
    _complex_close(changed['R'],covariance,'reference-invariant physical covariance')
    _close(changed['max_csm_difference'],0.,'reference covariance roundoff',atol=1e-13)


def _distributed_components(weight,reference,noise):
    a=np.array([1.,.5,2.,-.5]);w=np.asarray(weight);error=np.vdot(w,a)-a[reference]
    distortion=float(abs(error)**2);power=float(np.vdot(w,noise@w).real)
    return {'target_distortion':distortion,'noise_power':power,'total_mse':distortion+power,
        'reference_target_power':float(a[reference]**2),'normalized_mse':(distortion+power)/float(a[reference]**2),
        'target_noise_cross_term':0.,'cross_term_scope':'uncorrelated covariance model'}


def _distributed_supplement(number,report):
    data=report['results'];noise=np.eye(4);noise[0,2]=noise[2,0]=.2;noise[0,3]=noise[3,0]=.8
    if number==66:
        a=np.array([1.,.5,2.,-.5]);rx=np.outer(a,a)+noise
        for label,ref in (('central_node1',0),('central_node2',2)):
            w=np.linalg.solve(rx,a*a[ref])
            _close(data[label],_distributed_components(w,ref,noise),label)
        for run in data['trajectories'].values():
            for row in run['history']:
                for node,component in enumerate(row['cached_components']):
                    _close(component,_distributed_components(_z(row['cached_outputs'][node]),(0,2)[node],noise),'cached full component')
                for expected_node,solution in zip(row['active_nodes'],row['solutions']):
                    _close(solution['node'],expected_node);ref=(0,2)[expected_node]
                    _close(solution['reference'],ref)
                    _close(solution['missing_peers'],[]);_close(solution['receiver_peers'],[1-expected_node])
                    weight=_z(solution['weights']);raw=_z(solution['raw_projection'])
                    norms=np.sqrt(np.sum(abs(raw)**2,axis=1));unit=raw/norms[:,None]
                    _complex_close(solution['normalized_projection'],unit,'normalized broadcast coordinates')
                    _complex_close(solution['compressed_weights'],_z(solution['raw_compressed_weights'])*norms,'normalized receiver weights')
                    _complex_close(solution['proposal_local_weights'],weight[[0,1] if expected_node==0 else [2,3]],'proposal local coordinates')
                    _close(solution['condition_number'],float(np.linalg.cond(unit@rx@unit.conj().T)),'compressed condition')
                    _close(solution['components'],_distributed_components(weight,ref,noise),'proposal complete cost')
        return
    for role,row in data['rows'].items():
        if number==65:
            _close(row['analytic']['target_noise_cross_term'],0.)
            parts=[row['float']]
        else:
            c=row['float_components'];parts=list(c.values()) if role.startswith('packet_') else [c]
        for component in parts:
            # These are retained measurements, with signed terms whose sum is
            # independently checked. They are not arbitrary missing values.
            terms=[component[k] for k in ('target_distortion_power','noise_power','other_error_power',
                'target_noise_cross_power','target_other_cross_power','noise_other_cross_power')]
            _close(component['component_sum'],math.fsum(terms),'signed component sum')
            _close(component['decomposition_residual'],component['total_mse']-component['component_sum'],'decomposition residual')
            _close(component['samples'],800 if number==67 and role.startswith('packet_') and component is row['float_components']['packet'] else 28800)
            _close(component['reference_power'],.0068)
            if number==65:
                model=row['analytic']
                _close(component['target_distortion_power'],model['target_distortion'])
                _close(component['noise_power'],model['noise_power'])
                for key in ('other_error_power','target_noise_cross_power','target_other_cross_power','noise_other_cross_power'):
                    _close(component[key],0.,key,atol=2e-14)
    if number==67:
        _clock_report67(data)


def _clock_report67(data):
    """Finite declared-tone interpolation algebra, without the clock-state API."""
    phase=np.random.Generator(np.random.PCG64(1502026)).uniform(-np.pi,np.pi,6)
    positions=np.arange(1600,30400,dtype=float)
    def components(indices,rate=16000.):
        time=indices/rate;n=time*16000
        env=np.clip(np.minimum(n/1600,(31999-n)/1600),0,1)
        target=env*(.1*np.cos(2*np.pi*700*time+phase[0])+.06*np.cos(2*np.pi*1300*time+phase[1]))
        tones=env*math.sqrt(.0136)*np.cos(2*np.pi*np.array([1900,2300,2900,3500])[:,None]*time+phase[2:,None])
        return target,tones
    target,tones=components(positions)
    peer=4.25*target+2*tones[2]-.5*tones[3]
    raw_target,raw_tones=components(positions,16001.6)
    raw_peer=4.25*raw_target+2*raw_tones[2]-.5*raw_tones[3]
    point=positions*1.0001;lower=np.floor(point);fraction=point-lower
    tl,nl=components(lower,16001.6);tu,nu=components(lower+1,16001.6)
    corrected_target=(1-fraction)*tl+fraction*tu
    corrected_tones=(1-fraction)*nl+fraction*nu
    corrected_peer=4.25*corrected_target+2*corrected_tones[2]-.5*corrected_tones[3]
    energy=float(np.sum(peer**2))
    clock=data['clock']
    _close({k:v for k,v in clock.items() if k not in ('limitations','same_index_time_definition')},
        {'device_input_samples':32004,'reference_output_samples':32000,'invalid_reference_indices':[],
        'whole_vs_chunked_exact_equal':True,'uncorrected_peer_reconstruction_nmse':float(np.sum((raw_peer-peer)**2))/energy,
        'corrected_peer_reconstruction_nmse':float(np.sum((corrected_peer-peer)**2))/energy,
        'peer_reconstruction_reference_energy':energy,'same_index_drift_samples_last':31999*.0001/1.0001},
        'clock finite interpolation',rtol=3e-10,atol=3e-12)
    models={'compressed_white':(target,tones),'clock_misaligned_white':(raw_target,raw_tones),
            'clock_linear_corrected_white':(corrected_target,corrected_tones)}
    for role,(remote_target,remote_tones) in models.items():
        yt=(5/26)*target+(2/13)*4.25*remote_target
        yn=(2*tones[0]+tones[1])/13+(2/13)*(2*remote_tones[2]-.5*remote_tones[3])
        delta=yt-target
        expected={'samples':28800,'reference_power':float(np.mean(target**2)),
            'target_distortion_power':float(np.mean(delta**2)),'noise_power':float(np.mean(yn**2)),
            'other_error_power':0.,'target_noise_cross_power':float(2*np.mean(delta*yn)),
            'target_other_cross_power':0.,'noise_other_cross_power':0.,
            'total_mse':float(np.mean((delta+yn)**2)),
            'normalized_mse':float(np.mean((delta+yn)**2)/np.mean(target**2)),
            'component_sum':float(np.mean(delta**2)+np.mean(yn**2)+2*np.mean(delta*yn)),
            'decomposition_residual':0.}
        _close(data['rows'][role]['float_components'],expected,role+' finite components',rtol=3e-10,atol=3e-12)
    for role,distortion,noise in (('packet_zerofill_white',441/676,20/676),
                                ('packet_local_fallback_white',16/81,20/81)):
        for label,packet_fraction in (('packet',1.),('steady',800/28800)):
            d=.0068*((1-packet_fraction)*4/169+packet_fraction*distortion)
            n=.0068*((1-packet_fraction)*22/169+packet_fraction*noise)
            expected={'samples':800 if label=='packet' else 28800,'reference_power':.0068,
                'target_distortion_power':d,'noise_power':n,'other_error_power':0.,
                'target_noise_cross_power':0.,'target_other_cross_power':0.,'noise_other_cross_power':0.,
                'total_mse':d+n,'normalized_mse':(d+n)/.0068,'component_sum':d+n,'decomposition_residual':0.}
            _close(data['rows'][role]['float_components'][label],expected,role+'/'+label+' finite orthogonality',atol=2e-14)
    _close(_numeric_view(data['parameters']),{'sample_rate_hz':16000,'samples_per_channel':32000,
        'body_samples':32000,'complete_convolution_tail_samples':0,'seed':1502026,
        'target_frequencies_hz':[700,1300],'target_peak_amplitudes':[.1,.06],
        'noise_basis_frequencies_hz':[1900,2300,2900,3500],'noise_basis_peak_amplitude':math.sqrt(.0136),
        'target_steady_mean_square':.0068000000000000005,'a':[1.,.5,2.,-.5],
        'nodes_zero_based':[[0,1],[2,3]],'reference_microphones_zero_based':[0,2],
        'reference_target_factors':[1.,2.],
        'scoring_windows_samples':{'steady':[1600,30400],'packet':[16000,16800]},
        'steady_scoring_samples':28800,'packet_scoring_samples':800,'common_export_gain':1.,
        'clock':{'true_ppm':100.,'true_device_rate_hz':16001.6,'initial_offset_s':0.,'chunk_samples':257,
            'device_samples_for_right_endpoint_support':32004,'rate_is_known':True,'estimator_run':False},
        'packet':{'missing_interval_samples':[16000,16800],'known_flag':True,'packet_samples':160,
            'consecutive_missing_packets':5,'reference_timeline_preserved':True},
        'transport':{'remote_vector':[2.,-.5],'receiver_coefficient':2/13,
            'output_error_half_step_bound':(2/13)*(.5/32768)}},'figure67 numeric parameters')


def _report72(report):
    folder=ROOT/'codes/chapters/ch09/tracking_audio';manifest=_read(folder/'MANIFEST.json')
    quality._check_tracking_frame_arithmetic(manifest,folder/'array_noisy.wav')
    _close(report['input_sha256'],{name:hashlib.sha256(validate_parent_chain(folder/name).read_bytes()).hexdigest()
        for name in ('source.wav','array_noisy.wav','MANIFEST.json')},'lifecycle inputs')
    analysis=manifest['pcm_analysis'];frames=analysis['frames'];valid=[not 82<=i<=105 for i in range(197)]
    _close(frames['observation_valid'],valid)
    _close(report['analysis_scores'],analysis['scores'],'verified actual PCM frame scores')
    phases=[1 if i<2 else 2 if i<100 else 0 if i<106 else 3 if i<108 else 4 for i in range(197)]
    last=[320*i+511 if i<82 or i>105 else 26431 for i in range(197)]
    _close(report['plot_data'],{'state_time_s':[(320*i+511)/32000 for i in range(197)],
        'available_time_s':[(320*i+1024)/32000 for i in range(197)],
        'observation_angle_deg':frames['observation_angle_deg'],'filtered_angle_deg':frames['filtered_angle_deg'],
        'valid':valid,'last_valid_measurement_ticks':last,'state_age_s':[(320*i+511-last[i])/32000 for i in range(197)],
        'available_age_s':[(320*i+1024-last[i])/32000 for i in range(197)],'lifecycle_phase_code':phases},'full plotted clocks')
    for name,expiry,publications,clock in (('lifecycle',100,187,'available'),('state_age_comparison',102,189,'state')):
        section=report[name]
        parameters={'sample_rate_hz':16000,'window_samples':512,'hop_samples':160,'tick_denominator_hz':32000,
            'confirmation_frames':3,'maximum_age_ticks':6400,'maximum_age_s':.2,'age_clock':clock,
            'expiry_comparison':'strictly greater; equality retained','measurement_epoch':'receiver window centre',
            'event_publication_epoch':'complete receiver window availability',
            'expiry_order':'age old candidate/track before applying the incoming observation'}
        _close(section['parameters'],parameters)
        _close(section['counts'],{'frames':197,'allocated_local_ids':2,'confirmed_publications':publications,'expired_events':1})
        events=[(0,'candidate',1),(2,'confirmed',1),(expiry,'expired',1),(106,'candidate',2),(108,'confirmed',2)]
        expected_events=[]
        for i,event,identity in events:
            row={'frame':i,'event':event,'track_id':identity,'state_tick':320*i+511,'available_tick':320*i+1024,
                'state_time_s':(320*i+511)/32000,'publication_time_s':(320*i+1024)/32000}
            if event=='expired':row.update(age_ticks=6593 if clock=='available' else 6720,
                                           age_s=(6593 if clock=='available' else 6720)/32000)
            expected_events.append(row)
        _close(section['events'],expected_events,'full lifecycle event chronology')
        rows=[]
        for i in range(197):
            confirmed=2<=i<expiry or i>=108;absent=expiry<=i<106
            tick=None if absent else 320*i+511 if valid[i] else 26431
            age=None if tick is None else 320*i+(1024 if clock=='available' else 511)-tick
            rows.append({'frame':i,'observation_valid':valid[i],'state_tick':320*i+511,'available_tick':320*i+1024,
                'state_time_s':(320*i+511)/32000,'publication_time_s':(320*i+1024)/32000,
                'phase':'absent' if absent else 'confirmed' if confirmed else 'tentative',
                'track_id':None if absent else 1 if i<106 else 2,
                'consecutive_valid_frames':i+1 if i<82 else 0 if i<106 else i-105,
                'last_valid_measurement_tick':tick,'last_valid_measurement_time_s':None if tick is None else tick/32000,
                'age_ticks':age,'age_s':None if age is None else age/32000,'publish_confirmed':confirmed,
                'events':[event for frame,event,_ in events if frame==i]})
        _close(section['rows'],rows,'all197 lifecycle rows')
    _close(report['parameters'],report['lifecycle']['parameters'])


def _report73(report):
    folder=quality.CHANNEL_AUDIO_ROOT;manifest=_read(folder/'MANIFEST.json')
    pcm=quality._check_channel_pcm(folder,manifest)
    quality._check_channel_figure_report(report,manifest,pcm,folder)
    _channel_supplement(report,pcm)


_FILE_ORACLES={56:quality._check_tracking_information_report,57:quality._check_engineering_limits_report,
    59:quality._check_selection_scenarios_report,60:quality._check_weighted_noise_report,61:quality._check_equal_drr_report}


def validate_numeric_report(number,report,*,path=None):
    """Validate a parsed report; existing file oracles also inspect its original path."""
    _finite_tree(report)
    _check_numeric_schema(number, report)
    if number==66:
        _distributed_complex_shapes(report)
    header = 'figure' if number in {68,69,70,71,*range(76,88)} else 'schema_version'
    _close(report[header], number if header=='figure' else 1, header)
    direct={13:_report13,20:_report20,21:_report21,45:_report45,68:_report68,69:_report69,70:_report70,
        71:_report71,72:_report72,73:_report73,74:quality._check_selection_figure_report,
        76:_report76,77:_report77,78:_report78,79:_report79,80:_report80,81:_report81,82:_report82,83:_report83}
    if number in direct:
        direct[number](report)
    elif number in (84,85,86,87):
        _math_report(number,report)
    elif number in _FILE_ORACLES:
        if path is None:raise ValueError('a report path is required for this independent file oracle')
        _FILE_ORACLES[number](path)
    elif number in (62,63,64,75):
        if path is None:raise ValueError('an imaging report path is required')
        quality._check_imaging_figure_report(path,number)
        _imaging_supplement(number,report)
    elif number in (65,66,67):
        if path is None:raise ValueError('a distributed report path is required')
        quality._check_distributed_figure_report(path,number)
        _distributed_supplement(number,report)
    else:
        raise ValueError('no independent numeric oracle for figure'+str(number))


def check_figure_reports(errors):
    """Check fixed35 identities and numbers; PNG byte binding is a separate layer."""
    from scripts.figure_report_contracts import FIGURE_REPORTS
    expected={13,20,21,45,56,57,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,
              76,77,78,79,80,81,82,83,84,85,86,87}
    if {row.number for row in FIGURE_REPORTS}!=expected or len(FIGURE_REPORTS)!=35:
        errors.append('portable numeric report registry differs from the fixed35 oracles')
        return
    declared={ROOT/row.report_relative for row in FIGURE_REPORTS}
    actual={p for p in (ROOT/'codes/chapters').glob('*/reports/figure*.json')}
    if actual!=declared:
        errors.append('portable figure report file inventory differs from the fixed registry')
    for contract in FIGURE_REPORTS:
        path=ROOT/contract.report_relative
        try:
            report=_read(path)
            _identity(contract,report)
            validate_numeric_report(contract.number,report,path=path)
        except (OSError,ValueError,KeyError,TypeError,IndexError,OverflowError,EOFError,wave.Error,struct.error,np.linalg.LinAlgError) as exc:
            errors.append(f'figure{contract.number} numeric report: {exc}')


def _channel_supplement(report,pcm):
    _close(_numeric_view(report['parameters']),{'sample_rate_hz':16000,'samples_per_channel':32000,
        'target_frequency_hz':500.,'noise_basis_frequencies_hz':[1500.,2500.,3500.],'amplitude_each':.1,
        'fade_samples_each_end':320,'common_export_gain':1.,'tail_samples':0,'scoring_interval_samples':[2400,29600],
        'scoring_samples_per_channel':27200,'target_cycles_in_score':850,'noise_basis_cycles_in_score':[2550,4250,5950],
        'mixing_matrix_B':[[1.,0.,0.],[1.,1.,0.],[0.,1.,1.]],'steering':[1.,1.,1.],
        'known_failed_channel':0,'selected_channels':[1,2]})
    models={'reference':[[1,0,0,0]],'healthy_array':[[1,1,0,0],[1,1,1,0],[1,0,1,1]],
        'faulty_array':[[0,0,0,0],[1,1,1,0],[1,0,1,1]],'healthy_output':[[1,.5,0,.5]],
        'stale_output':[[0,-.5,0,.5]],'recomputed_output':[[1,.5,1,.5]]}
    for role,coefficients in models.items():
        power=[.005*sum(v*v for v in row) for row in coefficients]
        error=[.005*((row[0]-1)**2+sum(v*v for v in row[1:])) for row in coefficients]
        _close(report['samples'][role]['float_measurements'],{'scoring_interval_samples':[2400,29600],
            'samples_per_channel':27200,'mean_square_per_channel':power,
            'reference_error_mean_square_per_channel':error,'reference_NMSE_per_channel':[v/.005 for v in error],
            'target_frequency_phasor_real_per_channel':[0.]*len(coefficients),
            'target_frequency_phasor_imag_per_channel':[-.1*row[0] for row in coefficients],
            'phasor_convention':'2/N sum x[n] exp(-j*2*pi*500*n/fs); .1 sine -> -.1j; diagnostic only'},role+' float tone model')
        record=report['files']['channel_'+role+'.wav'];channels=len(coefficients)
        _close(_numeric_view({k:v for k,v in record.items() if k!='quantization_max_abs_error'}),
            {'channels':channels,'samples_per_channel':32000,'sample_rate_hz':16000,'sample_width_bytes':2,
             'peak':max(abs(x) for channel in pcm[role] for x in channel)/32768},role+' file format/actual peak')
        quant=record['quantization_max_abs_error']
        if type(quant)is not float or not 0<=quant<=.5/32768+1e-15:
            raise ValueError('retained channel quantization measurement outside half-step bound')


# Fixed published numeric field/type catalog. Values are checked by the oracles
# above; this independent finite schema rejects unreviewed new numeric fields,
# removed fields and Boolean substitutions. Array indices share one type path;
# actual lengths and every element are checked by the corresponding oracle.
_NUMERIC_SCHEMA_TEXT = {
    13: """
        b:/configuration/additive_noise
        f:/configuration/amplitude_decay_coefficient
        f:/configuration/column_t60_s/*
        i:/configuration/direct_delay_samples
        i:/configuration/fft_samples
        i:/configuration/frame_samples
        i:/configuration/frame_start_max_inclusive
        i:/configuration/frame_start_min
        f:/configuration/phat_additive_epsilon
        f:/configuration/row_drr_db/*
        i:/configuration/sample_rate_hz
        i:/configuration/search_lag_samples/*
        i:/configuration/signal_samples
        i:/configuration/source_lowpass_hz
        b:/configuration/source_peak_normalized
        i:/configuration/success_tolerance_samples
        f:/configuration/tail_duration_factor
        i:/configuration/true_gcc_lag_samples
        f:/configuration/wilson_z
        i:/schema_version
        i:/statistics/exceptions/*/*
        i:/statistics/failures/*/*
        f:/statistics/high/*/*
        f:/statistics/low/*/*
        f:/statistics/median/*/*
        f:/statistics/q1/*/*
        f:/statistics/q3/*/*
        f:/statistics/rate/*/*
        i:/statistics/seeds/*
        i:/statistics/successes/*/*
        f:/statistics/trial_contrasts/*/*/*
        b:/statistics/trial_correct/*/*/*
        i:/statistics/trials
    """,
    20: """
        i:/measurement/active_frames
        f:/measurement/active_threshold_fraction_of_reference_peak
        i:/measurement/display_frequency_range_hz/*
        i:/measurement/frame_count
        i:/measurement/frequency_bins
        i:/measurement/frequency_range_hz/*
        f:/measurement/last_frame_start_ms
        i:/measurement/other_frames
        i:/measurement/quiet_frames
        f:/measurement/quiet_threshold_fraction_of_reference_peak
        b:/parameters/center
        f:/parameters/decay_parameter_seconds
        i:/parameters/delay_frames
        i:/parameters/direct_delay_samples
        f:/parameters/drr_db
        f:/parameters/duration_seconds
        i:/parameters/fft_size
        i:/parameters/hop_samples
        i:/parameters/iterations
        b:/parameters/padding
        f:/parameters/relative_loading
        f:/parameters/relative_power_floor
        i:/parameters/rir_length_samples
        i:/parameters/sample_rate_hz
        i:/parameters/seed
        i:/parameters/taps
        f:/results/input/active_alignment/estimate_energy
        f:/results/input/active_alignment/reference_conjugate_dot_estimate/*
        f:/results/input/active_alignment/reference_energy
        f:/results/input/active_scale_aligned_nmse_db
        f:/results/input/quiet_energy
        f:/results/input/quiet_energy_ratio_db
        f:/results/input/total_energy
        f:/results/output/active_alignment/estimate_energy
        f:/results/output/active_alignment/reference_conjugate_dot_estimate/*
        f:/results/output/active_alignment/reference_energy
        f:/results/output/active_scale_aligned_nmse_db
        f:/results/output/quiet_energy
        f:/results/output/quiet_energy_ratio_db
        f:/results/output/total_energy
        i:/schema_version
    """,
    21: """
        i:/counts/missing
        i:/counts/observed
        i:/counts/observed_replacements
        i:/counts/reflection_at_least_5_percent
        i:/counts/replacement_draws
        i:/counts/resampling
        f:/frames/ess_before_resampling/*
        f:/frames/estimate_deg/*
        i:/frames/index/*
        b:/frames/missing/*
        f:/frames/observation_deg/*
        f:/frames/reflection_fraction/*
        b:/frames/replacement_draw/*
        b:/frames/resampled/*
        f:/frames/truth_deg/*
        i:/parameters/angle_bounds_deg/*
        i:/parameters/frames
        f:/parameters/likelihood_clutter_probability
        i:/parameters/missing_frame_interval_half_open/*
        f:/parameters/observation_replacement_probability
        i:/parameters/observation_std_deg
        i:/parameters/particles
        i:/parameters/process_std_deg_per_frame
        f:/parameters/resample_ess_below
        f:/parameters/velocity_process_std_deg_per_frame
        f:/phd_illustration/grid_deg/*
        f:/phd_illustration/integral_targets
        f:/phd_illustration/intensity_targets_per_degree/*
        b:/phd_illustration/not_a_filter_run
        i:/schema_version
        f:/scores/minimum_ess
        f:/scores/pf_missing_rmse_deg
        f:/scores/pf_observed_rmse_deg
        f:/scores/raw_observed_rmse_deg
        i:/seed
    """,
    45: """
        f:/errors_against_linear/fft4/*
        f:/errors_against_linear/fft8_keep4/*
        i:/frame_starts_samples/*
        i:/hop_samples
        f:/impulse_response/*
        f:/input/*
        f:/linear_full/*
        f:/periodic_hann/*
        i:/retained_output_samples
        i:/schema_version
        f:/synthesis_denominator/*
        f:/wola_outputs/fft4/*
        f:/wola_outputs/fft8_keep4/*
        f:/wola_outputs/identity/*
    """,
    56: """
        f:/bearing_deg/*
        i:/comparison_scale
        i:/initial_position_m/*
        i:/jacobian_rows_without_positive_denominators/*/*
        f:/positions_xy_m/*/*
        f:/same_scalar_state/independent_assumption_variance
        i:/same_scalar_state/measurement_variances/*
        i:/same_scalar_state/observations/*
        f:/same_scalar_state/posterior_mean/*
        f:/same_scalar_state/posterior_variance/*
        i:/same_scalar_state/prior_mean
        i:/same_scalar_state/prior_variance
        f:/same_scalar_state/rho/*
        i:/scale_null_direction/*
        i:/schema_version
        i:/stationary_observer_xy_m/*
        f:/time_s/*
        i:/velocity_m_s/*
    """,
    57: """
        f:/finite_tail/finite_db/*
        f:/finite_tail/four_point_endpoint_extrapolation_s
        f:/finite_tail/infinite_db/*
        f:/finite_tail/infinite_t60_s
        f:/finite_tail/reverse_energy/*
        i:/finite_tail/sample_rate_hz
        fi:/finite_tail/squared_impulse/*
        f:/finite_tail/time_s/*
        i:/nonpreemptive/A_deadline_ms
        i:/nonpreemptive/A_finish_ms
        i:/nonpreemptive/A_release_ms
        i:/nonpreemptive/A_response_ms
        i:/nonpreemptive/A_start_ms
        i:/nonpreemptive/B_finish_ms
        i:/nonpreemptive/B_start_ms
        i:/nonpreemptive/blocking_supremum_plus_service_ms
        f:/nonpreemptive/utilization
        i:/schema_version
        i:/soft_update/initial_noise
        i:/soft_update/initial_probability
        i:/soft_update/power/*
        i:/soft_update/rows/*/indicator
    """,
    59: """
        f:/analytic_decision/crossing_q_single
        f:/analytic_decision/endpoint_costs/fir3/*
        f:/analytic_decision/endpoint_costs/fir9/*
        f:/analytic_decision/interval_worst/fir3
        f:/analytic_decision/interval_worst/fir9
        f:/analytic_decision/q_single_weight/*
        f:/analytic_decision/scene_worst/fir3
        f:/analytic_decision/scene_worst/fir9
        f:/analytic_nmse/dual/fir3
        f:/analytic_nmse/dual/fir9
        f:/analytic_nmse/single/fir3
        f:/analytic_nmse/single/fir9
        i:/integer_pcm/dual/fir3/integer_error_squared_sum
        i:/integer_pcm/dual/fir3/integer_reference_squared_sum
        i:/integer_pcm/dual/fir3/output_window/*
        i:/integer_pcm/dual/fir3/scored_samples
        i:/integer_pcm/dual/fir3/source_window/*
        i:/integer_pcm/dual/fir9/integer_error_squared_sum
        i:/integer_pcm/dual/fir9/integer_reference_squared_sum
        i:/integer_pcm/dual/fir9/output_window/*
        i:/integer_pcm/dual/fir9/scored_samples
        i:/integer_pcm/dual/fir9/source_window/*
        i:/integer_pcm/single/fir3/integer_error_squared_sum
        i:/integer_pcm/single/fir3/integer_reference_squared_sum
        i:/integer_pcm/single/fir3/output_window/*
        i:/integer_pcm/single/fir3/scored_samples
        i:/integer_pcm/single/fir3/source_window/*
        i:/integer_pcm/single/fir9/integer_error_squared_sum
        i:/integer_pcm/single/fir9/integer_reference_squared_sum
        i:/integer_pcm/single/fir9/output_window/*
        i:/integer_pcm/single/fir9/scored_samples
        i:/integer_pcm/single/fir9/source_window/*
        f:/pcm_decision/crossing_q_single
        f:/pcm_decision/endpoint_costs/fir3/*
        f:/pcm_decision/endpoint_costs/fir9/*
        f:/pcm_decision/interval_worst/fir3
        f:/pcm_decision/interval_worst/fir9
        f:/pcm_decision/q_single_weight/*
        f:/pcm_decision/scene_worst/fir3
        f:/pcm_decision/scene_worst/fir9
        f:/pcm_nmse/dual/fir3
        f:/pcm_nmse/dual/fir9
        f:/pcm_nmse/single/fir3
        f:/pcm_nmse/single/fir9
        i:/schema_version
    """,
    60: """
        f:/analytic_mse/gls
        f:/analytic_mse/ols
        f:/analytic_mse/reversed
        i:/integer_pcm/weighted_gls.wav/integer_error_squared_sum
        i:/integer_pcm/weighted_gls.wav/integer_reference_squared_sum
        f:/integer_pcm/weighted_gls.wav/mse
        f:/integer_pcm/weighted_gls.wav/nmse
        i:/integer_pcm/weighted_gls.wav/scored_samples
        i:/integer_pcm/weighted_ols.wav/integer_error_squared_sum
        i:/integer_pcm/weighted_ols.wav/integer_reference_squared_sum
        f:/integer_pcm/weighted_ols.wav/mse
        f:/integer_pcm/weighted_ols.wav/nmse
        i:/integer_pcm/weighted_ols.wav/scored_samples
        i:/integer_pcm/weighted_reversed.wav/integer_error_squared_sum
        i:/integer_pcm/weighted_reversed.wav/integer_reference_squared_sum
        f:/integer_pcm/weighted_reversed.wav/mse
        f:/integer_pcm/weighted_reversed.wav/nmse
        i:/integer_pcm/weighted_reversed.wav/scored_samples
        i:/schema_version
        i:/score_window/*
        f:/weights/gls/*
        f:/weights/ols/*
        f:/weights/reversed/*
    """,
    61: """
        f:/analytic_output_reflection_nmse/a
        f:/analytic_output_reflection_nmse/b
        i:/frequencies_hz/*
        i:/integer_pcm/response_full_a.wav/integer_error_squared_sum
        i:/integer_pcm/response_full_a.wav/integer_reference_squared_sum
        f:/integer_pcm/response_full_a.wav/mse
        f:/integer_pcm/response_full_a.wav/nmse
        i:/integer_pcm/response_full_a.wav/scored_samples
        i:/integer_pcm/response_full_b.wav/integer_error_squared_sum
        i:/integer_pcm/response_full_b.wav/integer_reference_squared_sum
        f:/integer_pcm/response_full_b.wav/mse
        f:/integer_pcm/response_full_b.wav/nmse
        i:/integer_pcm/response_full_b.wav/scored_samples
        f:/reflection_squared_gain/a/*
        f:/reflection_squared_gain/b/*
        f:/rir_drr_db
        i:/schema_version
        i:/score_window/*
    """,
    62: """
        f:/results/spherical_scan/A/imag/*/*
        f:/results/spherical_scan/A/real/*/*
        f:/results/spherical_scan/R/imag/*/*
        f:/results/spherical_scan/R/real/*/*
        f:/results/spherical_scan/W/imag/*/*
        f:/results/spherical_scan/W/real/*/*
        f:/results/spherical_scan/b/*
        f:/results/spherical_scan/dirty_grid_sum
        f:/results/spherical_scan/frequency_hz
        f:/results/spherical_scan/grid_axis_m/*
        f:/results/spherical_scan/grid_m/*/*
        i:/results/spherical_scan/grid_shape/*
        i:/results/spherical_scan/local_peaks/*/peak_grid_index
        f:/results/spherical_scan/local_peaks/*/peak_position_m/*
        f:/results/spherical_scan/local_peaks/*/peak_power
        f:/results/spherical_scan/local_peaks/*/plane_peak_position_m/*
        f:/results/spherical_scan/local_peaks/*/plane_peak_power
        f:/results/spherical_scan/local_peaks/*/position_error_m
        i:/results/spherical_scan/local_peaks/*/source_index
        i:/results/spherical_scan/local_peaks/*/true_grid_index
        f:/results/spherical_scan/local_peaks/*/true_position_m/*
        f:/results/spherical_scan/microphones_m/*/*
        f:/results/spherical_scan/plane_mismatch_b/*
        f:/results/spherical_scan/psf_columns/*/*
        f:/results/spherical_scan/q/*
        f:/results/spherical_scan/reference_change/R/imag/*/*
        f:/results/spherical_scan/reference_change/R/real/*/*
        f:/results/spherical_scan/reference_change/max_csm_difference
        f:/results/spherical_scan/reference_change/q/*
        f:/results/spherical_scan/reference_change/reference_position_m/*
        f:/results/spherical_scan/reference_position_m/*
        f:/results/spherical_scan/sound_speed_m_s
        i:/results/spherical_scan/source_grid_indices/*
        f:/results/spherical_scan/source_positions_m/*/*
        f:/results/spherical_scan/source_power_sum
        f:/results/two_cell/A/imag/*/*
        f:/results/two_cell/A/real/*/*
        f:/results/two_cell/P/*/*
        f:/results/two_cell/W/imag/*/*
        f:/results/two_cell/W/real/*/*
        f:/results/two_cell/cases/coherent/R/imag/*/*
        f:/results/two_cell/cases/coherent/R/real/*/*
        f:/results/two_cell/cases/coherent/b/*
        f:/results/two_cell/cases/coherent/full_csm_residual/frobenius_squared
        f:/results/two_cell/cases/coherent/full_csm_residual/reference_frobenius_squared
        f:/results/two_cell/cases/coherent/full_csm_residual/relative_frobenius
        b:/results/two_cell/cases/coherent/full_csm_residual/squared_underflow/frobenius_squared
        b:/results/two_cell/cases/coherent/full_csm_residual/squared_underflow/reference_frobenius_squared
        b:/results/two_cell/cases/coherent/full_csm_residual/squared_underflow/upper_unweighted_squared
        f:/results/two_cell/cases/coherent/full_csm_residual/upper_unweighted_squared
        f:/results/two_cell/cases/coherent/gs/history/*/*
        i:/results/two_cell/cases/coherent/gs/iterations
        i:/results/two_cell/cases/coherent/gs/passes_per_iteration
        f:/results/two_cell/cases/coherent/gs/q/*
        f:/results/two_cell/cases/coherent/gs/relaxation
        f:/results/two_cell/cases/coherent/gs/scan_residual_history/*/*
        f:/results/two_cell/cases/coherent/q_inverse/*
        f:/results/two_cell/cases/independent/R/imag/*/*
        f:/results/two_cell/cases/independent/R/real/*/*
        f:/results/two_cell/cases/independent/b/*
        f:/results/two_cell/cases/independent/full_csm_residual/frobenius_squared
        f:/results/two_cell/cases/independent/full_csm_residual/reference_frobenius_squared
        f:/results/two_cell/cases/independent/full_csm_residual/relative_frobenius
        b:/results/two_cell/cases/independent/full_csm_residual/squared_underflow/frobenius_squared
        b:/results/two_cell/cases/independent/full_csm_residual/squared_underflow/reference_frobenius_squared
        b:/results/two_cell/cases/independent/full_csm_residual/squared_underflow/upper_unweighted_squared
        f:/results/two_cell/cases/independent/full_csm_residual/upper_unweighted_squared
        f:/results/two_cell/cases/independent/gs/history/*/*
        i:/results/two_cell/cases/independent/gs/iterations
        i:/results/two_cell/cases/independent/gs/passes_per_iteration
        f:/results/two_cell/cases/independent/gs/q/*
        f:/results/two_cell/cases/independent/gs/relaxation
        f:/results/two_cell/cases/independent/gs/scan_residual_history/*/*
        f:/results/two_cell/cases/independent/q_inverse/*
        f:/results/two_cell/coherence_gamma/*
        f:/results/two_cell/coherence_inverse_q/*/*
        f:/results/two_cell/q/*
        i:/schema_version
    """,
    63: """
        f:/results/forward/history/*/*
        i:/results/forward/iterations
        i:/results/forward/passes_per_iteration
        f:/results/forward/q/*
        f:/results/forward/relaxation
        f:/results/forward/scan_residual_history/*/*
        f:/results/forward_backward/history/*/*
        i:/results/forward_backward/iterations
        i:/results/forward_backward/passes_per_iteration
        f:/results/forward_backward/q/*
        f:/results/forward_backward/relaxation
        f:/results/forward_backward/scan_residual_history/*/*
        f:/results/log_plot_floor
        f:/results/rank_one_clean_sc/clean_map/*
        i:/results/rank_one_clean_sc/completed_iterations
        f:/results/rank_one_clean_sc/damping
        f:/results/rank_one_clean_sc/relative_peak_stop_factor
        f:/results/rank_one_clean_sc/residual_csm/imag/*/*
        f:/results/rank_one_clean_sc/residual_csm/real/*/*
        f:/results/rank_one_clean_sc/steps/*/allocated_power
        f:/results/rank_one_clean_sc/steps/*/component_csm/imag/*/*
        f:/results/rank_one_clean_sc/steps/*/component_csm/real/*/*
        f:/results/rank_one_clean_sc/steps/*/h/imag/*
        f:/results/rank_one_clean_sc/steps/*/h/real/*
        i:/results/rank_one_clean_sc/steps/*/iteration
        i:/results/rank_one_clean_sc/steps/*/peak_index
        f:/results/rank_one_clean_sc/steps/*/peak_power
        f:/results/rank_one_clean_sc/steps/*/residual_scan/*
        f:/results/rank_one_cumulative/*
        f:/results/two_cell/A/imag/*/*
        f:/results/two_cell/A/real/*/*
        f:/results/two_cell/P/*/*
        f:/results/two_cell/W/imag/*/*
        f:/results/two_cell/W/real/*/*
        f:/results/two_cell/cases/coherent/R/imag/*/*
        f:/results/two_cell/cases/coherent/R/real/*/*
        f:/results/two_cell/cases/coherent/b/*
        f:/results/two_cell/cases/coherent/full_csm_residual/frobenius_squared
        f:/results/two_cell/cases/coherent/full_csm_residual/reference_frobenius_squared
        f:/results/two_cell/cases/coherent/full_csm_residual/relative_frobenius
        b:/results/two_cell/cases/coherent/full_csm_residual/squared_underflow/frobenius_squared
        b:/results/two_cell/cases/coherent/full_csm_residual/squared_underflow/reference_frobenius_squared
        b:/results/two_cell/cases/coherent/full_csm_residual/squared_underflow/upper_unweighted_squared
        f:/results/two_cell/cases/coherent/full_csm_residual/upper_unweighted_squared
        f:/results/two_cell/cases/coherent/gs/history/*/*
        i:/results/two_cell/cases/coherent/gs/iterations
        i:/results/two_cell/cases/coherent/gs/passes_per_iteration
        f:/results/two_cell/cases/coherent/gs/q/*
        f:/results/two_cell/cases/coherent/gs/relaxation
        f:/results/two_cell/cases/coherent/gs/scan_residual_history/*/*
        f:/results/two_cell/cases/coherent/q_inverse/*
        f:/results/two_cell/cases/independent/R/imag/*/*
        f:/results/two_cell/cases/independent/R/real/*/*
        f:/results/two_cell/cases/independent/b/*
        f:/results/two_cell/cases/independent/full_csm_residual/frobenius_squared
        f:/results/two_cell/cases/independent/full_csm_residual/reference_frobenius_squared
        f:/results/two_cell/cases/independent/full_csm_residual/relative_frobenius
        b:/results/two_cell/cases/independent/full_csm_residual/squared_underflow/frobenius_squared
        b:/results/two_cell/cases/independent/full_csm_residual/squared_underflow/reference_frobenius_squared
        b:/results/two_cell/cases/independent/full_csm_residual/squared_underflow/upper_unweighted_squared
        f:/results/two_cell/cases/independent/full_csm_residual/upper_unweighted_squared
        f:/results/two_cell/cases/independent/gs/history/*/*
        i:/results/two_cell/cases/independent/gs/iterations
        i:/results/two_cell/cases/independent/gs/passes_per_iteration
        f:/results/two_cell/cases/independent/gs/q/*
        f:/results/two_cell/cases/independent/gs/relaxation
        f:/results/two_cell/cases/independent/gs/scan_residual_history/*/*
        f:/results/two_cell/cases/independent/q_inverse/*
        f:/results/two_cell/coherence_gamma/*
        f:/results/two_cell/coherence_inverse_q/*/*
        f:/results/two_cell/q/*
        i:/schema_version
    """,
    64: """
        f:/results/calibrated_scan/*
        f:/results/different_csm_objective_optima/*
        f:/results/intercept_prediction
        f:/results/region_db_relative_to_one/*
        f:/results/second_channel_gain/*
        f:/results/uncalibrated_scan/*
        f:/results/unequal_amplitude_diagonal_controls/*
        i:/schema_version
    """,
    65: """
        f:/results/rows/central_correlated/analytic/noise_power
        f:/results/rows/central_correlated/analytic/normalized_mse
        f:/results/rows/central_correlated/analytic/reference_target_power
        f:/results/rows/central_correlated/analytic/target_distortion
        f:/results/rows/central_correlated/analytic/target_noise_cross_term
        f:/results/rows/central_correlated/analytic/total_mse
        f:/results/rows/central_correlated/float/component_sum
        f:/results/rows/central_correlated/float/decomposition_residual
        f:/results/rows/central_correlated/float/noise_other_cross_power
        f:/results/rows/central_correlated/float/noise_power
        f:/results/rows/central_correlated/float/normalized_mse
        f:/results/rows/central_correlated/float/other_error_power
        f:/results/rows/central_correlated/float/reference_power
        i:/results/rows/central_correlated/float/samples
        f:/results/rows/central_correlated/float/target_distortion_power
        f:/results/rows/central_correlated/float/target_noise_cross_power
        f:/results/rows/central_correlated/float/target_other_cross_power
        f:/results/rows/central_correlated/float/total_mse
        i:/results/rows/central_correlated/pcm/D
        i:/results/rows/central_correlated/pcm/E
        f:/results/rows/central_correlated/pcm/nmse
        i:/results/rows/central_correlated/pcm/samples
        f:/results/rows/central_white/analytic/noise_power
        f:/results/rows/central_white/analytic/normalized_mse
        f:/results/rows/central_white/analytic/reference_target_power
        f:/results/rows/central_white/analytic/target_distortion
        f:/results/rows/central_white/analytic/target_noise_cross_term
        f:/results/rows/central_white/analytic/total_mse
        f:/results/rows/central_white/float/component_sum
        f:/results/rows/central_white/float/decomposition_residual
        f:/results/rows/central_white/float/noise_other_cross_power
        f:/results/rows/central_white/float/noise_power
        f:/results/rows/central_white/float/normalized_mse
        f:/results/rows/central_white/float/other_error_power
        f:/results/rows/central_white/float/reference_power
        i:/results/rows/central_white/float/samples
        f:/results/rows/central_white/float/target_distortion_power
        f:/results/rows/central_white/float/target_noise_cross_power
        f:/results/rows/central_white/float/target_other_cross_power
        f:/results/rows/central_white/float/total_mse
        i:/results/rows/central_white/pcm/D
        i:/results/rows/central_white/pcm/E
        f:/results/rows/central_white/pcm/nmse
        i:/results/rows/central_white/pcm/samples
        f:/results/rows/compressed_correlated/analytic/noise_power
        f:/results/rows/compressed_correlated/analytic/normalized_mse
        f:/results/rows/compressed_correlated/analytic/reference_target_power
        f:/results/rows/compressed_correlated/analytic/target_distortion
        f:/results/rows/compressed_correlated/analytic/target_noise_cross_term
        f:/results/rows/compressed_correlated/analytic/total_mse
        f:/results/rows/compressed_correlated/float/component_sum
        f:/results/rows/compressed_correlated/float/decomposition_residual
        f:/results/rows/compressed_correlated/float/noise_other_cross_power
        f:/results/rows/compressed_correlated/float/noise_power
        f:/results/rows/compressed_correlated/float/normalized_mse
        f:/results/rows/compressed_correlated/float/other_error_power
        f:/results/rows/compressed_correlated/float/reference_power
        i:/results/rows/compressed_correlated/float/samples
        f:/results/rows/compressed_correlated/float/target_distortion_power
        f:/results/rows/compressed_correlated/float/target_noise_cross_power
        f:/results/rows/compressed_correlated/float/target_other_cross_power
        f:/results/rows/compressed_correlated/float/total_mse
        i:/results/rows/compressed_correlated/pcm/D
        i:/results/rows/compressed_correlated/pcm/E
        f:/results/rows/compressed_correlated/pcm/nmse
        i:/results/rows/compressed_correlated/pcm/samples
        f:/results/rows/compressed_white/analytic/noise_power
        f:/results/rows/compressed_white/analytic/normalized_mse
        f:/results/rows/compressed_white/analytic/reference_target_power
        f:/results/rows/compressed_white/analytic/target_distortion
        f:/results/rows/compressed_white/analytic/target_noise_cross_term
        f:/results/rows/compressed_white/analytic/total_mse
        f:/results/rows/compressed_white/float/component_sum
        f:/results/rows/compressed_white/float/decomposition_residual
        f:/results/rows/compressed_white/float/noise_other_cross_power
        f:/results/rows/compressed_white/float/noise_power
        f:/results/rows/compressed_white/float/normalized_mse
        f:/results/rows/compressed_white/float/other_error_power
        f:/results/rows/compressed_white/float/reference_power
        i:/results/rows/compressed_white/float/samples
        f:/results/rows/compressed_white/float/target_distortion_power
        f:/results/rows/compressed_white/float/target_noise_cross_power
        f:/results/rows/compressed_white/float/target_other_cross_power
        f:/results/rows/compressed_white/float/total_mse
        i:/results/rows/compressed_white/pcm/D
        i:/results/rows/compressed_white/pcm/E
        f:/results/rows/compressed_white/pcm/nmse
        i:/results/rows/compressed_white/pcm/samples
        f:/results/rows/local_node1/analytic/noise_power
        f:/results/rows/local_node1/analytic/normalized_mse
        f:/results/rows/local_node1/analytic/reference_target_power
        f:/results/rows/local_node1/analytic/target_distortion
        f:/results/rows/local_node1/analytic/target_noise_cross_term
        f:/results/rows/local_node1/analytic/total_mse
        f:/results/rows/local_node1/float/component_sum
        f:/results/rows/local_node1/float/decomposition_residual
        f:/results/rows/local_node1/float/noise_other_cross_power
        f:/results/rows/local_node1/float/noise_power
        f:/results/rows/local_node1/float/normalized_mse
        f:/results/rows/local_node1/float/other_error_power
        f:/results/rows/local_node1/float/reference_power
        i:/results/rows/local_node1/float/samples
        f:/results/rows/local_node1/float/target_distortion_power
        f:/results/rows/local_node1/float/target_noise_cross_power
        f:/results/rows/local_node1/float/target_other_cross_power
        f:/results/rows/local_node1/float/total_mse
        i:/results/rows/local_node1/pcm/D
        i:/results/rows/local_node1/pcm/E
        f:/results/rows/local_node1/pcm/nmse
        i:/results/rows/local_node1/pcm/samples
        f:/results/rows/stale_correlated/analytic/noise_power
        f:/results/rows/stale_correlated/analytic/normalized_mse
        f:/results/rows/stale_correlated/analytic/reference_target_power
        f:/results/rows/stale_correlated/analytic/target_distortion
        f:/results/rows/stale_correlated/analytic/target_noise_cross_term
        f:/results/rows/stale_correlated/analytic/total_mse
        f:/results/rows/stale_correlated/float/component_sum
        f:/results/rows/stale_correlated/float/decomposition_residual
        f:/results/rows/stale_correlated/float/noise_other_cross_power
        f:/results/rows/stale_correlated/float/noise_power
        f:/results/rows/stale_correlated/float/normalized_mse
        f:/results/rows/stale_correlated/float/other_error_power
        f:/results/rows/stale_correlated/float/reference_power
        i:/results/rows/stale_correlated/float/samples
        f:/results/rows/stale_correlated/float/target_distortion_power
        f:/results/rows/stale_correlated/float/target_noise_cross_power
        f:/results/rows/stale_correlated/float/target_other_cross_power
        f:/results/rows/stale_correlated/float/total_mse
        i:/results/rows/stale_correlated/pcm/D
        i:/results/rows/stale_correlated/pcm/E
        f:/results/rows/stale_correlated/pcm/nmse
        i:/results/rows/stale_correlated/pcm/samples
        i:/results/window/*
        i:/schema_version
    """,
    66: """
        f:/results/central_node1/noise_power
        f:/results/central_node1/normalized_mse
        f:/results/central_node1/reference_target_power
        f:/results/central_node1/target_distortion
        f:/results/central_node1/target_noise_cross_term
        f:/results/central_node1/total_mse
        f:/results/central_node2/noise_power
        f:/results/central_node2/normalized_mse
        f:/results/central_node2/reference_target_power
        f:/results/central_node2/target_distortion
        f:/results/central_node2/target_noise_cross_term
        f:/results/central_node2/total_mse
        f:/results/model/Rs/*/*
        f:/results/model/a/*
        f:/results/model/correlated/*/*
        i:/results/model/nodes/*/*
        i:/results/model/references/*
        f:/results/model/white/*/*
        f:/results/trajectories/round_robin/cached_outputs/*/imag/*
        f:/results/trajectories/round_robin/cached_outputs/*/real/*
        f:/results/trajectories/round_robin/compressions/*/imag/*
        f:/results/trajectories/round_robin/compressions/*/real/*
        i:/results/trajectories/round_robin/history/*/active_nodes/*
        f:/results/trajectories/round_robin/history/*/cached_components/*/noise_power
        f:/results/trajectories/round_robin/history/*/cached_components/*/normalized_mse
        f:/results/trajectories/round_robin/history/*/cached_components/*/reference_target_power
        f:/results/trajectories/round_robin/history/*/cached_components/*/target_distortion
        f:/results/trajectories/round_robin/history/*/cached_components/*/target_noise_cross_term
        f:/results/trajectories/round_robin/history/*/cached_components/*/total_mse
        f:/results/trajectories/round_robin/history/*/cached_outputs/*/imag/*
        f:/results/trajectories/round_robin/history/*/cached_outputs/*/real/*
        f:/results/trajectories/round_robin/history/*/compressions_after/*/imag/*
        f:/results/trajectories/round_robin/history/*/compressions_after/*/real/*
        f:/results/trajectories/round_robin/history/*/compressions_before/*/imag/*
        f:/results/trajectories/round_robin/history/*/compressions_before/*/real/*
        f:/results/trajectories/round_robin/history/*/current_receiver_projections/*/imag/*/*
        f:/results/trajectories/round_robin/history/*/current_receiver_projections/*/real/*/*
        i:/results/trajectories/round_robin/history/*/epoch
        f:/results/trajectories/round_robin/history/*/normal_equation_relative_residuals/*
        i:/results/trajectories/round_robin/history/*/output_age_in_update_solves/*
        i:/results/trajectories/round_robin/history/*/output_last_update/*
        f:/results/trajectories/round_robin/history/*/outputs_at_last_solve/*/imag/*
        f:/results/trajectories/round_robin/history/*/outputs_at_last_solve/*/real/*
        f:/results/trajectories/round_robin/history/*/receiver_coefficients_raw_coordinates/*/imag/*
        f:/results/trajectories/round_robin/history/*/receiver_coefficients_raw_coordinates/*/real/*
        i:/results/trajectories/round_robin/history/*/receiver_peer_maps/*/*
        f:/results/trajectories/round_robin/history/*/solutions/*/components/noise_power
        f:/results/trajectories/round_robin/history/*/solutions/*/components/normalized_mse
        f:/results/trajectories/round_robin/history/*/solutions/*/components/reference_target_power
        f:/results/trajectories/round_robin/history/*/solutions/*/components/target_distortion
        f:/results/trajectories/round_robin/history/*/solutions/*/components/target_noise_cross_term
        f:/results/trajectories/round_robin/history/*/solutions/*/components/total_mse
        f:/results/trajectories/round_robin/history/*/solutions/*/compressed_weights/imag/*
        f:/results/trajectories/round_robin/history/*/solutions/*/compressed_weights/real/*
        f:/results/trajectories/round_robin/history/*/solutions/*/condition_number
        f:/results/trajectories/round_robin/history/*/solutions/*/effective_weights_after_broadcast/imag/*
        f:/results/trajectories/round_robin/history/*/solutions/*/effective_weights_after_broadcast/real/*
        f:/results/trajectories/round_robin/history/*/solutions/*/incoming_compressions/*/imag/*
        f:/results/trajectories/round_robin/history/*/solutions/*/incoming_compressions/*/real/*
        f:/results/trajectories/round_robin/history/*/solutions/*/new_compression/imag/*
        f:/results/trajectories/round_robin/history/*/solutions/*/new_compression/real/*
        i:/results/trajectories/round_robin/history/*/solutions/*/node
        f:/results/trajectories/round_robin/history/*/solutions/*/normalized_projection/imag/*/*
        f:/results/trajectories/round_robin/history/*/solutions/*/normalized_projection/real/*/*
        f:/results/trajectories/round_robin/history/*/solutions/*/proposal_local_weights/imag/*
        f:/results/trajectories/round_robin/history/*/solutions/*/proposal_local_weights/real/*
        f:/results/trajectories/round_robin/history/*/solutions/*/raw_compressed_weights/imag/*
        f:/results/trajectories/round_robin/history/*/solutions/*/raw_compressed_weights/real/*
        f:/results/trajectories/round_robin/history/*/solutions/*/raw_projection/imag/*/*
        f:/results/trajectories/round_robin/history/*/solutions/*/raw_projection/real/*/*
        i:/results/trajectories/round_robin/history/*/solutions/*/receiver_peers/*
        i:/results/trajectories/round_robin/history/*/solutions/*/reference
        f:/results/trajectories/round_robin/history/*/solutions/*/weights/imag/*
        f:/results/trajectories/round_robin/history/*/solutions/*/weights/real/*
        i:/results/trajectories/round_robin/history/*/total_solve_count
        i:/results/trajectories/round_robin/history/*/update_solve_count
        i:/results/trajectories/round_robin/initial_local_solve_count
        i:/results/trajectories/round_robin/max_update_solves
        f:/results/trajectories/round_robin/receiver_coefficients_raw_coordinates/*/imag/*
        f:/results/trajectories/round_robin/receiver_coefficients_raw_coordinates/*/real/*
        i:/results/trajectories/round_robin/receiver_peer_maps/*/*
        f:/results/trajectories/round_robin/relaxation
        f:/results/trajectories/round_robin/tolerance
        i:/results/trajectories/round_robin/total_solve_count
        i:/results/trajectories/round_robin/unused_update_budget
        i:/results/trajectories/round_robin/update_solve_count
        f:/results/trajectories/simultaneous/cached_outputs/*/imag/*
        f:/results/trajectories/simultaneous/cached_outputs/*/real/*
        f:/results/trajectories/simultaneous/compressions/*/imag/*
        f:/results/trajectories/simultaneous/compressions/*/real/*
        i:/results/trajectories/simultaneous/history/*/active_nodes/*
        f:/results/trajectories/simultaneous/history/*/cached_components/*/noise_power
        f:/results/trajectories/simultaneous/history/*/cached_components/*/normalized_mse
        f:/results/trajectories/simultaneous/history/*/cached_components/*/reference_target_power
        f:/results/trajectories/simultaneous/history/*/cached_components/*/target_distortion
        f:/results/trajectories/simultaneous/history/*/cached_components/*/target_noise_cross_term
        f:/results/trajectories/simultaneous/history/*/cached_components/*/total_mse
        f:/results/trajectories/simultaneous/history/*/cached_outputs/*/imag/*
        f:/results/trajectories/simultaneous/history/*/cached_outputs/*/real/*
        f:/results/trajectories/simultaneous/history/*/compressions_after/*/imag/*
        f:/results/trajectories/simultaneous/history/*/compressions_after/*/real/*
        f:/results/trajectories/simultaneous/history/*/compressions_before/*/imag/*
        f:/results/trajectories/simultaneous/history/*/compressions_before/*/real/*
        f:/results/trajectories/simultaneous/history/*/current_receiver_projections/*/imag/*/*
        f:/results/trajectories/simultaneous/history/*/current_receiver_projections/*/real/*/*
        i:/results/trajectories/simultaneous/history/*/epoch
        f:/results/trajectories/simultaneous/history/*/normal_equation_relative_residuals/*
        i:/results/trajectories/simultaneous/history/*/output_age_in_update_solves/*
        i:/results/trajectories/simultaneous/history/*/output_last_update/*
        f:/results/trajectories/simultaneous/history/*/outputs_at_last_solve/*/imag/*
        f:/results/trajectories/simultaneous/history/*/outputs_at_last_solve/*/real/*
        f:/results/trajectories/simultaneous/history/*/receiver_coefficients_raw_coordinates/*/imag/*
        f:/results/trajectories/simultaneous/history/*/receiver_coefficients_raw_coordinates/*/real/*
        i:/results/trajectories/simultaneous/history/*/receiver_peer_maps/*/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/components/noise_power
        f:/results/trajectories/simultaneous/history/*/solutions/*/components/normalized_mse
        f:/results/trajectories/simultaneous/history/*/solutions/*/components/reference_target_power
        f:/results/trajectories/simultaneous/history/*/solutions/*/components/target_distortion
        f:/results/trajectories/simultaneous/history/*/solutions/*/components/target_noise_cross_term
        f:/results/trajectories/simultaneous/history/*/solutions/*/components/total_mse
        f:/results/trajectories/simultaneous/history/*/solutions/*/compressed_weights/imag/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/compressed_weights/real/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/condition_number
        f:/results/trajectories/simultaneous/history/*/solutions/*/effective_weights_after_broadcast/imag/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/effective_weights_after_broadcast/real/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/incoming_compressions/*/imag/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/incoming_compressions/*/real/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/new_compression/imag/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/new_compression/real/*
        i:/results/trajectories/simultaneous/history/*/solutions/*/node
        f:/results/trajectories/simultaneous/history/*/solutions/*/normalized_projection/imag/*/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/normalized_projection/real/*/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/proposal_local_weights/imag/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/proposal_local_weights/real/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/raw_compressed_weights/imag/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/raw_compressed_weights/real/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/raw_projection/imag/*/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/raw_projection/real/*/*
        i:/results/trajectories/simultaneous/history/*/solutions/*/receiver_peers/*
        i:/results/trajectories/simultaneous/history/*/solutions/*/reference
        f:/results/trajectories/simultaneous/history/*/solutions/*/weights/imag/*
        f:/results/trajectories/simultaneous/history/*/solutions/*/weights/real/*
        i:/results/trajectories/simultaneous/history/*/total_solve_count
        i:/results/trajectories/simultaneous/history/*/update_solve_count
        i:/results/trajectories/simultaneous/initial_local_solve_count
        i:/results/trajectories/simultaneous/max_update_solves
        f:/results/trajectories/simultaneous/receiver_coefficients_raw_coordinates/*/imag/*
        f:/results/trajectories/simultaneous/receiver_coefficients_raw_coordinates/*/real/*
        i:/results/trajectories/simultaneous/receiver_peer_maps/*/*
        f:/results/trajectories/simultaneous/relaxation
        f:/results/trajectories/simultaneous/tolerance
        i:/results/trajectories/simultaneous/total_solve_count
        i:/results/trajectories/simultaneous/unused_update_budget
        i:/results/trajectories/simultaneous/update_solve_count
        f:/results/trajectories/simultaneous_half/cached_outputs/*/imag/*
        f:/results/trajectories/simultaneous_half/cached_outputs/*/real/*
        f:/results/trajectories/simultaneous_half/compressions/*/imag/*
        f:/results/trajectories/simultaneous_half/compressions/*/real/*
        i:/results/trajectories/simultaneous_half/history/*/active_nodes/*
        f:/results/trajectories/simultaneous_half/history/*/cached_components/*/noise_power
        f:/results/trajectories/simultaneous_half/history/*/cached_components/*/normalized_mse
        f:/results/trajectories/simultaneous_half/history/*/cached_components/*/reference_target_power
        f:/results/trajectories/simultaneous_half/history/*/cached_components/*/target_distortion
        f:/results/trajectories/simultaneous_half/history/*/cached_components/*/target_noise_cross_term
        f:/results/trajectories/simultaneous_half/history/*/cached_components/*/total_mse
        f:/results/trajectories/simultaneous_half/history/*/cached_outputs/*/imag/*
        f:/results/trajectories/simultaneous_half/history/*/cached_outputs/*/real/*
        f:/results/trajectories/simultaneous_half/history/*/compressions_after/*/imag/*
        f:/results/trajectories/simultaneous_half/history/*/compressions_after/*/real/*
        f:/results/trajectories/simultaneous_half/history/*/compressions_before/*/imag/*
        f:/results/trajectories/simultaneous_half/history/*/compressions_before/*/real/*
        f:/results/trajectories/simultaneous_half/history/*/current_receiver_projections/*/imag/*/*
        f:/results/trajectories/simultaneous_half/history/*/current_receiver_projections/*/real/*/*
        i:/results/trajectories/simultaneous_half/history/*/epoch
        f:/results/trajectories/simultaneous_half/history/*/normal_equation_relative_residuals/*
        i:/results/trajectories/simultaneous_half/history/*/output_age_in_update_solves/*
        i:/results/trajectories/simultaneous_half/history/*/output_last_update/*
        f:/results/trajectories/simultaneous_half/history/*/outputs_at_last_solve/*/imag/*
        f:/results/trajectories/simultaneous_half/history/*/outputs_at_last_solve/*/real/*
        f:/results/trajectories/simultaneous_half/history/*/receiver_coefficients_raw_coordinates/*/imag/*
        f:/results/trajectories/simultaneous_half/history/*/receiver_coefficients_raw_coordinates/*/real/*
        i:/results/trajectories/simultaneous_half/history/*/receiver_peer_maps/*/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/components/noise_power
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/components/normalized_mse
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/components/reference_target_power
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/components/target_distortion
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/components/target_noise_cross_term
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/components/total_mse
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/compressed_weights/imag/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/compressed_weights/real/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/condition_number
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/effective_weights_after_broadcast/imag/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/effective_weights_after_broadcast/real/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/incoming_compressions/*/imag/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/incoming_compressions/*/real/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/new_compression/imag/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/new_compression/real/*
        i:/results/trajectories/simultaneous_half/history/*/solutions/*/node
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/normalized_projection/imag/*/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/normalized_projection/real/*/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/proposal_local_weights/imag/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/proposal_local_weights/real/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/raw_compressed_weights/imag/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/raw_compressed_weights/real/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/raw_projection/imag/*/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/raw_projection/real/*/*
        i:/results/trajectories/simultaneous_half/history/*/solutions/*/receiver_peers/*
        i:/results/trajectories/simultaneous_half/history/*/solutions/*/reference
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/weights/imag/*
        f:/results/trajectories/simultaneous_half/history/*/solutions/*/weights/real/*
        i:/results/trajectories/simultaneous_half/history/*/total_solve_count
        i:/results/trajectories/simultaneous_half/history/*/update_solve_count
        i:/results/trajectories/simultaneous_half/initial_local_solve_count
        i:/results/trajectories/simultaneous_half/max_update_solves
        f:/results/trajectories/simultaneous_half/receiver_coefficients_raw_coordinates/*/imag/*
        f:/results/trajectories/simultaneous_half/receiver_coefficients_raw_coordinates/*/real/*
        i:/results/trajectories/simultaneous_half/receiver_peer_maps/*/*
        f:/results/trajectories/simultaneous_half/relaxation
        f:/results/trajectories/simultaneous_half/tolerance
        i:/results/trajectories/simultaneous_half/total_solve_count
        i:/results/trajectories/simultaneous_half/unused_update_budget
        i:/results/trajectories/simultaneous_half/update_solve_count
        i:/schema_version
    """,
    67: """
        f:/results/clock/corrected_peer_reconstruction_nmse
        i:/results/clock/device_input_samples
        f:/results/clock/peer_reconstruction_reference_energy
        i:/results/clock/reference_output_samples
        f:/results/clock/same_index_drift_samples_last
        f:/results/clock/uncorrected_peer_reconstruction_nmse
        b:/results/clock/whole_vs_chunked_exact_equal
        f:/results/cross_terms/*/*
        f:/results/parameters/a/*
        i:/results/parameters/body_samples
        i:/results/parameters/clock/chunk_samples
        i:/results/parameters/clock/device_samples_for_right_endpoint_support
        b:/results/parameters/clock/estimator_run
        f:/results/parameters/clock/initial_offset_s
        b:/results/parameters/clock/rate_is_known
        f:/results/parameters/clock/true_device_rate_hz
        f:/results/parameters/clock/true_ppm
        f:/results/parameters/common_export_gain
        i:/results/parameters/complete_convolution_tail_samples
        i:/results/parameters/nodes_zero_based/*/*
        i:/results/parameters/noise_basis_frequencies_hz/*
        f:/results/parameters/noise_basis_peak_amplitude
        i:/results/parameters/packet/consecutive_missing_packets
        b:/results/parameters/packet/known_flag
        i:/results/parameters/packet/missing_interval_samples/*
        i:/results/parameters/packet/packet_samples
        b:/results/parameters/packet/reference_timeline_preserved
        i:/results/parameters/packet_scoring_samples
        i:/results/parameters/reference_microphones_zero_based/*
        f:/results/parameters/reference_target_factors/*
        i:/results/parameters/sample_rate_hz
        i:/results/parameters/samples_per_channel
        i:/results/parameters/scoring_windows_samples/packet/*
        i:/results/parameters/scoring_windows_samples/steady/*
        i:/results/parameters/seed
        i:/results/parameters/steady_scoring_samples
        i:/results/parameters/target_frequencies_hz/*
        f:/results/parameters/target_peak_amplitudes/*
        f:/results/parameters/target_steady_mean_square
        f:/results/parameters/transport/output_error_half_step_bound
        f:/results/parameters/transport/receiver_coefficient
        f:/results/parameters/transport/remote_vector/*
        f:/results/rows/clock_linear_corrected_white/float_components/component_sum
        f:/results/rows/clock_linear_corrected_white/float_components/decomposition_residual
        f:/results/rows/clock_linear_corrected_white/float_components/noise_other_cross_power
        f:/results/rows/clock_linear_corrected_white/float_components/noise_power
        f:/results/rows/clock_linear_corrected_white/float_components/normalized_mse
        f:/results/rows/clock_linear_corrected_white/float_components/other_error_power
        f:/results/rows/clock_linear_corrected_white/float_components/reference_power
        i:/results/rows/clock_linear_corrected_white/float_components/samples
        f:/results/rows/clock_linear_corrected_white/float_components/target_distortion_power
        f:/results/rows/clock_linear_corrected_white/float_components/target_noise_cross_power
        f:/results/rows/clock_linear_corrected_white/float_components/target_other_cross_power
        f:/results/rows/clock_linear_corrected_white/float_components/total_mse
        i:/results/rows/clock_linear_corrected_white/pcm/D
        i:/results/rows/clock_linear_corrected_white/pcm/E
        f:/results/rows/clock_linear_corrected_white/pcm/nmse
        i:/results/rows/clock_linear_corrected_white/pcm/samples
        f:/results/rows/clock_misaligned_white/float_components/component_sum
        f:/results/rows/clock_misaligned_white/float_components/decomposition_residual
        f:/results/rows/clock_misaligned_white/float_components/noise_other_cross_power
        f:/results/rows/clock_misaligned_white/float_components/noise_power
        f:/results/rows/clock_misaligned_white/float_components/normalized_mse
        f:/results/rows/clock_misaligned_white/float_components/other_error_power
        f:/results/rows/clock_misaligned_white/float_components/reference_power
        i:/results/rows/clock_misaligned_white/float_components/samples
        f:/results/rows/clock_misaligned_white/float_components/target_distortion_power
        f:/results/rows/clock_misaligned_white/float_components/target_noise_cross_power
        f:/results/rows/clock_misaligned_white/float_components/target_other_cross_power
        f:/results/rows/clock_misaligned_white/float_components/total_mse
        i:/results/rows/clock_misaligned_white/pcm/D
        i:/results/rows/clock_misaligned_white/pcm/E
        f:/results/rows/clock_misaligned_white/pcm/nmse
        i:/results/rows/clock_misaligned_white/pcm/samples
        f:/results/rows/compressed_white/float_components/component_sum
        f:/results/rows/compressed_white/float_components/decomposition_residual
        f:/results/rows/compressed_white/float_components/noise_other_cross_power
        f:/results/rows/compressed_white/float_components/noise_power
        f:/results/rows/compressed_white/float_components/normalized_mse
        f:/results/rows/compressed_white/float_components/other_error_power
        f:/results/rows/compressed_white/float_components/reference_power
        i:/results/rows/compressed_white/float_components/samples
        f:/results/rows/compressed_white/float_components/target_distortion_power
        f:/results/rows/compressed_white/float_components/target_noise_cross_power
        f:/results/rows/compressed_white/float_components/target_other_cross_power
        f:/results/rows/compressed_white/float_components/total_mse
        i:/results/rows/compressed_white/pcm/D
        i:/results/rows/compressed_white/pcm/E
        f:/results/rows/compressed_white/pcm/nmse
        i:/results/rows/compressed_white/pcm/samples
        f:/results/rows/packet_local_fallback_white/float_components/packet/component_sum
        f:/results/rows/packet_local_fallback_white/float_components/packet/decomposition_residual
        f:/results/rows/packet_local_fallback_white/float_components/packet/noise_other_cross_power
        f:/results/rows/packet_local_fallback_white/float_components/packet/noise_power
        f:/results/rows/packet_local_fallback_white/float_components/packet/normalized_mse
        f:/results/rows/packet_local_fallback_white/float_components/packet/other_error_power
        f:/results/rows/packet_local_fallback_white/float_components/packet/reference_power
        i:/results/rows/packet_local_fallback_white/float_components/packet/samples
        f:/results/rows/packet_local_fallback_white/float_components/packet/target_distortion_power
        f:/results/rows/packet_local_fallback_white/float_components/packet/target_noise_cross_power
        f:/results/rows/packet_local_fallback_white/float_components/packet/target_other_cross_power
        f:/results/rows/packet_local_fallback_white/float_components/packet/total_mse
        f:/results/rows/packet_local_fallback_white/float_components/steady/component_sum
        f:/results/rows/packet_local_fallback_white/float_components/steady/decomposition_residual
        f:/results/rows/packet_local_fallback_white/float_components/steady/noise_other_cross_power
        f:/results/rows/packet_local_fallback_white/float_components/steady/noise_power
        f:/results/rows/packet_local_fallback_white/float_components/steady/normalized_mse
        f:/results/rows/packet_local_fallback_white/float_components/steady/other_error_power
        f:/results/rows/packet_local_fallback_white/float_components/steady/reference_power
        i:/results/rows/packet_local_fallback_white/float_components/steady/samples
        f:/results/rows/packet_local_fallback_white/float_components/steady/target_distortion_power
        f:/results/rows/packet_local_fallback_white/float_components/steady/target_noise_cross_power
        f:/results/rows/packet_local_fallback_white/float_components/steady/target_other_cross_power
        f:/results/rows/packet_local_fallback_white/float_components/steady/total_mse
        i:/results/rows/packet_local_fallback_white/pcm/packet/D
        i:/results/rows/packet_local_fallback_white/pcm/packet/E
        f:/results/rows/packet_local_fallback_white/pcm/packet/nmse
        i:/results/rows/packet_local_fallback_white/pcm/packet/samples
        i:/results/rows/packet_local_fallback_white/pcm/steady/D
        i:/results/rows/packet_local_fallback_white/pcm/steady/E
        f:/results/rows/packet_local_fallback_white/pcm/steady/nmse
        i:/results/rows/packet_local_fallback_white/pcm/steady/samples
        f:/results/rows/packet_zerofill_white/float_components/packet/component_sum
        f:/results/rows/packet_zerofill_white/float_components/packet/decomposition_residual
        f:/results/rows/packet_zerofill_white/float_components/packet/noise_other_cross_power
        f:/results/rows/packet_zerofill_white/float_components/packet/noise_power
        f:/results/rows/packet_zerofill_white/float_components/packet/normalized_mse
        f:/results/rows/packet_zerofill_white/float_components/packet/other_error_power
        f:/results/rows/packet_zerofill_white/float_components/packet/reference_power
        i:/results/rows/packet_zerofill_white/float_components/packet/samples
        f:/results/rows/packet_zerofill_white/float_components/packet/target_distortion_power
        f:/results/rows/packet_zerofill_white/float_components/packet/target_noise_cross_power
        f:/results/rows/packet_zerofill_white/float_components/packet/target_other_cross_power
        f:/results/rows/packet_zerofill_white/float_components/packet/total_mse
        f:/results/rows/packet_zerofill_white/float_components/steady/component_sum
        f:/results/rows/packet_zerofill_white/float_components/steady/decomposition_residual
        f:/results/rows/packet_zerofill_white/float_components/steady/noise_other_cross_power
        f:/results/rows/packet_zerofill_white/float_components/steady/noise_power
        f:/results/rows/packet_zerofill_white/float_components/steady/normalized_mse
        f:/results/rows/packet_zerofill_white/float_components/steady/other_error_power
        f:/results/rows/packet_zerofill_white/float_components/steady/reference_power
        i:/results/rows/packet_zerofill_white/float_components/steady/samples
        f:/results/rows/packet_zerofill_white/float_components/steady/target_distortion_power
        f:/results/rows/packet_zerofill_white/float_components/steady/target_noise_cross_power
        f:/results/rows/packet_zerofill_white/float_components/steady/target_other_cross_power
        f:/results/rows/packet_zerofill_white/float_components/steady/total_mse
        i:/results/rows/packet_zerofill_white/pcm/packet/D
        i:/results/rows/packet_zerofill_white/pcm/packet/E
        f:/results/rows/packet_zerofill_white/pcm/packet/nmse
        i:/results/rows/packet_zerofill_white/pcm/packet/samples
        i:/results/rows/packet_zerofill_white/pcm/steady/D
        i:/results/rows/packet_zerofill_white/pcm/steady/E
        f:/results/rows/packet_zerofill_white/pcm/steady/nmse
        i:/results/rows/packet_zerofill_white/pcm/steady/samples
        i:/schema_version
    """,
    68: """
        i:/figure
        i:/parameters/algorithm_delay_samples
        f:/parameters/base_weights_real_imag/*/*
        f:/parameters/common_export_gain
        i:/parameters/fade_samples
        f:/parameters/frequencies_hz/*
        f:/parameters/high_weight_multiplier_real_imag/flip/*
        f:/parameters/high_weight_multiplier_real_imag/quadrature/*
        f:/parameters/high_weight_multiplier_real_imag/reference/*
        i:/parameters/propagation_delay_samples
        i:/parameters/sample_rate_hz
        i:/parameters/samples_in_scoring_window
        i:/parameters/samples_per_channel
        i:/parameters/scoring_cycles_per_frequency/*
        i:/parameters/scoring_interval_samples/*
        f:/parameters/source_amplitude_per_frequency
        f:/parameters/source_duration_s
        b:/samples/flip/analytic/common_pure_delay_possible
        f:/samples/flip/analytic/mean_square
        f:/samples/flip/analytic/normalized_squared_reference_error
        f:/samples/flip/analytic/reference_mean_square
        f:/samples/flip/analytic/response_real_imag/*/*
        f:/samples/flip/analytic/total_reference_mse
        i:/samples/flip/float_measurements/channels
        f:/samples/flip/float_measurements/error_squared_sum
        f:/samples/flip/float_measurements/mean_square
        f:/samples/flip/float_measurements/normalized_squared_reference_error
        f:/samples/flip/float_measurements/peak
        f:/samples/flip/float_measurements/phasor_real_imag/*/*
        i:/samples/flip/float_measurements/real_projection_columns
        f:/samples/flip/float_measurements/reference_phasor_real_imag/*/*
        f:/samples/flip/float_measurements/reference_squared_sum
        f:/samples/flip/float_measurements/response_real_imag/*/*
        i:/samples/flip/float_measurements/samples_per_channel
        i:/samples/flip/float_measurements/scoring_interval_samples/*
        f:/samples/flip/float_measurements/squared_sum
        f:/samples/flip/float_measurements/total_reference_mse
        i:/samples/flip/pcm_integer_measurements/channels
        i:/samples/flip/pcm_integer_measurements/integer_denominator_D_all_channels
        i:/samples/flip/pcm_integer_measurements/integer_denominator_D_per_channel
        i:/samples/flip/pcm_integer_measurements/integer_error_squared_sum_per_channel/*
        i:/samples/flip/pcm_integer_measurements/integer_reference_squared_sum
        i:/samples/flip/pcm_integer_measurements/integer_squared_sum_E_all_channels
        i:/samples/flip/pcm_integer_measurements/integer_squared_sum_E_per_channel/*
        f:/samples/flip/pcm_integer_measurements/mean_square_all_channels
        f:/samples/flip/pcm_integer_measurements/mean_square_per_channel/*
        f:/samples/flip/pcm_integer_measurements/normalized_squared_reference_error_per_channel/*
        i:/samples/flip/pcm_integer_measurements/samples_per_channel
        i:/samples/flip/pcm_integer_measurements/scoring_interval_samples/*
        f:/samples/flip/pcm_integer_measurements/total_reference_mse_per_channel/*
        i:/samples/flip/pcm_measurements/channels
        f:/samples/flip/pcm_measurements/error_squared_sum
        f:/samples/flip/pcm_measurements/mean_square
        f:/samples/flip/pcm_measurements/normalized_squared_reference_error
        f:/samples/flip/pcm_measurements/peak
        f:/samples/flip/pcm_measurements/phasor_real_imag/*/*
        i:/samples/flip/pcm_measurements/real_projection_columns
        f:/samples/flip/pcm_measurements/reference_phasor_real_imag/*/*
        f:/samples/flip/pcm_measurements/reference_squared_sum
        f:/samples/flip/pcm_measurements/response_real_imag/*/*
        i:/samples/flip/pcm_measurements/samples_per_channel
        i:/samples/flip/pcm_measurements/scoring_interval_samples/*
        f:/samples/flip/pcm_measurements/squared_sum
        f:/samples/flip/pcm_measurements/total_reference_mse
        b:/samples/quadrature/analytic/common_pure_delay_possible
        f:/samples/quadrature/analytic/mean_square
        f:/samples/quadrature/analytic/normalized_squared_reference_error
        f:/samples/quadrature/analytic/reference_mean_square
        f:/samples/quadrature/analytic/response_real_imag/*/*
        f:/samples/quadrature/analytic/total_reference_mse
        i:/samples/quadrature/float_measurements/channels
        f:/samples/quadrature/float_measurements/error_squared_sum
        f:/samples/quadrature/float_measurements/mean_square
        f:/samples/quadrature/float_measurements/normalized_squared_reference_error
        f:/samples/quadrature/float_measurements/peak
        f:/samples/quadrature/float_measurements/phasor_real_imag/*/*
        i:/samples/quadrature/float_measurements/real_projection_columns
        f:/samples/quadrature/float_measurements/reference_phasor_real_imag/*/*
        f:/samples/quadrature/float_measurements/reference_squared_sum
        f:/samples/quadrature/float_measurements/response_real_imag/*/*
        i:/samples/quadrature/float_measurements/samples_per_channel
        i:/samples/quadrature/float_measurements/scoring_interval_samples/*
        f:/samples/quadrature/float_measurements/squared_sum
        f:/samples/quadrature/float_measurements/total_reference_mse
        i:/samples/quadrature/pcm_integer_measurements/channels
        i:/samples/quadrature/pcm_integer_measurements/integer_denominator_D_all_channels
        i:/samples/quadrature/pcm_integer_measurements/integer_denominator_D_per_channel
        i:/samples/quadrature/pcm_integer_measurements/integer_error_squared_sum_per_channel/*
        i:/samples/quadrature/pcm_integer_measurements/integer_reference_squared_sum
        i:/samples/quadrature/pcm_integer_measurements/integer_squared_sum_E_all_channels
        i:/samples/quadrature/pcm_integer_measurements/integer_squared_sum_E_per_channel/*
        f:/samples/quadrature/pcm_integer_measurements/mean_square_all_channels
        f:/samples/quadrature/pcm_integer_measurements/mean_square_per_channel/*
        f:/samples/quadrature/pcm_integer_measurements/normalized_squared_reference_error_per_channel/*
        i:/samples/quadrature/pcm_integer_measurements/samples_per_channel
        i:/samples/quadrature/pcm_integer_measurements/scoring_interval_samples/*
        f:/samples/quadrature/pcm_integer_measurements/total_reference_mse_per_channel/*
        i:/samples/quadrature/pcm_measurements/channels
        f:/samples/quadrature/pcm_measurements/error_squared_sum
        f:/samples/quadrature/pcm_measurements/mean_square
        f:/samples/quadrature/pcm_measurements/normalized_squared_reference_error
        f:/samples/quadrature/pcm_measurements/peak
        f:/samples/quadrature/pcm_measurements/phasor_real_imag/*/*
        i:/samples/quadrature/pcm_measurements/real_projection_columns
        f:/samples/quadrature/pcm_measurements/reference_phasor_real_imag/*/*
        f:/samples/quadrature/pcm_measurements/reference_squared_sum
        f:/samples/quadrature/pcm_measurements/response_real_imag/*/*
        i:/samples/quadrature/pcm_measurements/samples_per_channel
        i:/samples/quadrature/pcm_measurements/scoring_interval_samples/*
        f:/samples/quadrature/pcm_measurements/squared_sum
        f:/samples/quadrature/pcm_measurements/total_reference_mse
        b:/samples/reference/analytic/common_pure_delay_possible
        f:/samples/reference/analytic/mean_square
        f:/samples/reference/analytic/normalized_squared_reference_error
        f:/samples/reference/analytic/reference_mean_square
        f:/samples/reference/analytic/response_real_imag/*/*
        f:/samples/reference/analytic/total_reference_mse
        i:/samples/reference/float_measurements/channels
        f:/samples/reference/float_measurements/error_squared_sum
        f:/samples/reference/float_measurements/mean_square
        f:/samples/reference/float_measurements/normalized_squared_reference_error
        f:/samples/reference/float_measurements/peak
        f:/samples/reference/float_measurements/phasor_real_imag/*/*
        i:/samples/reference/float_measurements/real_projection_columns
        f:/samples/reference/float_measurements/reference_phasor_real_imag/*/*
        f:/samples/reference/float_measurements/reference_squared_sum
        f:/samples/reference/float_measurements/response_real_imag/*/*
        i:/samples/reference/float_measurements/samples_per_channel
        i:/samples/reference/float_measurements/scoring_interval_samples/*
        f:/samples/reference/float_measurements/squared_sum
        f:/samples/reference/float_measurements/total_reference_mse
        i:/samples/reference/pcm_integer_measurements/channels
        i:/samples/reference/pcm_integer_measurements/integer_denominator_D_all_channels
        i:/samples/reference/pcm_integer_measurements/integer_denominator_D_per_channel
        i:/samples/reference/pcm_integer_measurements/integer_error_squared_sum_per_channel/*
        i:/samples/reference/pcm_integer_measurements/integer_reference_squared_sum
        i:/samples/reference/pcm_integer_measurements/integer_squared_sum_E_all_channels
        i:/samples/reference/pcm_integer_measurements/integer_squared_sum_E_per_channel/*
        f:/samples/reference/pcm_integer_measurements/mean_square_all_channels
        f:/samples/reference/pcm_integer_measurements/mean_square_per_channel/*
        f:/samples/reference/pcm_integer_measurements/normalized_squared_reference_error_per_channel/*
        i:/samples/reference/pcm_integer_measurements/samples_per_channel
        i:/samples/reference/pcm_integer_measurements/scoring_interval_samples/*
        f:/samples/reference/pcm_integer_measurements/total_reference_mse_per_channel/*
        i:/samples/reference/pcm_measurements/channels
        f:/samples/reference/pcm_measurements/error_squared_sum
        f:/samples/reference/pcm_measurements/mean_square
        f:/samples/reference/pcm_measurements/normalized_squared_reference_error
        f:/samples/reference/pcm_measurements/peak
        f:/samples/reference/pcm_measurements/phasor_real_imag/*/*
        i:/samples/reference/pcm_measurements/real_projection_columns
        f:/samples/reference/pcm_measurements/reference_phasor_real_imag/*/*
        f:/samples/reference/pcm_measurements/reference_squared_sum
        f:/samples/reference/pcm_measurements/response_real_imag/*/*
        i:/samples/reference/pcm_measurements/samples_per_channel
        i:/samples/reference/pcm_measurements/scoring_interval_samples/*
        f:/samples/reference/pcm_measurements/squared_sum
        f:/samples/reference/pcm_measurements/total_reference_mse
    """,
    69: """
        i:/figure
        f:/parameters/common_export_gain
        i:/parameters/convolution_tail_samples
        i:/parameters/fade_samples
        i:/parameters/filter_length
        f:/parameters/frequencies_hz/*
        f:/parameters/ncc_activity_rms
        f:/parameters/ncc_coherence_threshold
        i:/parameters/ncc_frame_samples
        f:/parameters/nlms_epsilon
        f:/parameters/nlms_step_size
        fi:/parameters/playback_gain_intervals/*/*
        i:/parameters/sample_rate_hz
        i:/parameters/samples_per_channel
        i:/parameters/scoring_cycles_per_frequency/stable/*
        i:/parameters/scoring_cycles_per_frequency/step/*
        i:/parameters/scoring_intervals_samples/stable/*
        i:/parameters/scoring_intervals_samples/step/*
        i:/parameters/scoring_intervals_samples/tail/*
        f:/parameters/source_amplitudes/*
        i:/parameters/source_samples
        f:/parameters/true_path_coefficients/*
        i:/parameters/true_path_delays_samples/*
        f:/pcm_tail_activity/activity_rms
        i:/pcm_tail_activity/available_after_samples/*
        f:/pcm_tail_activity/coherence_threshold
        i:/pcm_tail_activity/frame_samples
        f:/pcm_tail_activity/tail_absolute_ncc/*
        i:/pcm_tail_activity/tail_frames
        i:/pcm_tail_activity/tail_interval_samples/*
        b:/pcm_tail_activity/tail_ncc_defined/*
        b:/pcm_tail_activity/tail_near_end_truth_present/*
        i:/pcm_tail_activity/tail_states/*
        f:/plotted_tail_centered_rms/*
        f:/samples/early/analytic/source_constant_envelope_mean_square
        f:/samples/early/analytic/stable/mean_square
        f:/samples/early/analytic/stable/signed_source_multiplier
        f:/samples/early/analytic/step/mean_square
        f:/samples/early/analytic/step/signed_source_multiplier
        b:/samples/early/analytic/tail/ideal_zero
        i:/samples/early/float_measurements/stable/channels
        f:/samples/early/float_measurements/stable/mean_square
        f:/samples/early/float_measurements/stable/peak
        i:/samples/early/float_measurements/stable/samples_per_channel
        i:/samples/early/float_measurements/stable/scoring_interval_samples/*
        i:/samples/early/float_measurements/step/channels
        f:/samples/early/float_measurements/step/mean_square
        f:/samples/early/float_measurements/step/peak
        i:/samples/early/float_measurements/step/samples_per_channel
        i:/samples/early/float_measurements/step/scoring_interval_samples/*
        i:/samples/early/float_measurements/tail/channels
        f:/samples/early/float_measurements/tail/mean_square
        f:/samples/early/float_measurements/tail/peak
        i:/samples/early/float_measurements/tail/samples_per_channel
        i:/samples/early/float_measurements/tail/scoring_interval_samples/*
        i:/samples/early/pcm_integer_measurements/stable/channels
        i:/samples/early/pcm_integer_measurements/stable/integer_denominator_D_all_channels
        i:/samples/early/pcm_integer_measurements/stable/integer_denominator_D_per_channel
        i:/samples/early/pcm_integer_measurements/stable/integer_squared_sum_E_all_channels
        i:/samples/early/pcm_integer_measurements/stable/integer_squared_sum_E_per_channel/*
        f:/samples/early/pcm_integer_measurements/stable/mean_square_all_channels
        f:/samples/early/pcm_integer_measurements/stable/mean_square_per_channel/*
        i:/samples/early/pcm_integer_measurements/stable/samples_per_channel
        i:/samples/early/pcm_integer_measurements/stable/scoring_interval_samples/*
        i:/samples/early/pcm_integer_measurements/step/channels
        i:/samples/early/pcm_integer_measurements/step/integer_denominator_D_all_channels
        i:/samples/early/pcm_integer_measurements/step/integer_denominator_D_per_channel
        i:/samples/early/pcm_integer_measurements/step/integer_squared_sum_E_all_channels
        i:/samples/early/pcm_integer_measurements/step/integer_squared_sum_E_per_channel/*
        f:/samples/early/pcm_integer_measurements/step/mean_square_all_channels
        f:/samples/early/pcm_integer_measurements/step/mean_square_per_channel/*
        i:/samples/early/pcm_integer_measurements/step/samples_per_channel
        i:/samples/early/pcm_integer_measurements/step/scoring_interval_samples/*
        i:/samples/early/pcm_integer_measurements/tail/channels
        i:/samples/early/pcm_integer_measurements/tail/integer_denominator_D_all_channels
        i:/samples/early/pcm_integer_measurements/tail/integer_denominator_D_per_channel
        i:/samples/early/pcm_integer_measurements/tail/integer_squared_sum_E_all_channels
        i:/samples/early/pcm_integer_measurements/tail/integer_squared_sum_E_per_channel/*
        f:/samples/early/pcm_integer_measurements/tail/mean_square_all_channels
        f:/samples/early/pcm_integer_measurements/tail/mean_square_per_channel/*
        i:/samples/early/pcm_integer_measurements/tail/samples_per_channel
        i:/samples/early/pcm_integer_measurements/tail/scoring_interval_samples/*
        i:/samples/early/pcm_measurements/stable/channels
        f:/samples/early/pcm_measurements/stable/mean_square
        f:/samples/early/pcm_measurements/stable/peak
        i:/samples/early/pcm_measurements/stable/samples_per_channel
        i:/samples/early/pcm_measurements/stable/scoring_interval_samples/*
        i:/samples/early/pcm_measurements/step/channels
        f:/samples/early/pcm_measurements/step/mean_square
        f:/samples/early/pcm_measurements/step/peak
        i:/samples/early/pcm_measurements/step/samples_per_channel
        i:/samples/early/pcm_measurements/step/scoring_interval_samples/*
        i:/samples/early/pcm_measurements/tail/channels
        f:/samples/early/pcm_measurements/tail/mean_square
        f:/samples/early/pcm_measurements/tail/peak
        i:/samples/early/pcm_measurements/tail/samples_per_channel
        i:/samples/early/pcm_measurements/tail/scoring_interval_samples/*
        f:/samples/early_residual/analytic/source_constant_envelope_mean_square
        f:/samples/early_residual/analytic/stable/mean_square
        f:/samples/early_residual/analytic/stable/signed_source_multiplier
        f:/samples/early_residual/analytic/step/mean_square
        f:/samples/early_residual/analytic/step/signed_source_multiplier
        b:/samples/early_residual/analytic/tail/ideal_zero
        i:/samples/early_residual/float_measurements/stable/channels
        f:/samples/early_residual/float_measurements/stable/mean_square
        f:/samples/early_residual/float_measurements/stable/peak
        i:/samples/early_residual/float_measurements/stable/samples_per_channel
        i:/samples/early_residual/float_measurements/stable/scoring_interval_samples/*
        i:/samples/early_residual/float_measurements/step/channels
        f:/samples/early_residual/float_measurements/step/mean_square
        f:/samples/early_residual/float_measurements/step/peak
        i:/samples/early_residual/float_measurements/step/samples_per_channel
        i:/samples/early_residual/float_measurements/step/scoring_interval_samples/*
        i:/samples/early_residual/float_measurements/tail/channels
        f:/samples/early_residual/float_measurements/tail/mean_square
        f:/samples/early_residual/float_measurements/tail/peak
        i:/samples/early_residual/float_measurements/tail/samples_per_channel
        i:/samples/early_residual/float_measurements/tail/scoring_interval_samples/*
        i:/samples/early_residual/pcm_integer_measurements/stable/channels
        i:/samples/early_residual/pcm_integer_measurements/stable/integer_denominator_D_all_channels
        i:/samples/early_residual/pcm_integer_measurements/stable/integer_denominator_D_per_channel
        i:/samples/early_residual/pcm_integer_measurements/stable/integer_squared_sum_E_all_channels
        i:/samples/early_residual/pcm_integer_measurements/stable/integer_squared_sum_E_per_channel/*
        f:/samples/early_residual/pcm_integer_measurements/stable/mean_square_all_channels
        f:/samples/early_residual/pcm_integer_measurements/stable/mean_square_per_channel/*
        i:/samples/early_residual/pcm_integer_measurements/stable/samples_per_channel
        i:/samples/early_residual/pcm_integer_measurements/stable/scoring_interval_samples/*
        i:/samples/early_residual/pcm_integer_measurements/step/channels
        i:/samples/early_residual/pcm_integer_measurements/step/integer_denominator_D_all_channels
        i:/samples/early_residual/pcm_integer_measurements/step/integer_denominator_D_per_channel
        i:/samples/early_residual/pcm_integer_measurements/step/integer_squared_sum_E_all_channels
        i:/samples/early_residual/pcm_integer_measurements/step/integer_squared_sum_E_per_channel/*
        f:/samples/early_residual/pcm_integer_measurements/step/mean_square_all_channels
        f:/samples/early_residual/pcm_integer_measurements/step/mean_square_per_channel/*
        i:/samples/early_residual/pcm_integer_measurements/step/samples_per_channel
        i:/samples/early_residual/pcm_integer_measurements/step/scoring_interval_samples/*
        i:/samples/early_residual/pcm_integer_measurements/tail/channels
        i:/samples/early_residual/pcm_integer_measurements/tail/integer_denominator_D_all_channels
        i:/samples/early_residual/pcm_integer_measurements/tail/integer_denominator_D_per_channel
        i:/samples/early_residual/pcm_integer_measurements/tail/integer_squared_sum_E_all_channels
        i:/samples/early_residual/pcm_integer_measurements/tail/integer_squared_sum_E_per_channel/*
        f:/samples/early_residual/pcm_integer_measurements/tail/mean_square_all_channels
        f:/samples/early_residual/pcm_integer_measurements/tail/mean_square_per_channel/*
        i:/samples/early_residual/pcm_integer_measurements/tail/samples_per_channel
        i:/samples/early_residual/pcm_integer_measurements/tail/scoring_interval_samples/*
        i:/samples/early_residual/pcm_measurements/stable/channels
        f:/samples/early_residual/pcm_measurements/stable/mean_square
        f:/samples/early_residual/pcm_measurements/stable/peak
        i:/samples/early_residual/pcm_measurements/stable/samples_per_channel
        i:/samples/early_residual/pcm_measurements/stable/scoring_interval_samples/*
        i:/samples/early_residual/pcm_measurements/step/channels
        f:/samples/early_residual/pcm_measurements/step/mean_square
        f:/samples/early_residual/pcm_measurements/step/peak
        i:/samples/early_residual/pcm_measurements/step/samples_per_channel
        i:/samples/early_residual/pcm_measurements/step/scoring_interval_samples/*
        i:/samples/early_residual/pcm_measurements/tail/channels
        f:/samples/early_residual/pcm_measurements/tail/mean_square
        f:/samples/early_residual/pcm_measurements/tail/peak
        i:/samples/early_residual/pcm_measurements/tail/samples_per_channel
        i:/samples/early_residual/pcm_measurements/tail/scoring_interval_samples/*
        f:/samples/echo/analytic/source_constant_envelope_mean_square
        f:/samples/echo/analytic/stable/mean_square
        f:/samples/echo/analytic/stable/signed_source_multiplier
        f:/samples/echo/analytic/step/mean_square
        f:/samples/echo/analytic/step/signed_source_multiplier
        b:/samples/echo/analytic/tail/ideal_zero
        i:/samples/echo/float_measurements/stable/channels
        f:/samples/echo/float_measurements/stable/mean_square
        f:/samples/echo/float_measurements/stable/peak
        i:/samples/echo/float_measurements/stable/samples_per_channel
        i:/samples/echo/float_measurements/stable/scoring_interval_samples/*
        i:/samples/echo/float_measurements/step/channels
        f:/samples/echo/float_measurements/step/mean_square
        f:/samples/echo/float_measurements/step/peak
        i:/samples/echo/float_measurements/step/samples_per_channel
        i:/samples/echo/float_measurements/step/scoring_interval_samples/*
        i:/samples/echo/float_measurements/tail/channels
        f:/samples/echo/float_measurements/tail/mean_square
        f:/samples/echo/float_measurements/tail/peak
        i:/samples/echo/float_measurements/tail/samples_per_channel
        i:/samples/echo/float_measurements/tail/scoring_interval_samples/*
        i:/samples/echo/pcm_integer_measurements/stable/channels
        i:/samples/echo/pcm_integer_measurements/stable/integer_denominator_D_all_channels
        i:/samples/echo/pcm_integer_measurements/stable/integer_denominator_D_per_channel
        i:/samples/echo/pcm_integer_measurements/stable/integer_squared_sum_E_all_channels
        i:/samples/echo/pcm_integer_measurements/stable/integer_squared_sum_E_per_channel/*
        f:/samples/echo/pcm_integer_measurements/stable/mean_square_all_channels
        f:/samples/echo/pcm_integer_measurements/stable/mean_square_per_channel/*
        i:/samples/echo/pcm_integer_measurements/stable/samples_per_channel
        i:/samples/echo/pcm_integer_measurements/stable/scoring_interval_samples/*
        i:/samples/echo/pcm_integer_measurements/step/channels
        i:/samples/echo/pcm_integer_measurements/step/integer_denominator_D_all_channels
        i:/samples/echo/pcm_integer_measurements/step/integer_denominator_D_per_channel
        i:/samples/echo/pcm_integer_measurements/step/integer_squared_sum_E_all_channels
        i:/samples/echo/pcm_integer_measurements/step/integer_squared_sum_E_per_channel/*
        f:/samples/echo/pcm_integer_measurements/step/mean_square_all_channels
        f:/samples/echo/pcm_integer_measurements/step/mean_square_per_channel/*
        i:/samples/echo/pcm_integer_measurements/step/samples_per_channel
        i:/samples/echo/pcm_integer_measurements/step/scoring_interval_samples/*
        i:/samples/echo/pcm_integer_measurements/tail/channels
        i:/samples/echo/pcm_integer_measurements/tail/integer_denominator_D_all_channels
        i:/samples/echo/pcm_integer_measurements/tail/integer_denominator_D_per_channel
        i:/samples/echo/pcm_integer_measurements/tail/integer_squared_sum_E_all_channels
        i:/samples/echo/pcm_integer_measurements/tail/integer_squared_sum_E_per_channel/*
        f:/samples/echo/pcm_integer_measurements/tail/mean_square_all_channels
        f:/samples/echo/pcm_integer_measurements/tail/mean_square_per_channel/*
        i:/samples/echo/pcm_integer_measurements/tail/samples_per_channel
        i:/samples/echo/pcm_integer_measurements/tail/scoring_interval_samples/*
        i:/samples/echo/pcm_measurements/stable/channels
        f:/samples/echo/pcm_measurements/stable/mean_square
        f:/samples/echo/pcm_measurements/stable/peak
        i:/samples/echo/pcm_measurements/stable/samples_per_channel
        i:/samples/echo/pcm_measurements/stable/scoring_interval_samples/*
        i:/samples/echo/pcm_measurements/step/channels
        f:/samples/echo/pcm_measurements/step/mean_square
        f:/samples/echo/pcm_measurements/step/peak
        i:/samples/echo/pcm_measurements/step/samples_per_channel
        i:/samples/echo/pcm_measurements/step/scoring_interval_samples/*
        i:/samples/echo/pcm_measurements/tail/channels
        f:/samples/echo/pcm_measurements/tail/mean_square
        f:/samples/echo/pcm_measurements/tail/peak
        i:/samples/echo/pcm_measurements/tail/samples_per_channel
        i:/samples/echo/pcm_measurements/tail/scoring_interval_samples/*
        f:/samples/late/analytic/source_constant_envelope_mean_square
        f:/samples/late/analytic/stable/mean_square
        f:/samples/late/analytic/stable/signed_source_multiplier
        f:/samples/late/analytic/step/mean_square
        f:/samples/late/analytic/step/signed_source_multiplier
        b:/samples/late/analytic/tail/ideal_zero
        i:/samples/late/float_measurements/stable/channels
        f:/samples/late/float_measurements/stable/mean_square
        f:/samples/late/float_measurements/stable/peak
        i:/samples/late/float_measurements/stable/samples_per_channel
        i:/samples/late/float_measurements/stable/scoring_interval_samples/*
        i:/samples/late/float_measurements/step/channels
        f:/samples/late/float_measurements/step/mean_square
        f:/samples/late/float_measurements/step/peak
        i:/samples/late/float_measurements/step/samples_per_channel
        i:/samples/late/float_measurements/step/scoring_interval_samples/*
        i:/samples/late/float_measurements/tail/channels
        f:/samples/late/float_measurements/tail/mean_square
        f:/samples/late/float_measurements/tail/peak
        i:/samples/late/float_measurements/tail/samples_per_channel
        i:/samples/late/float_measurements/tail/scoring_interval_samples/*
        i:/samples/late/pcm_integer_measurements/stable/channels
        i:/samples/late/pcm_integer_measurements/stable/integer_denominator_D_all_channels
        i:/samples/late/pcm_integer_measurements/stable/integer_denominator_D_per_channel
        i:/samples/late/pcm_integer_measurements/stable/integer_squared_sum_E_all_channels
        i:/samples/late/pcm_integer_measurements/stable/integer_squared_sum_E_per_channel/*
        f:/samples/late/pcm_integer_measurements/stable/mean_square_all_channels
        f:/samples/late/pcm_integer_measurements/stable/mean_square_per_channel/*
        i:/samples/late/pcm_integer_measurements/stable/samples_per_channel
        i:/samples/late/pcm_integer_measurements/stable/scoring_interval_samples/*
        i:/samples/late/pcm_integer_measurements/step/channels
        i:/samples/late/pcm_integer_measurements/step/integer_denominator_D_all_channels
        i:/samples/late/pcm_integer_measurements/step/integer_denominator_D_per_channel
        i:/samples/late/pcm_integer_measurements/step/integer_squared_sum_E_all_channels
        i:/samples/late/pcm_integer_measurements/step/integer_squared_sum_E_per_channel/*
        f:/samples/late/pcm_integer_measurements/step/mean_square_all_channels
        f:/samples/late/pcm_integer_measurements/step/mean_square_per_channel/*
        i:/samples/late/pcm_integer_measurements/step/samples_per_channel
        i:/samples/late/pcm_integer_measurements/step/scoring_interval_samples/*
        i:/samples/late/pcm_integer_measurements/tail/channels
        i:/samples/late/pcm_integer_measurements/tail/integer_denominator_D_all_channels
        i:/samples/late/pcm_integer_measurements/tail/integer_denominator_D_per_channel
        i:/samples/late/pcm_integer_measurements/tail/integer_squared_sum_E_all_channels
        i:/samples/late/pcm_integer_measurements/tail/integer_squared_sum_E_per_channel/*
        f:/samples/late/pcm_integer_measurements/tail/mean_square_all_channels
        f:/samples/late/pcm_integer_measurements/tail/mean_square_per_channel/*
        i:/samples/late/pcm_integer_measurements/tail/samples_per_channel
        i:/samples/late/pcm_integer_measurements/tail/scoring_interval_samples/*
        i:/samples/late/pcm_measurements/stable/channels
        f:/samples/late/pcm_measurements/stable/mean_square
        f:/samples/late/pcm_measurements/stable/peak
        i:/samples/late/pcm_measurements/stable/samples_per_channel
        i:/samples/late/pcm_measurements/stable/scoring_interval_samples/*
        i:/samples/late/pcm_measurements/step/channels
        f:/samples/late/pcm_measurements/step/mean_square
        f:/samples/late/pcm_measurements/step/peak
        i:/samples/late/pcm_measurements/step/samples_per_channel
        i:/samples/late/pcm_measurements/step/scoring_interval_samples/*
        i:/samples/late/pcm_measurements/tail/channels
        f:/samples/late/pcm_measurements/tail/mean_square
        f:/samples/late/pcm_measurements/tail/peak
        i:/samples/late/pcm_measurements/tail/samples_per_channel
        i:/samples/late/pcm_measurements/tail/scoring_interval_samples/*
        f:/samples/late_residual/analytic/source_constant_envelope_mean_square
        f:/samples/late_residual/analytic/stable/mean_square
        f:/samples/late_residual/analytic/stable/signed_source_multiplier
        f:/samples/late_residual/analytic/step/mean_square
        f:/samples/late_residual/analytic/step/signed_source_multiplier
        b:/samples/late_residual/analytic/tail/ideal_zero
        i:/samples/late_residual/float_measurements/stable/channels
        f:/samples/late_residual/float_measurements/stable/mean_square
        f:/samples/late_residual/float_measurements/stable/peak
        i:/samples/late_residual/float_measurements/stable/samples_per_channel
        i:/samples/late_residual/float_measurements/stable/scoring_interval_samples/*
        i:/samples/late_residual/float_measurements/step/channels
        f:/samples/late_residual/float_measurements/step/mean_square
        f:/samples/late_residual/float_measurements/step/peak
        i:/samples/late_residual/float_measurements/step/samples_per_channel
        i:/samples/late_residual/float_measurements/step/scoring_interval_samples/*
        i:/samples/late_residual/float_measurements/tail/channels
        f:/samples/late_residual/float_measurements/tail/mean_square
        f:/samples/late_residual/float_measurements/tail/peak
        i:/samples/late_residual/float_measurements/tail/samples_per_channel
        i:/samples/late_residual/float_measurements/tail/scoring_interval_samples/*
        i:/samples/late_residual/pcm_integer_measurements/stable/channels
        i:/samples/late_residual/pcm_integer_measurements/stable/integer_denominator_D_all_channels
        i:/samples/late_residual/pcm_integer_measurements/stable/integer_denominator_D_per_channel
        i:/samples/late_residual/pcm_integer_measurements/stable/integer_squared_sum_E_all_channels
        i:/samples/late_residual/pcm_integer_measurements/stable/integer_squared_sum_E_per_channel/*
        f:/samples/late_residual/pcm_integer_measurements/stable/mean_square_all_channels
        f:/samples/late_residual/pcm_integer_measurements/stable/mean_square_per_channel/*
        i:/samples/late_residual/pcm_integer_measurements/stable/samples_per_channel
        i:/samples/late_residual/pcm_integer_measurements/stable/scoring_interval_samples/*
        i:/samples/late_residual/pcm_integer_measurements/step/channels
        i:/samples/late_residual/pcm_integer_measurements/step/integer_denominator_D_all_channels
        i:/samples/late_residual/pcm_integer_measurements/step/integer_denominator_D_per_channel
        i:/samples/late_residual/pcm_integer_measurements/step/integer_squared_sum_E_all_channels
        i:/samples/late_residual/pcm_integer_measurements/step/integer_squared_sum_E_per_channel/*
        f:/samples/late_residual/pcm_integer_measurements/step/mean_square_all_channels
        f:/samples/late_residual/pcm_integer_measurements/step/mean_square_per_channel/*
        i:/samples/late_residual/pcm_integer_measurements/step/samples_per_channel
        i:/samples/late_residual/pcm_integer_measurements/step/scoring_interval_samples/*
        i:/samples/late_residual/pcm_integer_measurements/tail/channels
        i:/samples/late_residual/pcm_integer_measurements/tail/integer_denominator_D_all_channels
        i:/samples/late_residual/pcm_integer_measurements/tail/integer_denominator_D_per_channel
        i:/samples/late_residual/pcm_integer_measurements/tail/integer_squared_sum_E_all_channels
        i:/samples/late_residual/pcm_integer_measurements/tail/integer_squared_sum_E_per_channel/*
        f:/samples/late_residual/pcm_integer_measurements/tail/mean_square_all_channels
        f:/samples/late_residual/pcm_integer_measurements/tail/mean_square_per_channel/*
        i:/samples/late_residual/pcm_integer_measurements/tail/samples_per_channel
        i:/samples/late_residual/pcm_integer_measurements/tail/scoring_interval_samples/*
        i:/samples/late_residual/pcm_measurements/stable/channels
        f:/samples/late_residual/pcm_measurements/stable/mean_square
        f:/samples/late_residual/pcm_measurements/stable/peak
        i:/samples/late_residual/pcm_measurements/stable/samples_per_channel
        i:/samples/late_residual/pcm_measurements/stable/scoring_interval_samples/*
        i:/samples/late_residual/pcm_measurements/step/channels
        f:/samples/late_residual/pcm_measurements/step/mean_square
        f:/samples/late_residual/pcm_measurements/step/peak
        i:/samples/late_residual/pcm_measurements/step/samples_per_channel
        i:/samples/late_residual/pcm_measurements/step/scoring_interval_samples/*
        i:/samples/late_residual/pcm_measurements/tail/channels
        f:/samples/late_residual/pcm_measurements/tail/mean_square
        f:/samples/late_residual/pcm_measurements/tail/peak
        i:/samples/late_residual/pcm_measurements/tail/samples_per_channel
        i:/samples/late_residual/pcm_measurements/tail/scoring_interval_samples/*
        f:/samples/wrong_gain_residual/analytic/source_constant_envelope_mean_square
        f:/samples/wrong_gain_residual/analytic/stable/mean_square
        f:/samples/wrong_gain_residual/analytic/stable/signed_source_multiplier
        f:/samples/wrong_gain_residual/analytic/step/mean_square
        f:/samples/wrong_gain_residual/analytic/step/signed_source_multiplier
        b:/samples/wrong_gain_residual/analytic/tail/ideal_zero
        i:/samples/wrong_gain_residual/float_measurements/stable/channels
        f:/samples/wrong_gain_residual/float_measurements/stable/mean_square
        f:/samples/wrong_gain_residual/float_measurements/stable/peak
        i:/samples/wrong_gain_residual/float_measurements/stable/samples_per_channel
        i:/samples/wrong_gain_residual/float_measurements/stable/scoring_interval_samples/*
        i:/samples/wrong_gain_residual/float_measurements/step/channels
        f:/samples/wrong_gain_residual/float_measurements/step/mean_square
        f:/samples/wrong_gain_residual/float_measurements/step/peak
        i:/samples/wrong_gain_residual/float_measurements/step/samples_per_channel
        i:/samples/wrong_gain_residual/float_measurements/step/scoring_interval_samples/*
        i:/samples/wrong_gain_residual/float_measurements/tail/channels
        f:/samples/wrong_gain_residual/float_measurements/tail/mean_square
        f:/samples/wrong_gain_residual/float_measurements/tail/peak
        i:/samples/wrong_gain_residual/float_measurements/tail/samples_per_channel
        i:/samples/wrong_gain_residual/float_measurements/tail/scoring_interval_samples/*
        i:/samples/wrong_gain_residual/pcm_integer_measurements/stable/channels
        i:/samples/wrong_gain_residual/pcm_integer_measurements/stable/integer_denominator_D_all_channels
        i:/samples/wrong_gain_residual/pcm_integer_measurements/stable/integer_denominator_D_per_channel
        i:/samples/wrong_gain_residual/pcm_integer_measurements/stable/integer_squared_sum_E_all_channels
        i:/samples/wrong_gain_residual/pcm_integer_measurements/stable/integer_squared_sum_E_per_channel/*
        f:/samples/wrong_gain_residual/pcm_integer_measurements/stable/mean_square_all_channels
        f:/samples/wrong_gain_residual/pcm_integer_measurements/stable/mean_square_per_channel/*
        i:/samples/wrong_gain_residual/pcm_integer_measurements/stable/samples_per_channel
        i:/samples/wrong_gain_residual/pcm_integer_measurements/stable/scoring_interval_samples/*
        i:/samples/wrong_gain_residual/pcm_integer_measurements/step/channels
        i:/samples/wrong_gain_residual/pcm_integer_measurements/step/integer_denominator_D_all_channels
        i:/samples/wrong_gain_residual/pcm_integer_measurements/step/integer_denominator_D_per_channel
        i:/samples/wrong_gain_residual/pcm_integer_measurements/step/integer_squared_sum_E_all_channels
        i:/samples/wrong_gain_residual/pcm_integer_measurements/step/integer_squared_sum_E_per_channel/*
        f:/samples/wrong_gain_residual/pcm_integer_measurements/step/mean_square_all_channels
        f:/samples/wrong_gain_residual/pcm_integer_measurements/step/mean_square_per_channel/*
        i:/samples/wrong_gain_residual/pcm_integer_measurements/step/samples_per_channel
        i:/samples/wrong_gain_residual/pcm_integer_measurements/step/scoring_interval_samples/*
        i:/samples/wrong_gain_residual/pcm_integer_measurements/tail/channels
        i:/samples/wrong_gain_residual/pcm_integer_measurements/tail/integer_denominator_D_all_channels
        i:/samples/wrong_gain_residual/pcm_integer_measurements/tail/integer_denominator_D_per_channel
        i:/samples/wrong_gain_residual/pcm_integer_measurements/tail/integer_squared_sum_E_all_channels
        i:/samples/wrong_gain_residual/pcm_integer_measurements/tail/integer_squared_sum_E_per_channel/*
        f:/samples/wrong_gain_residual/pcm_integer_measurements/tail/mean_square_all_channels
        f:/samples/wrong_gain_residual/pcm_integer_measurements/tail/mean_square_per_channel/*
        i:/samples/wrong_gain_residual/pcm_integer_measurements/tail/samples_per_channel
        i:/samples/wrong_gain_residual/pcm_integer_measurements/tail/scoring_interval_samples/*
        i:/samples/wrong_gain_residual/pcm_measurements/stable/channels
        f:/samples/wrong_gain_residual/pcm_measurements/stable/mean_square
        f:/samples/wrong_gain_residual/pcm_measurements/stable/peak
        i:/samples/wrong_gain_residual/pcm_measurements/stable/samples_per_channel
        i:/samples/wrong_gain_residual/pcm_measurements/stable/scoring_interval_samples/*
        i:/samples/wrong_gain_residual/pcm_measurements/step/channels
        f:/samples/wrong_gain_residual/pcm_measurements/step/mean_square
        f:/samples/wrong_gain_residual/pcm_measurements/step/peak
        i:/samples/wrong_gain_residual/pcm_measurements/step/samples_per_channel
        i:/samples/wrong_gain_residual/pcm_measurements/step/scoring_interval_samples/*
        i:/samples/wrong_gain_residual/pcm_measurements/tail/channels
        f:/samples/wrong_gain_residual/pcm_measurements/tail/mean_square
        f:/samples/wrong_gain_residual/pcm_measurements/tail/peak
        i:/samples/wrong_gain_residual/pcm_measurements/tail/samples_per_channel
        i:/samples/wrong_gain_residual/pcm_measurements/tail/scoring_interval_samples/*
    """,
    70: """
        i:/figure
        i:/parameters/aligned_auxiliary_history_delay_samples
        f:/parameters/amplitude
        i:/parameters/common_auxiliary_history_delay_samples
        f:/parameters/common_export_gain
        i:/parameters/complete_retained_support_samples/*
        f:/parameters/frequency_hz
        i:/parameters/hop_samples
        i:/parameters/maximum_history_padding_samples
        i:/parameters/observation_delays_samples/auxiliary
        i:/parameters/observation_delays_samples/reference
        i:/parameters/physical_propagation_tail_samples
        i:/parameters/pulse_samples
        i:/parameters/pulse_start_frames_in_period/*
        i:/parameters/reference_protection_delay_samples
        i:/parameters/sample_rate_hz
        i:/parameters/samples_per_channel
        i:/parameters/scoring_interval_samples/*
        i:/parameters/scoring_samples_per_channel
        i:/parameters/scoring_source_periods
        i:/parameters/source_period_frames
        i:/parameters/source_period_samples
        i:/parameters/source_periods
        i:/parameters/source_samples
        i:/parameters/source_time_gaps_samples/aligned
        i:/parameters/source_time_gaps_samples/common
        i:/pcm_integer_errors/aligned_residual/D
        i:/pcm_integer_errors/aligned_residual/E
        f:/pcm_integer_errors/aligned_residual/nmse
        i:/pcm_integer_errors/common_residual/D
        i:/pcm_integer_errors/common_residual/E
        f:/pcm_integer_errors/common_residual/nmse
        f:/power_decay_control/T60_seconds
        f:/power_decay_control/hop_seconds
        i:/power_decay_control/hops/*
        f:/power_decay_control/relative_power/*
        f:/power_decay_control/wrong_amplitude_factor_used_as_power/*
        i:/scoring_interval_samples/*
        i:/source_time_control/aligned_history/*
        i:/source_time_control/common_history/*
        i:/source_time_control/receiver_indices/*
        i:/source_time_control/reference/*
    """,
    71: """
        i:/figure
        f:/overlap_fit/decoded_pcm_overlap_application_gain/*
        f:/overlap_fit/float_application_gain/*
        f:/overlap_fit/matching/absolute_centered_correlation/*/*
        f:/overlap_fit/matching/candidate_mean_scores_identity_swap/*
        i:/overlap_fit/matching/current_indices_for_previous/*
        f:/overlap_fit/matching/margin
        f:/overlap_fit/matching/thresholds/minimum_centered_rms
        f:/overlap_fit/matching/thresholds/minimum_correlation
        f:/overlap_fit/matching/thresholds/minimum_margin
        f:/parameters/amplitudes/*
        f:/parameters/common_export_gain
        i:/parameters/fade_samples_each_end
        i:/parameters/first_block_samples/*
        f:/parameters/frequencies_hz/*
        i:/parameters/overlap_cycles_per_channel/*
        i:/parameters/overlap_samples/*
        i:/parameters/overlap_samples_per_channel
        i:/parameters/post_overlap_cycles_per_channel/*
        i:/parameters/sample_rate_hz
        i:/parameters/samples_per_channel
        i:/parameters/scoring_intervals_samples/post_overlap/*
        i:/parameters/scoring_intervals_samples/primary/*
        i:/parameters/scoring_samples_per_channel/post_overlap
        i:/parameters/scoring_samples_per_channel/primary
        f:/parameters/second_block_multiplier
        i:/parameters/second_block_samples/*
        i:/pcm_integer_measurements/corrected/post_overlap/channels
        i:/pcm_integer_measurements/corrected/post_overlap/integer_denominator_D_all_channels
        i:/pcm_integer_measurements/corrected/post_overlap/integer_denominator_D_per_channel
        i:/pcm_integer_measurements/corrected/post_overlap/integer_output_reference_cross_sum_per_channel/*
        i:/pcm_integer_measurements/corrected/post_overlap/integer_reference_error_squared_sum_per_channel/*
        i:/pcm_integer_measurements/corrected/post_overlap/integer_reference_squared_sum_per_channel/*
        i:/pcm_integer_measurements/corrected/post_overlap/integer_squared_sum_E_all_channels
        i:/pcm_integer_measurements/corrected/post_overlap/integer_squared_sum_E_per_channel/*
        f:/pcm_integer_measurements/corrected/post_overlap/mean_square_all_channels
        f:/pcm_integer_measurements/corrected/post_overlap/mean_square_per_channel/*
        i:/pcm_integer_measurements/corrected/post_overlap/pcm_decode_divisor
        f:/pcm_integer_measurements/corrected/post_overlap/projection_gain_per_channel/*
        f:/pcm_integer_measurements/corrected/post_overlap/reference_NMSE_per_channel/*
        f:/pcm_integer_measurements/corrected/post_overlap/reference_error_mean_square_per_channel/*
        i:/pcm_integer_measurements/corrected/post_overlap/samples_per_channel
        i:/pcm_integer_measurements/corrected/post_overlap/scoring_interval_samples/*
        i:/pcm_integer_measurements/corrected/primary/channels
        i:/pcm_integer_measurements/corrected/primary/integer_denominator_D_all_channels
        i:/pcm_integer_measurements/corrected/primary/integer_denominator_D_per_channel
        i:/pcm_integer_measurements/corrected/primary/integer_output_reference_cross_sum_per_channel/*
        i:/pcm_integer_measurements/corrected/primary/integer_reference_error_squared_sum_per_channel/*
        i:/pcm_integer_measurements/corrected/primary/integer_reference_squared_sum_per_channel/*
        i:/pcm_integer_measurements/corrected/primary/integer_squared_sum_E_all_channels
        i:/pcm_integer_measurements/corrected/primary/integer_squared_sum_E_per_channel/*
        f:/pcm_integer_measurements/corrected/primary/mean_square_all_channels
        f:/pcm_integer_measurements/corrected/primary/mean_square_per_channel/*
        i:/pcm_integer_measurements/corrected/primary/pcm_decode_divisor
        f:/pcm_integer_measurements/corrected/primary/projection_gain_per_channel/*
        f:/pcm_integer_measurements/corrected/primary/reference_NMSE_per_channel/*
        f:/pcm_integer_measurements/corrected/primary/reference_error_mean_square_per_channel/*
        i:/pcm_integer_measurements/corrected/primary/samples_per_channel
        i:/pcm_integer_measurements/corrected/primary/scoring_interval_samples/*
        i:/pcm_integer_measurements/naive/post_overlap/channels
        i:/pcm_integer_measurements/naive/post_overlap/integer_denominator_D_all_channels
        i:/pcm_integer_measurements/naive/post_overlap/integer_denominator_D_per_channel
        i:/pcm_integer_measurements/naive/post_overlap/integer_output_reference_cross_sum_per_channel/*
        i:/pcm_integer_measurements/naive/post_overlap/integer_reference_error_squared_sum_per_channel/*
        i:/pcm_integer_measurements/naive/post_overlap/integer_reference_squared_sum_per_channel/*
        i:/pcm_integer_measurements/naive/post_overlap/integer_squared_sum_E_all_channels
        i:/pcm_integer_measurements/naive/post_overlap/integer_squared_sum_E_per_channel/*
        f:/pcm_integer_measurements/naive/post_overlap/mean_square_all_channels
        f:/pcm_integer_measurements/naive/post_overlap/mean_square_per_channel/*
        i:/pcm_integer_measurements/naive/post_overlap/pcm_decode_divisor
        f:/pcm_integer_measurements/naive/post_overlap/projection_gain_per_channel/*
        f:/pcm_integer_measurements/naive/post_overlap/reference_NMSE_per_channel/*
        f:/pcm_integer_measurements/naive/post_overlap/reference_error_mean_square_per_channel/*
        i:/pcm_integer_measurements/naive/post_overlap/samples_per_channel
        i:/pcm_integer_measurements/naive/post_overlap/scoring_interval_samples/*
        i:/pcm_integer_measurements/naive/primary/channels
        i:/pcm_integer_measurements/naive/primary/integer_denominator_D_all_channels
        i:/pcm_integer_measurements/naive/primary/integer_denominator_D_per_channel
        i:/pcm_integer_measurements/naive/primary/integer_output_reference_cross_sum_per_channel/*
        i:/pcm_integer_measurements/naive/primary/integer_reference_error_squared_sum_per_channel/*
        i:/pcm_integer_measurements/naive/primary/integer_reference_squared_sum_per_channel/*
        i:/pcm_integer_measurements/naive/primary/integer_squared_sum_E_all_channels
        i:/pcm_integer_measurements/naive/primary/integer_squared_sum_E_per_channel/*
        f:/pcm_integer_measurements/naive/primary/mean_square_all_channels
        f:/pcm_integer_measurements/naive/primary/mean_square_per_channel/*
        i:/pcm_integer_measurements/naive/primary/pcm_decode_divisor
        f:/pcm_integer_measurements/naive/primary/projection_gain_per_channel/*
        f:/pcm_integer_measurements/naive/primary/reference_NMSE_per_channel/*
        f:/pcm_integer_measurements/naive/primary/reference_error_mean_square_per_channel/*
        i:/pcm_integer_measurements/naive/primary/samples_per_channel
        i:/pcm_integer_measurements/naive/primary/scoring_interval_samples/*
        i:/pcm_integer_measurements/polarity/post_overlap/channels
        i:/pcm_integer_measurements/polarity/post_overlap/integer_denominator_D_all_channels
        i:/pcm_integer_measurements/polarity/post_overlap/integer_denominator_D_per_channel
        i:/pcm_integer_measurements/polarity/post_overlap/integer_output_reference_cross_sum_per_channel/*
        i:/pcm_integer_measurements/polarity/post_overlap/integer_reference_error_squared_sum_per_channel/*
        i:/pcm_integer_measurements/polarity/post_overlap/integer_reference_squared_sum_per_channel/*
        i:/pcm_integer_measurements/polarity/post_overlap/integer_squared_sum_E_all_channels
        i:/pcm_integer_measurements/polarity/post_overlap/integer_squared_sum_E_per_channel/*
        f:/pcm_integer_measurements/polarity/post_overlap/mean_square_all_channels
        f:/pcm_integer_measurements/polarity/post_overlap/mean_square_per_channel/*
        i:/pcm_integer_measurements/polarity/post_overlap/pcm_decode_divisor
        f:/pcm_integer_measurements/polarity/post_overlap/projection_gain_per_channel/*
        f:/pcm_integer_measurements/polarity/post_overlap/reference_NMSE_per_channel/*
        f:/pcm_integer_measurements/polarity/post_overlap/reference_error_mean_square_per_channel/*
        i:/pcm_integer_measurements/polarity/post_overlap/samples_per_channel
        i:/pcm_integer_measurements/polarity/post_overlap/scoring_interval_samples/*
        i:/pcm_integer_measurements/polarity/primary/channels
        i:/pcm_integer_measurements/polarity/primary/integer_denominator_D_all_channels
        i:/pcm_integer_measurements/polarity/primary/integer_denominator_D_per_channel
        i:/pcm_integer_measurements/polarity/primary/integer_output_reference_cross_sum_per_channel/*
        i:/pcm_integer_measurements/polarity/primary/integer_reference_error_squared_sum_per_channel/*
        i:/pcm_integer_measurements/polarity/primary/integer_reference_squared_sum_per_channel/*
        i:/pcm_integer_measurements/polarity/primary/integer_squared_sum_E_all_channels
        i:/pcm_integer_measurements/polarity/primary/integer_squared_sum_E_per_channel/*
        f:/pcm_integer_measurements/polarity/primary/mean_square_all_channels
        f:/pcm_integer_measurements/polarity/primary/mean_square_per_channel/*
        i:/pcm_integer_measurements/polarity/primary/pcm_decode_divisor
        f:/pcm_integer_measurements/polarity/primary/projection_gain_per_channel/*
        f:/pcm_integer_measurements/polarity/primary/reference_NMSE_per_channel/*
        f:/pcm_integer_measurements/polarity/primary/reference_error_mean_square_per_channel/*
        i:/pcm_integer_measurements/polarity/primary/samples_per_channel
        i:/pcm_integer_measurements/polarity/primary/scoring_interval_samples/*
        i:/pcm_integer_measurements/reference/post_overlap/channels
        i:/pcm_integer_measurements/reference/post_overlap/integer_denominator_D_all_channels
        i:/pcm_integer_measurements/reference/post_overlap/integer_denominator_D_per_channel
        i:/pcm_integer_measurements/reference/post_overlap/integer_output_reference_cross_sum_per_channel/*
        i:/pcm_integer_measurements/reference/post_overlap/integer_reference_error_squared_sum_per_channel/*
        i:/pcm_integer_measurements/reference/post_overlap/integer_reference_squared_sum_per_channel/*
        i:/pcm_integer_measurements/reference/post_overlap/integer_squared_sum_E_all_channels
        i:/pcm_integer_measurements/reference/post_overlap/integer_squared_sum_E_per_channel/*
        f:/pcm_integer_measurements/reference/post_overlap/mean_square_all_channels
        f:/pcm_integer_measurements/reference/post_overlap/mean_square_per_channel/*
        i:/pcm_integer_measurements/reference/post_overlap/pcm_decode_divisor
        f:/pcm_integer_measurements/reference/post_overlap/projection_gain_per_channel/*
        f:/pcm_integer_measurements/reference/post_overlap/reference_NMSE_per_channel/*
        f:/pcm_integer_measurements/reference/post_overlap/reference_error_mean_square_per_channel/*
        i:/pcm_integer_measurements/reference/post_overlap/samples_per_channel
        i:/pcm_integer_measurements/reference/post_overlap/scoring_interval_samples/*
        i:/pcm_integer_measurements/reference/primary/channels
        i:/pcm_integer_measurements/reference/primary/integer_denominator_D_all_channels
        i:/pcm_integer_measurements/reference/primary/integer_denominator_D_per_channel
        i:/pcm_integer_measurements/reference/primary/integer_output_reference_cross_sum_per_channel/*
        i:/pcm_integer_measurements/reference/primary/integer_reference_error_squared_sum_per_channel/*
        i:/pcm_integer_measurements/reference/primary/integer_reference_squared_sum_per_channel/*
        i:/pcm_integer_measurements/reference/primary/integer_squared_sum_E_all_channels
        i:/pcm_integer_measurements/reference/primary/integer_squared_sum_E_per_channel/*
        f:/pcm_integer_measurements/reference/primary/mean_square_all_channels
        f:/pcm_integer_measurements/reference/primary/mean_square_per_channel/*
        i:/pcm_integer_measurements/reference/primary/pcm_decode_divisor
        f:/pcm_integer_measurements/reference/primary/projection_gain_per_channel/*
        f:/pcm_integer_measurements/reference/primary/reference_NMSE_per_channel/*
        f:/pcm_integer_measurements/reference/primary/reference_error_mean_square_per_channel/*
        i:/pcm_integer_measurements/reference/primary/samples_per_channel
        i:/pcm_integer_measurements/reference/primary/scoring_interval_samples/*
        i:/waveform_display_interval_samples/*
    """,
    72: """
        f:/analysis_scores/filtered_missing_rmse_deg
        f:/analysis_scores/filtered_valid_rmse_deg
        i:/analysis_scores/frame_count
        i:/analysis_scores/initialized_missing_count
        i:/analysis_scores/missing_observation_count
        f:/analysis_scores/raw_valid_rmse_deg
        i:/analysis_scores/valid_observation_count
        i:/lifecycle/counts/allocated_local_ids
        i:/lifecycle/counts/confirmed_publications
        i:/lifecycle/counts/expired_events
        i:/lifecycle/counts/frames
        f:/lifecycle/events/*/age_s
        i:/lifecycle/events/*/age_ticks
        i:/lifecycle/events/*/available_tick
        i:/lifecycle/events/*/frame
        f:/lifecycle/events/*/publication_time_s
        i:/lifecycle/events/*/state_tick
        f:/lifecycle/events/*/state_time_s
        i:/lifecycle/events/*/track_id
        i:/lifecycle/parameters/confirmation_frames
        i:/lifecycle/parameters/hop_samples
        f:/lifecycle/parameters/maximum_age_s
        i:/lifecycle/parameters/maximum_age_ticks
        i:/lifecycle/parameters/sample_rate_hz
        i:/lifecycle/parameters/tick_denominator_hz
        i:/lifecycle/parameters/window_samples
        f:/lifecycle/rows/*/age_s
        i:/lifecycle/rows/*/age_ticks
        i:/lifecycle/rows/*/available_tick
        i:/lifecycle/rows/*/consecutive_valid_frames
        i:/lifecycle/rows/*/frame
        i:/lifecycle/rows/*/last_valid_measurement_tick
        f:/lifecycle/rows/*/last_valid_measurement_time_s
        b:/lifecycle/rows/*/observation_valid
        f:/lifecycle/rows/*/publication_time_s
        b:/lifecycle/rows/*/publish_confirmed
        i:/lifecycle/rows/*/state_tick
        f:/lifecycle/rows/*/state_time_s
        i:/lifecycle/rows/*/track_id
        i:/parameters/confirmation_frames
        i:/parameters/hop_samples
        f:/parameters/maximum_age_s
        i:/parameters/maximum_age_ticks
        i:/parameters/sample_rate_hz
        i:/parameters/tick_denominator_hz
        i:/parameters/window_samples
        f:/plot_data/available_age_s/*
        f:/plot_data/available_time_s/*
        f:/plot_data/filtered_angle_deg/*
        i:/plot_data/last_valid_measurement_ticks/*
        i:/plot_data/lifecycle_phase_code/*
        f:/plot_data/observation_angle_deg/*
        f:/plot_data/state_age_s/*
        f:/plot_data/state_time_s/*
        b:/plot_data/valid/*
        i:/schema_version
        i:/state_age_comparison/counts/allocated_local_ids
        i:/state_age_comparison/counts/confirmed_publications
        i:/state_age_comparison/counts/expired_events
        i:/state_age_comparison/counts/frames
        f:/state_age_comparison/events/*/age_s
        i:/state_age_comparison/events/*/age_ticks
        i:/state_age_comparison/events/*/available_tick
        i:/state_age_comparison/events/*/frame
        f:/state_age_comparison/events/*/publication_time_s
        i:/state_age_comparison/events/*/state_tick
        f:/state_age_comparison/events/*/state_time_s
        i:/state_age_comparison/events/*/track_id
        i:/state_age_comparison/parameters/confirmation_frames
        i:/state_age_comparison/parameters/hop_samples
        f:/state_age_comparison/parameters/maximum_age_s
        i:/state_age_comparison/parameters/maximum_age_ticks
        i:/state_age_comparison/parameters/sample_rate_hz
        i:/state_age_comparison/parameters/tick_denominator_hz
        i:/state_age_comparison/parameters/window_samples
        f:/state_age_comparison/rows/*/age_s
        i:/state_age_comparison/rows/*/age_ticks
        i:/state_age_comparison/rows/*/available_tick
        i:/state_age_comparison/rows/*/consecutive_valid_frames
        i:/state_age_comparison/rows/*/frame
        i:/state_age_comparison/rows/*/last_valid_measurement_tick
        f:/state_age_comparison/rows/*/last_valid_measurement_time_s
        b:/state_age_comparison/rows/*/observation_valid
        f:/state_age_comparison/rows/*/publication_time_s
        b:/state_age_comparison/rows/*/publish_confirmed
        i:/state_age_comparison/rows/*/state_tick
        f:/state_age_comparison/rows/*/state_time_s
        i:/state_age_comparison/rows/*/track_id
    """,
    73: """
        f:/analytic/healthy_output/noise_mean_square
        f:/analytic/healthy_output/output_mean_square
        f:/analytic/healthy_output/reference_NMSE
        f:/analytic/healthy_output/reference_error_mean_square
        f:/analytic/healthy_output/target_gain
        f:/analytic/healthy_output/target_mean_square
        f:/analytic/recomputed_output/noise_mean_square
        f:/analytic/recomputed_output/output_mean_square
        f:/analytic/recomputed_output/reference_NMSE
        f:/analytic/recomputed_output/reference_error_mean_square
        f:/analytic/recomputed_output/target_gain
        f:/analytic/recomputed_output/target_mean_square
        f:/analytic/stale_output/noise_mean_square
        f:/analytic/stale_output/output_mean_square
        f:/analytic/stale_output/reference_NMSE
        f:/analytic/stale_output/reference_error_mean_square
        f:/analytic/stale_output/target_gain
        f:/analytic/stale_output/target_mean_square
        i:/files/channel_faulty_array.wav/channels
        f:/files/channel_faulty_array.wav/peak
        f:/files/channel_faulty_array.wav/quantization_max_abs_error
        i:/files/channel_faulty_array.wav/sample_rate_hz
        i:/files/channel_faulty_array.wav/sample_width_bytes
        i:/files/channel_faulty_array.wav/samples_per_channel
        i:/files/channel_healthy_array.wav/channels
        f:/files/channel_healthy_array.wav/peak
        f:/files/channel_healthy_array.wav/quantization_max_abs_error
        i:/files/channel_healthy_array.wav/sample_rate_hz
        i:/files/channel_healthy_array.wav/sample_width_bytes
        i:/files/channel_healthy_array.wav/samples_per_channel
        i:/files/channel_healthy_output.wav/channels
        f:/files/channel_healthy_output.wav/peak
        f:/files/channel_healthy_output.wav/quantization_max_abs_error
        i:/files/channel_healthy_output.wav/sample_rate_hz
        i:/files/channel_healthy_output.wav/sample_width_bytes
        i:/files/channel_healthy_output.wav/samples_per_channel
        i:/files/channel_recomputed_output.wav/channels
        f:/files/channel_recomputed_output.wav/peak
        f:/files/channel_recomputed_output.wav/quantization_max_abs_error
        i:/files/channel_recomputed_output.wav/sample_rate_hz
        i:/files/channel_recomputed_output.wav/sample_width_bytes
        i:/files/channel_recomputed_output.wav/samples_per_channel
        i:/files/channel_reference.wav/channels
        f:/files/channel_reference.wav/peak
        f:/files/channel_reference.wav/quantization_max_abs_error
        i:/files/channel_reference.wav/sample_rate_hz
        i:/files/channel_reference.wav/sample_width_bytes
        i:/files/channel_reference.wav/samples_per_channel
        i:/files/channel_stale_output.wav/channels
        f:/files/channel_stale_output.wav/peak
        f:/files/channel_stale_output.wav/quantization_max_abs_error
        i:/files/channel_stale_output.wav/sample_rate_hz
        i:/files/channel_stale_output.wav/sample_width_bytes
        i:/files/channel_stale_output.wav/samples_per_channel
        f:/model/full_weights/*
        f:/model/noise_covariance/*/*
        f:/model/selected_covariance/*/*
        f:/model/selected_steering/*
        f:/model/selected_weights/*
        f:/model/selection_matrix/*/*
        f:/model/stale_remaining_weights/*
        f:/parameters/amplitude_each
        f:/parameters/common_export_gain
        i:/parameters/fade_samples_each_end
        i:/parameters/known_failed_channel
        f:/parameters/mixing_matrix_B/*/*
        i:/parameters/noise_basis_cycles_in_score/*
        f:/parameters/noise_basis_frequencies_hz/*
        i:/parameters/sample_rate_hz
        i:/parameters/samples_per_channel
        i:/parameters/scoring_interval_samples/*
        i:/parameters/scoring_samples_per_channel
        i:/parameters/selected_channels/*
        f:/parameters/steering/*
        i:/parameters/tail_samples
        i:/parameters/target_cycles_in_score
        f:/parameters/target_frequency_hz
        f:/plot_data/effective_weights_after_prescribed_fault/*/*
        f:/plot_data/projection_gain_pcm_diagnostic/*
        f:/plot_data/reference_NMSE_analytic/*
        f:/plot_data/reference_NMSE_pcm/*
        i:/plot_data/samples_per_channel
        i:/plot_data/scoring_interval_samples/*
        f:/plot_data/target_gain_analytic/*
        i:/plot_data/waveform_interval_samples/*
        f:/plot_data/waveform_pcm/healthy_output/*
        f:/plot_data/waveform_pcm/recomputed_output/*
        f:/plot_data/waveform_pcm/reference/*
        f:/plot_data/waveform_pcm/stale_output/*
        f:/plot_data/waveform_time_ms/*
        f:/samples/faulty_array/float_measurements/mean_square_per_channel/*
        f:/samples/faulty_array/float_measurements/reference_NMSE_per_channel/*
        f:/samples/faulty_array/float_measurements/reference_error_mean_square_per_channel/*
        i:/samples/faulty_array/float_measurements/samples_per_channel
        i:/samples/faulty_array/float_measurements/scoring_interval_samples/*
        f:/samples/faulty_array/float_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/faulty_array/float_measurements/target_frequency_phasor_real_per_channel/*
        i:/samples/faulty_array/pcm_integer_measurements/channels
        i:/samples/faulty_array/pcm_integer_measurements/integer_denominator_D_all_channels
        i:/samples/faulty_array/pcm_integer_measurements/integer_denominator_D_per_channel
        i:/samples/faulty_array/pcm_integer_measurements/integer_output_reference_cross_sum_per_channel/*
        i:/samples/faulty_array/pcm_integer_measurements/integer_reference_error_squared_sum_per_channel/*
        i:/samples/faulty_array/pcm_integer_measurements/integer_reference_squared_sum_per_channel/*
        i:/samples/faulty_array/pcm_integer_measurements/integer_squared_sum_E_all_channels
        i:/samples/faulty_array/pcm_integer_measurements/integer_squared_sum_E_per_channel/*
        f:/samples/faulty_array/pcm_integer_measurements/mean_square_all_channels
        f:/samples/faulty_array/pcm_integer_measurements/mean_square_per_channel/*
        i:/samples/faulty_array/pcm_integer_measurements/pcm_decode_divisor
        f:/samples/faulty_array/pcm_integer_measurements/projection_gain_per_channel/*
        f:/samples/faulty_array/pcm_integer_measurements/reference_NMSE_per_channel/*
        f:/samples/faulty_array/pcm_integer_measurements/reference_error_mean_square_per_channel/*
        i:/samples/faulty_array/pcm_integer_measurements/samples_per_channel
        i:/samples/faulty_array/pcm_integer_measurements/scoring_interval_samples/*
        f:/samples/faulty_array/pcm_measurements/mean_square_per_channel/*
        f:/samples/faulty_array/pcm_measurements/reference_NMSE_per_channel/*
        f:/samples/faulty_array/pcm_measurements/reference_error_mean_square_per_channel/*
        i:/samples/faulty_array/pcm_measurements/samples_per_channel
        i:/samples/faulty_array/pcm_measurements/scoring_interval_samples/*
        f:/samples/faulty_array/pcm_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/faulty_array/pcm_measurements/target_frequency_phasor_real_per_channel/*
        f:/samples/healthy_array/float_measurements/mean_square_per_channel/*
        f:/samples/healthy_array/float_measurements/reference_NMSE_per_channel/*
        f:/samples/healthy_array/float_measurements/reference_error_mean_square_per_channel/*
        i:/samples/healthy_array/float_measurements/samples_per_channel
        i:/samples/healthy_array/float_measurements/scoring_interval_samples/*
        f:/samples/healthy_array/float_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/healthy_array/float_measurements/target_frequency_phasor_real_per_channel/*
        i:/samples/healthy_array/pcm_integer_measurements/channels
        i:/samples/healthy_array/pcm_integer_measurements/integer_denominator_D_all_channels
        i:/samples/healthy_array/pcm_integer_measurements/integer_denominator_D_per_channel
        i:/samples/healthy_array/pcm_integer_measurements/integer_output_reference_cross_sum_per_channel/*
        i:/samples/healthy_array/pcm_integer_measurements/integer_reference_error_squared_sum_per_channel/*
        i:/samples/healthy_array/pcm_integer_measurements/integer_reference_squared_sum_per_channel/*
        i:/samples/healthy_array/pcm_integer_measurements/integer_squared_sum_E_all_channels
        i:/samples/healthy_array/pcm_integer_measurements/integer_squared_sum_E_per_channel/*
        f:/samples/healthy_array/pcm_integer_measurements/mean_square_all_channels
        f:/samples/healthy_array/pcm_integer_measurements/mean_square_per_channel/*
        i:/samples/healthy_array/pcm_integer_measurements/pcm_decode_divisor
        f:/samples/healthy_array/pcm_integer_measurements/projection_gain_per_channel/*
        f:/samples/healthy_array/pcm_integer_measurements/reference_NMSE_per_channel/*
        f:/samples/healthy_array/pcm_integer_measurements/reference_error_mean_square_per_channel/*
        i:/samples/healthy_array/pcm_integer_measurements/samples_per_channel
        i:/samples/healthy_array/pcm_integer_measurements/scoring_interval_samples/*
        f:/samples/healthy_array/pcm_measurements/mean_square_per_channel/*
        f:/samples/healthy_array/pcm_measurements/reference_NMSE_per_channel/*
        f:/samples/healthy_array/pcm_measurements/reference_error_mean_square_per_channel/*
        i:/samples/healthy_array/pcm_measurements/samples_per_channel
        i:/samples/healthy_array/pcm_measurements/scoring_interval_samples/*
        f:/samples/healthy_array/pcm_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/healthy_array/pcm_measurements/target_frequency_phasor_real_per_channel/*
        f:/samples/healthy_output/float_measurements/mean_square_per_channel/*
        f:/samples/healthy_output/float_measurements/reference_NMSE_per_channel/*
        f:/samples/healthy_output/float_measurements/reference_error_mean_square_per_channel/*
        i:/samples/healthy_output/float_measurements/samples_per_channel
        i:/samples/healthy_output/float_measurements/scoring_interval_samples/*
        f:/samples/healthy_output/float_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/healthy_output/float_measurements/target_frequency_phasor_real_per_channel/*
        i:/samples/healthy_output/pcm_integer_measurements/channels
        i:/samples/healthy_output/pcm_integer_measurements/integer_denominator_D_all_channels
        i:/samples/healthy_output/pcm_integer_measurements/integer_denominator_D_per_channel
        i:/samples/healthy_output/pcm_integer_measurements/integer_output_reference_cross_sum_per_channel/*
        i:/samples/healthy_output/pcm_integer_measurements/integer_reference_error_squared_sum_per_channel/*
        i:/samples/healthy_output/pcm_integer_measurements/integer_reference_squared_sum_per_channel/*
        i:/samples/healthy_output/pcm_integer_measurements/integer_squared_sum_E_all_channels
        i:/samples/healthy_output/pcm_integer_measurements/integer_squared_sum_E_per_channel/*
        f:/samples/healthy_output/pcm_integer_measurements/mean_square_all_channels
        f:/samples/healthy_output/pcm_integer_measurements/mean_square_per_channel/*
        i:/samples/healthy_output/pcm_integer_measurements/pcm_decode_divisor
        f:/samples/healthy_output/pcm_integer_measurements/projection_gain_per_channel/*
        f:/samples/healthy_output/pcm_integer_measurements/reference_NMSE_per_channel/*
        f:/samples/healthy_output/pcm_integer_measurements/reference_error_mean_square_per_channel/*
        i:/samples/healthy_output/pcm_integer_measurements/samples_per_channel
        i:/samples/healthy_output/pcm_integer_measurements/scoring_interval_samples/*
        f:/samples/healthy_output/pcm_measurements/mean_square_per_channel/*
        f:/samples/healthy_output/pcm_measurements/reference_NMSE_per_channel/*
        f:/samples/healthy_output/pcm_measurements/reference_error_mean_square_per_channel/*
        i:/samples/healthy_output/pcm_measurements/samples_per_channel
        i:/samples/healthy_output/pcm_measurements/scoring_interval_samples/*
        f:/samples/healthy_output/pcm_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/healthy_output/pcm_measurements/target_frequency_phasor_real_per_channel/*
        f:/samples/recomputed_output/float_measurements/mean_square_per_channel/*
        f:/samples/recomputed_output/float_measurements/reference_NMSE_per_channel/*
        f:/samples/recomputed_output/float_measurements/reference_error_mean_square_per_channel/*
        i:/samples/recomputed_output/float_measurements/samples_per_channel
        i:/samples/recomputed_output/float_measurements/scoring_interval_samples/*
        f:/samples/recomputed_output/float_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/recomputed_output/float_measurements/target_frequency_phasor_real_per_channel/*
        i:/samples/recomputed_output/pcm_integer_measurements/channels
        i:/samples/recomputed_output/pcm_integer_measurements/integer_denominator_D_all_channels
        i:/samples/recomputed_output/pcm_integer_measurements/integer_denominator_D_per_channel
        i:/samples/recomputed_output/pcm_integer_measurements/integer_output_reference_cross_sum_per_channel/*
        i:/samples/recomputed_output/pcm_integer_measurements/integer_reference_error_squared_sum_per_channel/*
        i:/samples/recomputed_output/pcm_integer_measurements/integer_reference_squared_sum_per_channel/*
        i:/samples/recomputed_output/pcm_integer_measurements/integer_squared_sum_E_all_channels
        i:/samples/recomputed_output/pcm_integer_measurements/integer_squared_sum_E_per_channel/*
        f:/samples/recomputed_output/pcm_integer_measurements/mean_square_all_channels
        f:/samples/recomputed_output/pcm_integer_measurements/mean_square_per_channel/*
        i:/samples/recomputed_output/pcm_integer_measurements/pcm_decode_divisor
        f:/samples/recomputed_output/pcm_integer_measurements/projection_gain_per_channel/*
        f:/samples/recomputed_output/pcm_integer_measurements/reference_NMSE_per_channel/*
        f:/samples/recomputed_output/pcm_integer_measurements/reference_error_mean_square_per_channel/*
        i:/samples/recomputed_output/pcm_integer_measurements/samples_per_channel
        i:/samples/recomputed_output/pcm_integer_measurements/scoring_interval_samples/*
        f:/samples/recomputed_output/pcm_measurements/mean_square_per_channel/*
        f:/samples/recomputed_output/pcm_measurements/reference_NMSE_per_channel/*
        f:/samples/recomputed_output/pcm_measurements/reference_error_mean_square_per_channel/*
        i:/samples/recomputed_output/pcm_measurements/samples_per_channel
        i:/samples/recomputed_output/pcm_measurements/scoring_interval_samples/*
        f:/samples/recomputed_output/pcm_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/recomputed_output/pcm_measurements/target_frequency_phasor_real_per_channel/*
        f:/samples/reference/float_measurements/mean_square_per_channel/*
        f:/samples/reference/float_measurements/reference_NMSE_per_channel/*
        f:/samples/reference/float_measurements/reference_error_mean_square_per_channel/*
        i:/samples/reference/float_measurements/samples_per_channel
        i:/samples/reference/float_measurements/scoring_interval_samples/*
        f:/samples/reference/float_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/reference/float_measurements/target_frequency_phasor_real_per_channel/*
        i:/samples/reference/pcm_integer_measurements/channels
        i:/samples/reference/pcm_integer_measurements/integer_denominator_D_all_channels
        i:/samples/reference/pcm_integer_measurements/integer_denominator_D_per_channel
        i:/samples/reference/pcm_integer_measurements/integer_output_reference_cross_sum_per_channel/*
        i:/samples/reference/pcm_integer_measurements/integer_reference_error_squared_sum_per_channel/*
        i:/samples/reference/pcm_integer_measurements/integer_reference_squared_sum_per_channel/*
        i:/samples/reference/pcm_integer_measurements/integer_squared_sum_E_all_channels
        i:/samples/reference/pcm_integer_measurements/integer_squared_sum_E_per_channel/*
        f:/samples/reference/pcm_integer_measurements/mean_square_all_channels
        f:/samples/reference/pcm_integer_measurements/mean_square_per_channel/*
        i:/samples/reference/pcm_integer_measurements/pcm_decode_divisor
        f:/samples/reference/pcm_integer_measurements/projection_gain_per_channel/*
        f:/samples/reference/pcm_integer_measurements/reference_NMSE_per_channel/*
        f:/samples/reference/pcm_integer_measurements/reference_error_mean_square_per_channel/*
        i:/samples/reference/pcm_integer_measurements/samples_per_channel
        i:/samples/reference/pcm_integer_measurements/scoring_interval_samples/*
        f:/samples/reference/pcm_measurements/mean_square_per_channel/*
        f:/samples/reference/pcm_measurements/reference_NMSE_per_channel/*
        f:/samples/reference/pcm_measurements/reference_error_mean_square_per_channel/*
        i:/samples/reference/pcm_measurements/samples_per_channel
        i:/samples/reference/pcm_measurements/scoring_interval_samples/*
        f:/samples/reference/pcm_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/reference/pcm_measurements/target_frequency_phasor_real_per_channel/*
        f:/samples/stale_output/float_measurements/mean_square_per_channel/*
        f:/samples/stale_output/float_measurements/reference_NMSE_per_channel/*
        f:/samples/stale_output/float_measurements/reference_error_mean_square_per_channel/*
        i:/samples/stale_output/float_measurements/samples_per_channel
        i:/samples/stale_output/float_measurements/scoring_interval_samples/*
        f:/samples/stale_output/float_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/stale_output/float_measurements/target_frequency_phasor_real_per_channel/*
        i:/samples/stale_output/pcm_integer_measurements/channels
        i:/samples/stale_output/pcm_integer_measurements/integer_denominator_D_all_channels
        i:/samples/stale_output/pcm_integer_measurements/integer_denominator_D_per_channel
        i:/samples/stale_output/pcm_integer_measurements/integer_output_reference_cross_sum_per_channel/*
        i:/samples/stale_output/pcm_integer_measurements/integer_reference_error_squared_sum_per_channel/*
        i:/samples/stale_output/pcm_integer_measurements/integer_reference_squared_sum_per_channel/*
        i:/samples/stale_output/pcm_integer_measurements/integer_squared_sum_E_all_channels
        i:/samples/stale_output/pcm_integer_measurements/integer_squared_sum_E_per_channel/*
        f:/samples/stale_output/pcm_integer_measurements/mean_square_all_channels
        f:/samples/stale_output/pcm_integer_measurements/mean_square_per_channel/*
        i:/samples/stale_output/pcm_integer_measurements/pcm_decode_divisor
        f:/samples/stale_output/pcm_integer_measurements/projection_gain_per_channel/*
        f:/samples/stale_output/pcm_integer_measurements/reference_NMSE_per_channel/*
        f:/samples/stale_output/pcm_integer_measurements/reference_error_mean_square_per_channel/*
        i:/samples/stale_output/pcm_integer_measurements/samples_per_channel
        i:/samples/stale_output/pcm_integer_measurements/scoring_interval_samples/*
        f:/samples/stale_output/pcm_measurements/mean_square_per_channel/*
        f:/samples/stale_output/pcm_measurements/reference_NMSE_per_channel/*
        f:/samples/stale_output/pcm_measurements/reference_error_mean_square_per_channel/*
        i:/samples/stale_output/pcm_measurements/samples_per_channel
        i:/samples/stale_output/pcm_measurements/scoring_interval_samples/*
        f:/samples/stale_output/pcm_measurements/target_frequency_phasor_imag_per_channel/*
        f:/samples/stale_output/pcm_measurements/target_frequency_phasor_real_per_channel/*
        i:/schema_version
    """,
    74: """
        f:/numerical_case/candidates/*/DI_dB
        f:/numerical_case/candidates/*/WNG_dB
        f:/numerical_case/candidates/*/absolute_diagonal_loading
        f:/numerical_case/candidates/*/actual_NMSE
        f:/numerical_case/candidates/*/actual_noise_power
        f:/numerical_case/candidates/*/actual_response_real_imag/*
        f:/numerical_case/candidates/*/nominal_NMSE
        f:/numerical_case/candidates/*/nominal_response_real_imag/*
        f:/numerical_case/candidates/*/relative_diagonal_loading
        f:/numerical_case/candidates/*/response_amplitude_error
        f:/numerical_case/candidates/*/weights_real_imag/*/*
        f:/numerical_case/diffuse_coherence/*/*
        f:/numerical_case/noise_covariance_real_imag/*/*/*
        f:/numerical_case/parameters/WNG_lower_limit_dB
        f:/numerical_case/parameters/actual_target_deg
        f:/numerical_case/parameters/frequency_hz
        f:/numerical_case/parameters/independent_sensor_noise_power
        f:/numerical_case/parameters/interference_deg
        f:/numerical_case/parameters/interference_power
        f:/numerical_case/parameters/latency_upper_limit_ms
        f:/numerical_case/parameters/nominal_target_deg
        f:/numerical_case/parameters/positions_m/*/*
        i:/numerical_case/parameters/reference_channel
        f:/numerical_case/parameters/response_amplitude_error_upper_limit
        f:/numerical_case/parameters/sound_speed_m_s
        f:/numerical_case/parameters/target_power
        f:/numerical_case/steering_vectors_real_imag/actual/*/*
        f:/numerical_case/steering_vectors_real_imag/interference/*/*
        f:/numerical_case/steering_vectors_real_imag/nominal/*/*
        f:/plot_data/DI_dB_separate_diffuse_model/*
        f:/plot_data/WNG_dB/*
        f:/plot_data/WNG_lower_limit_dB
        f:/plot_data/actual_NMSE/*
        f:/plot_data/noise_power_normalized_by_target/*
        f:/plot_data/response_amplitude_error/*
        f:/plot_data/response_amplitude_error_upper_limit
        i:/schema_version
    """,
    75: """
        f:/results/experiment/A/imag/*/*
        f:/results/experiment/A/real/*/*
        f:/results/experiment/D/*/*
        f:/results/experiment/G/*/*
        f:/results/experiment/P/*/*
        f:/results/experiment/R/imag/*/*
        f:/results/experiment/R/real/*/*
        f:/results/experiment/W/imag/*/*
        f:/results/experiment/W/real/*/*
        f:/results/experiment/b/*
        f:/results/experiment/csm_dictionary_real/*/*
        f:/results/experiment/csm_gradient_at_solution/*
        f:/results/experiment/csm_target_real/*
        f:/results/experiment/h/*
        f:/results/experiment/losses/CSM/csm_frobenius_squared
        f:/results/experiment/losses/CSM/scan_squared
        f:/results/experiment/losses/GS/csm_frobenius_squared
        f:/results/experiment/losses/GS/scan_squared
        f:/results/experiment/losses/scan/csm_frobenius_squared
        f:/results/experiment/losses/scan/scan_squared
        f:/results/experiment/nonuniform_control/A/imag/*/*
        f:/results/experiment/nonuniform_control/A/real/*/*
        f:/results/experiment/nonuniform_control/D/*/*
        f:/results/experiment/nonuniform_control/G/*/*
        f:/results/experiment/nonuniform_control/P/*/*
        f:/results/experiment/nonuniform_control/W/imag/*/*
        f:/results/experiment/nonuniform_control/W/real/*/*
        f:/results/experiment/nonuniform_control/b/*
        f:/results/experiment/nonuniform_control/h/*
        f:/results/experiment/objective_common_scale
        f:/results/experiment/q_CSM/*
        f:/results/experiment/q_GS/*
        f:/results/experiment/q_scan/*
        f:/results/experiment/scaled_csm_objectives/CSM
        f:/results/experiment/scaled_csm_objectives/GS
        f:/results/experiment/scaled_csm_objectives/scan
        f:/results/plotted/axis_range/*/*
        f:/results/plotted/csm_squared/*/*
        f:/results/plotted/q1/*
        f:/results/plotted/q2/*
        f:/results/plotted/scan_squared/*/*
        i:/schema_version
    """,
    76: """
        f:/anchor_points/in_phase_0deg/diffuseness
        f:/anchor_points/in_phase_0deg/energy_density
        f:/anchor_points/in_phase_0deg/intensity_norm
        f:/anchor_points/open_sphere_first_null_kr
        f:/anchor_points/out_of_phase_180deg/diffuseness
        f:/anchor_points/out_of_phase_180deg/energy_density
        f:/anchor_points/out_of_phase_180deg/intensity_norm
        f:/anchor_points/quadrature_90deg/diffuseness
        f:/anchor_points/quadrature_90deg/energy_density
        f:/anchor_points/quadrature_90deg/intensity_norm
        i:/figure
    """,
    77: """
        f:/anchor_3dir/ls_magnitudes/*
        f:/anchor_3dir/magls_magnitudes/*
        f:/circular_4khz_72dir/ls_magnitudes/*
        f:/circular_4khz_72dir/ls_weighted_mag_mse
        f:/circular_4khz_72dir/ls_weights/imag/*
        f:/circular_4khz_72dir/ls_weights/real/*
        f:/circular_4khz_72dir/magls_magnitudes/*
        f:/circular_4khz_72dir/magls_weighted_mag_mse
        f:/circular_4khz_72dir/magls_weights/imag/*
        f:/circular_4khz_72dir/magls_weights/real/*
        b:/circular_4khz_72dir/monotonically_non_increasing
        f:/circular_4khz_72dir/objective_history/*
        i:/circular_4khz_72dir/parameters/direction_count
        f:/circular_4khz_72dir/parameters/frequency_hz
        i:/circular_4khz_72dir/parameters/iterations
        f:/circular_4khz_72dir/parameters/radius_m
        f:/circular_4khz_72dir/parameters/regularization
        f:/circular_4khz_72dir/parameters/sound_speed_m_s
        f:/circular_4khz_72dir/target_magnitudes/*
        i:/figure
        f:/reference_mixing_points/beta_0.0/ic_out
        f:/reference_mixing_points/beta_0.0/snr_gain_db
        f:/reference_mixing_points/beta_0.2/ic_out
        f:/reference_mixing_points/beta_0.2/snr_gain_db
        f:/reference_mixing_points/beta_0.5/ic_out
        f:/reference_mixing_points/beta_0.5/snr_gain_db
        f:/reference_mixing_points/beta_1.0/ic_out
        f:/reference_mixing_points/beta_1.0/snr_gain_db
    """,
    78: """
        f:/analytic_reference_mse_ch0/wearable_air_wind_2ch
        f:/analytic_reference_mse_ch0/wearable_fused_output
        f:/analytic_reference_mse_ch0/wearable_inear_equalized
        f:/analytic_reference_mse_ch0/wearable_inear_occluded
        f:/anchor_250hz_sinr_db
        f:/anchor_4000hz_sinr_db
        i:/figure
        f:/pcm_reference_mse_ch0/wearable_air_wind_2ch
        f:/pcm_reference_mse_ch0/wearable_fused_output
        f:/pcm_reference_mse_ch0/wearable_inear_equalized
        f:/pcm_reference_mse_ch0/wearable_inear_occluded
    """,
    79: """
        f:/comb_notches_2ms_up_to_4500hz/*
        i:/figure
        f:/pcm_tone_magnitudes_clean_ref/*
        f:/pcm_tone_magnitudes_hearthrough_comb/*
    """,
    80: """
        f:/absent_loss_curves_db/0.0001/*
        f:/absent_loss_curves_db/0.001/*
        f:/absent_loss_curves_db/0.01/*
        b:/bypass_examples/dual_overlap_bypass
        f:/bypass_examples/dual_overlap_ratio_db
        b:/bypass_examples/single_fricative_bypass
        f:/bypass_examples/single_fricative_ratio_db
        i:/figure
        f:/film_input/*/*
        f:/film_output/*/*
        f:/leak_ratio_db_grid/*
        f:/present_neg_sisdr_db/*
        f:/present_waveform_nmse_db/*
        f:/scale_alphas/*
    """,
    81: """
        i:/figure
        f:/fresnel_m_vs_spacing/*
        f:/pcm_integer_nmse/bubble_d10cm
        f:/pcm_integer_nmse/bubble_d2cm_illcond
        f:/pcm_integer_nmse/farfield_dsb
        f:/radial_gain_db/bubble_d10cm/*
        f:/radial_gain_db/bubble_d2cm_reg5e4/*
        f:/radial_gain_db/dsb_d10cm/*
        f:/radii_m/*
        f:/spacings_m/*
        f:/wng_db_vs_spacing/*
    """,
    82: """
        i:/der_cases/imperfect_2c/correct_matched_frames
        f:/der_cases/imperfect_2c/der
        b:/der_cases/imperfect_2c/exceeds_unit_der
        i:/der_cases/imperfect_2c/false_alarm_frames
        f:/der_cases/imperfect_2c/false_alarm_s
        i:/der_cases/imperfect_2c/miss_frames
        f:/der_cases/imperfect_2c/miss_s
        i:/der_cases/imperfect_2c/optimal_permutation/*
        f:/der_cases/imperfect_2c/recording_duration_s
        i:/der_cases/imperfect_2c/ref_speaker_frames
        f:/der_cases/imperfect_2c/ref_speaker_time_s
        i:/der_cases/imperfect_2c/speaker_confusion_frames
        f:/der_cases/imperfect_2c/speaker_confusion_s
        i:/der_cases/over_100/correct_matched_frames
        f:/der_cases/over_100/der
        b:/der_cases/over_100/exceeds_unit_der
        i:/der_cases/over_100/false_alarm_frames
        f:/der_cases/over_100/false_alarm_s
        i:/der_cases/over_100/miss_frames
        f:/der_cases/over_100/miss_s
        i:/der_cases/over_100/optimal_permutation/*
        f:/der_cases/over_100/recording_duration_s
        i:/der_cases/over_100/ref_speaker_frames
        f:/der_cases/over_100/ref_speaker_time_s
        i:/der_cases/over_100/speaker_confusion_frames
        f:/der_cases/over_100/speaker_confusion_s
        i:/der_cases/single_bound/correct_matched_frames
        f:/der_cases/single_bound/der
        b:/der_cases/single_bound/exceeds_unit_der
        i:/der_cases/single_bound/false_alarm_frames
        f:/der_cases/single_bound/false_alarm_s
        i:/der_cases/single_bound/miss_frames
        f:/der_cases/single_bound/miss_s
        i:/der_cases/single_bound/optimal_permutation/*
        f:/der_cases/single_bound/recording_duration_s
        i:/der_cases/single_bound/ref_speaker_frames
        f:/der_cases/single_bound/ref_speaker_time_s
        i:/der_cases/single_bound/speaker_confusion_frames
        f:/der_cases/single_bound/speaker_confusion_s
        i:/der_cases/tsvad_joint/correct_matched_frames
        f:/der_cases/tsvad_joint/der
        b:/der_cases/tsvad_joint/exceeds_unit_der
        i:/der_cases/tsvad_joint/false_alarm_frames
        f:/der_cases/tsvad_joint/false_alarm_s
        i:/der_cases/tsvad_joint/miss_frames
        f:/der_cases/tsvad_joint/miss_s
        i:/der_cases/tsvad_joint/optimal_permutation/*
        f:/der_cases/tsvad_joint/recording_duration_s
        i:/der_cases/tsvad_joint/ref_speaker_frames
        f:/der_cases/tsvad_joint/ref_speaker_time_s
        i:/der_cases/tsvad_joint/speaker_confusion_frames
        f:/der_cases/tsvad_joint/speaker_confusion_s
        i:/figure
        f:/gss_interference_mask/*
        f:/gss_target_exclusive_mask/*
        f:/posteriors_joint/*/*
        f:/posteriors_spectral_only/*/*
        i:/ref_activity/*/*
        f:/single_label_bound/der_lower_bound
        i:/single_label_bound/min_miss_frames
        f:/single_label_bound/min_miss_s
        i:/single_label_bound/overlap_frames
        i:/single_label_bound/ref_speaker_frames
        f:/single_label_bound/ref_speaker_time_s
    """,
    83: """
        f:/architecture_chunk_latency_s/chunked_joint
        f:/architecture_chunk_latency_s/offline_gss
        f:/architecture_chunk_latency_s/streaming_tse
        f:/audio_decision_times_s/*
        f:/consonant_loss_percent/*
        i:/figure
        b:/lag120ms_valid_video_mask/*
        f:/lag120ms_video_arrival_times_s/*
        i:/lag120ms_video_indices/*
        f:/lags_ms/*
        f:/occlusion_alpha_v/*
        f:/occlusion_mask_av/*
        f:/occlusion_mask_fused/*
        b:/sync_valid_video_mask/*
        f:/sync_video_arrival_times_s/*
        i:/sync_video_indices/*
    """,
    84: """
        f:/anchors/circular_aliasing_l2_error
        f:/anchors/hann_50pct_max_abs_error
        f:/anchors/parseval_energy
        f:/anchors/quadrature_phasor_magnitude
        f:/anchors/quadrature_phasor_phase_deg
        f:/anchors/recovered_group_delay_us
        f:/anchors/squared_hann_50pct_ripple
        f:/anchors/zero_padded_max_abs_error
        i:/figure
    """,
    85: """
        f:/anchors/condition_number_2
        f:/anchors/norm_sq_b
        f:/anchors/norm_sq_orthogonal
        f:/anchors/norm_sq_parallel
        f:/anchors/schur_determinant
        f:/anchors/tikhonov_eps0p1_error_norm
        f:/anchors/unregularized_error_norm
        f:/anchors/woodbury_interferer_attenuation_ratio
        i:/figure
    """,
    86: """
        f:/anchors/khatri_rao_vec_max_abs_error
        f:/anchors/max_gev_sinr
        f:/anchors/r_nn_eigenvalues/*
        i:/anchors/sparse_unique_virtual_lags
        f:/anchors/unwhitened_matched_sinr
        f:/anchors/whitening_sinr_gain_db
        i:/figure
    """,
    87: """
        f:/anchors/dft_dc_pseudo_variance
        f:/anchors/dft_nyquist_pseudo_variance
        f:/anchors/improper_pseudo_variance_abs
        f:/anchors/kalman_step1_post_var0
        f:/anchors/kalman_step1_pred_var0
        f:/anchors/lcmv_output_power
        f:/anchors/mvdr_output_power
        f:/anchors/wiener_mmse
        i:/figure
        i:/seed
    """,
}
_NUMERIC_SCHEMAS = {number: {path: frozenset(types) for line in content.split()
    for types, path in [line.split(":", 1)]} for number, content in _NUMERIC_SCHEMA_TEXT.items()}
