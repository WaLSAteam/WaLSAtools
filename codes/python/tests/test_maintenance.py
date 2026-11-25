"""Regression cases for specific defects in the 1.0 series."""
import importlib

import numpy as np
import pytest
from scipy.signal import coherence, csd, welch

from WaLSAtools import WaLSAtools
from WaLSAtools.analysis_modules import WaLSA_speclizer as spectra
from WaLSAtools.analysis_modules.WaLSA_confidence import WaLSA_confidence
from WaLSAtools.analysis_modules.WaLSA_detrend_apod import WaLSA_detrend_apod
from WaLSAtools.analysis_modules.WaLSA_k_omega import averpower, ko_dist
from WaLSAtools.analysis_modules.WaLSA_pod import WaLSA_pod, spod
from WaLSAtools.analysis_modules.WaLSA_wavelet import significance, Morlet


@pytest.mark.parametrize('size', [9, 10])
def test_fft_labels_keep_original_scaling_and_raw_coefficients(size):
    time = 0.25 * np.arange(size)
    x = np.sin(2*np.pi*2*np.arange(size)/size)
    power, frequencies, _, coeff = spectra.getpowerFFT(
        x, time, amplitude=True, nodetrendapod=True, nosignificance=True, silent=True)
    expected = np.fft.fft(x)[1:size//2+1]
    np.testing.assert_array_equal(coeff, expected)
    np.testing.assert_allclose(frequencies, np.arange(1, size//2+1)/(size*.25))
    np.testing.assert_array_equal(power, 2*np.abs(expected)**2 / (1/(2*.25)/(size//2)))


def test_small_permutation_sample_does_not_select_minimum():
    np.testing.assert_array_equal(WaLSA_confidence(np.arange(10)[None,:], .05), [9])


def test_lomb_scargle_default_matches_explicit_confidence():
    t = np.arange(96.) + .1*np.sin(np.arange(96.))
    x = np.sin(2*np.pi*.13*t)
    kwargs = dict(nperm=20, nodetrendapod=True, apod=0, silent=True)
    np.random.seed(11)
    result = spectra.getpowerLS(x, t, **kwargs)
    np.random.seed(11)
    expected = spectra.getpowerLS(x, t, siglevel=.95, **kwargs)
    np.testing.assert_array_equal(result[2], expected[2])


def test_welch_short_signal_and_requested_randomisation():
    time = .5*np.arange(64)
    x = np.cos(2*np.pi*time/4)
    f, p = welch(x, fs=2, nperseg=64, noverlap=32)
    powers, frequencies, thresholds = spectra.welch_psd(x, time, nperm=5, silent=True)
    np.testing.assert_allclose(frequencies, f)
    np.testing.assert_allclose(powers, p)
    assert thresholds.shape == powers.shape
    assert spectra.welch_psd(x, time, nosignificance=True, silent=True)[2] is None


def test_cross_welch_forwards_window_and_overlap_to_all_estimators():
    t = .25*np.arange(256)
    x = np.cos(2*np.pi*t/8)
    y = np.cos(2*np.pi*t/8+np.pi/3)
    f, magnitude, phase, px, py, fc, c = WaLSAtools(
        time=t, data1=x, data2=y, method='welch', window='boxcar',
        nperseg=64, noverlap=11, nosignificance=True, silent=True)
    kwargs = dict(fs=4, window='boxcar', nperseg=64, noverlap=11)
    _, pxy = csd(x, y, **kwargs)
    np.testing.assert_allclose(f, welch(x, **kwargs)[0])
    np.testing.assert_allclose(magnitude, abs(pxy))
    np.testing.assert_allclose(phase, np.angle(pxy, deg=True))
    np.testing.assert_allclose(px, welch(x, **kwargs)[1])
    np.testing.assert_allclose(py, welch(y, **kwargs)[1])
    np.testing.assert_allclose(fc, coherence(x, y, **kwargs)[0])
    np.testing.assert_allclose(c, coherence(x, y, **kwargs)[1])


@pytest.mark.parametrize('flag', ['dominantfreq','averagedpower'])
def test_single_cube_flag_keeps_existing_four_outputs(flag):
    t = np.arange(16.)
    x = np.sin(2*np.pi*t/8)
    cube = x[:,None,None]*np.ones((1,2,2))
    d, p, f, per_pixel = WaLSAtools(
        cube, t, method='fft', nodetrendapod=True, silent=True, **{flag: True})
    np.testing.assert_allclose(d, .125)
    np.testing.assert_allclose(p, per_pixel.mean(axis=(0,1)))
    assert f[p.argmax()] == .125


def test_mean_only_detrending_keeps_ramp():
    x = np.arange(10.)
    np.testing.assert_allclose(WaLSA_detrend_apod(x, meantemporal=True, apod=0, silent=True), x-x.mean())


def test_eemd_forwards_explicit_count_without_changing_effective_default(monkeypatch):
    calls = []
    class Engine:
        def __init__(self, trials=100):
            calls.append(trials)
        def noise_seed(self, seed):
            pass
        def eemd(self, x, time):
            return x
    monkeypatch.setattr(spectra, 'EEMD', Engine)
    spectra.apply_eemd([1,2], [0,1])
    spectra.apply_eemd([1,2], [0,1], num_realizations=7)
    assert calls == [100,7]


def test_emd_no_preprocessing_keeps_public_outputs():
    time = .25*np.arange(64)
    signal = np.cos(2*np.pi*time/4)
    result = spectra.getEMD_HHT(signal, time, nodetrendapod=True, nperm=3, silent=True)
    assert len(result) == 8
    np.testing.assert_allclose(result[5], spectra.apply_emd(signal,time))


def test_wavelet_background_same_for_series_and_variance():
    scales = 2.**(np.arange(8)/4)
    x = np.arange(32.)
    dt, lag = .25, .7
    _, bg = significance(x, dt, scales, alpha=lag)
    _, scalar = significance(x.var(), dt, scales, alpha=lag)
    f = dt/(scales*Morlet().flambda())
    expected = x.var()*(1-lag**2)/(1+lag**2-2*lag*np.cos(2*np.pi*f))
    np.testing.assert_allclose(bg, expected)
    np.testing.assert_allclose(scalar, expected)
    first, _ = significance(x, dt, scales, sigma_test=1, alpha=lag, dof=32)
    second, _ = significance(x, dt/60, scales/60, sigma_test=1, alpha=lag, dof=32)
    np.testing.assert_allclose(first, second)


def test_square_odd_radial_grid_and_true_bin_mean():
    size = 9
    x, y = np.indices((size,size))
    time = np.arange(16)[:,None,None]
    cube = np.cos(2*np.pi*(3*time/16+2*x/size+y/size))
    radius = np.hypot(np.minimum(x,size-x),np.minimum(y,size-y))
    bins = np.rint(radius * (size//2+1)/(size/2)).astype(int)
    np.testing.assert_array_equal(ko_dist(size,size),bins)
    transform = np.fft.fftn(cube)
    radial = averpower(cube)
    shell = bins[2,1]
    expected = np.mean(abs(transform[3][bins==shell])**2)
    np.testing.assert_allclose(radial[shell,3],expected)


def test_komega_no_filter_reconstructs_with_or_without_maps():
    data = np.random.default_rng(1).normal(size=(32,8,8))
    time = .5*np.arange(len(data))
    for maps in [False, True]:
        result = WaLSAtools(data, time, method='k-omega', filtering=True,
                            no_temporal_filt=True, no_spatial_filt=True, processing_maps=maps)
        np.testing.assert_allclose(result[3],data,atol=1e-12)
        assert (result[8] is not None)==maps


def test_pod_retains_small_variation_and_uses_supplied_cadence():
    t = .1*np.arange(64)
    x = 1e9+np.cos(2*np.pi*1.25*t)
    data = np.broadcast_to(x[:,None,None],(64,2,2))
    result = WaLSA_pod(data,t,nperm=3,top_frequencies=[],silent=True)
    np.testing.assert_allclose(result['input_data'][:,0,0],x-x.mean(),atol=2e-7)
    np.testing.assert_allclose(result['cumulative_eigenvalues'][0],100)
    np.testing.assert_allclose(result['frequencies'][np.argmax(result['combined_welch_psd'])],1.25,atol=.01)


def test_spod_zero_width_on_integer_cube_matches_float_inputs():
    x = np.random.default_rng(3).integers(0,20,(8,3,3))
    first = spod(x, filter_size=0, silent=True)
    second = spod(x.astype(float), filter_size=0, silent=True)
    np.testing.assert_allclose(first[0],second[0])
    np.testing.assert_allclose(first[1],second[1])


def test_interactive_quickstart_still_calls_menu(monkeypatch):
    entry = importlib.import_module('WaLSAtools.WaLSAtools')
    calls = []
    monkeypatch.setattr(entry,'interactive',lambda: calls.append(True))
    assert repr(WaLSAtools) == ''
    assert WaLSAtools() is None
    assert entry.main([]) == 0
    assert calls == [True]*3


def test_wavelet_reconstruction_uses_matching_scale_step_and_handles_empty_band():
    from WaLSAtools.analysis_modules.WaLSA_detrend_apod import WaLSA_wave_recon, getWavelet

    t = .5*np.arange(256)
    x = np.cos(2*np.pi*t/16)
    _, period, coi, scales, wave = getWavelet(x,t,dj=1/32)
    keep = ((period[:,None] > 10) & (period[:,None] <= 20) & (period[:,None] < coi))
    expected = ((keep * wave.real) / np.sqrt(scales[:,None])).sum(axis=0)
    expected *= np.sqrt(.5)/32/(.776*np.pi**(-.25))
    np.testing.assert_allclose(WaLSA_wave_recon(x,.5,lo_cutoff=10,hi_cutoff=20), expected)
    np.testing.assert_allclose(WaLSA_wave_recon(x,.5,lo_cutoff=1000,hi_cutoff=2000), 0)


def test_emd_honours_explicit_preprocessing_and_can_select_no_modes(monkeypatch):
    t = .25*np.arange(64)
    x = 10 + np.cos(2*np.pi*t/4)
    seen = []
    real = spectra.apply_emd

    def capture(signal, time):
        seen.append(signal.copy())
        return real(signal, time)

    monkeypatch.setattr(spectra,'apply_emd',capture)
    monkeypatch.setattr(spectra,'test_imf_significance',lambda modes, draws: [.8]*len(modes))
    args = dict(nperm=3, silent=True, significant_imfs=True)
    result = spectra.getEMD_HHT(x,t,**args)
    np.testing.assert_array_equal(seen[0],x)
    assert not len(result[5])
    np.testing.assert_array_equal(result[0], np.zeros(100))
    spectra.getEMD_HHT(x,t,meantemporal=True,apod=0,**args)
    np.testing.assert_allclose(seen[1],x-x.mean())


@pytest.mark.parametrize('size, length', [(8,32),(9,31)])
@pytest.mark.parametrize('torus', [False, True])
def test_komega_passband_selects_even_and_odd_plane_waves(size,length,torus):
    time = .5*np.arange(length)
    a,b = np.indices((size,size))
    f = 3/(length*.5)
    temporal = f*time[:,None,None]
    wanted = np.cos(2*np.pi*(temporal + 2*a/size+b/size))
    radial_stop = np.cos(2*np.pi*(temporal + a/size))
    time_stop = np.cos(2*np.pi*(2*temporal + 2*a/size+b/size))
    result = WaLSAtools(wanted+radial_stop+time_stop,time,method='k-omega',
                        filtering=True,spatial_torus=torus,temporal_torus=torus,
                        k1=2*np.pi*2/size,k2=2*np.pi*2.5/size,f1=f,f2=f,
                        processing_maps=True)
    np.testing.assert_allclose(result[3],wanted,atol=1e-12)
    np.testing.assert_allclose(result[1],np.pi * np.arange(1,len(result[1])+1)/len(result[1]))
    np.testing.assert_allclose(result[2],np.arange(1,len(result[2])+1)/(length*.5))


@pytest.mark.parametrize('size, length', [(32,32), (33,31)])
def test_komega_legacy_smooth_tapers_with_paired_fft_bins(size, length):
    time = .5*np.arange(length)
    x = np.arange(size)[None,:,None]
    phase = 2*np.pi*4*time[:,None,None]/(length*.5)
    modes = [np.broadcast_to(np.cos(phase+2*np.pi*k*x/size), (length,size,size))
             for k in (2,4,6,7,10,11)]
    dk = 2*np.pi/size
    df = 1/(length*.5)
    result = WaLSAtools(sum(modes),time,method='k-omega',filtering=True,
                        k1=3*dk,k2=9*dk,f1=3*df,f2=9*df,processing_maps=True)

    # The legacy 3-to-9-bin torus contains nested rings of widths 0, 1 and 3.
    spatial_gain = [0, 1/3, 2/3, 1, 1/3, 0]
    # At bin 4 the half-open Gaussian slices provide weights at offsets
    # 1 and 5 in a six-bin kernel. A real wave uses their mean.
    temporal_gain = (np.exp(-2/18) + np.exp(-6/18))/2
    expected = temporal_gain * sum(gain*mode for gain, mode in zip(spatial_gain,modes))
    np.testing.assert_allclose(result[3],expected,atol=1e-12)
    for index in (5, 8):
        mask = result[index]['data'] if index == 5 else result[index]
        unshifted = np.fft.ifftshift(mask)
        reflected = np.ix_(*[(-np.arange(n)) % n for n in unshifted.shape])
        np.testing.assert_array_equal(unshifted,unshifted[reflected])


def test_komega_smooth_filter_keeps_single_dc_and_nyquist_weights():
    time = np.arange(32.)
    signal = 2 + .5*np.cos(np.pi*time)
    cube = np.broadcast_to(signal[:,None,None], (32,8,8))
    result = WaLSAtools(cube,time,method='k-omega',filtering=True,
                        no_spatial_filt=True,f1=.25/32,f2=.5)
    # A half-open 16-bin Gaussian with sigma=3 has an endpoint at 7.5 bins
    # from its centre and a discrete peak at 0.5 bins from centre.
    gain = np.exp(-(7.5**2 - .5**2)/18)
    np.testing.assert_allclose(result[3],gain*cube,atol=1e-12)


@pytest.mark.parametrize('raw', [False, True])
def test_odd_short_fft_with_enabled_randomisation(raw):
    t = np.arange(9.)
    x = np.sin(2*np.pi*2*t/9)
    power,freq,confidence,_ = spectra.getpowerFFT(
        x,t,nperm=5,nodetrendapod=raw,silent=True)
    assert freq[power.argmax()] == pytest.approx(2/9)
    assert np.isfinite(confidence).all()


def test_eemd_keeps_trial_settings_and_surrogate_count_distinct(monkeypatch):
    t = .25*np.arange(64)
    x = np.cos(2*np.pi*t/4)
    engines = []

    class Engine:
        def __init__(self,trials):
            self.trials = trials
            self.calls = 0
            engines.append(self)

        def noise_seed(self, seed):
            pass

        def eemd(self, values, times):
            self.calls += 1
            return np.asarray([values])

    monkeypatch.setattr(spectra,'EEMD',Engine)
    spectra.getEMD_HHT(x,t,EEMD=True,nperm=3,num_realizations=7,noise_std=.1,silent=True)
    assert [(e.trials,e.calls,e.noise_width) for e in engines] == [(7,1,.1),(7,3,.1)]


def test_fft_preserves_historical_power_for_float32_cadence():
    t = (np.arange(999)*0.01).astype(np.float32)
    x = np.sin(2*np.pi*50*np.arange(len(t))/len(t)).astype(np.float32)
    power,freq,_,_ = spectra.getpowerFFT(
        x,t,nodetrendapod=True,nosignificance=True,silent=True)
    cadence = np.median(np.diff(t))
    bins = np.arange(len(t)//2+1)[1:]
    old = 1./(cadence*2)*bins/(len(t)//2)
    expected = 2*abs(np.fft.fft(x)[1:len(bins)+1])**2/old[0]
    np.testing.assert_array_equal(power,expected)
    np.testing.assert_allclose(freq,bins/(len(t)*float(cadence)),rtol=1e-15)


def test_installed_komega_plot_without_example_colour_table(tmp_path,monkeypatch):
    from matplotlib import pyplot as plt
    from WaLSAtools import WaLSA_plot_k_omega

    monkeypatch.chdir(tmp_path)
    image = np.arange(9.).reshape(3,3)
    for colour in ['viridis','plasma']:
        ax = WaLSA_plot_k_omega(image,np.arange(3.),np.arange(3.),
                                 cmap=colour,dxpix=16,dypix=16)
        assert ax.images[0].cmap.name == colour
        ax.figure.canvas.draw()
        plt.close(ax.figure)


@pytest.mark.parametrize('bypass', [False, True])
def test_fft_randomisation_uses_observed_preprocessing_on_raw_permutations(bypass):
    t = np.arange(32.)
    signal = np.linspace(-1, 5, t.size) + np.sin(2 * np.pi * t / 8)
    options = dict(apod=.25, polyfit=1, cadence=1., silent=True)

    def process(samples):
        return samples if bypass else WaLSA_detrend_apod(samples, **options)

    def statistic(samples):
        return 2 * np.abs(np.fft.fft(samples)[1:t.size // 2 + 1]) ** 2

    np.random.seed(91)
    power, _, level, _ = spectra.getpowerFFT(
        signal, t, apod=.25, polyfit=1, nodetrendapod=bypass,
        nperm=20, siglevel=.95, silent=True)
    np.random.seed(91)
    shuffled = np.column_stack([
        statistic(process(np.random.permutation(signal))) for _ in range(20)
    ])
    # At the 5% upper tail, historical indexing selects the largest of 20 draws.
    legacy_df = 1 / t.size
    np.testing.assert_array_equal(power, statistic(process(signal)) / legacy_df)
    np.testing.assert_array_equal(level, shuffled.max(axis=1) / legacy_df)
