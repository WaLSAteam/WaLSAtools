# -----------------------------------------------------------------------------------------------------
# WaLSAtools - Wave analysis tools
# Copyright (C) 2025 WaLSA Team - Shahin Jafarzadeh et al.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# Note: If you use WaLSAtools for research, please consider citing:
# Jafarzadeh, S., Jess, D. B., Stangalini, M. et al. 2025, Nature Reviews Methods Primers, in press.
# -----------------------------------------------------------------------------------------------------

import numpy as np # type: ignore
from astropy.timeseries import LombScargle # type: ignore
from tqdm import tqdm # type: ignore
from .WaLSA_detrend_apod import WaLSA_detrend_apod # type: ignore
from .WaLSA_confidence import WaLSA_confidence # type: ignore
from .WaLSA_wavelet import cwt, significance # type: ignore
from scipy.signal import welch # type: ignore

# -------------------------------------- Main Function ----------------------------------------
def WaLSA_speclizer(signal=None, time=None, method=None, 
                    dominantfreq=False, averagedpower=False, **kwargs):
    """
    Main function to prepare data and call the specific spectral analysis method.
    
    Parameters:
        signal (array): The input signal (1D or 3D).
        time (array): The time array of the signal.
        method (str): The method to use for spectral analysis ('fft', 'lombscargle', etc.)
        **kwargs: Additional parameters for data preparation and analysis methods.
    
    Returns:
        Power spectrum, frequencies, and significance (if applicable).
    """
    if not dominantfreq or not averagedpower:
        if method == 'fft':
            return getpowerFFT(signal, time=time, **kwargs)
        elif method == 'lombscargle':
            return getpowerLS(signal, time=time, **kwargs)
        elif method == 'wavelet':
            return getpowerWavelet(signal, time=time, **kwargs)
        elif method == 'welch':
            return welch_psd(signal, time=time, **kwargs)
        elif method == 'emd':
            return getEMD_HHT(signal, time=time, **kwargs)
        else:
            raise ValueError(f"Unknown method: {method}")
    else:
        return get_dominant_averaged(signal, time=time, method=method, **kwargs)

# -------------------------------------- FFT ----------------------------------------
def getpowerFFT(signal, time, **kwargs):
    """
    Perform FFT analysis on the given signal.
    
    Parameters:
        signal (array): The input signal (1D).
        time (array): The time array corresponding to the signal.
        siglevel (float): Significance level for the confidence intervals. Default: 0.95.
        nperm (int): Number of permutations for significance testing. Default: 1000.
        nosignificance (bool): If True, skip significance calculation. Default: False.
        apod (float): Extent of apodization edges (of a Tukey window). Default: 0.1.
        pxdetrend (int): Subtract linear trend with time per pixel. Options: 1 (simple) or 2 (advanced). Default: 2.
        polyfit (int): Degree of polynomial fit for detrending the data. If set, a polynomial fit (instead of linear) is applied. Default: None.
        meantemporal (bool): If True, apply simple temporal detrending by subtracting the mean signal from the data, skipping fitting procedures. Default: False.
        meandetrend (bool): If True, subtract the linear trend with time for the image means (spatial detrending). Default: False.
        recon (bool): If True, perform Fourier reconstruction of the input time series. This does not preserve amplitudes but is useful for examining frequencies far from the low-frequency range.  Default: False.
        resample_original (bool): If True, and if recon set True, approximate values close to the original are returned for comparison. Default: False.
        nodetrendapod (bool): If True, neither detrending nor apodization is performed. Default: False.
        amplitude (bool): If True, return the amplitudes of the Fourier transform. Default: False.
        silent (bool): If True, suppress print statements. Default: False.
        **kwargs: Additional parameters for the analysis method.
    
    Returns:
        Power spectrum, frequencies, significance, amplitudes
    """
    # Define default values for the optional parameters
    defaults = {
        'siglevel': 0.95,
        'significance': None,
        'nperm': 1000,
        'nosignificance': False,
        'apod': 0.1,
        'pxdetrend': 2,
        'meandetrend': False,
        'polyfit': None,
        'meantemporal': False,
        'recon': False,
        'resample_original': False,
        'nodetrendapod': False,
        'amplitude': False,
        'silent': False
    }

    # Update defaults with any user-provided keyword arguments
    params = {**defaults, **kwargs}
    
    params['siglevel'] = 1 - params['siglevel'] # different convention

    tdiff = np.diff(time)
    cadence = np.median(tdiff)

    # Perform detrending and apodization
    if not params['nodetrendapod']:
        apocube = WaLSA_detrend_apod(
            signal, 
            apod=params['apod'], 
            meandetrend=params['meandetrend'], 
            pxdetrend=params['pxdetrend'], 
            polyfit=params['polyfit'], 
            meantemporal=params['meantemporal'], 
            recon=params['recon'], 
            cadence=cadence, 
            resample_original=params['resample_original'], 
            silent=params['silent']
        )
    else:
        apocube = signal
    
    nt = len(apocube)  # Length of the time series (1D)

    # Calculate the frequencies
    frequencies = 1. / (cadence * 2) * np.arange(int(np.floor(nt / 2)) + 1) / (nt // 2)
    frequencies = frequencies[1:]
    powermap = np.zeros(len(frequencies))
    signal = apocube
    spec = np.fft.fft(signal)
    power = 2 * np.abs(spec[1:len(frequencies)]) ** 2
    powermap[:] = power / frequencies[0]

    if params['amplitude']:
        amplitudes = np.zeros((len(signal), len(frequencies)), dtype=np.complex_)
        amplitudes = spec[1:len(frequencies)]
    else: 
        amplitudes = None

    # Calculate significance if requested
    if not params['nosignificance']:
        ps_perm = np.zeros((len(frequencies), params['nperm']))
        for ip in range(params['nperm']):
            perm_signal = np.random.permutation(signal)  # Permuting the original signal
            apocube = WaLSA_detrend_apod(
                perm_signal, 
                apod=params['apod'], 
                meandetrend=params['meandetrend'], 
                pxdetrend=params['pxdetrend'], 
                polyfit=params['polyfit'], 
                meantemporal=params['meantemporal'], 
                recon=params['recon'], 
                cadence=cadence, 
                resample_original=params['resample_original'], 
                silent=True
            )
            perm_spec = np.fft.fft(perm_signal)
            perm_power = 2 * np.abs(perm_spec[1:len(frequencies)]) ** 2
            ps_perm[:, ip] = perm_power

        # Correct call to WaLSA_confidence
        significance = WaLSA_confidence(ps_perm, siglevel=params['siglevel'], nf=len(frequencies))
        significance = significance / frequencies[0]
    else:
        significance = None

    if not params['silent']:
        print("FFT processed.")

    return powermap, frequencies, significance, amplitudes

# -------------------------------------- Lomb-Scargle ----------------------------------------
def getpowerLS(signal, time, **kwargs):
    """
    Perform Lomb-Scargle analysis on the given signal.
    
    Parameters:
        signal (array): The input signal (1D).
        time (array): The time array corresponding to the signal.
        siglevel (float): Significance level for the confidence intervals. Default: 0.95.
        nperm (int): Number of permutations for significance testing. Default: 1000.
        dy (array): Errors or observational uncertainties associated with the time series.
        fit_mean (bool): If True, include a constant offset as part of the model at each frequency. This improves accuracy, especially for incomplete phase coverage.
        center_data (bool): If True, pre-center the data by subtracting the weighted mean of the input data. This is especially important if fit_mean=False.
        nterms (int): Number of terms to use in the Fourier fit. Default: 1.
        normalization (str): The normalization method for the periodogram. Options: 'standard', 'model', 'log', 'psd'. Default: 'standard'.
        nosignificance (bool): If True, skip significance calculation. Default: False.
        apod (float): Extent of apodization edges (of a Tukey window). Default: 0.1.
        pxdetrend (int): Subtract linear trend with time per pixel. Options: 1 (simple) or 2 (advanced). Default: 2.
        polyfit (int): Degree of polynomial fit for detrending the data. If set, a polynomial fit (instead of linear) is applied. Default: None.
        meantemporal (bool): If True, apply simple temporal detrending by subtracting the mean signal from the data, skipping fitting procedures. Default: False.
        meandetrend (bool): If True, subtract the linear trend with time for the image means (spatial detrending). Default: False.
        recon (bool): If True, perform Fourier reconstruction of the input time series. This does not preserve amplitudes but is useful for examining frequencies far from the low-frequency range.  Default: False.
        resample_original (bool): If True, and if recon set True, approximate values close to the original are returned for comparison. Default: False.
        nodetrendapod (bool): If True, neither detrending nor apodization is performed. Default: False..
        silent (bool): If True, suppress print statements. Default: False.
        **kwargs: Additional parameters for the analysis method.
    
    Returns:
        Power spectrum, frequencies, and significance (if applicable).
    """
    # Define default values for the optional parameters
    defaults = {
        'siglevel': 0.05,
        'nperm': 1000,
        'nosignificance': False,
        'apod': 0.1,
        'pxdetrend': 2,
        'meandetrend': False,
        'polyfit': None,
        'meantemporal': False,
        'recon': False,
        'resample_original': False,
        'silent': False,  # Ensure 'silent' is included in defaults
        'nodetrendapod': False,
        'dy': None, # Error or sequence of observational errors associated with times t
        'fit_mean': False, # If True include a constant offset as part of the model at each frequency. This can lead to more accurate results, especially in the case of incomplete phase coverage
        'center_data': True, # If True pre-center the data by subtracting the weighted mean of the input data. This is especially important if fit_mean = False
        'nterms': 1, # Number of terms to use in the Fourier fit
        'normalization': 'standard' # The normalization to use for the periodogram. Options are 'standard', 'model', 'log', 'psd'
    }

    # Update defaults with any user-provided keyword arguments
    params = {**defaults, **kwargs}
    
    params['siglevel'] = 1 - params['siglevel'] # different convention
    
    tdiff = np.diff(time)
    cadence = np.median(tdiff)

    # Perform detrending and apodization
    if not params['nodetrendapod']:
        apocube = WaLSA_detrend_apod(
            signal, 
            apod=params['apod'], 
            meandetrend=params['meandetrend'], 
            pxdetrend=params['pxdetrend'], 
            polyfit=params['polyfit'], 
            meantemporal=params['meantemporal'], 
            recon=params['recon'], 
            cadence=cadence, 
            resample_original=params['resample_original'], 
            silent=params['silent']
        )
    else:
        apocube = signal
    
    frequencies, power = LombScargle(time, apocube, dy=params['dy'], fit_mean=params['fit_mean'], 
                                     center_data=params['center_data'], nterms=params['nterms'], 
                                     normalization=params['normalization']).autopower()

    # Calculate significance if needed
    if not params['nosignificance']:
        ps_perm = np.zeros((len(frequencies), params['nperm']))
        for ip in range(params['nperm']):
            perm_signal = np.random.permutation(signal)  # Permuting the original signal
            apocubep = WaLSA_detrend_apod(
                perm_signal, 
                apod=params['apod'], 
                meandetrend=params['meandetrend'], 
                pxdetrend=params['pxdetrend'], 
                polyfit=params['polyfit'], 
                meantemporal=params['meantemporal'], 
                recon=params['recon'], 
                cadence=cadence, 
                resample_original=params['resample_original'], 
                silent=True
            )
            frequencies, perm_power = LombScargle(time, apocubep, dy=params['dy'], fit_mean=params['fit_mean'], 
                                     center_data=params['center_data'], nterms=params['nterms'], 
                                     normalization=params['normalization']).autopower()
            ps_perm[:, ip] = perm_power
        significance = WaLSA_confidence(ps_perm, siglevel=params['siglevel'], nf=len(frequencies))
    else:
        significance = None

    if not params['silent']:
        print("Lomb-Scargle processed.")

    return power, frequencies, significance

# -------------------------------------- Wavelet ----------------------------------------
def getpowerWavelet(signal, time, **kwargs):
    """
    Perform wavelet analysis using the pycwt package.
    
    Parameters:
        signal (array): The input signal (1D).
        time (array): The time array corresponding to the signal.
        siglevel (float): Significance level for the confidence intervals. Default: 0.95.
        nperm (int): Number of permutations for significance testing. Default: 1000.
        mother (str): The mother wavelet function to use. Default: 'morlet'.
        GWS (bool): If True, calculate the Global Wavelet Spectrum. Default: False.
        RGWS (bool): If True, calculate the Refined Global Wavelet Spectrum (time-integrated power, excluding COI and insignificant areas). Default: False.
        dj (float): Scale spacing. Smaller values result in better scale resolution but slower calculations. Default: 0.025.
        s0 (float): Initial (smallest) scale of the wavelet. Default: 2 * dt.
        J (int): Number of scales minus one. Scales range from s0 up to s0 * 2**(J * dj), giving a total of (J + 1) scales. Default: (log2(N * dt / s0)) / dj.
        lag1 (float): Lag-1 autocorrelation. Default: 0.0.
        apod (float): Extent of apodization edges (of a Tukey window). Default: 0.1.
        pxdetrend (int): Subtract linear trend with time per pixel. Options: 1 (simple) or 2 (advanced). Default: 2.
        polyfit (int): Degree of polynomial fit for detrending the data. If set, a polynomial fit (instead of linear) is applied. Default: None.
        meantemporal (bool): If True, apply simple temporal detrending by subtracting the mean signal from the data, skipping fitting procedures. Default: False.
        meandetrend (bool): If True, subtract the linear trend with time for the image means (spatial detrending). Default: False.
        recon (bool): If True, perform Fourier reconstruction of the input time series. This does not preserve amplitudes but is useful for examining frequencies far from the low-frequency range.  Default: False.
        resample_original (bool): If True, and if recon set True, approximate values close to the original are returned for comparison. Default: False.
        nodetrendapod (bool): If True, neither detrending nor apodization is performed. Default: False.
        silent (bool): If True, suppress print statements. Default: False.
        **kwargs: Additional parameters for the analysis method.
    
    Returns:
        power: The wavelet power spectrum.
        periods: Corresponding periods.
        sig_slevel: The significance levels.
        coi: The cone of influence.
        Optionally, if global_power=True:
        global_power: Global wavelet power spectrum.
        global_conf: Confidence levels for the global wavelet spectrum.
        Optionally, if RGWS=True:
        rgws_periods: Periods for the refined global wavelet spectrum.
        rgws_power: Refined global wavelet power spectrum.
    """
    # Define default values for the optional parameters similar to IDL
    defaults = {
        'siglevel': 0.95,
        'mother': 'morlet',  # Morlet wavelet as the mother function
        'dj': 0.025,  # Scale spacing
        's0': -1,  # Initial scale
        'J': -1,  # Number of scales
        'lag1': 0.0,  # Lag-1 autocorrelation
        'apod': 0.1,  # Tukey window apodization
        'silent': False,
        'pxdetrend': 2,
        'meandetrend': False,
        'polyfit': None,
        'meantemporal': False,
        'recon': False,
        'resample_original': False,
        'GWS': False,  # If True, calculate global wavelet spectrum
        'RGWS': False,  # If True, calculate refined global wavelet spectrum (excluding COI)
        'nperm': 1000,   # Number of permutations for significance calculation
        'nodetrendapod': False
    }

    # Update defaults with any user-provided keyword arguments
    params = {**defaults, **kwargs}

    tdiff = np.diff(time)
    cadence = np.median(tdiff)

    # Perform detrending and apodization
    if not params['nodetrendapod']:
        apocube = WaLSA_detrend_apod(
            signal, 
            apod=params['apod'], 
            meandetrend=params['meandetrend'], 
            pxdetrend=params['pxdetrend'], 
            polyfit=params['polyfit'], 
            meantemporal=params['meantemporal'], 
            recon=params['recon'], 
            cadence=cadence, 
            resample_original=params['resample_original'], 
            silent=params['silent']
        )
    else:
        apocube = signal

    n = len(apocube)
    dt = cadence

    # Standardize the signal before the wavelet transform
    std_signal = apocube.std()
    norm_signal = apocube / std_signal

    # Determine the initial scale s0 if not provided
    if params['s0'] == -1:
        params['s0'] = 2 * dt

    # Determine the number of scales J if not provided
    if params['J'] == -1:
        params['J'] = int((np.log(float(n) * dt / params['s0']) / np.log(2)) / params['dj'])

    # Perform wavelet transform
    W, scales, frequencies, coi, _, _ = cwt(
        norm_signal,
        dt,
        dj=params['dj'],
        s0=params['s0'],
        J=params['J'],
        wavelet=params['mother']
    )

    power = np.abs(W) ** 2  # Wavelet power spectrum
    periods = 1 / frequencies  # Convert frequencies to periods

    # Calculate the significance levels
    signif, _ = significance(
        1.0,  # Normalized variance
        dt,
        scales,
        0,  # Ignored for now, used for background noise (e.g., red noise)
        params['lag1'],
        significance_level=params['siglevel'],
        wavelet=params['mother']
    )

    # Calculate the significance level for each power value
    sig_matrix = np.ones([len(scales), n]) * signif[:, None]
    sig_slevel = power / sig_matrix

    # Calculate global power if requested
    if params['GWS']:
        global_power = np.mean(power, axis=1)  # Average along the time axis

        dof = n - scales   # the -scale corrects for padding at edges
        global_conf, _ = significance(
            norm_signal,                   # time series data
            dt,              # time step between values
            scales,          # scale vector
            1,               # sigtest = 1 for "time-average" test
            0.0,
            significance_level=params['siglevel'],
            dof=dof,
            wavelet=params['mother']
        )
    else:
        global_power = None
        global_conf = None

    # Calculate refined global wavelet spectrum (RGWS) if requested
    if params['RGWS']:
        isig = sig_slevel < 1.0
        power[isig] = np.nan
        ipower = np.full_like(power, np.nan)
        for i in range(len(coi)):
            pcol = power[:, i]
            valid_idx = periods < coi[i]
            ipower[valid_idx, i] = pcol[valid_idx]
        rgws_power = np.nansum(ipower, axis=1)
    else:
        rgws_power = None

    if not params['silent']:
        print("Wavelet (" + params['mother'] + ") processed.")

    return power, periods, sig_slevel, coi, global_power, global_conf, rgws_power

# -------------------------------------- Welch ----------------------------------------
def welch_psd(signal, time, **kwargs):
    """
    Calculate Welch Power Spectral Density (PSD) and significance levels.
    
    Parameters:
        signal (array): The 1D time series signal.
        time (array): The time array corresponding to the signal.
        nperseg (int, optional): Length of each segment for analysis. Default: 256.
        noverlap (int, optional): Number of points to overlap between segments. Default: 128.
        window (str, optional): Type of window function used in the Welch method. Default: 'hann'.
        siglevel (float, optional): Significance level for confidence intervals. Default: 0.95.
        nperm (int, optional): Number of permutations for significance testing. Default: 1000.
        silent (bool, optional): If True, suppress print statements. Default: False.
        **kwargs: Additional parameters for the analysis method.
    
    Returns:
        frequencies: Frequencies at which the PSD is estimated.
        psd: Power Spectral Density values.
        significance: Significance levels for the PSD.
    """
    # Define default values for the optional parameters similar to FFT
    defaults = {
        'siglevel': 0.95,
        'nperm': 1000,  # Number of permutations for significance calculation
        'window': 'hann',  # Window type for Welch method
        'nperseg': 256,  # Length of each segment
        'silent': False,
        'noverlap': 128  # Number of points to overlap between segments
    }

    # Update defaults with any user-provided keyword arguments
    params = {**defaults, **kwargs}

    tdiff = np.diff(time)
    cadence = np.median(tdiff)

    # Calculate Welch PSD
    # DC no option for detrend, it might be updated later on
    frequencies, psd = welch(signal, fs=1.0 / cadence, window=params['window'], nperseg=params['nperseg'], noverlap=params['noverlap'])

    # Calculate significance levels using permutation
    ps_perm = np.zeros((len(frequencies), params['nperm']))
    for ip in range(params['nperm']):
        perm_signal = np.random.permutation(signal)  # Permuting the original signal
        _, perm_psd = welch(perm_signal, fs=1.0 / cadence, window=params['window'], nperseg=params['nperseg'], noverlap=params['noverlap'])
        ps_perm[:, ip] = perm_psd

    # Calculate significance level for Welch PSD
    significance = np.percentile(ps_perm, params['siglevel'] * 100, axis=1)

    if not params['silent']:
        print("Welch processed.")

    return psd, frequencies, significance

# -------------------------------------- EMD / EEMD ----------------------------------------
import numpy as np # type: ignore
from .PyEMD import EMD, EEMD # type: ignore
from scipy.stats import norm # type: ignore
from scipy.signal import hilbert, welch # type: ignore
from scipy.signal import find_peaks # type: ignore
from scipy.fft import fft, fftfreq # type: ignore

# Function to apply EMD and return Intrinsic Mode Functions (IMFs)
def apply_emd(signal, time):
    emd = EMD()
    emd.FIXE_H = 5
    imfs = emd.emd(signal, time) # max_imfs=7, emd.FIXE_H = 10
    return imfs

# Function to apply EEMD and return Intrinsic Mode Functions (IMFs)
def apply_eemd(signal, time, noise_std=0.2, num_realizations=1000):
    eemd = EEMD()
    eemd.FIXE_H = 5
    eemd.noise_seed(12345)
    eemd.noise_width = noise_std
    eemd.ensemble_size = num_realizations
    imfs = eemd.eemd(signal, time)
    return imfs

# Function to generate white noise IMFs for significance testing
def generate_white_noise_imfs(signal_length, time, num_realizations, use_eemd=None, noise_std=0.2):
    white_noise_imfs = []
    if use_eemd:
        eemd = EEMD()
        eemd.FIXE_H = 5
        eemd.noise_seed(12345)
        eemd.noise_width = noise_std
        eemd.ensemble_size = num_realizations
        for _ in range(num_realizations):
            white_noise = np.random.normal(size=signal_length)
            imfs = eemd.eemd(white_noise, time)
            white_noise_imfs.append(imfs)
    else:
        emd = EMD()
        emd.FIXE_H = 5
        for _ in range(num_realizations):
            white_noise = np.random.normal(size=signal_length)
            imfs = emd.emd(white_noise, time)
            white_noise_imfs.append(imfs)
    return white_noise_imfs

# Function to calculate the energy of each IMF
def calculate_imf_energy(imfs):
    return [np.sum(imf**2) for imf in imfs]

# Function to determine the significance of the IMFs
def test_imf_significance(imfs, white_noise_imfs):
    imf_energies = calculate_imf_energy(imfs)
    num_imfs = len(imfs)
    
    white_noise_energies = [calculate_imf_energy(imf_set) for imf_set in white_noise_imfs]
    
    significance_levels = []
    for i in range(num_imfs):
        # Collect the i-th IMF energy from each white noise realization
        white_noise_energy_dist = [imf_energies[i] for imf_energies in white_noise_energies if i < len(imf_energies)]
        mean_energy = np.mean(white_noise_energy_dist)
        std_energy = np.std(white_noise_energy_dist)
        
        z_score = (imf_energies[i] - mean_energy) / std_energy
        significance_level = 1 - norm.cdf(z_score)
        significance_levels.append(significance_level)
    
    return significance_levels

# Function to compute instantaneous frequency of IMFs using Hilbert transform
def compute_instantaneous_frequency(imfs, time):
    instantaneous_frequencies = []
    for imf in imfs:
        analytic_signal = hilbert(imf)
        instantaneous_phase = np.unwrap(np.angle(analytic_signal))
        instantaneous_frequency = np.diff(instantaneous_phase) / (2.0 * np.pi * np.diff(time))
        instantaneous_frequency = np.abs(instantaneous_frequency)  # Ensure frequencies are positive
        instantaneous_frequencies.append(instantaneous_frequency)
    return instantaneous_frequencies

# Function to compute HHT power spectrum
def compute_hht_power_spectrum(imfs, instantaneous_frequencies, freq_bins):
    power_spectrum = np.zeros_like(freq_bins)
    
    for i in range(len(imfs)):
        amplitude = np.abs(hilbert(imfs[i]))
        power = amplitude**2
        for j in range(len(instantaneous_frequencies[i])):
            freq_idx = np.searchsorted(freq_bins, instantaneous_frequencies[i][j])
            if freq_idx < len(freq_bins):
                power_spectrum[freq_idx] += power[j]
    
    return power_spectrum

# Function to compute significance level for HHT power spectrum
def compute_significance_level(white_noise_imfs, freq_bins, time):
    all_power_spectra = []
    for imfs in white_noise_imfs:
        instantaneous_frequencies = compute_instantaneous_frequency(imfs, time)
        power_spectrum = compute_hht_power_spectrum(imfs, instantaneous_frequencies, freq_bins)
        all_power_spectra.append(power_spectrum)
    
    # Compute the 95th percentile power spectrum as the significance level
    significance_level = np.percentile(all_power_spectra, 95, axis=0)
    
    return significance_level

## Custom rounding function
def custom_round(freq):
    if freq < 1:
        return round(freq, 1)
    else:
        return round(freq)

def smooth_power_spectrum(power_spectrum, window_size=5):
    return np.convolve(power_spectrum, np.ones(window_size)/window_size, mode='same')

# Function to identify and print significant peaks
def identify_significant_peaks(freq_bins, power_spectrum, significance_level):
    peaks, _ = find_peaks(power_spectrum, height=significance_level)
    significant_frequencies = freq_bins[peaks]
    rounded_frequencies = [custom_round(freq) for freq in significant_frequencies]
    rounded_frequencies_two = [round(freq, 1) for freq in significant_frequencies]
    print("Significant Frequencies (Hz):", rounded_frequencies_two)
    print("Significant Frequencies (Hz):", rounded_frequencies)
    return significant_frequencies

# Function to generate randomized signals
def generate_randomized_signals(signal, num_realizations):
    randomized_signals = []
    for _ in range(num_realizations):
        randomized_signal = np.random.permutation(signal)
        randomized_signals.append(randomized_signal)
    return randomized_signals

# Function to calculate the FFT power spectrum for each IMF
def calculate_fft_psd_spectra(imfs, time):
    psd_spectra = []
    for imf in imfs:
        N = len(imf)
        T = time[1] - time[0]  # Assuming uniform sampling
        yf = fft(imf)
        xf = fftfreq(N, T)[:N // 2]
        psd = (2.0 / N) * (np.abs(yf[:N // 2]) ** 2) / (N * T)
        psd_spectra.append((xf, psd))
    return psd_spectra

# Function to calculate the FFT power spectrum for randomized signals
def calculate_fft_psd_spectra_randomized(imfs, time, num_realizations):
    all_psd_spectra = []
    for imf in imfs:
        randomized_signals = generate_randomized_signals(imf, num_realizations)
        for signal in randomized_signals:
            N = len(signal)
            T = time[1] - time[0]  # Assuming uniform sampling
            yf = fft(signal)
            psd = (2.0 / N) * (np.abs(yf[:N // 2]) ** 2) / (N * T)
            all_psd_spectra.append(psd)
    return np.array(all_psd_spectra)

# Function to calculate the 95th percentile confidence level
def calculate_confidence_level(imfs, time, num_realizations):
    confidence_levels = []
    for imf in imfs:
        all_psd_spectra = calculate_fft_psd_spectra_randomized([imf], time, num_realizations)
        confidence_level = np.percentile(all_psd_spectra, 95, axis=0)
        confidence_levels.append(confidence_level)
    return confidence_levels

# Function to calculate the Welch PSD for each IMF
def calculate_welch_psd_spectra(imfs, fs):
    psd_spectra = []
    for imf in imfs:
        f, psd = welch(imf, fs=fs)
        psd_spectra.append((f, psd))
    return psd_spectra

# Function to calculate the Welch PSD for randomized signals
def calculate_welch_psd_spectra_randomized(imfs, fs, num_realizations):
    all_psd_spectra = []
    for imf in imfs:
        randomized_signals = generate_randomized_signals(imf, num_realizations)
        for signal in randomized_signals:
            f, psd = welch(signal, fs=fs)
            all_psd_spectra.append(psd)
    return np.array(all_psd_spectra)

# Function to calculate the 95th percentile confidence level
def calculate_confidence_level_welch(imfs, fs, num_realizations):
    confidence_levels = []
    for imf in imfs:
        all_psd_spectra = calculate_welch_psd_spectra_randomized([imf], fs, num_realizations)
        confidence_level = np.percentile(all_psd_spectra, 95, axis=0)
        confidence_levels.append(confidence_level)
    return confidence_levels

# ************** Main EMD/EEMD routine **************
def getEMD_HHT(signal, time, **kwargs):
    """
    Calculate EMD/EEMD and HHT
    
    Parameters:
        signal (array): The input signal (1D).
        time (array): The time array of the signal.
        siglevel (float): Significance level for the confidence intervals. Default: 0.95.
        nperm (int): Number of permutations for significance testing. Default: 1000.
        EEMD (bool): If True, use Ensemble Empirical Mode Decomposition (EEMD) instead of Empirical Mode Decomposition (EMD). Default: False.
        Welch_psd (bool): If True, calculate Welch PSD spectra instead of FFT PSD spectra (for the psd_spectra and psd_confidence_levels). Default: False.
        apod (float): Extent of apodization edges (of a Tukey window). Default: 0.1.
        pxdetrend (int): Subtract linear trend with time per pixel. Options: 1 (simple) or 2 (advanced). Default: 2.
        polyfit (int): Degree of polynomial fit for detrending the data. If set, a polynomial fit (instead of linear) is applied. Default: None.
        meantemporal (bool): If True, apply simple temporal detrending by subtracting the mean signal from the data, skipping fitting procedures. Default: False.
        meandetrend (bool): If True, subtract the linear trend with time for the image means (spatial detrending). Default: False.
        recon (bool): If True, perform Fourier reconstruction of the input time series. This does not preserve amplitudes but is useful for examining frequencies far from the low-frequency range. Default: False.
        resample_original (bool): If True, and if recon is set to True, approximate values close to the original are returned for comparison. Default: False.
        nodetrendapod (bool): If True, neither detrending nor apodization is performed. Default: False.
        silent (bool): If True, suppress print statements. Default: False.
        **kwargs: Additional parameters for the analysis method.
    
    Returns:
        frequencies: Frequencies at which the PSD is estimated.
        psd: Power Spectral Density values.
        significance: Significance levels for the PSD.
        ‘marginal’ spectrum of Hilbert-Huang Transform (HHT)
    """
    # Define default values for the optional parameters similar to FFT
    defaults = {
        'siglevel': 0.95,
        'nperm': 1000,  # Number of permutations for significance calculation
        'apod': 0.1,  # Tukey window apodization
        'silent': False,
        'pxdetrend': 2,
        'meandetrend': False,
        'polyfit': None,
        'meantemporal': False,
        'recon': False,
        'resample_original': False,
        'EEMD': False,  # If True, calculate Ensemble Empirical Mode Decomposition (EEMD)
        'nodetrendapod': False,
        'Welch_psd': False,  # If True, calculate Welch PSD spectra
        'significant_imfs': False  # If True, return only significant IMFs (and for associated calculations)
    }

    # Update defaults with any user-provided keyword arguments
    params = {**defaults, **kwargs}

    tdiff = np.diff(time)
    cadence = np.median(tdiff)

    # Perform detrending and apodization
    if not params['nodetrendapod']:
        apocube = WaLSA_detrend_apod(
            signal, 
            apod=params['apod'], 
            meandetrend=params['meandetrend'], 
            pxdetrend=params['pxdetrend'], 
            polyfit=params['polyfit'], 
            meantemporal=params['meantemporal'], 
            recon=params['recon'], 
            cadence=cadence, 
            resample_original=params['resample_original'], 
            silent=params['silent']
        )

    n = len(apocube)
    dt = cadence
    fs = 1 / dt  # Sampling frequency

    if params['EEMD']:
        imfs = apply_eemd(signal, time)
        use_eemd = True
    else:
        imfs = apply_emd(signal, time)
        use_eemd = False

    white_noise_imfs = generate_white_noise_imfs(len(signal), time, params['nperm'], use_eemd=use_eemd)

    IMF_significance_levels = test_imf_significance(imfs, white_noise_imfs)

    # Determine significant IMFs
    significant_imfs_val = [imf for imf, sig_level in zip(imfs, IMF_significance_levels) if sig_level < params['siglevel']]
    if params['significant_imfs']:
        imfs = significant_imfs_val

    # Compute instantaneous frequencies of IMFs
    instantaneous_frequencies = compute_instantaneous_frequency(imfs, time)

    if params['Welch_psd']:
        # Calculate and plot Welch PSD spectra of (significant) IMFs
        psd_spectra = calculate_welch_psd_spectra(imfs, fs)
        # Calculate 95th percentile confidence levels for Welch PSD spectra
        psd_confidence_levels = calculate_confidence_level_welch(imfs, fs, params['nperm'])
    else:
        # Calculate and plot FFT PSD spectra of (significant) IMFs
        psd_spectra = calculate_fft_psd_spectra(imfs, time)
        # Calculate 95th percentile confidence levels for FFT PSD spectra
        psd_confidence_levels = calculate_confidence_level(imfs, time, params['nperm'])

    # Define frequency bins for HHT power spectrum
    max_freq = max([freq.max() for freq in instantaneous_frequencies])
    HHT_freq_bins = np.linspace(0, max_freq, 100)

    # Compute HHT power spectrum of (significant) IMFs
    HHT_power_spectrum = compute_hht_power_spectrum(imfs, instantaneous_frequencies, HHT_freq_bins)
    # Compute significance level for HHT power spectrum
    HHT_significance_level = compute_significance_level(white_noise_imfs, HHT_freq_bins, time)    

    if not params['silent']:
        if params['EEMD']:
            print("EEMD processed.")
        else:
            print("EMD processed.")

    return HHT_power_spectrum, HHT_significance_level, HHT_freq_bins, psd_spectra, psd_confidence_levels, imfs, IMF_significance_levels, instantaneous_frequencies


# ----------------------- Dominant Frequency & Averaged Power -------------------------------
def get_dominant_averaged(cube, time, **kwargs):
    """
    Analyze a 3D data cube to compute the dominant frequency and averaged power.

    Parameters:
        cube (3D array): Input data cube, expected in either 'txy' or 'xyt' format.
        time (array): Time array of the data cube.
        method (str): Analysis method ('fft' or 'wavelet').
        format (str): Format of the data cube ('txy' or 'xyt'). Default is 'txy'.
        **kwargs: Additional parameters specific to the analysis method.

    Returns:
        dominantfreq (float): Dominant frequency of the data cube.
        averagedpower (float): Averaged power of the data cube.
    """
    defaults = {
        'method': 'fft',
        'format': 'txy',
        'silent': False,
        'GWS': False,  # If True, calculate global wavelet spectrum
        'RGWS': True,  # If True, calculate refined global wavelet spectrum
    }

    # Update defaults with any user-provided keyword arguments
    params = {**defaults, **kwargs}

    # Check and adjust format if necessary
    if params['format'] == 'xyt':
        cube = np.transpose(cube, (2, 0, 1))  # Convert 'xyt' to 'txy'
    elif params['format'] != 'txy':
        raise ValueError("Unsupported format. Choose 'txy' or 'xyt'.")

    # Initialize arrays for storing results across spatial coordinates
    nt, nx, ny = cube.shape
    dominantfreq = np.zeros((nx, ny))

    if params['method'] == 'fft':
        method_name = 'FFT'
    elif params['method'] == 'lombscargle':
        method_name = 'Lomb-Scargle'
    elif params['method'] == 'wavelet':
        method_name = 'Wavelet'
    elif params['method'] == 'welch':
        method_name = 'Welch'

    if not params['silent']:
            print(f"Processing {method_name} for a 3D cube with format '{params['format']}' and shape {cube.shape}.")
            print(f"Calculating Dominant frequencies and/or averaged power spectrum ({method_name}) ....")

    # Iterate over spatial coordinates and apply the chosen analysis method
    if params['method'] == 'fft':
        for x in tqdm(range(nx), desc="Processing x", leave=True):
            for y in range(ny):
                signal = cube[:, x, y]
                fft_power, fft_freqs, _, _ = getpowerFFT(signal, time, silent=True, nosignificance=True, **kwargs)
                if x == 0 and y == 0:
                    powermap = np.zeros((nx, ny, len(fft_freqs)))
                powermap[x, y, :] = fft_power
                # Determine the dominant frequency for this pixel
                dominantfreq[x, y] = fft_freqs[np.argmax(fft_power)]

        # Calculate the averaged power over all pixels
        averagedpower = np.mean(powermap, axis=(0, 1))
        frequencies = fft_freqs
        print("\nAnalysis completed.")

    elif params['method'] == 'lombscargle':
        for x in tqdm(range(nx), desc="Processing x", leave=True):
            for y in range(ny):
                signal = cube[:, x, y]
                ls_power, ls_freqs, _ = getpowerLS(signal, time, silent=True, **kwargs)
                if x == 0 and y == 0:
                    powermap = np.zeros((nx, ny, len(ls_freqs)))
                powermap[x, y, :] = ls_power
                # Determine the dominant frequency for this pixel
                dominantfreq[x, y] = ls_freqs[np.argmax(ls_power)]

        # Calculate the averaged power over all pixels
        averagedpower = np.mean(powermap, axis=(0, 1))
        frequencies = ls_freqs
        print("\nAnalysis completed.")

    elif params['method'] == 'wavelet':
        for x in tqdm(range(nx), desc="Processing x", leave=True):
            for y in range(ny):
                signal = cube[:, x, y]
                if params['GWS']:
                    _, wavelet_periods, _, _, wavelet_power, _, _ = getpowerWavelet(signal, time, silent=True, **kwargs)
                elif params['RGWS']:
                    _, wavelet_periods, _, _, _, _, wavelet_power = getpowerWavelet(signal, time, silent=True, **kwargs)
                if x == 0 and y == 0:
                    powermap = np.zeros((nx, ny, len(wavelet_periods)))
                wavelet_freq = 1. / wavelet_periods
                powermap[x, y, :] = wavelet_power
                # Determine the dominant frequency for this pixel
                dominantfreq[x, y] = wavelet_freq[np.argmax(wavelet_power)]

        # Calculate the averaged power over all pixels
        averagedpower = np.mean(powermap, axis=(0, 1))
        frequencies = wavelet_freq
        print("\nAnalysis completed.")

    elif params['method'] == 'welch':
        for x in tqdm(range(nx), desc="Processing x", leave=True):
            for y in range(ny):
                signal = cube[:, x, y]
                welch_power, welch_freqs, _ = welch_psd(signal, time, silent=True, **kwargs)
                if x == 0 and y == 0:
                    powermap = np.zeros((nx, ny, len(welch_freqs)))
                powermap[x, y, :] = welch_power
                # Determine the dominant frequency for this pixel
                dominantfreq[x, y] = welch_freqs[np.argmax(welch_power)]

        # Calculate the averaged power over all pixels
        averagedpower = np.mean(powermap, axis=(0, 1))
        frequencies = welch_freqs
        print("\nAnalysis completed.")

    else:
        raise ValueError("Unsupported method. Choose 'fft', 'lombscargle', 'wavelet', or 'welch'.")
    
    return dominantfreq, averagedpower, frequencies, powermap
