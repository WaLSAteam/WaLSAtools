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
from scipy.signal import coherence, csd # type: ignore
from .WaLSA_speclizer import WaLSA_speclizer # type: ignore
from .WaLSA_wavelet import xwt, wct # type: ignore
from .WaLSA_detrend_apod import WaLSA_detrend_apod # type: ignore
from .WaLSA_wavelet_confidence import WaLSA_wavelet_confidence # type: ignore

# -------------------------------------- Main Function ----------------------------------------
def WaLSA_cross_spectra(signal=None, time=None, method=None, **kwargs):
    """
    Compute the cross-spectra between two time series.
    
    Parameters:
        data1 (array): The first time series (1D).
        data2 (array): The first time series (1D).
        time (array): The time array of the signals.
        method (str): The method to use for the analysis. Options are 'welch' and 'wavelet'.
    """
    if method == 'welch':
        return getcross_spectrum_Welch(signal, time=time, **kwargs)
    elif method == 'wavelet':
        return getcross_spectrum_Wavelet(signal, time=time, **kwargs)
    else:
        raise ValueError(f"Unknown method: {method}")


# -------------------------------------- Welch ----------------------------------------
def getcross_spectrum_Welch(signal, time, **kwargs):
    """
    Calculate cross-spectral relationships of two time series whose amplitudes and powers
    are computed using the WaLSA_speclizer routine.
    The cross-spectrum is complex valued, thus its magnitude is returned as the
    co-spectrum. The phase lags between the two time series are are estimated 
    from the imaginary and real arguments of the complex cross spectrum.
    The coherence is calculated from the normalized square of the amplitude 
    of the complex cross-spectrum
    
    Parameters:
        data1 (array): The first 1D time series signal.
        data2 (array): The second 1D time series signal.
        time (array): The time array corresponding to the signals.
        nperseg (int, optional): Length of each segment for analysis. Default: 256.
        noverlap (int, optional): Number of points to overlap between segments. Default: 128.
        window (str, optional): Type of window function used in the Welch method. Default: 'hann'.
        silent (bool, optional): If True, suppress print statements. Default: False.
        **kwargs: Additional parameters for the analysis method.
    
    Returns:
        frequencies, cospectrum, phase_angle, power_data1, power_data2, freq_coh, coh
    """
    # Define default values for the optional parameters
    defaults = {
        'data1': None,        # First data array
        'data1': None,        # Second data array
        'nperseg': 256,      # Number of points per segments to average
        'noverlap': None,    # Number of points to overlap between segments
        'nfft': None,        # Number of points to use in the FFT
        'window': 'hann',  # Window type for Welch method
        # 'nodetrendapod': None,  # No detrending or apodization applied
        # 'pxdetrend': 2,  # Detrend parameter
        # 'meandetrend': None,  # Detrend parameter
        # 'polyfit': None,    # Detrend parameter
        # 'meantemporal': None,  # Detrend parameter
        # 'recon': None      # Detrend parameter
    }

    # Update defaults with any user-provided keyword arguments
    params = {**defaults, **kwargs}

    data1 = params['data1']
    data2 = params['data2']
    # DC in this way you can give params as input for the welch and csd function
    del params['data1']
    del params['data2']
    

    dummy = signal+2

    tdiff = np.diff(time)
    cadence = np.median(tdiff)
    
    # Power spectrum for data1
    power_data1, frequencies, _ = WaLSA_speclizer(signal=data1, time=time, method='welch', 
                                               amplitude=True, nosignificance=True, silent=kwargs.pop('silent', True), **params)

    # Power spectrum for data2
    power_data2, _, _ = WaLSA_speclizer(signal=data2, time=time, method='welch', 
                                               amplitude=True, nosignificance=True, silent=kwargs.pop('silent', False), **params)

    # Calculate cross-spectrum
    _, crosspower = csd(data1, data2, fs=1.0/cadence, **params)

    cospectrum = np.abs(crosspower)
    
    # Calculate phase lag
    phase_angle = np.angle(crosspower, deg=True)

    # Calculate coherence
    freq_coh, coh = coherence(data1, data2, 1.0/cadence, window=params['window'], 
                              nperseg=params['nperseg'], noverlap=params['noverlap'], nfft=params['nfft'])
    
    return frequencies, cospectrum, phase_angle, power_data1, power_data2, freq_coh, coh


# -------------------------------------- Wavelet ----------------------------------------
def getcross_spectrum_Wavelet(signal, time, **kwargs):
    """
    Calculate the cross-power spectrum, coherence spectrum, and phase angles 
    between two signals using Wavelet Transform.

    Parameters:
        data1 (array): The first 1D time series signal.
        data2 (array): The second 1D time series signal.
        time (array): The time array corresponding to the signals.
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

    Returns:
    cross_power : np.ndarray
        Cross power spectrum between `data1` and `data2`.
        
    cross_periods : np.ndarray
        Periods corresponding to the cross-power spectrum.
        
    cross_sig : np.ndarray
        Significance of the cross-power spectrum.
        
    cross_coi : np.ndarray
        Cone of influence for the cross-power spectrum.
        
    coherence : np.ndarray
        Coherence spectrum between `data1` and `data2`.
        
    coh_periods : np.ndarray
        Periods corresponding to the coherence spectrum.
        
    coh_sig : np.ndarray
        Significance of the coherence spectrum.
        
    corr_coi : np.ndarray
        Cone of influence for the coherence spectrum.
        
    phase_angle : np.ndarray
        2D array containing x- and y-components of the phase direction 
        arrows for each frequency and time point.
        
    Notes
    -----
    - Cross-power, coherence, and phase angles are calculated using 
      **cross-wavelet transform (XWT)** and **wavelet coherence transform (WCT)**.
    - Arrows for phase direction are computed such that:
        - Arrows pointing downwards indicate anti-phase.
        - Arrows pointing upwards indicate in-phase.
        - Arrows pointing right indicate `data1` leading `data2`.
        - Arrows pointing left indicate `data2` leading `data1`.

    Examples
    --------
    >>> cross_power, cross_periods, cross_sig, cross_coi, coherence, coh_periods, coh_sig, corr_coi, phase_angle, dt = \
    >>> getcross_spectrum_Wavelet(signal, cadence, data1=signal1, data2=signal2)
    """

    # Define default values for optional parameters
    defaults = {
        'data1': None,               # First data array
        'data2': None,                # Second data array
        'mother': 'morlet',          # Type of mother wavelet
        'siglevel': 0.95,  # Significance level for cross-spectral analysis
        'dj': 0.025,                 # Spacing between scales
        's0': -1,  # Initial scale
        'J': -1,  # Number of scales
        'cache': False,               # Cache results to avoid recomputation
        'apod': 0.1,
        'silent': False,
        'pxdetrend': 2,
        'meandetrend': False,
        'polyfit': None,
        'meantemporal': False,
        'recon': False,
        'resample_original': False,
        'nperm': 1000               # Number of permutations for significance testing
    }

    # Update defaults with any user-provided keyword arguments
    params = {**defaults, **kwargs}

    data1 = params['data1']
    data2 = params['data2']
    dummy = signal+2

    tdiff = np.diff(time)
    cadence = np.median(tdiff)

    data1orig = data1.copy()
    data2orig = data2.copy()

    data1 = WaLSA_detrend_apod(
        data1, 
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

    data2 = WaLSA_detrend_apod(
        data2, 
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

    # Calculate signal properties
    nt = len(params['data1'])  # Number of time points

    # Determine the initial scale s0 if not provided
    if params['s0'] == -1:
        params['s0'] = 2 * cadence

    # Determine the number of scales J if not provided
    if params['J'] == -1:
        params['J'] = int((np.log(float(nt) * cadence / params['s0']) / np.log(2)) / params['dj'])

    # ----------- CROSS-WAVELET TRANSFORM (XWT) -----------
    W12, cross_coi, freq, signif = xwt(
        data1, data2, cadence, dj=params['dj'], s0=params['s0'], J=params['J'], 
        no_default_signif=True, wavelet=params['mother'], normalize=False
    )
    print('Wavelet cross-power spectrum calculated.')
    # Calculate the cross-power spectrum
    cross_power = np.abs(W12) ** 2
    cross_periods = 1 / freq  # Periods corresponding to the frequency axis
    
    #----------------------------------------------------------------------
    # Calculate significance levels using Monte Carlo randomization method
    #----------------------------------------------------------------------
    nxx, nyy = cross_power.shape
    cross_perm = np.zeros((nxx, nyy, params['nperm']))
    print('\nCalculating wavelet cross-power significance:')
    total_iterations = params['nperm']  # Total number of (x, y) combinations
    iteration_count = 0  # To keep track of progress
    for ip in range(params['nperm']):
        iteration_count += 1
        progress = (iteration_count / total_iterations) * 100
        print(f"\rProgress: {progress:.2f}%", end="")

        y_perm1 = np.random.permutation(data1orig)  # Permuting the original signal
        apocube1 = WaLSA_detrend_apod(
            y_perm1, 
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
        y_perm2 = np.random.permutation(data2orig)  # Permuting the original signal
        apocube2 = WaLSA_detrend_apod(
            y_perm2, 
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
        W12s, _, _, _ = xwt(
            apocube1, apocube2, cadence, dj=params['dj'], s0=params['s0'], J=params['J'], 
            no_default_signif=True, wavelet=params['mother'], normalize=False
        )
        cross_power_sig = np.abs(W12s) ** 2
        cross_perm[:, :, ip] = cross_power_sig

    signifn = WaLSA_wavelet_confidence(cross_perm, siglevel=params['siglevel'])

    cross_sig = cross_power / signifn

    # ----------- WAVELET COHERENCE TRANSFORM (WCT) -----------
    WCT, aWCT, coh_coi, freq, sig = wct(
        data1, data2, cadence, dj=params['dj'], s0=params['s0'], J=params['J'], 
        sig=False, wavelet=params['mother'], normalize=False, cache=params['cache']
    )
    print('Wavelet coherence calculated.')
    # Calculate the coherence spectrum
    coh_periods = 1 / freq  # Periods corresponding to the frequency axis
    coherence = WCT

    #----------------------------------------------------------------------
    # Calculate significance levels using Monte Carlo randomization method
    #----------------------------------------------------------------------
    nxx, nyy = coherence.shape
    coh_perm = np.zeros((nxx, nyy, params['nperm']))
    print('\nCalculating wavelet coherence significance:')
    total_iterations = params['nperm']  # Total number of permutations
    iteration_count = 0  # To keep track of progress
    for ip in range(params['nperm']):
        iteration_count += 1
        progress = (iteration_count / total_iterations) * 100
        print(f"\rProgress: {progress:.2f}%", end="")
        
        y_perm1 = np.random.permutation(data1orig)  # Permuting the original signal
        apocube1 = WaLSA_detrend_apod(
            y_perm1, 
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
        y_perm2 = np.random.permutation(data2orig)  # Permuting the original signal
        apocube2 = WaLSA_detrend_apod(
            y_perm2, 
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
        WCTs, _, _, _, _ = wct(
            apocube1, apocube2, cadence, dj=params['dj'], s0=params['s0'], J=params['J'], 
            sig=False, wavelet=params['mother'], normalize=False, cache=params['cache']
        )
        coh_perm[:, :, ip] = WCTs

    sig_coh = WaLSA_wavelet_confidence(coh_perm, siglevel=params['siglevel'])

    coh_sig = coherence / sig_coh  # Ratio > 1 means coherence is significant

    # --------------- PHASE ANGLES ---------------
    phase_angle = aWCT

    return (
        cross_power, cross_periods, cross_sig, cross_coi, 
        coherence, coh_periods, coh_sig, coh_coi, 
        phase_angle
    )
