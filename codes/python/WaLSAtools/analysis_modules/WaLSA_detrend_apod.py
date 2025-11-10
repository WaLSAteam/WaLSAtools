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
# Jafarzadeh, S., Jess, D. B., Stangalini, M. et al. 2025, Nature Reviews Methods Primers, 5, 21
# -----------------------------------------------------------------------------------------------------

from ast import If
import numpy as np # type: ignore
from scipy.optimize import curve_fit # type: ignore
from .WaLSA_wavelet import cwt, significance # type: ignore

# -------------------------------------- Wavelet ----------------------------------------
def getWavelet(signal, time, **kwargs):
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
        'dj': 1/32. ,  # Scale spacing
        's0': -1,  # Initial scale
        'J': -1,  # Number of scales
        'lag1': 0.0,  # Lag-1 autocorrelation
        'silent': False,
        'nperm': 1000,   # Number of permutations for significance calculation
    }

    # Update defaults with any user-provided keyword arguments
    params = {**defaults, **kwargs}

    tdiff = np.diff(time)
    cadence = np.median(tdiff)

    n = len(signal)
    dt = cadence

    # Standardize the signal before the wavelet transform
    std_signal = signal.std()
    norm_signal = signal / std_signal

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

    return power, periods, coi, scales, W

# -----------------------------------------------------------------------------------------------------

# Linear detrending function for curve fitting
def linear(x, a, b):
    return a + b * x

# Custom Tukey window implementation
def custom_tukey(nt, apod=0.1):
    apodrim = int(apod * nt)
    apodt = np.ones(nt)  # Initialize with ones
    
    # Apply sine-squared taper to the first 'apodrim' points
    taper = (np.sin(np.pi / 2. * np.arange(apodrim) / apodrim)) ** 2
    
    # Apply taper symmetrically at both ends
    apodt[:apodrim] = taper
    apodt[-apodrim:] = taper[::-1]  # Reverse taper for the last points
    
    return apodt

# Main detrending and apodization function
# apod=0: The Tukey window becomes a rectangular window (no tapering).
# apod=1: The Tukey window becomes a Hann window (fully tapered with a cosine function).
def WaLSA_detrend_apod(
    cube,
    apod=0.1,
    meandetrend=False,
    pxdetrend=2, # removed!
    polyfit=None,
    meantemporal=False,
    recon=False,
    cadence=None,
    resample_original=False,
    min_resample=None,
    max_resample=None,
    silent=False,
    dj=32,
    lo_cutoff=None,
    hi_cutoff=None,
    upper=False,
):
    """
    1D preprocessing pipeline:
    Detrend (mean / linear / polynomial) -> (optional) wavelet recon/filter ->
    (optional) rescale -> Apodize (Tukey) last.
    """
    
    # Ensure float array and local working copy
    apocube = np.asarray(cube, dtype=float).copy()
    nt = apocube.size

    # Time axis for detrending: use physical seconds if cadence provided
    if cadence is not None:
        t = np.arange(nt, dtype=float) * float(cadence)
    else:
        t = np.arange(nt, dtype=float)

    # ---------------- DETREND ----------------
    if meantemporal or meandetrend:
        # DC removal
        apocube -= apocube.mean()
    if polyfit is not None:
        # Polynomial detrend of given degree
        deg = int(polyfit)
        coeffs = np.polyfit(t, apocube, deg)
        trend = np.polyval(coeffs, t)
    else:
        # Linear detrend (slope + intercept)
        slope, intercept = np.polyfit(t, apocube, 1)
        trend = slope * t + intercept
        # mean_val = np.mean(apocube)
        # popt, _ = curve_fit(linear, t, apocube, p0=[mean_val, 0])
        # trend = linear(t, *popt)
    apocube -= trend

    # -------- OPTIONAL wavelet reconstruction / filtering --------
    if recon and cadence is not None:
        apocube = WaLSA_wave_recon(
            apocube, float(cadence), dj=dj,
            lo_cutoff=lo_cutoff, hi_cutoff=hi_cutoff, upper=upper
        )

    # -------- OPTIONAL rescale to preserve amplitudes --------
    if resample_original:
        if min_resample is None:
            min_resample = float(np.min(apocube))
        if max_resample is None:
            max_resample = float(np.max(apocube))
        apocube = np.interp(apocube, (np.min(apocube), np.max(apocube)), (min_resample, max_resample))

    # ---------------- APODIZATION (Tukey) ----------------
    if apod > 0:
        tukey_window = custom_tukey(nt, apod)
        apocube *= tukey_window

    if not silent:
        print("Detrending and apodization complete.")

    return apocube

# Wavelet-based reconstruction function (optional)
def WaLSA_wave_recon(ts, delt, dj=32, lo_cutoff=None, hi_cutoff=None, upper=False):
    """
    Reconstructs the wavelet-filtered time series based on given frequency cutoffs.
    """

    # Define duration based on the time series length
    dur = (len(ts) - 1) * delt

    # Assign default values if lo_cutoff or hi_cutoff is None
    if lo_cutoff is None:
        lo_cutoff = 0.0  # Default to 0 as in IDL

    if hi_cutoff is None:
        hi_cutoff = dur / (3.0 * np.sqrt(2)) 

    mother = 'morlet'
    num_points = len(ts)
    time_array = np.linspace(0, (num_points - 1) * delt, num_points)

    _, period, coi, scales, wave= getWavelet(signal=ts, time=time_array, method='wavelet', siglevel=0.99, apod=0.1, mother=mother)

    # Ensure good_per_idx and bad_per_idx are properly assigned
    if upper:
        good_per_idx = np.where(period > hi_cutoff)[0][0]
        bad_per_idx = len(period)
    else:
        good_per_idx = np.where(period > lo_cutoff)[0][0]
        bad_per_idx = np.where(period > hi_cutoff)[0][0]

    # set the power inside the CoI equal to zero 
	# (i.e., exclude points inside the CoI -- subject to edge effect)
    iampl = np.zeros((len(ts), len(period)), dtype=float)
    for i in range(len(ts)):
        pcol = np.real(wave[:, i])  # Extract real part of wavelet transform for this time index
        ii = np.where(period < coi[i])[0]  # Find indices where period is less than COI at time i
        if ii.size > 0:
            iampl[i, ii] = pcol[ii]  # Assign values where condition is met

    # Initialize reconstructed signal array
    recon_sum = np.zeros(len(ts))
    # Summation over valid period indices
    for i in range(good_per_idx, bad_per_idx):
        recon_sum += iampl[:, i] / np.sqrt(scales[i])

    # Apply normalization factor
    recon_all = dj * np.sqrt(delt) * recon_sum / (0.766 * (np.pi ** -0.25))

    return recon_all
