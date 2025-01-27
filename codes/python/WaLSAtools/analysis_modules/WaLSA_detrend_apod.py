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


import numpy as np
from scipy.optimize import curve_fit
from .WaLSA_wavelet import cwt, Morlet

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
def WaLSA_detrend_apod(cube, apod=0.1, meandetrend=False, pxdetrend=2, polyfit=None, meantemporal=False,
                       recon=False, cadence=None, resample_original=False, min_resample=None, 
                       max_resample=None, silent=False, dj=32, lo_cutoff=None, hi_cutoff=None, upper=False):
    
    nt = len(cube)  # Assume input is a 1D signal
    apocube = np.copy(cube)  # Create a copy of the input signal
    t = np.arange(nt)  # Time array

    # Apply Tukey window (apodization)
    if apod > 0:
        tukey_window = custom_tukey(nt, apod)
        apocube = apocube * tukey_window  # Apodize the signal
    
    # Mean detrend (optional)
    # DC this looks like removing the mean from the signal, why is there a fit?
    if meandetrend:
        avg_signal = np.mean(apocube)
        # time = np.arange(nt)
        # mean_fit_params, _ = curve_fit(linear, time, avg_signal)
        # mean_trend = linear(time, *mean_fit_params)
        apocube -= avg_signal#mean_trend
    
    # Pixel-based detrending (temporal detrend)
    if pxdetrend > 0:
        mean_val = np.mean(apocube)
        if meantemporal:
            # Simple temporal detrend by subtracting the mean
            apocube -= mean_val
        else:
            # Advanced detrend (linear or polynomial fit)
            if polyfit is not None:
                poly_coeffs = np.polyfit(t, apocube, polyfit)
                trend = np.polyval(poly_coeffs, t)
            else:
                popt, _ = curve_fit(linear, t, apocube, p0=[mean_val, 0])
                trend = linear(t, *popt)
            apocube -= trend

    # Wavelet-based Fourier reconstruction (optional)
    if recon and cadence:
        apocube = WaLSA_wave_recon(apocube, cadence, dj=dj, lo_cutoff=lo_cutoff, hi_cutoff=hi_cutoff, upper=upper)
    
    # Resampling to preserve amplitudes (optional)
    if resample_original:
        if min_resample is None:
            min_resample = np.min(apocube)
        if max_resample is None:
            max_resample = np.max(apocube)
        apocube = np.interp(apocube, (np.min(apocube), np.max(apocube)), (min_resample, max_resample))

    if not silent:
        print("Detrending and apodization complete.")
    
    return apocube

# Wavelet-based reconstruction function (optional)
def WaLSA_wave_recon(signal, cadence, dj=32, lo_cutoff=None, hi_cutoff=None, upper=False):
    mother = Morlet(6)
    dt = cadence
    n = len(signal)
    
    # Perform continuous wavelet transform (CWT) using the Morlet wavelet
    wave, scales, freqs, coi, fft, fftfreqs = cwt(signal, dt, dj, 2 * dt, -1, mother)
    period = 1 / freqs

    # Set default high cutoff if not provided
    if hi_cutoff is None:
        hi_cutoff = n * dt / (3. * np.sqrt(2))
    
    # Filter periods based on the cutoff values
    if upper:
        good_per = np.where(period > hi_cutoff)[0]
    else:
        good_per = np.where((period > lo_cutoff) & (period < hi_cutoff))[0]
    
    print(f"Filtering out frequencies below {1000. / hi_cutoff} mHz")
    
    # Filter out regions inside the cone of influence
    filtered_wave = np.zeros_like(wave)
    for i in range(n):
        filtered_wave[i, good_per] = np.real(wave[i, good_per])

    # Reconstruct the signal from filtered wavelet coefficients
    recon_signal = np.sum(filtered_wave / np.sqrt(scales), axis=1)
    recon_signal *= dj * np.sqrt(dt) / (0.766 * (np.pi ** -0.25))
    
    return recon_signal
