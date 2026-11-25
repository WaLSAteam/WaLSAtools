# --------------------------------------------------------------------------------------------------------------
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
# --------------------------------------------------------------------------------------------------------------
# The following codes are baed on those originally written by Rob Rutten, David B. Jess, and Samuel D. T. Grant
# --------------------------------------------------------------------------------------------------------------

import numpy as np # type: ignore
from scipy.optimize import curve_fit # type: ignore
from scipy.signal import convolve # type: ignore
# --------------------------------------------------------------------------------------------

def gaussian_function(sigma, width=None):
    """
    Create a Gaussian kernel that closely matches the IDL implementation.

    Parameters:
        sigma (float or list of floats): Standard deviation(s) of the Gaussian.
        width (int or list of ints, optional): Width of the Gaussian kernel.

    Returns:
        np.ndarray: The Gaussian kernel.
    """

    # sigma = np.array(sigma, dtype=np.float64)
    sigma = np.atleast_1d(np.array(sigma, dtype=np.float64))  # Ensure sigma is at least 1D
    if np.any(sigma <= 0):
        raise ValueError("Sigma must be greater than 0.")
    
    # width = np.array(width)
    width = np.atleast_1d(np.array(width))  # Ensure width is always at least 1D
    width = np.maximum(width.astype(np.int64), 1)  # Convert to integers and ensure > 0
    if np.any(width <= 0):
        raise ValueError("Width must be greater than 0.")
    
    nSigma = np.size(sigma)
    if nSigma > 8:
        raise ValueError('Sigma can have no more than 8 elements')
    nWidth = np.size(width)
    if nWidth > nSigma:
        raise ValueError('Incorrect width specification')
    if (nWidth == 1) and (nSigma > 1):
        width = np.full(nSigma, width[0])

    # kernel = np.zeros(tuple(width.astype(int)), dtype=np.float64)
    kernel = np.zeros(tuple(width[:nSigma].astype(int)), dtype=np.float64)  # Match kernel size to nSigma

    temp = np.zeros(8, dtype=np.int64)
    temp[:nSigma] = width[:nSigma]
    width = temp

    if nSigma == 2:
        b = a = 0
        if nSigma >= 1:
            a = np.linspace(0, width[0] - 1, width[0]) - width[0] // 2 + (0 if width[0] % 2 else 0.5)
        if nSigma >= 2:
            b = np.linspace(0, width[1] - 1, width[1]) - width[1] // 2 + (0 if width[1] % 2 else 0.5)

        a1 = b1 = 0
        # Nested loop for kernel computation (to be completed for larger nSigma ....)
        for bb in range(width[1]): 
            b1 = (b[bb] ** 2) / (2 * sigma[1] ** 2) if nSigma >= 2 else 0
            for aa in range(width[0]): 
                a1 = (a[aa] ** 2) / (2 * sigma[0] ** 2) if nSigma >= 1 else 0
                kernel[aa, bb] = np.exp(
                    -np.clip((a1 + b1), -1e3, 1e3)
                )
    elif nSigma == 1:
        a = 0
        if nSigma >= 1:
            a = np.linspace(0, width[0] - 1, width[0]) - width[0] // 2 + (0 if width[0] % 2 else 0.5)

        a1 = 0

        # Nested loop for kernel computation
        for aa in range(width[0]): 
            a1 = (a[aa] ** 2) / (2 * sigma[0] ** 2) if nSigma >= 1 else 0
            kernel[aa] = np.exp(
                -np.clip(a1, -1e3, 1e3)
            )

    kernel = np.nan_to_num(kernel, nan=0.0, posinf=0.0, neginf=0.0)

    if np.sum(kernel) == 0:
        raise ValueError("Generated Gaussian kernel is invalid (all zeros).")
    
    return kernel


def walsa_radial_distances(shape):
    """
    Compute the radial distance array for a given shape.

    Parameters:
        shape (tuple): Shape of the array, typically (ny, nx).

    Returns:
        numpy.ndarray: Array of radial distances.
    """
    if not (isinstance(shape, tuple) and len(shape) == 2):
        raise ValueError("Shape must be a tuple with two elements, e.g., (ny, nx).")

    y, x = np.indices(shape)
    cy, cx = (np.array(shape) - 1) / 2
    return np.sqrt((x - cx) ** 2 + (y - cy) ** 2)

def avgstd(array):
    """
    Calculate the average and standard deviation of an array.
    """
    avrg = np.sum(array) / array.size
    stdev = np.sqrt(np.sum((array - avrg) ** 2) / array.size)
    return avrg, stdev

def linear(x, *p):
    """
    Compute a linear model y = p[0] + x * p[1].
    (used in temporal detrending )
    """
    if len(p) < 2:
        raise ValueError("Parameters p[0] and p[1] must be provided.")
    ymod = p[0] + x * p[1]
    return ymod


def gradient(xy, *p):
    """
    Gradient function for fitting spatial trends.
    Parameters:
        xy: Tuple of grid coordinates (x, y).
        p: Coefficients [offset, slope_x, slope_y].
    """
    x, y = xy  # Unpack the tuple
    if len(p) < 3:
        raise ValueError("Parameters p[0], p[1], and p[2] must be provided.")
    return p[0] + x * p[1] + y * p[2]


def apod3dcube(cube, apod):
    """
    Apodizes a 3D cube in all three coordinates, with detrending.

    Parameters:
        cube : Input 3D data cube with dimensions (nx, ny, nt).
        apod (float): Apodization factor (0 means no apodization).

    Returns:
        Apodized 3D cube.
    """
    # Get cube dimensions
    nt, nx, ny = cube.shape
    apocube = np.zeros_like(cube, dtype=np.float32)

    # Define temporal apodization
    apodt = np.ones(nt, dtype=np.float32)
    if apod != 0:
        apodrimt = nt * apod
        apodrimt = int(apodrimt)  # Ensure apodrimt is an integer
        apodt[:apodrimt] = (np.sin(np.pi / 2. * np.arange(apodrimt) / apodrimt)) ** 2
        apodt = apodt * np.roll(np.flip(apodt), 1)  
        # Apply symmetrical apodization

    # Temporal detrending (mean-image trend, not per pixel)
    ttrend = np.zeros(nt, dtype=np.float32)
    tf = np.arange(nt) + 1.0
    for it in range(nt):
        img = cube[it, :, :]
        # ttrend[it], _ = avgstd(img)
        ttrend[it] = np.mean(img)

    # Fit the trend with a linear model
    fitp, _ = curve_fit(linear, tf, ttrend, p0=[1000.0, 0.0])
    # fit = fitp[0] + tf * fitp[1]
    fit = linear(tf, *fitp)

    # Temporal apodization per (x, y) column
    for it in range(nt):
        img = cube[it, :, :]
        apocube[it, :, :] = (img - fit[it]) * apodt[it]

    # Define spatial apodization
    apodx = np.ones(nx, dtype=np.float32)
    apody = np.ones(ny, dtype=np.float32)
    if apod != 0:
        apodrimx = apod * nx
        apodrimy = apod * ny
        apodrimx = int(apodrimx)  # Ensure apodrimx is an integer
        apodrimy = int(apodrimy)  # Ensure apodrimy is an integer
        apodx[:apodrimx] = (np.sin(np.pi / 2. * np.arange(int(apodrimx)) / apodrimx)) ** 2
        apody[:apodrimy] = (np.sin(np.pi / 2. * np.arange(int(apodrimy)) / apodrimy)) ** 2
        apodx = apodx * np.roll(np.flip(apodx), 1)
        apody = apody * np.roll(np.flip(apody), 1)
        apodxy = np.outer(apodx, apody)
    else:
        apodxy = np.outer(apodx, apody)

    # The original fit used constant spatial predictors. Its least-squares
    # model is the image mean. Keep that model without a rank-deficient fit.
    for it in range(nt):
        img = apocube[it, :, :]
        avg = np.mean(img)
        apocube[it, :, :] = (img - avg) * apodxy + avg

    return apocube


def ko_dist(sx, sy, double=False):
    """
    Set up Pythagorean distance array from the origin.
    """
    # The spatial transform is unshifted. DC is at index (0, 0).
    kx = np.rint(np.fft.fftfreq(sx) * sx)
    ky = np.rint(np.fft.fftfreq(sy) * sy)
    # Preserve the original Nyquist-normalised ring spacing.
    afstanden = np.hypot(kx[:, None] / (sx / 2), ky[None, :] / (sy / 2))
    afstanden *= min(sx, sy) // 2 + 1

    # Convert to integers if 'double' is False
    if not double:
        afstanden = np.round(afstanden).astype(int)
    
    return afstanden


def averpower(cube):
    """
    Compute 2D (k_h, f) power array by circular averaging over k_x, k_y.

    Parameters:
    cube : Input 3D data cube of dimensions (nx, ny, nt).

    Returns:
    avpow : 2D array of average power over distances, dimensions (maxdist+1, nt/2+1).
    """
    # Get cube dimensions
    nt, nx, ny = cube.shape

    # Perform FFT in all three dimensions (first in time direction)
    fftcube = np.fft.fft(np.fft.fft(np.fft.fft(cube, axis=0)[:nt//2+1, :, :], axis=1), axis=2)

    # Set up distances
    afstanden = ko_dist(nx, ny)  # Integer-rounded Pythagoras array

    maxdist = min(nx, ny) // 2 + 1  # Largest quarter circle
    
    # Initialize average power array
    avpow = np.zeros((maxdist + 1, nt // 2 + 1), dtype=np.float64)

    # Compute average power over all k_h distances, building power(k_h, f)
    for i in range(maxdist + 1):
        where_indices = np.where(afstanden == i)
        for j in range(nt // 2 + 1):
            w1 = fftcube[j, :, :][where_indices]
            avpow[i, j] = np.mean(np.abs(w1) ** 2) if w1.size else 0.0

    return avpow


def walsa_kopower_funct(datacube, **kwargs):
    """
    Calculate k-omega diagram
    (Fourier power at temporal frequency f against horizontal spatial wavenumber k_h)

    Origonally written in IDL by Rob Rutten (RR) assembly of Alfred de Wijn's routines (2010)
    - Translated into Pythn by Shahin Jafarzadeh (2024)

    Parameters:
        datacube: Input data cube [t, x, y].
        arcsecpx (float): Spatial sampling in arcsec/pixel.
        cadence (float): Temporal sampling in seconds.
        apod: fractional extent of apodization edges; default 0.1
        kmax: maximum k_h axis as fraction of Nyquist value, default 0.2
        fmax: maximum f axis as fraction of Nyquist value, default 0.5
        minpower: minimum of power plot range, default maxpower-5
        maxpower: maximum of power plot range, default alog10(max(power))

    Returns:
        k-omega power map.
    """

    # Set default parameter values
    defaults = {
        'apod': 0.1,
        'kmax': 1.0,
        'fmax': 1.0,
        'minpower': None,
        'maxpower': None
    }

    # Update defaults with user-provided values
    params = {**defaults, **kwargs}

    # Apodize the cube
    apocube = apod3dcube(datacube, params['apod'])

    # Compute radially-averaged power
    avpow = averpower(apocube)

    return avpow

# -------------------------------------- Main Function ----------------------------------------
def WaLSA_k_omega(signal, time=None, **kwargs):
    """
    NAME: WaLSA_k_omega
        part of -- WaLSAtools --
    
    * Main function to calculate and plot k-omega diagram.

    ORIGINAL CODE: QUEEns Fourier Filtering (QUEEFF) code
    WRITTEN, ANNOTATED, TESTED AND UPDATED in IDL BY:
    (1) Dr. David B. Jess
    (2) Dr. Samuel D. T. Grant
    The original code along with its manual can be downloaded at: https://bit.ly/37mx9ic

    WaLSA_k_omega (in IDL): A lightly modified version of the original code (i.e., a few additional keywords added) by Dr. Shahin Jafarzadeh
    - Translated into Pythn by Shahin Jafarzadeh (2024)

    Parameters:
        signal (array): Input datacube, normally in the form of [x, y, t] or [t, x, y]. Note that the input datacube must have identical x and y dimensions. If not, the datacube will be cropped accordingly.
        time (array): Time array corresponding to the input datacube.
        pixelsize (float): Spatial sampling of the input datacube. If not given, it is plotted in units of 'pixel'.
        filtering (bool): If True, filtering is applied, and the filtered datacube (filtered_cube) is returned. Otherwise, None is returned. Default: False.
        f1 (float): Optional lower (temporal) frequency to filter, in Hz.
        f2 (float): Optional upper (temporal) frequency to filter, in Hz.
        k1 (float): Optional lower (spatial) wavenumber to filter, in units of pixelsize^-1 (k = (2 * π) / wavelength).
        k2 (float): Optional upper (spatial) wavenumber to filter, in units of pixelsize^-1.
        spatial_torus (bool): If True, makes the annulus used for spatial filtering have a Gaussian-shaped profile, useful for preventing aliasing. Default: True.
        temporal_torus (bool): If True, makes the temporal filter have a Gaussian-shaped profile, useful for preventing aliasing. Default: True.
        no_spatial_filt (bool): If True, ensures no spatial filtering is performed on the dataset (i.e., only temporal filtering is applied).
        no_temporal_filt (bool): If True, ensures no temporal filtering is performed on the dataset (i.e., only spatial filtering is applied).
        silent (bool): If True, suppresses the k-ω diagram plot.
        mode (int): Output power mode: 0 = log10(power) (default), 1 = linear power, 2 = sqrt(power) = amplitude.
        processing_maps (bool): If True, the function returns the processing maps (spatial_fft_map, torus_map, spatial_fft_filtered_map, temporal_fft, temporal_filter, temporal_frequencies, spatial_frequencies). Otherwise, they are all returned as None. Default: False.

    OUTPUTS:
        power : 2D array of power (see mode for the scale).
        frequencies : 1D array of frequencies (in mHz).
        wavenumber : 1D array of wavenumber (in pixelsize^-1).
        filtered_cube : 3D array of filtered datacube (if filtering is set).
        processing_maps (if set to True)

    IF YOU USE THIS CODE, THEN PLEASE CITE THE ORIGINAL PUBLICATION WHERE IT WAS USED:
    Jess et al. 2017, ApJ, 842, 59 (http://adsabs.harvard.edu/abs/2017ApJ...842...59J)
    """

    # Set default parameter values
    defaults = {
        'pixelsize': 1,
        'format': 'txy',
        'filtering': False,
        'f1': None,
        'f2': None,
        'k1': None,
        'k2': None,
        'spatial_torus': True,
        'temporal_torus': True,
        'no_spatial_filt': False,
        'no_temporal_filt': False,
        'silent': False,
        'xlog': False,
        'ylog': False,
        'xrange': None,
        'yrange': None,
        'nox2': False,
        'noy2': False,
        'mode': 0,
        'xtitle': 'Wavenumber',
        'xtitle_units': '(pixel⁻¹)',
        'ytitle': 'Frequency',
        'yttitle_units': '(Hz)',
        'x2ndaxistitle': 'Spatial size',
        'y2ndaxistitle': 'Period',
        'x2ndaxistitle_units': '(pixel)',
        'y2ndaxistitle_units': '(s)',
        'processing_maps': False
    }

    # Update defaults with user-provided values
    params = {**defaults, **kwargs}

    tdiff = np.diff(time)
    cadence = np.median(tdiff)

    pixelsize = params['pixelsize']
    filtered_cube = None
    spatial_fft_map = None
    torus_map = None
    spatial_fft_filtered_map = None
    temporal_fft = None
    temporal_frequencies = None
    spatial_frequencies = None
    temporal_filter = None

    # Check and adjust format if necessary
    if params['format'] == 'xyt':
        signal = np.transpose(signal, (2, 0, 1))  # Convert 'xyt' to 'txy'
    elif params['format'] != 'txy':
        raise ValueError("Unsupported format. Choose 'txy' or 'xyt'.")
    if not params['silent']:
            print(f"Processing k-ω analysis for a 3D cube with format '{params['format']}' and shape {signal.shape}.")

    # Input dimensions
    nt, nx, ny = signal.shape
    if nx != ny:
        min_dim = min(nx, ny)
        signal = signal[:,:min_dim, :min_dim]
    nt, nx, ny = signal.shape

    # Calculating the Nyquist frequencies
    spatial_nyquist = (2 * np.pi) / (pixelsize * 2)
    temporal_nyquist = 1 / (cadence * 2)

    print("")
    print(f"Input datacube size (t,x,y): {signal.shape}")
    print("")
    print("Spatially, the important values are:")
    print(f"    2-pixel size = {pixelsize * 2:.2f} {params['x2ndaxistitle_units']}")
    print(f"    Nyquist wavenumber = {spatial_nyquist:.2f} {params['xtitle_units']}")
    if params['no_spatial_filt']:
        print("*** NO SPATIAL FILTERING WILL BE PERFORMED ***")
    print("")
    print("Temporally, the important values are:")
    print(f"    2-element duration (Nyquist period) = {cadence * 2:.2f} {params['y2ndaxistitle_units']}")
    print(f"    Time series duration = {cadence * signal.shape[0]:.2f} {params['y2ndaxistitle_units']}")
    temporal_nyquist = 1 / (cadence * 2)
    print(f"    Nyquist frequency = {temporal_nyquist:.2f} {params['yttitle_units']}")
    if params['no_temporal_filt']:
        print("***NO TEMPORAL FILTERING WILL BE PERFORMED***")
    print("")

    # Generate k-omega power map
    print("Constructing a k-ω diagram of the input datacube..........")
    print("")
    # Make the k-omega diagram using the proven method of Rob Rutten
    kopower = walsa_kopower_funct(signal)

    # Scales
    xsize_kopower = kopower.shape[0]
    dxsize_kopower = spatial_nyquist / float(xsize_kopower - 1)
    kopower_xscale = np.arange(xsize_kopower) * dxsize_kopower  # radians per spatial unit

    ysize_kopower = kopower.shape[1]
    dysize_kopower = 1.0 / (nt * cadence)
    kopower_yscale = (np.arange(ysize_kopower) * dysize_kopower) # in Hz

    # Generate Gaussian Kernel
    Gaussian_kernel = gaussian_function(sigma=[0.65, 0.65], width=3)
    Gaussian_kernel_norm = np.nansum(Gaussian_kernel)  # Normalize kernel sum

    # Copy kopower to kopower_plot
    kopower_plot = kopower.copy()

    # Convolve kopower (ignoring zero-th element, starting from index 1)
    kopower_plot[:, 1:] = convolve(
        kopower[:, 1:], 
        Gaussian_kernel, 
        mode='same'
    ) / Gaussian_kernel_norm

    # Normalize to frequency resolution (in mHz)
    freq = kopower_yscale[1:]
    if freq[0] == 0:
        freq0 = freq[1]
    else:
        freq0 = freq[0]
    kopower_plot /= freq0

    # Apply logarithmic or square root transformation based on mode
    if params['mode'] == 0:  # Logarithmic scaling
        kopower_plot = np.log10(kopower_plot)
    elif params['mode'] == 2:  # Square root scaling
        kopower_plot = np.sqrt(kopower_plot)

    # Normalize the power - preferred for plotting
    komegamap = np.clip(
        kopower_plot[1:, 1:],
        np.nanmin(kopower_plot[1:, 1:]),
        np.nanmax(kopower_plot[1:, 1:])
    )

    kopower_zscale = kopower_plot[1:,1:]

    # Rotate kopower counterclockwise by 90 degrees
    komegamap = np.rot90(komegamap, k=1)
    # Flip vertically (y-axis)
    komegamap = np.flip(komegamap, axis=0) 

    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Filtering implementation
    if params['filtering']:
        # Extract parameters from the dictionary
        k1 = params['k1']
        k2 = params['k2']
        f1 = params['f1']
        f2 = params['f2']

        if params['no_spatial_filt']:
            k1 = kopower_xscale[1]  # Default lower wavenumber
            k2 = np.nanmax(kopower_xscale)  # Default upper wavenumber
        # Ensure k1 and k2 are within valid bounds
        if k1 is None or k1 <= 0.0:
            k1 = kopower_xscale[1]
        if k2 is None or k2 > np.nanmax(kopower_xscale):
            k2 = np.nanmax(kopower_xscale)

        if params['no_temporal_filt']:
            f1 = kopower_yscale[1]
            f2 = np.nanmax(kopower_yscale)
        # Ensure f1 and f2 are within valid bounds
        if f1 is None or f1 <= 0.0:
            f1 = kopower_yscale[1] 
        if f2 is None or f2 > np.nanmax(kopower_yscale):
            f2 = np.nanmax(kopower_yscale) 

        print("Start filtering (in k-ω space) ......")
        print("")
        print(f"The preserved wavenumbers are [{k1:.3f}, {k2:.3f}] {params['xtitle_units']}")
        print(f"The preserved spatial sizes are [{(2 * np.pi) / k2:.3f}, {(2 * np.pi) / k1:.3f}] {params['x2ndaxistitle_units']}")
        print("")
        print(f"The preserved frequencies are [{f1:.3f}, {f2:.3f}] {params['yttitle_units']}")
        print(f"The preserved periods are [{int(1 / (f2))}, {int(1 / (f1))}] {params['y2ndaxistitle_units']}")
        print("")

        # Perform the 3D Fourier transform
        print("Making a 3D Fourier transform of the input datacube ..........")
        threedft = np.fft.fftshift(np.fft.fftn(signal))

        # Read both axes directly from the shifted DFT grid.
        dk = 2 * np.pi / (nx * pixelsize)
        df = 1.0 / (nt * cadence)
        spatial_frequencies = np.fft.fftshift(np.fft.fftfreq(nx, pixelsize)) * (2 * np.pi)
        temporal_frequencies = np.fft.fftshift(np.fft.fftfreq(nt, cadence))
        spatial_indices = np.rint(spatial_frequencies / dk)
        spatial_bins = np.hypot(spatial_indices[:, None], spatial_indices[None, :])
        absolute_time_bin = np.rint(np.abs(temporal_frequencies) / df).astype(int)
        if k1 > k2 or f1 > f2:
            raise ValueError("Passband bounds must be in increasing order.")

        if params['no_spatial_filt']:
            spatial_ring_filter = np.ones((nx, ny))
        else:
            low, high = k1 / dk, k2 / dk
            tolerance = 8 * np.finfo(float).eps * max(1.0, high)
            annulus = (spatial_bins >= low - tolerance) & (spatial_bins <= high + tolerance)
            spatial_ring_filter = annulus.astype(float)
            if params['spatial_torus']:
                # Preserve the legacy nested-ring taper in FFT-bin units.
                # The shifted FFT origin replaces the image's midpoint.
                positive_k = np.arange(nx // 2 + 1) * dk
                first_k = np.argmin(np.abs(positive_k - k1))
                last_k = np.argmin(np.abs(positive_k - k2))
                half_depth = (last_k - first_k) // 2
                centre = (first_k + last_k) // 2
                if half_depth:
                    spatial_ring_filter.fill(0)
                    for width in range(half_depth + 1):
                        # The outer ring overwrote the preceding central
                        # slice pair in the original torus.
                        if width != half_depth - 1:
                            spatial_ring_filter += (
                                (spatial_bins > centre - width)
                                & (spatial_bins <= centre + width + 1)
                            )
                    peak = spatial_ring_filter.max()
                    if peak:
                        spatial_ring_filter /= peak

        # Keep the original nearest-bin selection and Gaussian kernel width.
        positive_f = np.arange(nt // 2 + 1) * df
        first_bin = np.argmin(np.abs(positive_f - f1))
        last_bin = np.argmin(np.abs(positive_f - f2))
        temporal_pass = (absolute_time_bin >= first_bin) & (absolute_time_bin <= last_bin)
        if params['no_temporal_filt']:
            temporal_filter = np.ones(nt, dtype=float)
        elif params['temporal_torus'] and first_bin < last_bin:
            filter_width = last_bin - first_bin

            # Determine sigma based on filter width
            if filter_width < 25:
                sigma = 3
            if filter_width >= 25 and filter_width < 30:
                sigma = 4
            if filter_width >= 30 and filter_width < 40:
                sigma = 5
            if filter_width >= 40 and filter_width < 45:
                sigma = 6
            if filter_width >= 45 and filter_width < 50:
                sigma = 7
            if filter_width >= 50 and filter_width < 55:
                sigma = 8
            if filter_width >= 55 and filter_width < 60:
                sigma = 9
            if filter_width >= 60 and filter_width < 65:
                sigma = 10
            if filter_width >= 65 and filter_width < 70:
                sigma = 11
            if filter_width >= 70 and filter_width < 80:
                sigma = 12
            if filter_width >= 80 and filter_width < 90:
                sigma = 13
            if filter_width >= 90 and filter_width < 100:
                sigma = 14
            if filter_width >= 100 and filter_width < 110:
                sigma = 15
            if filter_width >= 110 and filter_width < 130:
                sigma = 16
            if filter_width >= 130:
                sigma = 17

            temporal_gaussian = gaussian_function(sigma=sigma, width=filter_width)
            temporal_gaussian /= temporal_gaussian.max()
            # The original signed slices were half-open. Pair their weights
            # to preserve their contribution to the real inverse transform.
            gains = np.r_[temporal_gaussian, 0.0]
            temporal_filter = np.zeros(nt, dtype=float)
            indices = absolute_time_bin[temporal_pass] - first_bin
            temporal_filter[temporal_pass] = (gains[indices] + gains[filter_width - indices]) / 2
            # DC and the even-length Nyquist bin are each their own partner.
            if first_bin == 0:
                temporal_filter[absolute_time_bin == 0] = temporal_gaussian[0]
            if nt % 2 == 0 and last_bin == nt // 2:
                temporal_filter[absolute_time_bin == last_bin] = temporal_gaussian[0]
        else:
            temporal_filter = temporal_pass.astype(float)

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # Create useful variables for plotting (if needed), demonstrating the filtering process
        if params['processing_maps']:
            # Define the spatial frequency step for plotting
            spatial_dx = spatial_frequencies[1] - spatial_frequencies[0]
            spatial_dy = spatial_frequencies[1] - spatial_frequencies[0]

            # Torus map
            torus_map = {
                'data': spatial_ring_filter,
                'dx': spatial_dx,
                'dy': spatial_dy,
                'xc': 0,
                'yc': 0,
                'time': '',
                'units': 'pixels'
            }

            # Compute the total spatial FFT
            spatial_fft = np.nansum(threedft, axis=0)

            # Spatial FFT map
            spatial_fft_map = {
                'data': np.log10(spatial_fft),
                'dx': spatial_dx,
                'dy': spatial_dy,
                'xc': 0,
                'yc': 0,
                'time': '',
                'units': 'pixels'
            }

            # Spatial FFT filtered and its map
            spatial_fft_filtered = spatial_fft * spatial_ring_filter
            spatial_fft_filtered_map = {
                'data': np.log10(np.maximum(spatial_fft_filtered, 1e-15)),
                'dx': spatial_dx,
                'dy': spatial_dy,
                'xc': 0,
                'yc': 0,
                'time': '',
                'units': 'pixels'
            }

            # Compute the total temporal FFT
            temporal_fft = np.nansum(np.nansum(threedft, axis=2), axis=1)

        else:
            spatial_fft_map = None
            torus_map = None
            spatial_fft_filtered_map = None
            temporal_fft = None

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        # Apply the frequency-domain masks
        for i in range(nt): 
            threedft[i, :, :] *= spatial_ring_filter

        for x in range(nx):  
            for y in range(ny):  
                threedft[:, x, y] *= temporal_filter

        # Inverse FFT to get the filtered cube
        # filtered_cube = np.real(np.fft.ifftn(threedft, norm='ortho'))
        filtered_cube = np.real(np.fft.ifftn(np.fft.ifftshift(threedft)))
        if not params['processing_maps']:
            temporal_filter, temporal_frequencies, spatial_frequencies = None, None, None

    else:
        filtered_cube = None

        print("Filtered datacube generated.")

    return komegamap, kopower_xscale[1:], kopower_yscale[1:], filtered_cube, spatial_fft_map, torus_map, spatial_fft_filtered_map, temporal_fft, temporal_filter, temporal_frequencies, spatial_frequencies

