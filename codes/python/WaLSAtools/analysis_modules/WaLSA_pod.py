# --------------------------------------------------------------------------------------------------------
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
# --------------------------------------------------------------------------------------------------------
# The following codes are based on those originally written by Jonathan E. Higham & Luiz A. C. A. Schiavo
# --------------------------------------------------------------------------------------------------------

import numpy as np # type: ignore
from scipy.signal import welch, find_peaks # type: ignore
from scipy.linalg import svd # type: ignore
from scipy.optimize import curve_fit # type: ignore
import math

def print_pod_results(results):
    """
    Print a summary of the results from WaLSA_pod, including parameter descriptions, types, and shapes.
    
    Parameters:
    -----------
    results : dict
        The dictionary containing all POD results and relevant outputs.
    """
    descriptions = {
        'input_data': 'Original input data, mean subtracted (Shape: (Nt, Ny, Nx))',
        'spatial_mode': 'Reshaped spatial modes matching the dimensions of the input data (Shape: (Nmodes, Ny, Nx))',
        'temporal_coefficient': 'Temporal coefficients associated with each spatial mode (Shape: (Nmodes, Nt))',
        'eigenvalue': 'Eigenvalues corresponding to singular values squared (Shape: (Nmodes))',
        'eigenvalue_contribution': 'Eigenvalue contribution of each mode (Shape: (Nmodes))',
        'cumulative_eigenvalues': 'Cumulative percentage of eigenvalues for the first "num_cumulative_modes" modes (Shape: (num_cumulative_modes))',
        'combined_welch_psd': 'Combined Welch power spectral density for the temporal coefficients of the firts "num_modes" modes (Shape: (Nf))',
        'frequencies': 'Frequencies identified in the Welch spectrum (Shape: (Nf))',
        'combined_welch_significance': 'Significance threshold of the combined Welch spectrum (Shape: (Nf,))',
        'reconstructed': 'Reconstructed frame at the specified timestep (or for the entire time series) using the top "num_modes" modes (Shape: (Ny, Nx))',
        'sorted_frequencies': 'Frequencies identified in the Welch combined power spectrum (Shape: (Nfrequencies))',
        'frequency_filtered_modes': 'Frequency-filtered spatial POD modes for the first "num_top_frequencies" frequencies (Shape: (Nt, Ny, Nx, num_top_frequencies))',
        'frequency_filtered_modes_frequencies': 'Frequencies corresponding to the frequency-filtered modes (Shape: (num_top_frequencies))',
        'SPOD_spatial_modes': 'SPOD spatial modes if SPOD is used (Shape: (Nspod_modes, Ny, Nx))',
        'SPOD_temporal_coefficients': 'SPOD temporal coefficients if SPOD is used (Shape: (Nspod_modes, Nt))',
        'p': 'Left singular vectors (spatial modes) from SVD (Shape: (Nx, Nmodes))',
        's': 'Singular values from SVD (Shape: (Nmodes))',
        'a': 'Right singular vectors (temporal coefficients) from SVD (Shape: (Nmodes, Nt))'
    }
    
    print("\n---- POD/SPOD Results Summary ----\n")
    for key, value in results.items():
        desc = descriptions.get(key, 'No description available')
        shape = np.shape(value) if value is not None else 'None'
        dtype = type(value).__name__
        print(f"{key} ({dtype}, Shape: {shape}): {desc}")
    print("\n----------------------------------")

def spod(data, **kwargs):
    """
    Perform Spectral Proper Orthogonal Decomposition (SPOD) analysis on input data.
    Steps:
    1. Load the data.
    2. Compute a time average.
    3. Build the correlation matrix to compute the POD using the snapshot method.
    4. Compute eigenvalue decomposition for the correlation matrix for the fluctuation field. 
    The eigenvalues and eigenvectors can be computed by any eigenvalue function or by an SVD procedure, since the correlation matrix is symmetric positive-semidefinite.
    5. After obtaining the eigenvalues and eigenvectors, compute the temporal and spatial modes according to the snapshot method.
    
    Parameters:
    -----------
    data : np.ndarray
        3D data array with shape (time, x, y) or similar. 
    **kwargs : dict, optional
        Additional keyword arguments to configure the analysis.
    
    Returns:
    --------
    SPOD_spatial_modes : np.ndarray
    SPOD_temporal_coefficients : np.ndarray
    """
    # Set default parameter values
    defaults = {
        'silent': False,
        'num_modes': None,
        'filter_size': None,
        'periodic_matrix': True
    }
    
    # Update defaults with user-provided values
    params = {**defaults, **kwargs}

    if params['num_modes'] is None:
        params['num_modes'] = data.shape[0]
    if params['filter_size'] is None:
        params['filter_size'] = params['num_modes']

    if not params['silent']:
        print("Starting SPOD analysis ....")

    nsnapshots, ny, nx = data.shape

    # Center the data by subtracting the mean
    time_average = data.mean(axis=0)
    fluctuation_field = np.zeros(data.shape, dtype=np.result_type(data.dtype, np.float64))
    for n in range(nsnapshots):
        fluctuation_field[n,:,:] = data[n,:,:] - time_average  # Subtract time-average from each snapshot

    # Build the correlation matrix (snapshot method)
    correlation_matrix = np.zeros((nsnapshots, nsnapshots))
    for i in range(nsnapshots):
        for j in range(i, nsnapshots):
            correlation_matrix[i, j] = (fluctuation_field[i, :, :] * fluctuation_field[j, :, :]).sum() / (nsnapshots * nx * ny)
            correlation_matrix[j, i] = correlation_matrix[i, j]

    # SPOD correlation matrix with periodic boundary conditions
    nmatrix = np.zeros((3 * nsnapshots, 3 * nsnapshots))
    nmatrix[nsnapshots:2 * nsnapshots, nsnapshots:2 * nsnapshots] = correlation_matrix[0:nsnapshots, 0:nsnapshots]

    if params['periodic_matrix']:
        for i in range(3):
            for j in range(3):
                xs = i * nsnapshots  # x-offset for periodic positioning
                ys = j * nsnapshots  # y-offset for periodic positioning
                nmatrix[0 + xs:nsnapshots + xs, 0 + ys:nsnapshots + ys] = correlation_matrix[0:nsnapshots, 0:nsnapshots]

    # Apply Gaussian filter
    gk = np.zeros((2 * params['filter_size'] + 1))  # Create 1D Gaussian kernel
    esp = 8.0  # Exponential parameter controlling the spread of the filter
    sumgk = 0.0  # Sum of the Gaussian kernel for normalization

    for n in range(2 * params['filter_size'] + 1):
        k = -params['filter_size'] + n  # Offset from the center of the kernel
        gk[n] = math.exp(-esp * (k / max(1, params['filter_size']))**2)  # Gaussian filter formula
        sumgk += gk[n]  # Sum of kernel values for normalization

    # Filter the extended correlation matrix without feeding filtered values back in.
    unfiltered = nmatrix.copy()
    for j in range(nsnapshots, 2 * nsnapshots):
        for i in range(nsnapshots, 2 * nsnapshots):
            aux = 0.0  # Initialize variable for the weighted sum
            for n in range(2 * params['filter_size'] + 1):
                k = -params['filter_size'] + n  # Offset from the current index
                aux += unfiltered[i + k, j + k] * gk[n]  # Apply Gaussian weighting
            nmatrix[i, j] = aux / sumgk  # Normalize by the sum of the Gaussian kernel
    
    # Extract the SPOD correlation matrix from the central part of the filtered matrix
    spectral_matrix = nmatrix[nsnapshots:2 * nsnapshots, nsnapshots:2 * nsnapshots]

    # Perform SVD to compute SPOD
    UU, SS, VV = svd(spectral_matrix, full_matrices=True, compute_uv=True)

    # Temporal coefficients
    SPOD_temporal_coefficients = np.zeros((nsnapshots, nsnapshots))
    for k in range(nsnapshots):
        SPOD_temporal_coefficients[:, k] = np.sqrt(SS[k] * nsnapshots) * UU[:, k]

    # Extract spatial SPOD modes
    SPOD_spatial_modes = np.zeros((params['num_modes'], ny, nx))
    for m in range(params['num_modes']):
        if SS[m] == 0:
            continue
        for t in range(nsnapshots):
            SPOD_spatial_modes[m, :, :] += fluctuation_field[t, :, :] * SPOD_temporal_coefficients[t, m] / (SS[m] * nsnapshots)

    if not params['silent']:
        print(f"SPOD analysis completed.")
    
    return SPOD_spatial_modes, SPOD_temporal_coefficients


# -------------------------------------- Main Function ----------------------------------------
def WaLSA_pod(signal, time, **kwargs):
    """
    Perform Proper Orthogonal Decomposition (POD) analysis on input data.
    
    Parameters:
        signal (array): 3D data cube with shape (time, x, y) or similar.
        time (array): 1D array representing the time points for each time step in the data.
        num_modes (int, optional): Number of top modes to compute. Default is None (all modes).
        num_top_frequencies (int, optional): Number of top frequencies to consider. Default is None (all frequencies).
        top_frequencies (list, optional): List of top frequencies to consider. Default is None.
        num_cumulative_modes (int, optional): Number of cumulative modes to consider. Default is None (all modes).
        welch_nperseg (int, optional): Number of samples per segment for Welch's method. Default is 150.
        welch_noverlap (int, optional): Number of overlapping samples for Welch's method. Default is 25.
        welch_nfft (int, optional): Number of points for the FFT. Default is 2^14.
        welch_fs (int, optional): Sampling frequency for the data. By default inferred from time.
        nperm (int, optional): Number of permutations for significance testing. Default is 1000.
        siglevel (float, optional): Significance level for the Welch spectrum. Default is 0.95.
        timestep_to_reconstruct (int, optional): Timestep of the datacube to reconstruct using the top modes. Default is 0.
        num_modes_reconstruct (int, optional): Number of modes to use for reconstruction. Default is None (all modes).
        reconstruct_all (bool, optional): If True, reconstruct the entire time series using the top modes. Default is False.
        spod (bool, optional): If True, perform Spectral Proper Orthogonal Decomposition (SPOD) analysis. Default is False.
        spod_filter_size (int, optional): Filter size for SPOD analysis. Default is None.
        spod_num_modes (int, optional): Number of SPOD modes to compute. Default is None.
        print_results (bool, optional): If True, print a summary of results. Default is True.

    **kwargs : Additional keyword arguments to configure the analysis.
    
    Returns:
    results : dict
        A dictionary containing all computed POD results and relevant outputs. 
        See 'descriptions' in the 'print_pod_results' function on top of this page.
    """
    # Set default parameter values
    defaults = {
        'silent': False,
        'num_modes': None,
        'num_top_frequencies': None,
        'top_frequencies': None,
        'num_cumulative_modes': None,
        'welch_nperseg': 150,
        'welch_noverlap': 25,
        'welch_nfft': 2**14,
        'welch_fs': 1.0 / np.median(np.diff(time)),
        'nperm': 1000,
        'siglevel': 0.95,
        'timestep_to_reconstruct': 0,
        'reconstruct_all': False,
        'num_modes_reconstruct': None,
        'spod': False,
        'spod_filter_size': None,
        'spod_num_modes': None,
        'print_results': True  # Print summary of results by default
    }
    
    # Update defaults with user-provided values
    params = {**defaults, **kwargs}

    data = signal

    if params['num_modes'] is None:
        params['num_modes'] = min(data.shape[0], data.shape[1] * data.shape[2])

    if params['num_top_frequencies'] is None:
        params['num_top_frequencies'] = min(params['num_modes'], 10)
    
    if params['num_cumulative_modes'] is None:
        params['num_cumulative_modes'] = min(params['num_modes'], 10)

    if params['num_modes_reconstruct'] is None:
        params['num_modes_reconstruct'] = min(params['num_modes'], 10)

    if not params['silent']:
        print("Starting POD analysis ....")
        print(f"Processing a 3D cube with shape {data.shape}.")

    # The first step is to read in the data and then perform the SVD (or POD). Before we do this, 
    # we need to reshape the matrix such that it is an N x T matrix where N is the column vectorized set of each spatial image and T is the temporal domain. 
    # We also need to ensure that the mean is subtracted from the data; this will ensure we are looking only at the variance of the data, and mode 1 will not be contaminated with the mean image.
    # Reshape the 3D data into 2D array where each row is a vector from the original 3D data
    inp = np.reshape(data, (data.shape[0], data.shape[1] * data.shape[2])).astype(np.float64)
    inp = inp.T  # Transpose the matrix to have spatial vectors as columns

    # Center the data by subtracting the mean
    mean_per_row = np.nanmean(inp, axis=1, keepdims=True)
    mean_replicated = np.tile(mean_per_row, (1, data.shape[0]))
    inp_centered = inp - mean_replicated

    # Input data, mean subtracted
    input_dat = np.zeros((data.shape[0], data.shape[1], data.shape[2]))
    for im in range(data.shape[0]):
        input_dat[im, :, :] = np.reshape(inp_centered[:, im], (data.shape[1], data.shape[2]))

    # Perform SVD to compute POD
    p, s, a = svd(inp_centered, full_matrices=False)
    sorg = s.copy()  # Store original singular values
    eigenvalue = s**2  # Convert singular values to eigenvalues

    # Reshape spatial modes to have the same shape as the input data
    num_modes = p.shape[1]
    spatial_mode = np.zeros((num_modes, data.shape[1], data.shape[2]))
    for m in range(num_modes):
        spatial_mode[m, :, :] = np.reshape(p[:, m], (data.shape[1], data.shape[2]))

    # Extract temporal coefficients
    temporal_coefficient = a[:num_modes, :]

    # Calculate eigenvalue contributions
    eigenvalue_contribution = eigenvalue / np.sum(eigenvalue)

    # Calculate cumulative eigenvalues
    cumulative_eigenvalues = []
    for m in range(params['num_cumulative_modes']): 
        contm = 100 * eigenvalue[:m + 1] / np.sum(eigenvalue)
        cumulative_eigenvalues.append(np.sum(contm)) 

    if params['reconstruct_all']:
        # Reconstruct the entire time series using the top 'num_modes_reconstruct' modes
        reconstructed = np.zeros((data.shape[0], data.shape[1], data.shape[2]))
        for tindex in range(data.shape[0]):
            reconim = np.zeros((data.shape[1], data.shape[2]))
            for i in range(params['num_modes_reconstruct']):
                reconim=reconim+np.reshape(p[:, i], (data.shape[1], data.shape[2]))*a[i, tindex]*sorg[i]
            reconstructed[tindex, :, :] = reconim
    else:
        # Reconstruct the specified timestep using the top 'num_modes_reconstruct' modes
        reconstructed = np.zeros((data.shape[1], data.shape[2]))
        for i in range(params['num_modes_reconstruct']):
            reconstructed=reconstructed+np.reshape(p[:, i], (data.shape[1], data.shape[2]))*a[i, params['timestep_to_reconstruct']]*sorg[i]

    #---------------------------------------------------------------------------------
    # Combined Welch power spectrum and its significance
    #---------------------------------------------------------------------------------
    # Compute Welch power spectrum for each mode and combine them for 'num_modes' modes
    combined_welch_psd = []
    for m in range(params['num_modes']):
        frequencies, px = welch(a[m, :] - np.mean(a[m, :]), 
                      nperseg=params['welch_nperseg'], 
                      noverlap=params['welch_noverlap'], 
                      nfft=params['welch_nfft'], 
                      fs=params['welch_fs'])
        if m == 0:
            combined_welch_psd = np.zeros((len(frequencies),))
        combined_welch_psd += eigenvalue_contribution[m] * (px / np.sum(px))

    # Generate resampled peaks to compute significance threshold
    resampled_peaks = np.zeros((params['nperm'], len(frequencies)))
    for i in range(params['nperm']):
        resampled_data = np.random.randn(*a.shape)
        resampled_peak = np.zeros((len(frequencies),))
        for m in range(params['num_modes']):
            _, px = welch(resampled_data[m, :] - np.mean(resampled_data[m, :]), 
                          nperseg=params['welch_nperseg'], 
                          noverlap=params['welch_noverlap'], 
                          nfft=params['welch_nfft'], 
                          fs=params['welch_fs'])
            resampled_peak += eigenvalue_contribution[m] * (px / np.sum(px))
        resampled_peak /= (np.max(resampled_peak) + 1e-30) # a small epsilon added to avoid division by zero
        resampled_peaks[i, :] = resampled_peak
    # Calculate significance threshold
    combined_welch_significance = np.percentile(resampled_peaks, params['siglevel'] * 100, axis=0)
    #---------------------------------------------------------------------------------

    # Find peaks in the combined spectrum and sort them in descending order
    normalized_peak = combined_welch_psd / (np.max(combined_welch_psd) + 1e-30)
    peak_indices, _ = find_peaks(normalized_peak)
    sorted_indices = np.argsort(normalized_peak[peak_indices])[::-1]
    sorted_frequencies = frequencies[peak_indices][sorted_indices]

    # Generate a table of the top N frequencies
    # Cleaner single-line list comprehension
    table_data = [
        [float(np.round(freq, 2)), float(np.round(pwr, 2))] 
        for freq, pwr in zip(sorted_frequencies[:params['num_top_frequencies']], 
                            normalized_peak[peak_indices][sorted_indices][:params['num_top_frequencies']])
    ]

    # The frequencies that the POD is able to capture have now been identified (in the top 'num_top_frequencies' modes). 
    # These frequencies can be fitted to the temporal coefficients of the top 'num_top_frequencies' modes, 
    # allowing for a representation of the data described solely by these "pure" frequencies. 
    # This approach enables the reconstruction of the original data using the identified dominant frequencies,
    # resulting in frequency-filtered spatial POD modes.
    if params['top_frequencies'] is None:
        top_frequencies = sorted_frequencies[:params['num_top_frequencies']]
    else:    
        top_frequencies = params['top_frequencies']
    params['num_top_frequencies'] = len(top_frequencies)
    clean_data = np.zeros((inp.shape[0], inp.shape[1], len(top_frequencies)))

    for i in range(params['num_top_frequencies']):
        def model_fun(t, amplitude, phase):
            """
            Generate a sinusoidal model function.

            Parameters:
            t (array-like): The time variable.
            amplitude (float): The amplitude of the sinusoidal function.
            phase (float): The phase shift of the sinusoidal function.

            Returns:
            array-like: The computed sinusoidal values at each time point `t`.
            """
            return amplitude * np.sin(2 * np.pi * top_frequencies[i] * t + phase)
        
        for j in range(params['num_modes_reconstruct']):
            csignal = a[j,:]
            # Initial Guesses for Parameters: [Amplitude, Phase]
            initial_guess = [1, 0]
            # Nonlinear Fit
            fit_params, _ = curve_fit(model_fun, time, csignal, p0=initial_guess)
            # Clean Signal from Fit
            if j == 0:
                clean_signal=np.zeros((params['num_modes_reconstruct'],len(csignal)))
            clean_signal[j,:] = model_fun(time, *fit_params)
        # forming a set of clean data (reconstructed data at the fitted frequencies; frequency filtered reconstructed data)
        clean_data[:,:,i] = p[:,0:params['num_modes_reconstruct']]@np.diag(sorg[0:params['num_modes_reconstruct']])@clean_signal

    frequency_filtered_modes = np.zeros((data.shape[0], data.shape[1], data.shape[2], params['num_top_frequencies']))
    for jj in range(params['num_top_frequencies']):
        for frame_index in range(data.shape[0]):
            frequency_filtered_modes[frame_index, :, :, jj] = np.reshape(clean_data[:, frame_index, jj], (data.shape[1], data.shape[2]))

    if params['spod']:
        SPOD_spatial_modes, SPOD_temporal_coefficients = spod(data, num_modes=params['spod_num_modes'], filter_size=params['spod_filter_size'])
    else:
        SPOD_spatial_modes = None
        SPOD_temporal_coefficients = None

    results = {
        'input_data': input_dat, # Original input data, mean subtracted (Shape: (Nt, Ny, Nx))
        'spatial_mode': spatial_mode,  # POD spatial modes matching the dimensions of the input data (Shape: (Nmodes, Ny, Nx))
        'temporal_coefficient': temporal_coefficient,  # POD temporal coefficients associated with each spatial mode (Shape: (Nmodes, Nt))
        'eigenvalue': eigenvalue,  # Eigenvalues corresponding to singular values squared (Shape: (Nmodes))
        'eigenvalue_contribution': eigenvalue_contribution,  # Eigenvalue contribution of each mode (Shape: (Nmodes))
        'cumulative_eigenvalues': cumulative_eigenvalues,  # Cumulative percentage of eigenvalues for the first 'num_cumulative_modes' modes (Shape: (Ncumulative_modes))
        'combined_welch_psd': combined_welch_psd,  # Combined Welch power spectral density for the temporal coefficients (Shape: (Nf))
        'frequencies': frequencies,  # Frequencies identified in the Welch spectrum (Shape: (Nf))
        'combined_welch_significance': combined_welch_significance,  # Significance threshold of the combined Welch spectrum (Shape: (Nf))
        'reconstructed': reconstructed,  # Reconstructed frame at the specified timestep (or for the entire time series) using the top modes (Shape: (Ny, Nx))
        'sorted_frequencies': sorted_frequencies,  # Frequencies identified in the Welch spectrum (Shape: (Nfrequencies,))
        'frequency_filtered_modes': frequency_filtered_modes,  # Frequency-filtered spatial POD modes (Shape: (Nt, Ny, Nx, Ntop_frequencies))
        'frequency_filtered_modes_frequencies': top_frequencies,  # Frequencies corresponding to the frequency-filtered modes (Shape: (Ntop_frequencies))
        'SPOD_spatial_modes': SPOD_spatial_modes,  # SPOD spatial modes if SPOD is used (Shape: (Nspod_modes, Ny, Nx))
        'SPOD_temporal_coefficients': SPOD_temporal_coefficients,  # SPOD temporal coefficients if SPOD is used (Shape: (Nspod_modes, Nt))
        'p': p,  # Left singular vectors (spatial modes) from SVD (Shape: (Nx, Nmodes))
        's': sorg,  # Singular values from SVD (Shape: (Nmodes,))
        'a': a  # Right singular vectors (temporal coefficients) from SVD (Shape: (Nmodes, Nt))
    }

    if not params['silent']:
        print(f"POD analysis completed.")
        print(f"Top {params['num_top_frequencies']} frequencies and normalized power values:\n{table_data}")
        print(f"Total variance contribution of the first {params['num_modes']} modes: {np.sum(100 * eigenvalue[:params['num_modes']] / np.sum(eigenvalue)):.2f}%")

    if params['print_results']:
        print_pod_results(results)

    return results
