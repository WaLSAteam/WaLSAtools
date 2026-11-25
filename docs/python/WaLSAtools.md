---
template: main.html
---

# WaLSAtools

**WaLSAtools** is designed for ease of use and accessibility. Its interactive interface guides you through the analysis process, providing clear instructions and helpful information at each step. This section demonstrates how to use WaLSAtools and highlights its key features.

Before diving into the interactive demonstration, we recommend familiarizing yourself with the various analysis methods available in WaLSAtools. You can find detailed descriptions of these methods in the [Introduction](introduction.md) section. Additionally, this page provides several **Worked Examples** of different analysis techniques applied to synthetic datasets (see the left menu). To learn more about its capabilities and how to apply it to your research, we encourage you to explore the WaLSAtools documentation, the associated *[Nature Reviews Methods Primers](https://www.nature.com/articles/s43586-025-00392-0)* article ([full-text access to a view-only version](https://WaLSA.tools/nrmp) and its [Supplementary Information](https://WaLSA.tools/nrmp-si)), and the provided examples. If you use WaLSAtools in your work, please remember to cite it appropriately (see [Citation][9]).

The ["Under the Hood"](routines.md) section provides details on the individual routines used for wave analysis within the **WaLSAtools** package, for those interested in exploring the underlying code. However, we strongly encourage all users to perform their analyses by running WaLSAtools directly, as this ensures the correct execution of the analysis workflow and provides a more user-friendly experience.

  [9]: ../citation.md

!!! walsa-gear "Interactive Demonstration"

    **WaLSAtools** provides an interactive interface that simplifies wave analysis. To launch the interface, simply import the `WaLSAtools` package and run the `WaLSAtools` command in a Python terminal or Jupyter notebook.

    The interface will guide you through the following steps:

    1.  **Select a category of analysis:** Choose from single time series analysis or cross-correlation analysis.
    2.  **Choose the data type:** Specify the type of data you are working with (e.g., 1D signal, 3D datacube).
    3.  **Pick a specific analysis method:** Select the method most suitable for your data and research question.

    The interface will then provide information on the selected method, including its calling sequence, input parameters, and expected outputs.

    **Here's an example of the execution of WaLSAtools in a Jupyter notebook:**

<style>
    .dropdown-container {
        margin-left: 30px;
        margin-top: 20px;
        font-size: 0.9em;
        line-height: 2;
        display: flex;
        flex-direction: column;
        gap: 10px;
    }
    .dropdown-row {
        display: flex;
        align-items: center;
        gap: 10px;
    }   
    select {
        width: 270px !important;
        height: 33px !important;
        padding: 5px 10px;
        font-size: 1em;
        color: #333;
        border: 1px solid #ccc;
        border-radius: 2px;
        background-color: #f9f9f9;
        appearance: none;
        -webkit-appearance: none;
        -moz-appearance: none;
        background-image: url('data:image/svg+xml,%3Csvg xmlns="http://www.w3.org/2000/svg" width="12" height="12" viewBox="0 0 12 12"%3E%3Cpath d="M1,4 L6,9 L11,4" fill="none" stroke="%23000" stroke-width="2.0" /%3E%3C/svg%3E');
        background-repeat: no-repeat;
        background-position: right 10px center;
        background-size: 12px 12px;
        cursor: pointer;
    }
    select:focus {
        outline: none;
        border-color: #4caf50;
        box-shadow: 0 0 5px rgba(76, 175, 80, 0.6);
    }
    select:disabled {
        background-color: #eaeaea;
        cursor: not-allowed;
    }   
    .output-container {
        margin-left: 30px;
        margin-top: 5px;
        padding: 0;
        display: none; /* Hidden by default */
    }
    .parameters-table {
        border-collapse: collapse;
        border: 1px solid #222;
        width: calc(100% - 30px);
        box-sizing: border-box;
        table-layout: auto;
        margin-top: 20px;
        font-size: 0.9em;
    }
    .parameters-table td {
        border: 1px solid #222;
        padding: 8px;
        text-align: left;
    }
    .parameters-table th {
        padding: 8px;
        text-align: left;
    }   
    .code-container {
        font-family: monospace;
        position: relative;
        background-color: #f7f7f7;
        border: 1px solid #ddd;
        margin-left: 30px; 
        padding: 10px;
        padding-left: 40px; /* For line numbers */
        font-size: 14px;
        line-height: 1.6;
        display: inline-block;
        width: calc(100% - 30px);
        box-sizing: border-box;
    }
    .line-numbers {
        position: absolute;
        top: 10px;
        left: 10px;
        color: #888;
        text-align: right;
        line-height: 1.6;
        font-size: 14px;
    }
    .execute-btn {
        background-color: #4caf50;
        color: white;
        border: none;
        border-radius: 50%;
        width: 20px;
        height: 20px;
        display: flex;
        justify-content: center;
        align-items: center;
        font-size: 16px;
        font-weight: bold;
        cursor: pointer;
        position: absolute;
        left: -30px; /* Place the button outside the box */
        top: 2px;
    }
    .execute-btn:hover {
        background-color: #45a049;
    }
    .index-number {
        position: absolute;
        bottom: 0px;
        left: -30px;
        color: #888;
        font-size: 12px;
    }
    .python-label {
        position: absolute;
        bottom: 2px;
        right: 6px;
        font-size: 12px;
        color: #888;
    }
    .walsa-gear-code {
        border-left: 4px solid rgb(189, 26, 31); 
        background-color: rgba(255, 255, 255, 0.1);
        padding: 16px; 
        border-radius: 1px; 
        margin-bottom: 16px; 
        box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.1);
    }
</style>
<div class="walsa-gear-code">
<div style="font-family: Arial, sans-serif; margin: 20px; font-size: 0.9em; line-height: 0.9em;">
<div class="code-container">
    <!-- Line numbers -->
    <div class="line-numbers">
        1<br>
        2
    </div>
    <!-- Code area -->
    <div>
        <span style="color:DodgerBlue">from</span> WaLSAtools <span style="color:DodgerBlue">import</span> WaLSAtools<br>
        WaLSAtools
    </div>
    <!-- Execute button -->
    <!-- <button class="execute-btn">▶</button> -->
    <!-- Index number -->
    <div class="index-number">[1]</div>  
    <!-- Python label -->
    <div class="python-label">Python</div>
</div>
<!-- Logo -->
<div style="margin-left: 30px; margin-top: 20px;">
    <img src="https://walsa.team/images/walsatools_black.png" style="width: 300px; height: auto;">
</div>
<!-- Credits -->
<div style="margin-left: 30px; margin-top: 20px; font-size: 0.9em;">
    <p>© WaLSA Team (<a href="https://www.WaLSA.team" target="_blank" style="color: #4169E1; text-decoration: none;">www.WaLSA.team</a>)</p>
    <hr style="width: 70%; margin: 0; border: 0.98px solid #888; margin-bottom: 10px;">
    <p><strong>WaLSAtools</strong> v1.0.1 - Wave analysis tools</p>
    <p>Documentation: <a href="https://www.WaLSA.tools" target="_blank" style="color: #4169E1; text-decoration: none;">www.WaLSA.tools</a></p>
    <p>GitHub repository: <a href="https://www.github.com/WaLSAteam/WaLSAtools" target="_blank" style="color: #4169E1; text-decoration: none;">www.github.com/WaLSAteam/WaLSAtools</a></p>
    <hr style="width: 70%; margin: 0; border: 0.98px solid #888; margin-bottom: 10px;">
    <p>If you use <strong>WaLSAtools</strong> in your research, please cite:</p>
    <p>Jafarzadeh, S., Jess, D. B., Stangalini, M. et al. 2025, <em>Nature Reviews Methods Primers</em>, 5, 21</p>
    <hr style="width: 70%; margin: 0; border: 0.98px solid #888; margin-bottom: 15px;">
    <p>Choose a category, data type, and analysis method from the list below,</p>
    <p>to get hints on the calling sequence and parameters:</p>
</div>
<!-- Dropdown Menus -->
<div class="dropdown-container">
    <div class="dropdown-row">
        <label for="category" style="width: 90px !important; text-align: right !important;">Category:</label>
        <select id="category">
            <option value="">Select Category</option>
            <option value="a">Single Time Series Analysis</option>
            <option value="b">Cross-Correlation Between Two Time Series</option>
        </select>
    </div>
    <div class="dropdown-row">
        <label for="datatype" style="width: 90px !important; text-align: right !important;">Data Type:</label>
        <select id="datatype" disabled>
            <option value="">Select Data Type</option>
        </select>
    </div>
    <div class="dropdown-row">
        <label for="analysisMethod" style="width: 90px !important; text-align: right !important;">Method:</label>
        <select id="analysisMethod" disabled>
            <option value="">Select Method</option>
        </select>
    </div>
    <div class="dropdown-row">
        <label for="subMethod" id="subMethodLabel" style="width: 90px !important; text-align: right !important; display:none;">Sub-method:</label>
        <select id="subMethod" style="display:none;">
            <option value="">Select Sub-method</option>
            <option value="fft">FFT</option>
            <option value="wavelet">Wavelet</option>
            <option value="lombscargle">Lomb-Scargle</option>
            <option value="welch">Welch</option>
        </select>
    </div>
</div>
<div id="dropdownMessage" style="margin-left: 30px; margin-top: 15px; font-size: 0.9em; display: none;">
    Please select appropriate options from all dropdown menus.
</div>
<div class="output-container" id="outputContainer">
    <p style="font-size: 1.0em;">Calling Sequence:</p>
    <span id="callingSequence" style="font-size: 0.95em !important; margin: 0 !important; padding: 0 !important;"></span>
    <div id="warningNote" style="display: none; margin-left: 30px; margin-bottom: 15px;">
        <div style="border: 2px solid #27AE60; padding: 10px; background-color: #E9F7EF; border-radius: 4px;">
            <strong style="color: #196F3D;">Note: </strong>
            <span id="warningMessage" style="color: #000;"></span>
        </div>
    </div>
    <table class="parameters-table">
        <thead>
            <tr style="background-color: #fff;"><th colspan="3" style="text-align: left; color: #000; font-size: 110%;">Parameters (**kwargs)</th></tr>
            <tr style="background-color: #222;">
                <th style="color: #fff; border-right: 1px solid #ccc; text-align: left; white-space: nowrap;">Parameter</th>
                <th style="color: #fff; border-right: 1px solid #ccc; text-align: left; white-space: nowrap;">Type</th>
                <th style="color: #fff; text-align: left; width: 100%;">Description</th>
            </tr>
        </thead>
        <tbody id="parameterTableBody">
        </tbody>
    </table>
</div>
</div>
</div>
<script>
    const parameters = {
        single_series: {
            fft: {
                returnValues: 'power, frequency, significance, amplitude',
                parameters: {
                    signal: { type: 'array', description: 'The input signal (1D).' },
                    time: { type: 'array', description: 'The time array corresponding to the signal.' },
                    siglevel: { type: 'float', description: 'Significance level for the confidence intervals. Default: 0.95.' },
                    nperm: { type: 'int', description: 'Number of permutations for significance testing. Default: 1000.' },
                    nosignificance: { type: 'bool', description: 'If True, skip significance calculation. Default: False.' },
                    apod: { type: 'float', description: 'Extent of apodization edges (of a Tukey window). Default: 0.1.' },
                    pxdetrend: { type: 'int', description: 'Subtract linear trend with time per pixel. Options: 1 (simple) or 2 (advanced). Default: 2.' },
                    polyfit: { type: 'int', description: 'Degree of polynomial fit for detrending the data. If set, a polynomial fit (instead of linear) is applied. Default: None' },
                    meantemporal: { type: 'bool', description: 'If True, apply simple temporal detrending by subtracting the mean signal from the data, skipping fitting procedures. Default: False.' },
                    meandetrend: { type: 'bool', description: 'If True, subtract the linear trend with time for the image means (spatial detrending). Default: False.' },
                    recon: { type: 'bool', description: 'If True, perform Fourier reconstruction of the input time series. This does not preserve amplitudes but is useful for examining frequencies far from the low-frequency range. Default: False.' },
                    resample_original: { type: 'bool', description: 'If True, and if recon is set to True, approximate values close to the original are returned for comparison. Default: False.' },
                    nodetrendapod: { type: 'bool', description: 'If True, neither detrending nor apodization is performed. Default: False.' },
                    amplitude: { type: 'bool', description: ' If True, return the amplitudes of the Fourier transform. Default: False.' },
                    silent: { type: 'bool', description: 'If True, suppress print statements. Default: False.' }
                }
            },
            wavelet: {
                returnValues: 'power, period, significance, coi, gws_power, gws_significance, rgws_power',
                parameters: {
                    signal: { type: 'array', description: 'The input signal (1D).' },
                    time: { type: 'array', description: 'The time array corresponding to the signal.' },
                    siglevel: { type: 'float', description: 'Significance level for the confidence intervals. Default: 0.95.' },
                    nperm: { type: 'int', description: 'Number of permutations for significance testing. Default: 1000.' },
                    mother: { type: 'str', description: 'The mother wavelet function to use. Default: "morlet".' },
                    GWS: { type: 'bool', description: 'If True, calculate the Global Wavelet Spectrum. Default: False.' },
                    RGWS: { type: 'bool', description: 'If True, calculate the Refined Global Wavelet Spectrum (time-integrated power, excluding COI and insignificant areas). Default: False.' },
                    dj: { type: 'float', description: 'Scale spacing. Smaller values result in better scale resolution but slower calculations. Default: 0.025.' },
                    s0: { type: 'float', description: 'Initial (smallest) scale of the wavelet. Default: 2 * dt.' },
                    J: { type: 'int', description: 'Number of scales minus one. Scales range from s0 up to s0 * 2**(J * dj), giving a total of (J + 1) scales. Default: (log2(N * dt / s0)) / dj.' },
                    lag1: { type: 'float', description: 'Lag-1 autocorrelation. Default: 0.0.' },
                    apod: { type: 'float', description: 'Extent of apodization edges (of a Tukey window). Default: 0.1.' },
                    pxdetrend: { type: 'int', description: 'Subtract linear trend with time per pixel. Options: 1 (simple) or 2 (advanced). Default: 2.' },
                    polyfit: { type: 'int', description: 'Degree of polynomial fit for detrending the data. If set, a polynomial fit (instead of linear) is applied. Default: None' },
                    meantemporal: { type: 'bool', description: 'If True, apply simple temporal detrending by subtracting the mean signal from the data, skipping fitting procedures. Default: False.' },
                    meandetrend: { type: 'bool', description: 'If True, subtract the linear trend with time for the image means (spatial detrending). Default: False.' },
                    recon: { type: 'bool', description: 'If True, perform Fourier reconstruction of the input time series. This does not preserve amplitudes but is useful for examining frequencies far from the low-frequency range. Default: False.' },
                    resample_original: { type: 'bool', description: 'If True, and if recon is set to True, approximate values close to the original are returned for comparison. Default: False.' },
                    nodetrendapod: { type: 'bool', description: 'If True, neither detrending nor apodization is performed. Default: False.' },
                    silent: { type: 'bool', description: 'If True, suppress print statements. Default: False.' }
                }
            },
            lombscargle: {
                returnValues: 'power, frequency, significance',
                parameters: {
                    signal: { type: 'array', description: 'The input signal (1D).' },
                    time: { type: 'array', description: 'The time array corresponding to the signal.' },
                    siglevel: { type: 'float', description: 'Significance level for the confidence intervals. Default: 0.95.' },
                    nperm: { type: 'int', description: 'Number of permutations for significance testing. Default: 1000.' },
                    dy: { type: 'array', description: 'Errors or observational uncertainties associated with the time series.' },
                    fit_mean: { type: 'bool', description: 'If True, include a constant offset as part of the model at each frequency. This improves accuracy, especially for incomplete phase coverage.' },
                    center_data: { type: 'bool', description: 'If True, pre-center the data by subtracting the weighted mean of the input data. This is especially important if fit_mean=False.' },
                    nterms: { type: 'int', description: 'Number of terms to use in the Fourier fit. Default: 1.' },
                    normalization: { type: 'str', description: 'The normalization method for the periodogram. Options: "standard", "model", "log", "psd". Default: "standard".' },
                    nosignificance: { type: 'bool', description: 'If True, skip significance calculation. Default: False.' },
                    apod: { type: 'float', description: 'Extent of apodization edges (of a Tukey window). Default: 0.1.' },
                    pxdetrend: { type: 'int', description: 'Subtract linear trend with time per pixel. Options: 1 (simple) or 2 (advanced). Default: 2.' },
                    polyfit: { type: 'int', description: 'Degree of polynomial fit for detrending the data. If set, a polynomial fit (instead of linear) is applied. Default: None' },
                    meantemporal: { type: 'bool', description: 'If True, apply simple temporal detrending by subtracting the mean signal from the data, skipping fitting procedures. Default: False.' },
                    meandetrend: { type: 'bool', description: 'If True, subtract the linear trend with time for the image means (spatial detrending). Default: False.' },
                    recon: { type: 'bool', description: 'If True, perform Fourier reconstruction of the input time series. This does not preserve amplitudes but is useful for examining frequencies far from the low-frequency range. Default: False.' },
                    resample_original: { type: 'bool', description: 'If True, and if recon is set to True, approximate values close to the original are returned for comparison. Default: False.' },
                    nodetrendapod: { type: 'bool', description: 'If True, neither detrending nor apodization is performed. Default: False.' },
                    silent: { type: 'bool', description: 'If True, suppress print statements. Default: False.' }
                }
            },
            welch: {
                returnValues: 'power, frequency, significance',
                parameters: {
                    signal: { type: 'array', description: 'The 1D time series signal.' },
                    time: { type: 'array', description: 'The time array corresponding to the signal.' },
                    nperseg: { type: 'int', description: 'Length of each segment for analysis. Default: 256.' },
                    noverlap: { type: 'int', description: 'Number of points to overlap between segments. Default: 128.' },
                    window: { type: 'str', description: 'Type of window function used in the Welch method. Default: "hann".' },
                    siglevel: { type: 'float', description: 'Significance level for confidence intervals. Default: 0.95.' },
                    nperm: { type: 'int', description: 'Number of permutations for significance testing. Default: 1000.' },
                    silent: { type: 'bool', description: 'If True, suppress print statements. Default: False.' }
                }
            },
            emd: {
                returnValues: 'HHT_power, HHT_significance, HHT_frequency, psd_spectra, psd_significance, IMFs, IMF_significance, instantaneous_frequency',
                parameters: {
                    signal: { type: 'array', description: 'The input signal (1D).' },
                    time: { type: 'array', description: 'The time array of the signal.' },
                    siglevel: { type: 'float', description: 'Significance level for the confidence intervals. Default: 0.95.' },
                    nperm: { type: 'int', description: 'Number of permutations for significance testing. Default: 1000.' },
                    EEMD: { type: 'bool', description: 'If True, use Ensemble Empirical Mode Decomposition (EEMD) instead of Empirical Mode Decomposition (EMD). Default: False.' },
                    Welch_psd: { type: 'bool', description: 'If True, calculate Welch PSD spectra instead of FFT PSD spectra (for the psd_spectra and psd_confidence_levels). Default: False.' },
                    apod: { type: 'float', description: 'Extent of apodization edges (of a Tukey window). Default: 0.1.' },
                    pxdetrend: { type: 'int', description: 'Subtract linear trend with time per pixel. Options: 1 (simple) or 2 (advanced). Default: 2.' },
                    polyfit: { type: 'int', description: 'Degree of polynomial fit for detrending the data. If set, a polynomial fit (instead of linear) is applied. Default: None.' },
                    meantemporal: { type: 'bool', description: 'If True, apply simple temporal detrending by subtracting the mean signal from the data, skipping fitting procedures. Default: False.' },
                    meandetrend: { type: 'bool', description: 'If True, subtract the linear trend with time for the image means (spatial detrending). Default: False.' },
                    recon: { type: 'bool', description: 'If True, perform Fourier reconstruction of the input time series. This does not preserve amplitudes but is useful for examining frequencies far from the low-frequency range. Default: False.' },
                    resample_original: { type: 'bool', description: 'If True, and if recon is set to True, approximate values close to the original are returned for comparison. Default: False.' },
                    nodetrendapod: { type: 'bool', description: 'If True, neither detrending nor apodization is performed. Default: False.' },
                    silent: { type: 'bool', description: 'If True, suppress print statements. Default: False.' }
                }
            },
            komega: {
                returnValues: 'power, wavenumber, frequency, filtered_cube, spatial_fft_map, torus_map, spatial_fft_filtered_map, temporal_fft, temporal_filter, temporal_frequencies, spatial_frequencies',
                parameters: {
                    signal: { type: 'array', description: 'Input datacube, normally in the form of [x, y, t] or [t, x, y]. Note that the input datacube must have identical x and y dimensions. If not, the datacube will be cropped accordingly.' },
                    time: { type: 'array', description: 'Time array corresponding to the input datacube.' },
                    pixelsize: { type: 'float', description: 'Spatial sampling of the input datacube. If not given, it is plotted in units of "pixel".' },
                    filtering: { type: 'bool', description: 'If True, filtering is applied, and the filtered datacube (filtered_cube) is returned. Otherwise, None is returned. Default: False.' },
                    f1: { type: 'float', description: 'Optional lower (temporal) frequency to filter, in Hz.' },
                    f2: { type: 'float', description: 'Optional upper (temporal) frequency to filter, in Hz.' },
                    k1: { type: 'float', description: 'Optional lower (spatial) wavenumber to filter, in units of pixelsize^-1 (k = (2 * π) / wavelength).' },
                    k2: { type: 'float', description: 'Optional upper (spatial) wavenumber to filter, in units of pixelsize^-1.' },
                    spatial_torus: { type: 'bool', description: 'If True, makes the annulus used for spatial filtering have a Gaussian-shaped profile, useful for preventing aliasing. Default: True.' },
                    temporal_torus: { type: 'bool', description: 'If True, makes the temporal filter have a Gaussian-shaped profile, useful for preventing aliasing. Default: True.' },
                    no_spatial_filt: { type: 'bool', description: 'If True, ensures no spatial filtering is performed on the dataset (i.e., only temporal filtering is applied).' },
                    no_temporal_filt: { type: 'bool', description: 'If True, ensures no temporal filtering is performed on the dataset (i.e., only spatial filtering is applied).' },
                    silent: { type: 'bool', description: 'If True, suppresses the k-ω diagram plot.' },
                    smooth: { type: 'bool', description: 'If True, power is smoothed. Default: True.' },
                    mode: { type: 'int', description: 'Output power mode: 0 = log10(power) (default), 1 = linear power, 2 = sqrt(power) = amplitude.' },
                    processing_maps: { type: 'bool', description: 'If True, the function returns the processing maps (spatial_fft_map, torus_map, spatial_fft_filtered_map, temporal_fft, temporal_filter, temporal_frequencies, spatial_frequencies). Otherwise, they are all returned as None. Default: False.' }
                }
            },
            pod: {
                returnValues: 'pod_results',
                parameters: {
                    signal: { type: 'array', description: '3D data cube with shape (time, x, y) or similar.' },
                    time: { type: 'array', description: '1D array representing the time points for each time step in the data.' },
                    num_modes: { type: 'int, optional', description: 'Number of top modes to compute. Default is None (all modes).' },
                    num_top_frequencies: { type: 'int, optional', description: 'Number of top frequencies to consider. Default is None (all frequencies).' },
                    top_frequencies: { type: 'list, optional', description: 'List of top frequencies to consider. Default is None.' },
                    num_cumulative_modes: { type: 'int, optional', description: 'Number of cumulative modes to consider. Default is None (all modes).' },
                    welch_nperseg: { type: 'int, optional', description: "Number of samples per segment for Welch's method. Default is 150." },
                    welch_noverlap: { type: 'int, optional', description: "Number of overlapping samples for Welch's method. Default is 25." },
                    welch_nfft: { type: 'int, optional', description: 'Number of points for the FFT. Default is 2^14.' },
                    welch_fs: { type: 'int, optional', description: 'Sampling frequency for the data. Default is 2.' },
                    nperm: { type: 'int, optional', description: 'Number of permutations for significance testing. Default is 1000.' },
                    siglevel: { type: 'float, optional', description: 'Significance level for the Welch spectrum. Default is 0.95.' },
                    timestep_to_reconstruct: { type: 'int, optional', description: 'Timestep of the datacube to reconstruct using the top modes. Default is 0.' },
                    num_modes_reconstruct: { type: 'int, optional', description: 'Number of modes to use for reconstruction. Default is None (all modes).' },
                    spod: { type: 'bool, optional', description: 'If True, perform Spectral Proper Orthogonal Decomposition (SPOD) analysis. Default is False.' },
                    spod_filter_size: { type: 'int, optional', description: 'Filter size for SPOD analysis. Default is None.' },
                    spod_num_modes: { type: 'int, optional', description: 'Number of SPOD modes to compute. Default is None.' },
                    print_results: { type: 'bool, optional', description: 'If True, print a summary of results. Default is True.' }
                }
            },
            dominantfreq: {
                returnValues: 'power, frequency, significance',
                parameters: {
                    signal: { type: 'array', description: 'Input signal array (1D or 2D).' },
                    time: { type: 'array', description: 'Time array corresponding to the signal.' },
                    method: { type: 'string', description: 'Analysis method (e.g., fft, wavelet, etc.).' },
                    kwargs: { type: 'object', description: 'Additional optional parameters for customization.' }
                }
            }
        },
        cross_correlation: {
            wavelet: {
                returnValues: 'cross_power, cross_period, cross_sig, cross_coi, coherence, coh_period, coh_sig, coh_coi, phase_angle',
                parameters: {
                    data1: { type: 'array', description: 'The first 1D time series signal.' },
                    data2: { type: 'array', description: 'The second 1D time series signal.' },
                    time: { type: 'array', description: 'The time array corresponding to the signals.' },
                    siglevel: { type: 'float', description: 'Significance level for the confidence intervals. Default: 0.95.' },
                    nperm: { type: 'int', description: 'Number of permutations for significance testing. Default: 1000.' },
                    mother: { type: 'str', description: 'The mother wavelet function to use. Default: "morlet".' },
                    GWS: { type: 'bool', description: 'If True, calculate the Global Wavelet Spectrum. Default: False.' },
                    RGWS: { type: 'bool', description: 'If True, calculate the Refined Global Wavelet Spectrum (time-integrated power, excluding COI and insignificant areas). Default: False.' },
                    dj: { type: 'float', description: 'Scale spacing. Smaller values result in better scale resolution but slower calculations. Default: 0.025.' },
                    s0: { type: 'float', description: 'Initial (smallest) scale of the wavelet. Default: 2 * dt.' },
                    J: { type: 'int', description: 'Number of scales minus one. Scales range from s0 up to s0 * 2**(J * dj), giving a total of (J + 1) scales. Default: (log2(N * dt / s0)) / dj.' },
                    lag1: { type: 'float', description: 'Lag-1 autocorrelation. Default: 0.0.' },
                    apod: { type: 'float', description: 'Extent of apodization edges (of a Tukey window). Default: 0.1.' },
                    pxdetrend: { type: 'int', description: 'Subtract linear trend with time per pixel. Options: 1 (simple) or 2 (advanced). Default: 2.' },
                    polyfit: { type: 'int', description: 'Degree of polynomial fit for detrending the data. If set, a polynomial fit (instead of linear) is applied. Default: None.' },
                    meantemporal: { type: 'bool', description: 'If True, apply simple temporal detrending by subtracting the mean signal from the data, skipping fitting procedures. Default: False.' },
                    meandetrend: { type: 'bool', description: 'If True, subtract the linear trend with time for the image means (spatial detrending). Default: False.' },
                    recon: { type: 'bool', description: 'If True, perform Fourier reconstruction of the input time series. This does not preserve amplitudes but is useful for examining frequencies far from the low-frequency range. Default: False.' },
                    resample_original: { type: 'bool', description: 'If True, and if recon is set to True, approximate values close to the original are returned for comparison. Default: False.' },
                    nodetrendapod: { type: 'bool', description: 'If True, neither detrending nor apodization is performed. Default: False.' },
                    silent: { type: 'bool', description: 'If True, suppress print statements. Default: False.' }
                }
            },
            welch: {
                returnValues: 'frequency, cospectrum, phase_angle, power_data1, power_data2, frequency_coherence, coherence',
                parameters: {
                    data1: { type: 'array', description: 'The first 1D time series signal.' },
                    data2: { type: 'array', description: 'The second 1D time series signal.' },
                    time: { type: 'array', description: 'The time array corresponding to the signals.' },
                    nperseg: { type: 'int', description: 'Length of each segment for analysis. Default: 256.' },
                    noverlap: { type: 'int', description: 'Number of points to overlap between segments. Default: 128.' },
                    window: { type: 'str', description: 'Type of window function used in the Welch method. Default: "hann".' },
                    siglevel: { type: 'float', description: 'Significance level for confidence intervals. Default: 0.95.' },
                    nperm: { type: 'int', description: 'Number of permutations for significance testing. Default: 1000.' },
                    silent: { type: 'bool', description: 'If True, suppress print statements. Default: False.' }
                }
            },
            fft: {
                returnValues: 'frequency, cospectrum, phase_angle, power_data1, power_data2, frequency_coherence, coherence',
                warning: 'Selecting FFT for cross-correlation automatically runs Welch instead, to improve noise handling by averaging over segments. Consider adjusting "nperseg" to tune resolution.',
                parameters: {
                    data1: { type: 'array', description: 'The first 1D time series signal.' },
                    data2: { type: 'array', description: 'The second 1D time series signal.' },
                    time: { type: 'array', description: 'The time array corresponding to the signals.' },
                    nperseg: { type: 'int', description: 'Length of each segment for analysis. Default: 256.' },
                    noverlap: { type: 'int', description: 'Number of points to overlap between segments. Default: 128.' },
                    window: { type: 'str', description: 'Type of window function used in the Welch method. Default: "hann".' },
                    siglevel: { type: 'float', description: 'Significance level for confidence intervals. Default: 0.95.' },
                    nperm: { type: 'int', description: 'Number of permutations for significance testing. Default: 1000.' },
                    silent: { type: 'bool', description: 'If True, suppress print statements. Default: False.' }
                }
            }
        }
    };
    const categoryDropdown = document.getElementById('category');
    const datatypeDropdown = document.getElementById('datatype');
    const analysisMethodDropdown = document.getElementById('analysisMethod');
    const subMethodDropdown = document.getElementById('subMethod');
    const subMethodLabel = document.getElementById('subMethodLabel');
    const outputContainer = document.getElementById('outputContainer');
    const callingSequence = document.getElementById('callingSequence');
    const parameterTableBody = document.getElementById('parameterTableBody');
    function resetDropdown(dropdown, placeholder = "Select ...") {
        dropdown.innerHTML = `<option value="">${placeholder}</option>`;
        dropdown.disabled = true;
    }
    function hideOutput() {
        outputContainer.style.display = 'none';
    }
    function clearOutput() {
        callingSequence.innerHTML = '';
        parameterTableBody.innerHTML = '';
        hideOutput();
    }
    document.addEventListener('DOMContentLoaded', () => {
        updateOutput();
        // Attach event listeners
        categoryDropdown.addEventListener('change', () => {
            const category = categoryDropdown.value;
            resetDropdown(datatypeDropdown, "Select Data Type");
            resetDropdown(analysisMethodDropdown, "Select Method");
            resetDropdown(subMethodDropdown, "Select Sub-method");
            subMethodDropdown.style.display = 'none';
            subMethodLabel.style.display = 'none';
            clearOutput();
            if (category) {
                datatypeDropdown.disabled = false;
                if (category === 'a') {
                    datatypeDropdown.innerHTML += `
                        <option value="1">1D Signal</option>
                        <option value="2">3D Datacube</option>`;
                } else if (category === 'b') {
                    datatypeDropdown.innerHTML += `<option value="1">1D Signal</option>`;
                }
            }
            updateOutput();
        });
        datatypeDropdown.addEventListener('change', () => {
            const category = categoryDropdown.value;
            const datatype = datatypeDropdown.value;
            resetDropdown(analysisMethodDropdown, "Select Method");
            resetDropdown(subMethodDropdown, "Select Sub-method");
            subMethodDropdown.style.display = 'none';
            subMethodLabel.style.display = 'none';
            clearOutput();
            if (datatype) {
                analysisMethodDropdown.disabled = false;
                if (category === 'a' && datatype === '1') {
                    analysisMethodDropdown.innerHTML += `
                        <option value="fft">FFT</option>
                        <option value="wavelet">Wavelet</option>
                        <option value="lombscargle">Lomb-Scargle</option>
                        <option value="welch">Welch</option>
                        <option value="emd">EMD</option>`;
                } else if (category === 'a' && datatype === '2') {
                    analysisMethodDropdown.innerHTML += `
                        <option value="komega">k-omega</option>
                        <option value="pod">POD</option>
                        <option value="dominantfreq">Dominant Freq / Mean Power Spectrum</option>`;
                } else if (category === 'b') {
                    analysisMethodDropdown.innerHTML += `
                        <option value="fft">FFT</option>
                        <option value="wavelet">Wavelet</option>
                        <option value="welch">Welch</option>`;
                }
            }
            updateOutput();
        });
        analysisMethodDropdown.addEventListener('change', () => {
            const category = categoryDropdown.value;
            const datatype = datatypeDropdown.value;
            const analysisMethod = analysisMethodDropdown.value;
            resetDropdown(subMethodDropdown, "Select Sub-method");
            subMethodDropdown.style.display = 'none';
            subMethodLabel.style.display = 'none';
            subMethodDropdown.disabled = true;
            clearOutput();
            if (
                category === 'a' &&
                datatype === '2' &&
                analysisMethod === 'dominantfreq'
            ) {
                subMethodDropdown.style.display = 'inline-block';
                subMethodLabel.style.display = 'inline-block';
                subMethodDropdown.disabled = false;
                subMethodDropdown.innerHTML = `
                    <option value="">Select Sub-method</option>
                    <option value="fft">FFT</option>
                    <option value="wavelet">Wavelet</option>
                    <option value="lombscargle">Lomb-Scargle</option>
                    <option value="welch">Welch</option>`;
            }
            updateOutput();
        });
        subMethodDropdown.addEventListener('change', updateOutput);
        updateOutput();
    });
    // Update Output Container
    function updateOutput() {
        const category = categoryDropdown.value;
        const datatype = datatypeDropdown.value;
        const analysisMethod = analysisMethodDropdown.value;
        const subMethod = subMethodDropdown.value;
        const message = document.getElementById('dropdownMessage');      
        // Hide output and show a message if selections are incomplete
        if (!category || !datatype || !analysisMethod || (subMethodDropdown.style.display === 'inline-block' && !subMethod)) {
            message.style.display = 'block';
            hideOutput();
            return;
        }      
        // Hide the message when all selections are made
        message.style.display = 'none';
        // Retrieve returnValues based on category and analysisMethod
        let returnValues = '';
        if (category === 'a') {
            returnValues = parameters.single_series[analysisMethod]?.returnValues || 'No return values available';
        } else if (category === 'b') {
            returnValues = parameters.cross_correlation[analysisMethod]?.returnValues || 'No return values available';
        }
        // Construct the command string
        let command1 = '';
        if (category === 'a' && datatype === '1') {
            command1 = `${returnValues} = WaLSAtools(signal=INPUT_DATA, time=TIME_ARRAY, method='${analysisMethod}', **kwargs)`;
        } else if (category === 'b') {
            command1 = `${returnValues} = WaLSAtools(data1=INPUT_DATA1, data2=INPUT_DATA2, time=TIME_ARRAY, method='${analysisMethod}', **kwargs)`;
        } else if (category === 'a' && datatype === '2') {
            if (analysisMethod === 'dominantfreq') {
                retValues = 'dominant_frequency, mean_power, frequency, power_map'
                command1 = `${retValues} = WaLSAtools(signal=INPUT_DATA, time=TIME_ARRAY, averagedpower=True, dominantfreq=True, method='${subMethod}', **kwargs)`;
            } else {
                if (analysisMethod === 'komega') {
                    analysisMethodout = 'k-omega';
                } else {
                    analysisMethodout = analysisMethod;
                }
                command1 = `${returnValues} = WaLSAtools(signal=INPUT_DATA, time=TIME_ARRAY, method='${analysisMethodout}', **kwargs)`;
            }
        }
        const command = `
            <div style="display: flex; align-items: baseline;">
                <span style="color: #222; min-width: 4ch; margin: 0 !important; line-height: 1.5;">>>> </span>
                <pre style="
                    white-space: pre-wrap; 
                    word-wrap: break-word;  
                    color: #01016D; 
                    margin: 0 !important;
                    line-height: 1.5;
                ">${command1}</pre>
            </div>
        `;
        // Update calling sequence and parameter table
        callingSequence.innerHTML = command;
        const warningBox = document.getElementById('warningNote');
        const warningMessage = document.getElementById('warningMessage');
        if (category === 'b' && analysisMethod === 'fft') {
            warningMessage.innerText = parameters.cross_correlation.fft.warning;
            warningBox.style.display = 'block';
        } else {
            warningBox.style.display = 'none';
            warningMessage.innerText = '';
        }
        if (analysisMethod === 'dominantfreq') {
            updateParameterTable(subMethod);
        } else {
            updateParameterTable(analysisMethod);
        }
        outputContainer.style.display = 'block';
    }
    function updateParameterTable(datatype) {
        parameterTableBody.innerHTML = '';
        const paramData =
            parameters.single_series[datatype]?.parameters ||
            parameters.cross_correlation[datatype]?.parameters;
        if (!paramData) {
            parameterTableBody.innerHTML = `
                <tr>
                    <td colspan="3" style="text-align: center;">No parameters available.</td>
                </tr>`;
            return;
        }
        Object.entries(paramData).forEach(([key, value]) => {
            parameterTableBody.innerHTML += `
                <tr>
                    <td>${key}</td>
                    <td>${value.type}</td>
                    <td>${value.description}</td>
                </tr>`;
        });
    }    
    document.addEventListener('DOMContentLoaded', () => {
        updateOutput();
    });
</script>

??? source-code "Source code"
    ``` python linenums="1"
    --8<-- "codes/python/WaLSAtools/WaLSAtools.py"
    ```
