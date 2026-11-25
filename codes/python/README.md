<div align="center">
  <a href="https://WaLSA.tools" target="_blank">
    <img align="left" src="https://walsa.team/images/WaLSAtools_logo.svg" alt="WaLSAtools" width="350" height="auto">
  </a>
</div>

<br><br><br>

# WaLSAtools &ndash; Wave Analysis Tools

<p align="left">
  <a href="https://pypi.org/project/WaLSAtools/"><img src="https://img.shields.io/pypi/v/WaLSAtools.svg" alt="PyPI"></a>
  <a href="https://walsa.team" target="_blank"><img src="https://img.shields.io/badge/powered%20by-WaLSA%20Team-000d1a"></a>
  <a href="https://walsa.tools/license" target="_blank"><img src="https://img.shields.io/badge/license-Apache%202.0-green"></a>
  <a href="https://doi.org/10.5281/zenodo.14978610" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14978610.svg" alt="DOI"></a>
</p>

**WaLSAtools** is an open-source library for fundamental and advanced wave analysis in 1D signals, images, and multi-dimensional datasets. Designed as a one-stop resource, it empowers researchers across disciplines to extract insights from oscillatory phenomena. **WaLSAtools** is envisioned as a community-driven starting point — inviting collaboration to refine techniques, develop new methods, and advance discoveries wherever waves and fluctuations play a fundamental role. It encourages transparent, reproducible science by providing validated techniques, detailed documentation, and worked examples.

It currently provides a comprehensive set of analysis techniques, including:

- Fast Fourier Transform (FFT)
- Wavelet Analysis
- Lomb-Scargle Periodogram
- Welch Power Spectral Density
- Empirical Mode Decomposition (EMD)
- Hilbert and Hilbert-Huang Transforms
- k-ω (k-omega) Analysis
- Proper Orthogonal Decomposition (POD)
- Cross-Spectral Analysis

**WaLSAtools** features an interactive interface for both Python terminals and Jupyter notebooks, simplifying the wave analysis workflow.

For full documentation, see **[https://WaLSA.tools](https://WaLSA.tools)**.

# Installation

We strongly recommend installing **WaLSAtools** inside a new virtual environment to avoid dependency conflicts, especially with pre-installed packages.

For a full setup guide, see the [Beginner-Friendly Guide](https://walsa.tools/python/beginner-friendly-guide/).

## Install via PyPI

```bash
pip install --upgrade pip
pip install WaLSAtools
```

## Install from Source (GitHub)

```bash
git clone https://github.com/WaLSAteam/WaLSAtools.git
cd WaLSAtools/codes/python/
pip install .
```

For detailed instructions and troubleshooting tips, refer to:  
- [Installation Guide](https://walsa.tools/python/installation/)  
- [Troubleshooting](https://walsa.tools/python/troubleshooting/)

# Interactive Usage

After installing, launch the interactive interface by running:

```python
from WaLSAtools import WaLSAtools
WaLSAtools
```

This will open an interactive menu where you can:

- Select an analysis category.
- Choose your data type (1D or 3D).
- Pick an analysis method (FFT, wavelet, EMD, k-omega, POD, etc.).
- View customised calling sequences and parameter hints.

Interactive usage is supported both in standard Python terminals and Jupyter notebooks.

# Requirements

WaLSAtools requires Python ≥ 3.11 and automatically installs the following core dependencies:

`astropy`, `IPython`, `ipywidgets`, `matplotlib`, `numba`, `numpy`, `scipy`, `setuptools`, `scikit-image`, `tqdm`.

pip selects the required dependencies for your Python version. On Python 3.15 and later, pyFFTW is optional and wavelet analysis uses the existing SciPy FFT backend when pyFFTW is absent.

# License

WaLSAtools is licensed under the [Apache License, Version 2.0](https://walsa.tools/license/).

# Citation

If you use WaLSAtools in your research, please cite:

> Jafarzadeh, S., Jess, D. B., Stangalini, M., et al. 2025, *Nature Reviews Methods Primers*, 5, 21   
> [https://www.nature.com/articles/s43586-025-00392-0](https://www.nature.com/articles/s43586-025-00392-0)

Zenodo DOI: [10.5281/zenodo.14978610](https://doi.org/10.5281/zenodo.14978610)

# Contributing

We welcome contributions, suggestions, and feedback.

Please check the [Contribution Guide](https://walsa.tools/contribution/), or join our discussions via [GitHub Discussions](https://github.com/WaLSAteam/WaLSAtools/discussions).
