---
template: main.html
title: Troubleshooting
---

# Troubleshooting :material-bug-outline:

This page provides detailed solutions to common issues encountered during the installation and usage of **WaLSAtools**.

If you experience a problem not listed here, we recommend:

- First, carefully reviewing the [Beginner-Friendly Guide](https://walsa.tools/python/beginner-friendly-guide/) and this Troubleshooting page. 
- If the issue remains unresolved, posting your question in our [GitHub Discussions](https://github.com/WaLSAteam/WaLSAtools/discussions) area. 
- If the problem appears to be a confirmed bug, submitting a [GitHub Issue](https://github.com/WaLSAteam/WaLSAtools/issues) (see the [Contribution](https://walsa.tools/contribution/) guideline for details). 

!!! walsa-issues "New here?"
    Many installation and setup problems are already covered step-by-step in the [Beginner-Friendly Guide](https://walsa.tools/python/beginner-friendly-guide/). We highly recommend checking it if you haven't already.

-----

## Installation Issues :material-download-network:

### 1. `pip install WaLSAtools` fails

**Possible reasons and solutions:**

- **Outdated `pip`**:  
  Upgrade `pip` to the latest version:

    ```bash
    pip install --upgrade pip
    ```

- **Python version mismatch**:  
  Use `python --version` to check that you are running Python 3.11 or later. If `python setup.py install` reports missing `setuptools`, install from the `codes/python` directory using `python -m pip install .`.

- **Virtual environment not used**:  
  To avoid conflicts (especially with libraries like **NumPy 2.0+**), we strongly recommend creating a clean [virtual environment](https://walsa.tools/python/beginner-friendly-guide/#step-2-create-a-virtual-environment).

	!!! walsa-warning "Conda Users"
		If you are using [`Anaconda`](https://www.anaconda.com) or [`Miniconda`](https://www.anaconda.com/docs/getting-started/miniconda/main), be aware that some pre-installed packages (like numpy, scipy, matplotlib) may not match the latest stable versions or may be compiled differently, leading to compatibility problems.
		Solution:

		- Create a fresh `Conda` environment without preinstalled packages, or
		- Prefer a clean virtual environment created with `pyenv` or `venv` when working with **WaLSAtools**.
		- Always manually upgrade critical packages like `pip`, `setuptools`, and `wheel` after creating a new environment.

- **Firewall or proxy blocking installation**:  
  If you are behind a firewall or corporate proxy, configure pip accordingly:

    ```bash
    pip install WaLSAtools --proxy="your_proxy_address"
    ```

- **Temporary server issue**:  
  If installation fails from PyPI, wait a few minutes and try again.

### 2. Problems related to NumPy or compiled dependencies

**Specific Issue:**  
Some users reported installation failures due to **NumPy ABI (Application Binary Interface)** conflicts, especially when using **NumPy 2.0.0+**. See the [GitHub Discussion here](https://github.com/WaLSAteam/WaLSAtools/discussions/4#discussioncomment-12944031) for more details.

**Solution:**  
If you encounter NumPy-related errors:

- Ensure you are installing inside a **fresh virtual environment**.
- Manually install a compatible NumPy version:

    ```bash
    pip install "numpy<2.0"
    ```

- Then reinstall WaLSAtools:

    ```bash
    pip install WaLSAtools
    ```

!!! walsa-important "Important"
    Avoid using pre-installed system Python or OS-level packages (e.g., `/usr/bin/python3`) for WaLSAtools. Always prefer user-managed virtual environments.

### 3. Installation from source fails

**Checklist:**

- **Clone the repository properly**:

    ```bash
    git clone https://github.com/WaLSAteam/WaLSAtools.git
    cd WaLSAtools/codes/python
    ```

- **Install using pip inside the directory**:

    ```bash
    pip install .
    ```

- **Ensure build tools are installed**:
  - macOS: Install Xcode and command-line tools:

    ```bash
    xcode-select --install
    ```

  - Linux: Install build essentials:

    ```bash
    sudo apt install build-essential
    ```

  - Windows: Use pre-built binaries and ensure `pip` is up to date.

- **Permission errors**:  
  Never use `sudo pip install`. Always use virtual environments to avoid permission problems.

-----

## Usage Issues :material-cog-outline:

### 1. WaLSAtools interactive interface does not launch

**Solutions:**

- **Confirm environment activation**:  
  Ensure you are working inside the environment where WaLSAtools was installed.

- **Correct import syntax**:  
  Start a Python session and enter:

    ```python
    from WaLSAtools import WaLSAtools
    WaLSAtools
    ```

- **Notebook environment issue**:  
  Inside a Jupyter notebook, make sure you select the correct Python **kernel** linked to your WaLSAtools environment.

- **Conflicting libraries**:  
  Conflicts with old versions of packages like `matplotlib`, `scipy`, or `numpy` can prevent proper functioning. Use a clean environment.

### 2. Problems related to data input or unit handling

**Specific Issue:**  
One user reported crashes when passing data arrays with attached **units** (e.g., astropy `Quantity` arrays for time or signal). See the [specific issue and solution here](https://github.com/WaLSAteam/WaLSAtools/issues/3#issuecomment-2829730034) for more details.

**Solution:**  
**WaLSAtools expects raw `numpy.ndarray` inputs** without units.

- Remove units before passing:

    ```python
    data = data_with_units.value
    ```

- If using `astropy.Time` objects for time axes, convert them to seconds manually.

!!! walsa-note "Future Improvements"
    Support for native `Quantity` and `Time` inputs is planned for a future WaLSAtools release.

### 3. WaLSAtools behaves differently across environments

Different behavior (e.g., function outputs, warnings) can occur due to:

- Different dependency versions (e.g., `numpy`, `scipy`, `matplotlib`).
- Mixed installations (system vs. virtual environment).

**Solution:**  
Check installed versions and compare them with those listed in [`requirements.txt`](https://github.com/WaLSAteam/WaLSAtools/blob/main/codes/python/requirements.txt) under `WaLSAtools/codes/python`:
```python
import astropy
import IPython
import ipywidgets
import matplotlib
import numba
import numpy
import pyfftw
import scipy
import setuptools
import skimage
import tqdm

modules = {
    'astropy': astropy,
    'ipython': IPython,
    'ipywidgets': ipywidgets,
    'matplotlib': matplotlib,
    'numba': numba,
    'numpy': numpy,
    'pyFFTW': pyfftw,
    'scipy': scipy,
    'setuptools': setuptools,
    'scikit-image': skimage,
    'tqdm': tqdm
}

for name, module in modules.items():
    print(f"{name}: {module.__version__}")
```

Always install WaLSAtools inside a fresh environment with the recommended package versions (that is done automatically when installing WaLSAtools).

-----

## Additional Resources :material-book-open-outline:

- [Beginner-Friendly Guide](https://walsa.tools/python/beginner-friendly-guide/)
- [Installation Guide](https://walsa.tools/python/installation/)
- [Contribution Guidelines](https://walsa.tools/contribution/)
- [GitHub Discussions](https://github.com/WaLSAteam/WaLSAtools/discussions)

!!! walsa-issues "Still stuck?"
    If the solutions above don't work, search GitHub Issues or Discussions.  
    If your issue is not already reported, feel free to post a new discussion or bug report!

-----

## Final Notes :material-lightbulb-outline:

- Always read error messages carefully — they usually point directly to the problem.
- Keep your `pip`, `setuptools`, and `wheel` packages updated.
- Prefer working inside clean virtual environments.
- If you fix an issue not yet listed here, please consider contributing to this page!

<br>
