---
template: main.html
title: Installation
---

# Installation

## Install WaLSAtools :material-package-variant-closed:

This page provides a direct and simple installation guide for **WaLSAtools**.    
For a complete beginner-friendly setup from scratch, see the [Beginner's Guide](https://walsa.tools/python/beginner-friendly-guide/).

---

## Prerequisites :material-clipboard-text:

Before installing WaLSAtools, ensure you have the following:

- **Python 3.11 or later** installed. pip selects dependencies for your Python version.
  👉 You can download Python from [python.org](https://www.python.org/downloads/) if needed.

- **pip** (Python package installer) available.  
  👉 You can check by running:
    ```bash
    pip --version
    ```
    If you do not have pip, follow the [official installation guide](https://pip.pypa.io/en/stable/installation/).

    !!! walsa-info "Version used in the examples"
        The beginner’s guide uses **Python 3.12.8**.

- **Virtual Environment** (Recommended):
    To avoid conflicts with existing packages (especially due to major changes in packages like NumPy 2.0+),
    we strongly recommend installing WaLSAtools inside a fresh virtual environment.
    This keeps WaLSAtools and its dependencies isolated from your system-wide Python or other projects.
    If you are new to virtual environments, you can refer to our [Beginner's Guide](https://walsa.tools/python/beginner-friendly-guide/) for a full setup tutorial.

## Installation via pip (Recommended) :material-cloud-download-outline:

To install the latest stable release from PyPI (Python Package Index), simply run:
```bash
pip install WaLSAtools
```
This will automatically download and install WaLSAtools and its required dependencies into your environment.

On Python 3.15 and later, wavelet analysis uses the SciPy FFT backend if pyFFTW is absent. Wavelet coherence values can differ between the two backends. For direct calls to `cwt` with a non-power-of-two input length, keep the default `pad=True` when using SciPy.

!!! walsa-tip "Tip: Upgrade pip first"
    Before installing any package, it is good practice to upgrade pip:
    ```
    pip install --upgrade pip
    ```

## Installation from Source (Optional) :material-git:

If you prefer installing the development version directly from GitHub:

1. Clone (or Download) the WaLSAtools GitHub repository:
    ```bash
    git clone https://github.com/WaLSAteam/WaLSAtools.git
    ```

2. Navigate to the Python codes directory:
    ```bash
    cd WaLSAtools/codes/python/
    ```

3. Install using pip:
    ```bash
    python -m pip install .
    ```
    pip installs the build tools declared in `pyproject.toml` in a temporary environment.

This method ensures you have access to the latest updates and examples.

!!! walsa-info "Cloning the GitHub repository"
    Even if you install WaLSAtools via PyPI, we strongly recommend also cloning the WaLSAtools GitHub repository. This gives you full access to the source code, worked examples, and documentation files, and ensures you can easily explore, customize, or [contribute](https://walsa.tools/contribution/) to the project.

## Verifying the Installation :material-check-decagram:

After installation, verify that WaLSAtools is properly installed.    
You can check inside the terminal by launching Python and typing:  
```python
from WaLSAtools import WaLSAtools
WaLSAtools
```
If installed correctly, WaLSAtools’ interactive welcome menu will appear 🎉

## Need Help? :material-lifebuoy:

If you encounter any issues during installation, first check our [Troubleshooting](https://walsa.tools/python/troubleshooting/) and/or [Beginner's Guide](https://walsa.tools/python/beginner-friendly-guide/) pages.
If you still need help, browse the topics in our [GitHub discussion](https://github.com/WaLSAteam/WaLSAtools/discussions) — and feel free to post your question if your issue has not yet been addressed.
 
 <br>
 