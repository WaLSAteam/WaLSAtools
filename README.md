<div style="display: block; margin-bottom: 20px;">
  <a href="https://WaLSA.tools" target="_blank">
    <img align="left" src="docs/images/WaLSAtools_logo.svg" alt="WaLSAtools Documentation" width="400" height="auto">
  </a>
</div>

<br><br><br>

# WaLSAtools &ndash; Wave Analysis Tools

<p align="left">
    <a href="#"><img src="https://img.shields.io/badge/WaLSAtools-v1.0.0-0066cc"></a> 
    <a href="https://walsa.team" target="_blank"><img src="https://img.shields.io/badge/powered%20by-WaLSA%20Team-000d1a"></a>
    <a href="https://walsa.tools/license"><img src="https://img.shields.io/badge/license-Apache%202.0-green"></a>
    <a href="https://doi.org/10.5281/zenodo.14978610"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14978610.svg" alt="DOI"></a>
    <a href="https://github.com/WaLSAteam/WaLSAtools/blob/main/CONTRIBUTORS.md"><img src="https://img.shields.io/badge/contributors-28-gold" alt="Contributors"></a>
    <a href="https://pypi.org/project/WaLSAtools/"><img src="https://img.shields.io/pypi/v/WaLSAtools.svg" alt="PyPI Latest Release"></a>
    <a href="https://github.com/WaLSAteam/WaLSAtools/actions/workflows/ci.yml"><img src="https://github.com/WaLSAteam/WaLSAtools/workflows/docs/badge.svg"></a>
</p>

**WaLSAtools** is an open-source library for analysing a wide variety of wave phenomena in time series data &ndash; including 1D signals, images, and multi-dimensional datasets. It provides tools to extract meaningful insights from complex datasets and is applicable across diverse fields, including astrophysics, engineering, life, physical and environmental sciences, and biomedical studies, among others. The library is continuously expanding with new features and functionalities, ensuring it remains a valuable resource for wave analysis.

The core of WaLSAtools is built upon [Python](https://www.python.org). This ensures accessibility and ease of use for a broad audience. We are actively developing versions in other popular languages to further enhance accessibility, enabling researchers from various backgrounds to leverage the power of WaLSAtools for their wave analysis needs. Currently, WaLSAtools is partially implemented in IDL, with plans to expand its functionality and extend to other programming languages in the future.

WaLSAtools provides a suite of both fundamental and advanced tools, but it remains the user's responsibility to choose the method that best fits the nature of their dataset and the scientific questions being addressed. Selecting the appropriate analysis method is essential for ensuring reliable and scientifically valid results. The use of unsuitable or overly simplified techniques &ndash; without consideration of the data's properties or the research goals &ndash; can lead to incomplete or incorrect conclusions, and potentially to misinterpretation. This principle is central to our accompanying [Primer](https://www.nature.com/articles/s43586-025-00392-0), which emphasises the importance of methodological awareness in wave analysis across all disciplines.

Developed by the [WaLSA Team](https://WaLSA.team), WaLSAtools was initially inspired by the intricate wave dynamics observed in the Sun's atmosphere. However, its applications extend far beyond solar physics, offering a versatile toolkit for anyone working with oscillatory signals.

WaLSAtools promotes reproducibility and transparency in wave analysis. Its robust implementations of validated techniques ensure consistent and trustworthy results, empowering researchers to delve deeper into the complexities of their data. Through its interactive interface, WaLSAtools guides users through the analysis process, providing the necessary information and tools to perform various types of wave analysis with ease.

This repository is associated with a primer article titled *"Wave analysis tools"* in **[Nature Reviews Methods Primers](https://www.nature.com/articles/s43586-025-00392-0)** (NRMP), showcasing its capabilities through detailed analyses of synthetic datasets. The `examples/Worked_examples__NRMP` directories (for both Python and IDL) contain reproducible codes for generating all figures presented in the NRMP article, serving as a practical guide for applying WaLSAtools to real-world analyses.

**Free read-only access to the Primer**: [https://WaLSA.tools/nrmp](https://WaLSA.tools/nrmp)

**Supplementary Information**: [https://WaLSA.tools/nrmp-si](https://WaLSA.tools/nrmp-si)

---

⭐ **If you find WaLSAtools useful**, please consider giving it a star (click the ⭐ at the top right of this page).  
It helps others discover the project and supports the growth of an open, cross-disciplinary wave analysis community.

🔔 To stay updated on new releases and features, click **“Watch → Custom → Releases”** at the top of the page.

---

## **Key Features**

* **Wide Range of Wave Analysis Techniques:**  From foundational methods like FFT and wavelet analysis to advanced techniques such as EMD, k-ω, and POD analysis.
* **Cross-Disciplinary Applicability:**  Suitable for signal processing, oscillation studies, and multi-dimensional analysis in various fields.
* **Interactive Interfaces:** Simplified workflows through interactive menus for both Python and IDL.
* **Open Science Principles:** Promotes reproducibility and transparency in data analysis. 


## **Documentation**

<a href="https://WaLSA.tools" target="_blank"><img align="right" src="docs/images/misc/WaLSAtools_documentation_screenshot.png" alt="WaLSAtools Documentation" width="485" height="auto" /></a>

Complete documentation, including installation guides, method descriptions, and usage examples, is available online:

**[WaLSAtools Documentation](https://WaLSA.tools)**

The documentation includes:
- Step-by-step installation instructions.
- Descriptions of implemented methods.
- Examples applied to synthetic datasets.


## **Repository Structure**

<pre>
WaLSAtools/
├── codes/
│   ├── python/                      # Python implementation of WaLSAtools
│   │   ├── WaLSAtools/              # Core library
│   │   ├── setup.py                 # Setup script for Python
│   │   └── README.md                # Python-specific README
│   ├── idl/                         # IDL implementation of WaLSAtools
│   │   ├── WaLSAtools/              # Core library
│   │   ├── setup.pro                # Setup script for IDL
│   │   └── README.md                # IDL-specific README
├── docs/                            # Documentation for WaLSAtools
├── examples/                        # Worked examples directory
│   ├── python/                      # Python-specific examples
│   │   └── Worked_examples__NRMP/
│   ├── idl/                         # IDL-specific examples
│   │   └── Worked_examples__NRMP/
├── LICENSE                          # License information
├── NOTICE                           # Copyright Notice
└── README.md                        # Main repository README
</pre>

## **Installation**

Refer to the `README.md` files in the `codes/python` and `codes/idl` directories for language-specific installation instructions. Further details are in the [online documentation](https://WaLSA.tools).

**Quick Start:**

* **Python (via PyPI):** 
  ```bash
  pip install WaLSAtools
  ```
* **Python (from source):** 
  ```bash
  git clone https://github.com/WaLSAteam/WaLSAtools.git
  cd WaLSAtools/codes/python/
  pip install .
  ```
* **IDL:** 
  ```bash
  git clone https://github.com/WaLSAteam/WaLSAtools.git
  cd WaLSAtools/codes/idl/
  .run setup.pro
  ```

## **Interactive Usage**

* **Python:** Type `from WaLSAtools import WaLSAtools; WaLSAtools` in a terminal (in Python) or in a Jupyter notebook.
* **IDL:** Run the command `WaLSAtools` (in IDL).

The interface guides you through selecting an analysis category, data type, and method, providing hints on calling sequences and parameters.


## **Citing WaLSAtools**

If you use WaLSAtools in your research, please consider citing:

> [Jafarzadeh, S., Jess, D. B., Stangalini, M. et al. 2025, Nature Reviews Methods Primers, 5, 21](https://www.nature.com/articles/s43586-025-00392-0)

> [Jafarzadeh, S., Jess, D. B., Stangalini, M., et al. 2025, WaLSAtools, WaLSA Team (GitHub) / Zenodo, doi: 10.5281/zenodo.14978610](https://doi.org/10.5281/zenodo.14978610)


## **Contributing**

We welcome contributions from researchers and developers across all disciplines.   
To learn how to get involved, report issues, or propose improvements, please refer to our [Contribution Guide](https://walsa.tools/contribution/).

You can also explore the list of [Contributors](https://github.com/WaLSAteam/WaLSAtools/blob/main/CONTRIBUTORS.md) to date.


## **License**

WaLSAtools is licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).
See the LICENSE file for details.
