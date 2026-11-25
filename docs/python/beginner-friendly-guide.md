---
template: main.html
title: Getting Started
---

# Getting Started with WaLSAtools

## A Beginner-Friendly Guide from Setting Up Python to using WaLSAtools :material-shoe-print:

This guide will walk you through the full setup process for running **WaLSAtools**, even if you are completely new to Python. We recommend using a lightweight and flexible method using [`pyenv`](https://github.com/pyenv/pyenv), which lets you easily install and manage different versions of Python &mdash; without needing admin access.

!!! walsa-info "Note"
    You do **not** need to have Python pre-installed. However, even if Python is already installed on your system, we still recommend the steps below to ensure compatibility and avoid conflicts.

## Step 1: Install Python (via `pyenv`) :material-language-python:

!!! walsa-hint "Operating Systems"

    === "macOS"
        1. Install `pyenv` on macOS using [Homebrew](https://brew.sh):
            ```
            brew update
            brew install pyenv pyenv-virtualenv openssl readline sqlite3 xz zlib ncurses tcl-tk
            ```

        2. Update your shell configuration:   

            === "bash (default)"
                ```bash
                echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
                echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
                echo 'eval "$(pyenv init --path)"' >> ~/.bashrc
                echo 'eval "$(pyenv init -)"' >> ~/.bashrc
                echo 'eval "$(pyenv virtualenv-init -)"' >> ~/.bashrc
                source ~/.bashrc
                ```
            === "zsh (macOS default since Catalina)"
                ```bash
                echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.zshrc
                echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.zshrc
                echo 'eval "$(pyenv init --path)"' >> ~/.zshrc
                echo 'eval "$(pyenv init -)"' >> ~/.zshrc
                echo 'eval "$(pyenv virtualenv-init -)"' >> ~/.zshrc
                source ~/.zshrc
                ```

        3. Install a Python verion:
            ```
            pyenv install 3.12.8
            ```

    === "Linux (Ubuntu/Debian)"
        1. install system dependencies:
            ```bash
            sudo apt update && sudo apt install -y make build-essential libssl-dev zlib1g-dev \
            libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncursesw5-dev \
            xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev
            ```
        2. Install `pyenv` using the recommended installer:
            ```bash
            curl https://pyenv.run | bash
            ```
            Follow the post-installation instructions printed in the terminal to activate `pyenv`.

        3. Update your shell configuration:  

            === "bash"
                ```bash
                echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
                echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
                echo 'eval "$(pyenv init --path)"' >> ~/.bashrc
                echo 'eval "$(pyenv init -)"' >> ~/.bashrc
                echo 'eval "$(pyenv virtualenv-init -)"' >> ~/.bashrc
                source ~/.bashrc
                ```
            === "zsh"
                ```bash
                echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.zshrc
                echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.zshrc
                echo 'eval "$(pyenv init --path)"' >> ~/.zshrc
                echo 'eval "$(pyenv init -)"' >> ~/.zshrc
                echo 'eval "$(pyenv virtualenv-init -)"' >> ~/.zshrc
                source ~/.zshrc
                ```
               
        4. Install a Python verion:
            ```
            pyenv install 3.12.8
            ```

    === "Windows"
        We recommend using the official `pyenv` port for Windows:

        1. Install [pyenv-win](https://github.com/pyenv-win/pyenv-win#installation)

        2. Ensure `pyenv`, `pyenv-win` and Python executables are in your PATH
        
        3. Install a Python verion:
            ```
            pyenv install 3.12.8
            ```

    !!! walsa-info "Why Python 3.12.8?"
        Python 3.12.8 is the version used in this example. WaLSAtools requires Python 3.11 or later, and pip selects dependencies for your Python version. You can install multiple Python versions if you wish and create a separate virtual environment for each.

## Step 2: Create a Virtual Environment :material-cube-outline:

Create a clean environment for your installed Python version (e.g., 3.12.8; hence `walsa_env_3_12_8`) for WaLSAtools:

```bash
pyenv virtualenv 3.12.8 walsa_env_3_12_8
pyenv activate walsa_env_3_12_8
```

!!! walsa-info "If `pyenv` is not initialized properly"
    Double-check your .bashrc (or .bash_profile) / .zshrc entries and reload (`source`) them.

## Step 3: Install WaLSAtools <img src="https://walsa.team/images/WaLSAtools-black.svg" style="height: 1.2em; vertical-align: middle; margin-left: 0.2em; padding-bottom: 0.2em">

First, upgrade `pip`:
```bash
pip install --upgrade pip
```

Then install WaLSAtools from PyPI:
```bash
pip install WaLSAtools
```

!!! walsa-info "What if pip is not installed?"
    If you get an error like `command not found: pip`, you can install it with:
    ```bash
    python -m ensurepip --upgrade
    ```

## Step 4: Verify the Installation :material-code-tags-check:

Start Python inside your **terminal** (within the virtual enviroment):
```
python
```
This will open the Python interactive prompt (you will see something like >>>).   
Then, at the prompt, type:  
```python
from WaLSAtools import WaLSAtools
WaLSAtools
```
If installed correctly, WaLSAtools’ interactive welcome menu will appear 🎉       
Alternatively, you can verify the installation inside a Jupyter notebook (preferred; see below).

!!! walsa-code "Jupyter notebook is highly recommended"
    Although WaLSAtools runs from the terminal, we recommend using **Jupyter notebooks**. All our tutorials and worked examples are written in notebook format for easier use.

## Step 5: Clone the GitHub repository :octicons-mark-github-16:

Although WaLSAtools has already been installed via pip, we strongly recommend also cloning (downloading) the WaLSAtools GitHub repository. By cloning the repository, you will gain access to:

- The complete source codes,
- Worked examples demonstrating practical analyses,
- The full set of documentation files used to build this website.

Having the local repository allows you to explore the materials at your own pace, customize them if needed, and even contribute improvements by submitting a pull request (see the [Contribution](https://walsa.tools/contribution/) page for details).

To clone the repository, simply run:
```
git clone https://github.com/WaLSAteam/WaLSAtools.git
```
This will create a local folder named WaLSAtools containing everything you need.

## Step 6: Install and Use Jupyter Notebooks in VS Code (Highly Recommended) :material-microsoft-visual-studio-code:

We recommend using [Visual Studio Code (VS Code)](https://code.visualstudio.com) to work with **WaLSAtools**. It provides an easy, modern, and fully integrated environment to edit and run Python scripts and Jupyter notebooks — all inside one application!

VS Code supports interactive notebooks (.ipynb), standard Python scripts (.py), and includes many helpful tools (such as syntax highlighting, variable explorers, and inline output).

### Step-by-Step Guide

1. **Download and install VS Code**:  
   🔗 [https://code.visualstudio.com](https://code.visualstudio.com)

2. **Install the Python Extension**:  
    - Open VS Code.
    - Click the **Extensions** icon (left sidebar or press `Ctrl+Shift+X` / `Cmd+Shift+X`).
    - Search for `Python`.
    - Install the official **Python extension by Microsoft**.   
        (It includes full Jupyter notebook support.)

3. **Open your Project Folder**
    - Launch VS Code.
    - Open the folder where your notebooks or scripts are located (or to be located).    
        (Example: If you cloned WaLSAtools GitHub repo, open the WaLSAtools/codes/python/ folder.)

4.	**Create or open a notebook (.ipynb)**:    
	- You can open any example notebook provided in the WaLSAtools GitHub repository.
	- Or create a new one via:   
       File → New File → Jupyter Notebook (name it something like example.ipynb)

5.	**Select your Python environment (kernel)**:
    - When you open a notebook for the first time, VS Code may ask you to select a kernel.      
       💡 If it doesn’t ask, click “Select Kernel” or the Python version shown at the top right of the notebook interface to manually select it. Alternatively, you can press `Ctrl+Shift+P` (`Cmd+Shift+P` on Mac), search for “Python: Select Kernel”, and manually pick it.
    - Choose the Python interpreter corresponding to your WaLSAtools virtual environment (e.g., walsa_env_3_12_8).

6.	**Test WaLSAtools Inside the Notebook**:
    Inside a notebook's code cell, enter the following:
    ``` Python
    from WaLSAtools import WaLSAtools
    WaLSAtools
    ```
    - Click the :material-play-outline: Execute Cell button (on the left side of the code cell), or   
    - Click :material-play-outline: Run All at the top to execute the whole notebook.  

    If installed correctly, you should see the interactive welcome menu for WaLSAtools 🎉

## You are All Set :material-check:

You have now installed Python, set up a clean environment, and verified WaLSAtools.

Explore the [Analysis Tools](https://walsa.tools/python/WaLSAtools/) page and the worked examples.

!!! walsa-issues "Still stuck?"
    Check the [Troubleshooting](https://walsa.tools/python/troubleshooting/) page or post your question in our [GitHub discussion](https://github.com/WaLSAteam/WaLSAtools/discussions).

<br>
