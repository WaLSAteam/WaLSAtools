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

import numpy as np # type: ignore
from .analysis_modules.WaLSA_speclizer import WaLSA_speclizer
from .analysis_modules.WaLSA_k_omega import WaLSA_k_omega
from .analysis_modules.WaLSA_pod import WaLSA_pod
from .analysis_modules.WaLSA_cross_spectra import WaLSA_cross_spectra
from .analysis_modules.WaLSA_interactive import interactive

class WaLSAtools:
    """
    Main class to perform spectral analysis using different methods.
    
    It runs interactively if called without (or empty
    and acts as a function if called with parentheses and parameters.
    """
    
    def __call__(self, signal=None, time=None, method=None, **kwargs):
        """
        Main function to perform spectral analysis using different methods.
        
        Parameters:
            signal (array): The input signal (1D, or 3D).
            time (array): The time array of the signal.
            method (str): The method to use for spectral analysis ('fft', 'lombscargle', 'k-omega', etc.)
            **kwargs: Additional parameters for data preparation and analysis methods.
        
        Returns:
            Results of the selected analysis method.
        """
        if signal is None and time is None and method is None:
            interactive()  # Run interactive (help) mode if no inputs are provided
            return None
        
        method_map = {
            'fft': WaLSA_speclizer,
            'lombscargle': WaLSA_speclizer,
            'wavelet': WaLSA_speclizer,
            'welch': WaLSA_speclizer,
            'emd': WaLSA_speclizer,
            'k-omega': WaLSA_k_omega,
            'pod': WaLSA_pod
        }
        
        if method not in method_map:
            raise ValueError(f"Unknown method '{method}'. Please choose from {list(method_map.keys())}.")
        
        func = method_map[method]
        
        if method in ['fft', 'lombscargle', 'wavelet', 'welch', 'emd'] and 'data1' in kwargs and 'data2' in kwargs:
            return WaLSA_cross_spectra(signal=np.ones(10), time=time, method=method, **kwargs)
        else:
            return func(signal=signal, time=time, method=method, **kwargs)
    
    def __repr__(self):
        """
        Custom representation to allow interactive mode when 'WaLSAtools' is typed without parentheses.
        """
        interactive()  # Call the interactive (help) function
        return ''  # Return an empty string so that no additional output is shown

# Instantiate the object so that it works like a function and runs interactively
WaLSAtools = WaLSAtools()

def main(argv=None):
    """Open the installed command-line guide, or print command usage."""
    import argparse

    parser = argparse.ArgumentParser(
        prog='WaLSAtools', description='Open the WaLSAtools interactive guide.'
    )
    parser.add_argument('--interactive', action='store_true',
                        help='open the guide (also the default)')
    parser.parse_args(argv)
    try:
        interactive()
    except EOFError:
        return 0
    except KeyboardInterrupt:
        return 130
    return 0
