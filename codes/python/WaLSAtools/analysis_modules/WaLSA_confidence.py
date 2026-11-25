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

def WaLSA_confidence(ps_perm, siglevel=0.05, nf=None):
    """
    Find the confidence levels (significance levels) for the given power spectrum permutations.

    Parameters:
    -----------
    ps_perm : np.ndarray
        2D array where each row represents a power spectrum, and each column represents a permutation.
    siglevel : float
        Significance level (default 0.05 for 95% confidence).
    nf : int
        Number of frequencies (length of the power spectra). If None, inferred from the shape of ps_perm.

    Returns:
    --------
    signif : np.ndarray
        Significance levels for each frequency.
    """
    if nf is None:
        nf = ps_perm.shape[0]  # Number of frequencies inferred from the shape of ps_perm

    signif = np.zeros((nf))

    # Loop through each frequency
    for iv in range(nf):
        # Extract the permutation values for this frequency
        tmp = np.sort(ps_perm[iv, :])  # Sort the power spectrum permutation values

        # Calculate the number of permutations
        ntmp = len(tmp)

        # Find the significance threshold
        nsig = int(round(siglevel * ntmp))  # Number of permutations to cut off for the significance level

        # Set the confidence level for this frequency
        signif[iv] = tmp[-max(1, nsig)]  # Index zero would select the opposite tail.

    return signif
