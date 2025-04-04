# -----------------------------------------------------------------------------------------------------
# WaLSAtools: - Wave analysis tools
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

# The following code is a modified version of the PyCWT package.
# The original code and licence:
# PyCWT is released under the open source 3-Clause BSD license:
#
# Copyright (c) 2023 Sebastian Krieger, Nabil Freij, and contributors. All rights
# reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#----------------------------------------------------------------------------------
# Modifications by Shahin Jafarzadeh, 2024
#----------------------------------------------------------------------------------

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np # type: ignore
from tqdm import tqdm # type: ignore
from scipy.stats import chi2 # type: ignore

# Try to import the Python wrapper for FFTW.
try:
    import pyfftw.interfaces.scipy_fftpack as fft # type: ignore
    from multiprocessing import cpu_count

    # Fast planning, use all available threads.
    _FFTW_KWARGS_DEFAULT = {'planner_effort': 'FFTW_ESTIMATE',
                            'threads': cpu_count()}

    def fft_kwargs(signal, **kwargs):
        """Return optimized keyword arguments for FFTW"""
        kwargs.update(_FFTW_KWARGS_DEFAULT)
        kwargs['n'] = len(signal)  # do not pad
        return kwargs

# Otherwise, fall back to 2 ** n padded scipy FFTPACK
except ImportError:
    import scipy.fftpack as fft # type: ignore
    # Can be turned off, e.g. for MKL optimizations
    _FFT_NEXT_POW2 = True

    def fft_kwargs(signal, **kwargs):
        """Return next higher power of 2 for given signal to speed up FFT"""
        if _FFT_NEXT_POW2:
            return {'n': int(2 ** np.ceil(np.log2(len(signal))))}

from scipy.signal import lfilter # type: ignore
from os import makedirs
from os.path import exists, expanduser


def find(condition):
    """Returns the indices where ravel(condition) is true."""
    res, = np.nonzero(np.ravel(condition))
    return res


def ar1(x):
    """
    Allen and Smith autoregressive lag-1 autocorrelation coefficient.
    In an AR(1) model

        x(t) - <x> = gamma(x(t-1) - <x>) + \alpha z(t) ,

    where <x> is the process mean, gamma and \alpha are process
    parameters and z(t) is a Gaussian unit-variance white noise.

    Parameters
    ----------
    x : numpy.ndarray, list
        Univariate time series

    Returns
    -------
    g : float
        Estimate of the lag-one autocorrelation.
    a : float
        Estimate of the noise variance [var(x) ~= a**2/(1-g**2)]
    mu2 : float
        Estimated square on the mean of a finite segment of AR(1)
        noise, mormalized by the process variance.

    References
    ----------
    [1] Allen, M. R. and Smith, L. A. Monte Carlo SSA: detecting
        irregular oscillations in the presence of colored noise.
        *Journal of Climate*, **1996**, 9(12), 3373-3404.
        <http://dx.doi.org/10.1175/1520-0442(1996)009<3373:MCSDIO>2.0.CO;2>
    [2] http://www.madsci.org/posts/archives/may97/864012045.Eg.r.html

    """
    x = np.asarray(x)
    N = x.size
    xm = x.mean()
    x = x - xm

    # Estimates the lag zero and one covariance
    c0 = x.transpose().dot(x) / N
    c1 = x[0:N-1].transpose().dot(x[1:N]) / (N - 1)

    # According to A. Grinsteds' substitutions
    B = -c1 * N - c0 * N**2 - 2 * c0 + 2 * c1 - c1 * N**2 + c0 * N
    A = c0 * N**2
    C = N * (c0 + c1 * N - c1)
    D = B**2 - 4 * A * C

    if D > 0:
        g = (-B - D**0.5) / (2 * A)
    else:
        raise Warning('Cannot place an upperbound on the unbiased AR(1). '
                      'Series is too short or trend is to large.')

    # According to Allen & Smith (1996), footnote 4
    mu2 = -1 / N + (2 / N**2) * ((N - g**N) / (1 - g) -
                                 g * (1 - g**(N - 1)) / (1 - g)**2)
    c0t = c0 / (1 - mu2)
    a = ((1 - g**2) * c0t) ** 0.5

    return g, a, mu2


def ar1_spectrum(freqs, ar1=0.):
    """
    Lag-1 autoregressive theoretical power spectrum.

    Parameters
    ----------
    freqs : numpy.ndarray, list
        Frequencies at which to calculate the theoretical power
        spectrum.
    ar1 : float
        Autoregressive lag-1 correlation coefficient.

    Returns
    -------
    Pk : numpy.ndarray
        Theoretical discrete Fourier power spectrum of noise signal.

    References
    ----------
    [1] http://www.madsci.org/posts/archives/may97/864012045.Eg.r.html

    """
    # According to a post from the MadSci Network available at
    # http://www.madsci.org/posts/archives/may97/864012045.Eg.r.html,
    # the time-series spectrum for an auto-regressive model can be
    # represented as
    #
    # P_k = \frac{E}{\left|1- \sum\limits_{k=1}^{K} a_k \, e^{2 i \pi
    #   \frac{k f}{f_s} } \right|^2}
    #
    # which for an AR1 model reduces to
    #
    freqs = np.asarray(freqs)
    Pk = (1 - ar1 ** 2) / np.abs(1 - ar1 * np.exp(-2 * np.pi * 1j * freqs)) \
        ** 2

    return Pk


def rednoise(N, g, a=1.):
    """
    Red noise generator using filter.

    Parameters
    ----------
    N : int
        Length of the desired time series.
    g : float
        Lag-1 autocorrelation coefficient.
    a : float, optional
        Noise innovation variance parameter.

    Returns
    -------
    y : numpy.ndarray
        Red noise time series.

    """
    if g == 0:
        yr = np.randn(N, 1) * a
    else:
        # Twice the decorrelation time.
        tau = int(np.ceil(-2 / np.log(np.abs(g))))
        yr = lfilter([1, 0], [1, -g], np.random.randn(N + tau, 1) * a)
        yr = yr[tau:]

    return yr.flatten()


def rect(x, normalize=False):
    """TODO: describe what I do."""
    if type(x) in [int, float]:
        shape = [x, ]
    elif type(x) in [list, dict]:
        shape = x
    elif type(x) in [np.ndarray, np.ma.core.MaskedArray]:
        shape = x.shape
    X = np.zeros(shape)
    X[0] = X[-1] = 0.5
    X[1:-1] = 1

    if normalize:
        X /= X.sum()

    return X


def boxpdf(x):
    """
    Forces the probability density function of the input data to have
    a boxed distribution.

    Parameters
    ----------
    x (array like) :
        Input data

    Returns
    -------
    X (array like) :
        Boxed data varying between zero and one.
    Bx, By (array like) :
        Data lookup table.

    """
    import numpy as np
    x = np.asarray(x)
    n = x.size

    # Kind of 'unique'
    i = np.argsort(x)
    d = (np.diff(x[i]) != 0)
    j = find(np.concatenate([d, [True]]))
    X = x[i][j]

    j = np.concatenate([[0], j + 1])
    Y = 0.5 * (j[0:-1] + j[1:]) / n
    bX = np.interp(x, X, Y)

    return bX, X, Y


def get_cache_dir():
    """Returns the location of the cache directory."""
    # Sets cache directory according to user home path.
    cache_dir = '{}/.cache/pycwt/'.format(expanduser('~'))
    # Creates cache directory if not existant.
    if not exists(cache_dir):
        makedirs(cache_dir)
    # Returns cache directory.
    return cache_dir


import numpy as np
from scipy.special import gamma
from scipy.signal import convolve2d
from scipy.special.orthogonal import hermitenorm


class Morlet(object):
    """Implements the Morlet wavelet class.

    Note that the input parameters f and f0 are angular frequencies.
    f0 should be more than 0.8 for this function to be correct, its
    default value is f0 = 6.

    """

    def __init__(self, f0=6):
        self._set_f0(f0)
        self.name = 'Morlet'

    def psi_ft(self, f):
        """Fourier transform of the approximate Morlet wavelet."""
        return (np.pi ** -0.25) * np.exp(-0.5 * (f - self.f0) ** 2)

    def psi(self, t):
        """Morlet wavelet as described in Torrence and Compo (1998)."""
        return (np.pi ** -0.25) * np.exp(1j * self.f0 * t - t ** 2 / 2)

    def flambda(self):
        """Fourier wavelength as of Torrence and Compo (1998)."""
        return (4 * np.pi) / (self.f0 + np.sqrt(2 + self.f0 ** 2))

    def coi(self):
        """e-Folding Time as of Torrence and Compo (1998)."""
        return 1. / np.sqrt(2)

    def sup(self):
        """Wavelet support defined by the e-Folding time."""
        return 1. / self.coi

    def _set_f0(self, f0):
        # Sets the Morlet wave number, the degrees of freedom and the
        # empirically derived factors for the wavelet bases C_{\delta},
        # gamma, \delta j_0 (Torrence and Compo, 1998, Table 2)
        self.f0 = f0             # Wave number
        self.dofmin = 2          # Minimum degrees of freedom
        if self.f0 == 6:
            self.cdelta = 0.776  # Reconstruction factor
            self.gamma = 2.32    # Decorrelation factor for time averaging
            self.deltaj0 = 0.60  # Factor for scale averaging
        else:
            self.cdelta = -1
            self.gamma = -1
            self.deltaj0 = -1

    def smooth(self, W, dt, dj, scales):
        """Smoothing function used in coherence analysis.

        Parameters
        ----------
        W :
        dt :
        dj :
        scales :

        Returns
        -------
        T :

        """
        # The smoothing is performed by using a filter given by the absolute
        # value of the wavelet function at each scale, normalized to have a
        # total weight of unity, according to suggestions by Torrence &
        # Webster (1999) and by Grinsted et al. (2004).
        m, n = W.shape

        # Filter in time.
        k = 2 * np.pi * fft.fftfreq(fft_kwargs(W[0, :])['n'])
        k2 = k ** 2
        snorm = scales / dt
        # Smoothing by Gaussian window (absolute value of wavelet function)
        # using the convolution theorem: multiplication by Gaussian curve in
        # Fourier domain for each scale, outer product of scale and frequency
        F = np.exp(-0.5 * (snorm[:, np.newaxis] ** 2) * k2)  # Outer product
        smooth = fft.ifft(F * fft.fft(W, axis=1, **fft_kwargs(W[0, :])),
                          axis=1,  # Along Fourier frequencies
                          **fft_kwargs(W[0, :], overwrite_x=True))
        T = smooth[:, :n]  # Remove possibly padded region due to FFT

        if np.isreal(W).all():
            T = T.real

        # Filter in scale. For the Morlet wavelet it's simply a boxcar with
        # 0.6 width.
        wsize = self.deltaj0 / dj * 2
        win = rect(int(np.round(wsize)), normalize=True)
        T = convolve2d(T, win[:, np.newaxis], 'same')  # Scales are "vertical"

        return T


class Paul(object):
    """Implements the Paul wavelet class.

    Note that the input parameter f is the angular frequency and that
    the default order for this wavelet is m=4.

    """
    def __init__(self, m=4):
        self._set_m(m)
        self.name = 'Paul'

    # def psi_ft(self, f):
    #     """Fourier transform of the Paul wavelet."""
    #     return (2 ** self.m /
    #             np.sqrt(self.m * np.prod(range(2, 2 * self.m))) *
    #             f ** self.m * np.exp(-f) * (f > 0))
    
    def psi_ft(self, f): # modified by SJ
        """Fourier transform of the Paul wavelet with limits to prevent underflow."""
        expnt = -f
        expnt[expnt < -100] = -100  # Apply the threshold to avoid extreme values
        return (2 ** self.m /
                np.sqrt(self.m * np.prod(range(2, 2 * self.m))) *
                f ** self.m * np.exp(expnt) * (f > 0))

    def psi(self, t):
        """Paul wavelet as described in Torrence and Compo (1998)."""
        return (2 ** self.m * 1j ** self.m * np.prod(range(2, self.m - 1)) /
                np.sqrt(np.pi * np.prod(range(2, 2 * self.m + 1))) *
                (1 - 1j * t) ** (-(self.m + 1)))

    def flambda(self):
        """Fourier wavelength as of Torrence and Compo (1998)."""
        return 4 * np.pi / (2 * self.m + 1)

    def coi(self):
        """e-Folding Time as of Torrence and Compo (1998)."""
        return np.sqrt(2)

    def sup(self):
        """Wavelet support defined by the e-Folding time."""
        return 1 / self.coi

    def _set_m(self, m):
        # Sets the m derivative of a Gaussian, the degrees of freedom and the
        # empirically derived factors for the wavelet bases C_{\delta},
        # gamma, \delta j_0 (Torrence and Compo, 1998, Table 2)
        self.m = m               # Wavelet order
        self.dofmin = 2          # Minimum degrees of freedom
        if self.m == 4:
            self.cdelta = 1.132  # Reconstruction factor
            self.gamma = 1.17    # Decorrelation factor for time averaging
            self.deltaj0 = 1.50  # Factor for scale averaging
        else:
            self.cdelta = -1
            self.gamma = -1
            self.deltaj0 = -1


class DOG(object):
    """Implements the derivative of a Guassian wavelet class.

    Note that the input parameter f is the angular frequency and that
    for m=2 the DOG becomes the Mexican hat wavelet, which is then
    default.

    """
    def __init__(self, m=2):
        self._set_m(m)
        self.name = 'DOG'

    def psi_ft(self, f):
        """Fourier transform of the DOG wavelet."""
        return (- 1j ** self.m / np.sqrt(gamma(self.m + 0.5)) * f ** self.m *
                np.exp(- 0.5 * f ** 2))

    def psi(self, t):
        """DOG wavelet as described in Torrence and Compo (1998).

        The derivative of a Gaussian of order `n` can be determined using
        the probabilistic Hermite polynomials. They are explicitly
        written as:
            Hn(x) = 2 ** (-n / s) * n! * sum ((-1) ** m) *
                    (2 ** 0.5 * x) ** (n - 2 * m) / (m! * (n - 2*m)!)
        or in the recursive form:
            Hn(x) = x * Hn(x) - nHn-1(x)

        Source: http://www.ask.com/wiki/Hermite_polynomials

        """
        p = hermitenorm(self.m)
        return ((-1) ** (self.m + 1) * np.polyval(p, t) *
                np.exp(-t ** 2 / 2) / np.sqrt(gamma(self.m + 0.5)))

    def flambda(self):
        """Fourier wavelength as of Torrence and Compo (1998)."""
        return (2 * np.pi / np.sqrt(self.m + 0.5))

    def coi(self):
        """e-Folding Time as of Torrence and Compo (1998)."""
        return 1 / np.sqrt(2)

    def sup(self):
        """Wavelet support defined by the e-Folding time."""
        return 1 / self.coi

    def _set_m(self, m):
        # Sets the m derivative of a Gaussian, the degrees of freedom and the
        # empirically derived factors for the wavelet bases C_{\delta},
        # gamma, \delta j_0 (Torrence and Compo, 1998, Table 2).
        self.m = m               # m-derivative
        self.dofmin = 1          # Minimum degrees of freedom
        if self.m == 2:
            self.cdelta = 3.541  # Reconstruction factor
            self.gamma = 1.43    # Decorrelation factor for time averaging
            self.deltaj0 = 1.40  # Factor for scale averaging
        elif self.m == 6:
            self.cdelta = 1.966
            self.gamma = 1.37
            self.deltaj0 = 0.97
        else:
            self.cdelta = -1
            self.gamma = -1
            self.deltaj0 = -1


class MexicanHat(DOG):
    """Implements the Mexican hat wavelet class.

    This class inherits the DOG class using m=2.

    """
    def __init__(self):
        self.name = 'Mexican Hat'
        self._set_m(2)




def cwt(signal, dt, dj=1/12, s0=-1, J=-1, wavelet='morlet', freqs=None, pad=True):
    """Continuous wavelet transform of the signal at specified scales.

    Parameters
    ----------
    signal : numpy.ndarray, list
        Input signal array.
    dt : float
        Sampling interval.
    dj : float, optional
        Spacing between discrete scales. Default value is 1/12.
        Smaller values will result in better scale resolution, but
        slower calculation and plot.
    s0 : float, optional
        Smallest scale of the wavelet. Default value is 2*dt.
    J : float, optional
        Number of scales less one. Scales range from s0 up to
        s0 * 2**(J * dj), which gives a total of (J + 1) scales.
        Default is J = (log2(N * dt / so)) / dj.
    wavelet : instance of Wavelet class, or string
        Mother wavelet class. Default is Morlet wavelet.
    freqs : numpy.ndarray, optional
        Custom frequencies to use instead of the ones corresponding
        to the scales described above. Corresponding scales are
        calculated using the wavelet Fourier wavelength.
    pad : optional. Default is True. if set, then pad the time series 
        with enough zeroes to get N up to the next higher power of 2. 
        This prevents wraparound from the end of the time series to 
        the beginning, and also speeds up the FFT's used to do the 
        wavelet transform. This will not eliminate all edge effects.
        (added by SJ)

    Returns
    -------
    W : numpy.ndarray
        Wavelet transform according to the selected mother wavelet.
        Has (J+1) x N dimensions.
    sj : numpy.ndarray
        Vector of scale indices given by sj = s0 * 2**(j * dj),
        j={0, 1, ..., J}.
    freqs : array like
        Vector of Fourier frequencies (in 1 / time units) that
        corresponds to the wavelet scales.
    coi : numpy.ndarray
        Returns the cone of influence, which is a vector of N
        points containing the maximum Fourier period of useful
        information at that particular time. Periods greater than
        those are subject to edge effects.
    fft : numpy.ndarray
        Normalized fast Fourier transform of the input signal.
    fftfreqs : numpy.ndarray
        Fourier frequencies (in 1/time units) for the calculated
        FFT spectrum.

    Example
    -------
    >> mother = wavelet.Morlet(6.)
    >> wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(signal,
           0.25, 0.25, 0.5, 28, mother)

    """
    wavelet = _check_parameter_wavelet(wavelet)

    # Original signal length
    n0 = len(signal)
    # If no custom frequencies are set, then set default frequencies
    # according to input parameters `dj`, `s0` and `J`. Otherwise, set wavelet
    # scales according to Fourier equivalent frequencies.
    
    # Zero-padding if enabled (added + the entire code further modified by SJ)
    if pad:
        next_power_of_two = int(2 ** np.ceil(np.log2(n0)))
        signal_padded = np.zeros(next_power_of_two)
        signal_padded[:n0] = signal - np.mean(signal)  # Remove the mean and pad with zeros
    else:
        signal_padded = signal - np.mean(signal)  # Remove the mean without padding
    
    N = len(signal_padded)  # Length of the padded signal

    # Calculate scales and frequencies if not provided
    if freqs is None:
        # Smallest resolvable scale
        # if s0 == -1:
        #     s0 = 2 * dt / wavelet.flambda()
        # # Number of scales
        # if J == -1:
        #     J = int(np.round(np.log2(N * dt / s0) / dj))
        if s0 == -1:
            s0 = 2 * dt
        if J == -1:
            J = int((np.log(float(N) * dt / s0) / np.log(2)) / dj)
        # The scales as of Mallat 1999
        sj = s0 * 2 ** (np.arange(0, J + 1) * dj)
        # Fourier equivalent frequencies
        freqs = 1 / (wavelet.flambda() * sj)
    else:
        # The wavelet scales using custom frequencies.
        sj = 1 / (wavelet.flambda() * freqs)

    # Fourier transform of the (padded) signal
    signal_ft = fft.fft(signal_padded, **fft_kwargs(signal_padded))
    # Fourier angular frequencies
    ftfreqs = 2 * np.pi * fft.fftfreq(N, dt)

    # Creates wavelet transform matrix as outer product of scaled transformed
    # wavelets and transformed signal according to the convolution theorem.
    # (i)   Transform scales to column vector for outer product;
    # (ii)  Calculate 2D matrix [s, f] for each scale s and Fourier angular
    #       frequency f;
    # (iii) Calculate wavelet transform;
    sj_col = sj[:, np.newaxis]
    psi_ft_bar = ((sj_col * ftfreqs[1] * N) ** .5 *
                  np.conjugate(wavelet.psi_ft(sj_col * ftfreqs)))
    W = fft.ifft(signal_ft * psi_ft_bar, axis=1,
                 **fft_kwargs(signal_ft, overwrite_x=True))

    # Trim the wavelet transform to original signal length if padded
    if pad:
        W = W[:, :n0]  # Trim to the original signal length
    
    # Checks for NaN in transform results and removes them from the scales if
    # needed, frequencies and wavelet transform. Trims wavelet transform at
    # length `n0`.
    sel = np.invert(np.isnan(W).all(axis=1))
    if np.any(sel):
        sj = sj[sel]
        freqs = freqs[sel]
        W = W[sel, :]

    # Determines the cone-of-influence. Note that it is returned as a function
    # of time in Fourier periods. Uses triangualr Bartlett window with
    # non-zero end-points.
    coi = (n0 / 2 - np.abs(np.arange(0, n0) - (n0 - 1) / 2))
    coi = wavelet.flambda() * wavelet.coi() * dt * coi

    return (W[:, :n0], sj, freqs, coi, signal_ft[1:N//2] / N ** 0.5,
            ftfreqs[1:N//2] / (2 * np.pi))


def icwt(W, sj, dt, dj=1/12, wavelet='morlet'):
    """Inverse continuous wavelet transform.

    Parameters
    ----------
    W : numpy.ndarray
        Wavelet transform, the result of the `cwt` function.
    sj : numpy.ndarray
        Vector of scale indices as returned by the `cwt` function.
    dt : float
        Sample spacing.
    dj : float, optional
        Spacing between discrete scales as used in the `cwt`
        function. Default value is 0.25.
    wavelet : instance of Wavelet class, or string
        Mother wavelet class. Default is Morlet

    Returns
    -------
    iW : numpy.ndarray
        Inverse wavelet transform.

    Example
    -------
    >> mother = wavelet.Morlet()
    >> wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(var,
           0.25, 0.25, 0.5, 28, mother)
    >> iwave = wavelet.icwt(wave, scales, 0.25, 0.25, mother)

    """
    wavelet = _check_parameter_wavelet(wavelet)

    a, b = W.shape
    c = sj.size
    if a == c:
        sj = (np.ones([b, 1]) * sj).transpose()
    elif b == c:
        sj = np.ones([a, 1]) * sj
    else:
        raise Warning('Input array dimensions do not match.')

    # As of Torrence and Compo (1998), eq. (11)
    iW = (dj * np.sqrt(dt) / (wavelet.cdelta * wavelet.psi(0)) *
          (np.real(W) / np.sqrt(sj)).sum(axis=0))
    return iW


def significance(signal, dt, scales, sigma_test=0, alpha=None,
                 significance_level=0.95, dof=-1, wavelet='morlet'):
    """Significance test for the one dimensional wavelet transform.

    Parameters
    ----------
    signal : array like, float
        Input signal array. If a float number is given, then the
        variance is assumed to have this value. If an array is
        given, then its variance is automatically computed.
    dt : float
        Sample spacing.
    scales : array like
        Vector of scale indices given returned by `cwt` function.
    sigma_test : int, optional
        Sets the type of significance test to be performed.
        Accepted values are 0 (default), 1 or 2. See notes below for
        further details.
    alpha : float, optional
        Lag-1 autocorrelation, used for the significance levels.
        Default is 0.0.
    significance_level : float, optional
        Significance level to use. Default is 0.95.
    dof : variant, optional
        Degrees of freedom for significance test to be set
        according to the type set in sigma_test.
    wavelet : instance of Wavelet class, or string
        Mother wavelet class. Default is Morlet

    Returns
    -------
    signif : array like
        Significance levels as a function of scale.
    fft_theor (array like):
        Theoretical red-noise spectrum as a function of period.

    Notes
    -----
    If sigma_test is set to 0, performs a regular chi-square test,
    according to Torrence and Compo (1998) equation 18.

    If set to 1, performs a time-average test (equation 23). In
    this case, dof should be set to the number of local wavelet
    spectra that where averaged together. For the global
    wavelet spectra it would be dof=N, the number of points in
    the time-series.

    If set to 2, performs a scale-average test (equations 25 to
    28). In this case dof should be set to a two element vector
    [s1, s2], which gives the scale range that were averaged
    together. If, for example, the average between scales 2 and
    8 was taken, then dof=[2, 8].

    """
    wavelet = _check_parameter_wavelet(wavelet)

    try:
        n0 = len(signal)
    except TypeError:
        n0 = 1
    J = len(scales) - 1
    dj = np.log2(scales[1] / scales[0])

    if n0 == 1:
        variance = signal
    else:
        variance = signal.std() ** 2

    if alpha is None:
        alpha, _, _ = ar1(signal)

    period = scales * wavelet.flambda()  # Fourier equivalent periods
    freq = dt / period                   # Normalized frequency
    dofmin = wavelet.dofmin              # Degrees of freedom with no smoothing
    Cdelta = wavelet.cdelta              # Reconstruction factor
    gamma_fac = wavelet.gamma            # Time-decorrelation factor
    dj0 = wavelet.deltaj0                # Scale-decorrelation factor

    # Theoretical discrete Fourier power spectrum of the noise signal
    # following Gilman et al. (1963) and Torrence and Compo (1998),
    # equation 16.
    def pk(k, a, N):
        return (1 - a ** 2) / (1 + a ** 2 - 2 * a * np.cos(2 * np.pi * k / N))
    fft_theor = pk(freq, alpha, n0)
    fft_theor = variance * fft_theor     # Including time-series variance
    signif = fft_theor

    try:
        if dof == -1:
            dof = dofmin
    except ValueError:
        pass

    if sigma_test == 0:  # No smoothing, dof=dofmin, TC98 sec. 4
        dof = dofmin
        # As in Torrence and Compo (1998), equation 18.
        chisquare = chi2.ppf(significance_level, dof) / dof
        signif = fft_theor * chisquare
    elif sigma_test == 1:  # Time-averaged significance
        if len(dof) == 1:
            dof = np.zeros(1, J+1) + dof
        sel = find(dof < 1)
        dof[sel] = 1
        # As in Torrence and Compo (1998), equation 23:
        dof = dofmin * (1 + (dof * dt / gamma_fac / scales) ** 2) ** 0.5
        sel = find(dof < dofmin)
        dof[sel] = dofmin  # Minimum dof is dofmin
        for n, d in enumerate(dof):
            chisquare = chi2.ppf(significance_level, d) / d
            signif[n] = fft_theor[n] * chisquare
    elif sigma_test == 2:  # Time-averaged significance
        if len(dof) != 2:
            raise Exception('DOF must be set to [s1, s2], '
                            'the range of scale-averages')
        if Cdelta == -1:
            raise ValueError('Cdelta and dj0 not defined '
                             'for {} with f0={}'.format(wavelet.name,
                                                        wavelet.f0))
        s1, s2 = dof
        sel = find((scales >= s1) & (scales <= s2))
        navg = sel.size
        if navg == 0:
            raise ValueError('No valid scales between {} and {}.'.format(s1,
                                                                         s2))
        # As in Torrence and Compo (1998), equation 25.
        Savg = 1 / sum(1. / scales[sel])
        # Power-of-two mid point:
        Smid = np.exp((np.log(s1) + np.log(s2)) / 2.)
        # As in Torrence and Compo (1998), equation 28.
        dof = (dofmin * navg * Savg / Smid) * \
              ((1 + (navg * dj / dj0) ** 2) ** 0.5)
        # As in Torrence and Compo (1998), equation 27.
        fft_theor = Savg * sum(fft_theor[sel] / scales[sel])
        chisquare = chi2.ppf(significance_level, dof) / dof
        # As in Torrence and Compo (1998), equation 26.
        signif = (dj * dt / Cdelta / Savg) * fft_theor * chisquare
    else:
        raise ValueError('sigma_test must be either 0, 1, or 2.')

    return signif, fft_theor


def xwt(y1, y2, dt, dj=1/12, s0=-1, J=-1, significance_level=0.95,
        wavelet='morlet', normalize=False, no_default_signif=False):
    """Cross wavelet transform (XWT) of two signals.

    The XWT finds regions in time frequency space where the time series
    show high common power.

    Parameters
    ----------
    y1, y2 : numpy.ndarray, list
        Input signal array to calculate cross wavelet transform.
    dt : float
        Sample spacing.
    dj : float, optional
        Spacing between discrete scales. Default value is 1/12.
        Smaller values will result in better scale resolution, but
        slower calculation and plot.
    s0 : float, optional
        Smallest scale of the wavelet. Default value is 2*dt.
    J : float, optional
        Number of scales less one. Scales range from s0 up to
        s0 * 2**(J * dj), which gives a total of (J + 1) scales.
        Default is J = (log2(N*dt/so))/dj.
    wavelet : instance of a wavelet class, optional
        Mother wavelet class. Default is Morlet wavelet.
    significance_level : float, optional
        Significance level to use. Default is 0.95.
    normalize : bool, optional
        If set to true, normalizes CWT by the standard deviation of
        the signals.

    Returns
    -------
    xwt (array like):
        Cross wavelet transform according to the selected mother
        wavelet.
    x (array like):
        Intersected independent variable.
    coi (array like):
        Cone of influence, which is a vector of N points containing
        the maximum Fourier period of useful information at that
        particular time. Periods greater than those are subject to
        edge effects.
    freqs (array like):
        Vector of Fourier equivalent frequencies (in 1 / time units)
        that correspond to the wavelet scales.
    signif (array like):
        Significance levels as a function of scale.

    Notes
    -----
    Torrence and Compo (1998) state that the percent point function
    (PPF) -- inverse of the cumulative distribution function -- of a
    chi-square distribution at 95% confidence and two degrees of
    freedom is Z2(95%)=3.999. However, calculating the PPF using
    chi2.ppf gives Z2(95%)=5.991. To ensure similar significance
    intervals as in Grinsted et al. (2004), one has to use confidence
    of 86.46%.

    """
    wavelet = _check_parameter_wavelet(wavelet)

    # Makes sure input signal are numpy arrays.
    y1 = np.asarray(y1)
    y2 = np.asarray(y2)
    # Calculates the standard deviation of both input signals.
    std1 = y1.std()
    std2 = y2.std()
    # Normalizes both signals, if appropriate.
    if normalize:
        y1_normal = (y1 - y1.mean()) / std1
        y2_normal = (y2 - y2.mean()) / std2
    else:
        y1_normal = y1
        y2_normal = y2

    # Calculates the CWT of the time-series making sure the same parameters
    # are used in both calculations.
    _kwargs = dict(dj=dj, s0=s0, J=J, wavelet=wavelet)
    W1, sj, freq, coi, _, _ = cwt(y1_normal, dt, **_kwargs)
    W2, sj, freq, coi, _, _ = cwt(y2_normal, dt, **_kwargs)

    # Now the wavelet transform coherence
    # W12ini = W1 * W2.conj()
    # scales = np.ones([1, y1.size]) * sj[:, None]
    # # -- Normalization by Scale and Smoothing
    # W12 = wavelet.smooth(W12ini / scales, dt, dj, sj)

    # Calculates the cross CWT of y1 and y2.
    W12 = W1 * W2.conj()

    # And the significance tests. Note that the confidence level is calculated
    # using the percent point function (PPF) of the chi-squared cumulative
    # distribution function (CDF) instead of using Z1(95%) = 2.182 and
    # Z2(95%)=3.999 as suggested by Torrence & Compo (1998) and Grinsted et
    # al. (2004). If the CWT has been normalized, then std1 and std2 should
    # be reset to unity, otherwise the standard deviation of both series have
    # to be calculated.
    if normalize:
        std1 = std2 = 1.
    a1, _, _ = ar1(y1)
    a2, _, _ = ar1(y2)
    Pk1 = ar1_spectrum(freq * dt, a1)
    Pk2 = ar1_spectrum(freq * dt, a2)
    dof = wavelet.dofmin
    if not no_default_signif:
        PPF = chi2.ppf(significance_level, dof)
        signif = (std1 * std2 * (Pk1 * Pk2) ** 0.5 * PPF / dof)
    else:
        signif = np.asarray([0])

    # The resuts:
    return W12, coi, freq, signif


def wct(y1, y2, dt, dj=1/12, s0=-1, J=-1, sig=False,
        significance_level=0.95, wavelet='morlet', normalize=False, **kwargs):
    """Wavelet coherence transform (WCT).

    The WCT finds regions in time frequency space where the two time
    series co-vary, but do not necessarily have high power.

    Parameters
    ----------
    y1, y2 : numpy.ndarray, list
        Input signals.
    dt : float
        Sample spacing.
    dj : float, optional
        Spacing between discrete scales. Default value is 1/12.
        Smaller values will result in better scale resolution, but
        slower calculation and plot.
    s0 : float, optional
        Smallest scale of the wavelet. Default value is 2*dt.
    J : float, optional
        Number of scales less one. Scales range from s0 up to
        s0 * 2**(J * dj), which gives a total of (J + 1) scales.
        Default is J = (log2(N*dt/so))/dj.
    sig : bool 
        set to compute signficance, default is True
    significance_level (float, optional) :
        Significance level to use. Default is 0.95.
    normalize (boolean, optional) :
        If set to true, normalizes CWT by the standard deviation of
        the signals.

    Returns
    -------
    WCT : magnitude of coherence
    aWCT : phase angle of coherence
    coi (array like):
        Cone of influence, which is a vector of N points containing
        the maximum Fourier period of useful information at that
        particular time. Periods greater than those are subject to
        edge effects.
    freq (array like):
        Vector of Fourier equivalent frequencies (in 1 / time units)    coi :  
    sig :  Significance levels as a function of scale 
       if sig=True when called, otherwise zero.

    See also
    --------
    cwt, xwt

    """
    wavelet = _check_parameter_wavelet(wavelet)

    nt = len(y1)

    # Checking some input parameters
    if s0 == -1:
        # s0 = 2 * dt / wavelet.flambda()
        s0 = 2 * dt
    if J == -1:
        # Number of scales
        # J = int(np.round(np.log2(y1.size * dt / s0) / dj))
        J = int((np.log(float(nt) * dt / s0) / np.log(2)) / dj)

    # Makes sure input signals are numpy arrays.
    y1 = np.asarray(y1)
    y2 = np.asarray(y2)
    # Calculates the standard deviation of both input signals.
    std1 = y1.std()
    std2 = y2.std()
    # Normalizes both signals, if appropriate.
    if normalize:
        y1_normal = (y1 - y1.mean()) / std1
        y2_normal = (y2 - y2.mean()) / std2
    else:
        y1_normal = y1
        y2_normal = y2

    # Calculates the CWT of the time-series making sure the same parameters
    # are used in both calculations.
    _kwargs = dict(dj=dj, s0=s0, J=J, wavelet=wavelet)
    W1, sj, freq, coi, _, _ = cwt(y1_normal, dt, **_kwargs)
    W2, sj, freq, coi, _, _ = cwt(y2_normal, dt, **_kwargs)

    scales1 = np.ones([1, y1.size]) * sj[:, None]
    scales2 = np.ones([1, y2.size]) * sj[:, None]

    # Smooth the wavelet spectra before truncating -- Time Smoothing
    S1 = wavelet.smooth(np.abs(W1) ** 2 / scales1, dt, dj, sj)
    S2 = wavelet.smooth(np.abs(W2) ** 2 / scales2, dt, dj, sj)

    # Now the wavelet transform coherence
    W12 = W1 * W2.conj()
    scales = np.ones([1, y1.size]) * sj[:, None]
    # -- Normalization by Scale and Scale Smoothing
    S12 = wavelet.smooth(W12 / scales, dt, dj, sj)
    WCT = np.abs(S12) ** 2 / (S1 * S2)
    aWCT = np.angle(W12)

    # Calculates the significance using Monte Carlo simulations with 95%
    # confidence as a function of scale.
    if sig:
        a1, b1, c1 = ar1(y1)
        a2, b2, c2 = ar1(y2)

        sig = wct_significance(a1, a2, dt=dt, dj=dj, s0=s0, J=J,
                               significance_level=significance_level,
                               wavelet=wavelet, **kwargs)
    else:
        sig = np.asarray([0])

    return WCT, aWCT, coi, freq, sig


def wct_significance(al1, al2, dt, dj, s0, J, significance_level=0.95,
                     wavelet='morlet', mc_count=50, progress=True,
                     cache=True):
    """Wavelet coherence transform significance.

    Calculates WCT significance using Monte Carlo simulations with
    95% confidence.

    Parameters
    ----------
    al1, al2: float
        Lag-1 autoregressive coeficients of both time series.
    dt : float
        Sample spacing.
    dj : float, optional
        Spacing between discrete scales. Default value is 1/12.
        Smaller values will result in better scale resolution, but
        slower calculation and plot.
    s0 : float, optional
        Smallest scale of the wavelet. Default value is 2*dt.
    J : float, optional
        Number of scales less one. Scales range from s0 up to
        s0 * 2**(J * dj), which gives a total of (J + 1) scales.
        Default is J = (log2(N*dt/so))/dj.
    significance_level : float, optional
        Significance level to use. Default is 0.95.
    wavelet : instance of a wavelet class, optional
        Mother wavelet class. Default is Morlet wavelet.
    mc_count : integer, optional
        Number of Monte Carlo simulations. Default is 300.
    progress : bool, optional
        If `True` (default), shows progress bar on screen.
    cache : bool, optional
        If `True` (default) saves cache to file.

    Returns
    -------
    TODO

    """

    if cache:
        # Load cache if previously calculated. It is assumed that wavelet
        # analysis is performed using the wavelet's default parameters.
        aa = np.round(np.arctanh(np.array([al1, al2]) * 4))
        aa = np.abs(aa) + 0.5 * (aa < 0)
        cache_file = 'wct_sig_{:0.5f}_{:0.5f}_{:0.5f}_{:0.5f}_{:d}_{}'\
            .format(aa[0], aa[1], dj, s0 / dt, J, wavelet.name)
        cache_dir = get_cache_dir()
        try:
            dat = np.loadtxt('{}/{}.gz'.format(cache_dir, cache_file),
                             unpack=True)
            print('NOTE: WCT significance loaded from cache.\n')
            return dat
        except IOError:
            pass

    # Some output to the screen
    print('Calculating wavelet coherence significance')

    # Choose N so that largest scale has at least some part outside the COI
    ms = s0 * (2 ** (J * dj)) / dt
    N = int(np.ceil(ms * 6))
    noise1 = rednoise(N, al1, 1)
    nW1, sj, freq, coi, _, _ = cwt(noise1, dt=dt, dj=dj, s0=s0, J=J,
                                   wavelet=wavelet)

    period = np.ones([1, N]) / freq[:, None]
    coi = np.ones([J + 1, 1]) * coi[None, :]
    outsidecoi = (period <= coi)
    scales = np.ones([1, N]) * sj[:, None]
    sig95 = np.zeros(J + 1)
    maxscale = find(outsidecoi.any(axis=1))[-1]
    sig95[outsidecoi.any(axis=1)] = np.nan

    nbins = 1000
    wlc = np.ma.zeros([J + 1, nbins])
    # Displays progress bar with tqdm
    for _ in tqdm(range(mc_count), disable=not progress):
        # Generates two red-noise signals with lag-1 autoregressive
        # coefficients given by al1 and al2
        noise1 = rednoise(N, al1, 1)
        noise2 = rednoise(N, al2, 1)
        # Calculate the cross wavelet transform of both red-noise signals
        kwargs = dict(dt=dt, dj=dj, s0=s0, J=J, wavelet=wavelet)
        nW1, sj, freq, coi, _, _ = cwt(noise1, **kwargs)
        nW2, sj, freq, coi, _, _ = cwt(noise2, **kwargs)
        nW12 = nW1 * nW2.conj()
        # Smooth wavelet wavelet transforms and calculate wavelet coherence
        # between both signals.
        S1 = wavelet.smooth(np.abs(nW1) ** 2 / scales, dt, dj, sj)
        S2 = wavelet.smooth(np.abs(nW2) ** 2 / scales, dt, dj, sj)
        S12 = wavelet.smooth(nW12 / scales, dt, dj, sj)
        R2 = np.ma.array(np.abs(S12) ** 2 / (S1 * S2), mask=~outsidecoi)
        # Walks through each scale outside the cone of influence and builds a
        # coherence coefficient counter.
        for s in range(maxscale):
            cd = np.floor(R2[s, :] * nbins)
            for j, t in enumerate(cd[~cd.mask]):
                wlc[s, int(t)] += 1

    # After many, many, many Monte Carlo simulations, determine the
    # significance using the coherence coefficient counter percentile.
    wlc.mask = (wlc.data == 0.)
    R2y = (np.arange(nbins) + 0.5) / nbins
    for s in range(maxscale):
        sel = ~wlc[s, :].mask
        P = wlc[s, sel].data.cumsum()
        P = (P - 0.5) / P[-1]
        sig95[s] = np.interp(significance_level, P, R2y[sel])

    if cache:
        # Save the results on cache to avoid to many computations in the future
        np.savetxt('{}/{}.gz'.format(cache_dir, cache_file), sig95)

    # And returns the results
    return sig95


def _check_parameter_wavelet(wavelet):
    mothers = {'morlet': Morlet, 'paul': Paul, 'dog': DOG,
               'mexicanhat': MexicanHat}
    # Checks if input parameter is a string. For backwards
    # compatibility with Python 2 we check if instance is a
    # `str`.
    try:
        if isinstance(wavelet, str):
            return mothers[wavelet]()
    except NameError:
        if isinstance(wavelet, str):
            return mothers[wavelet]()
    # Otherwise, return itself.
    return wavelet
