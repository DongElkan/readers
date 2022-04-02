"""
Analysis of chromatographic peaks
"""
import numpy as np
from scipy import optimize, special


_ = np.seterr(all='warn', over='raise', invalid='ignore')


def _emg(x, h, sigma, mu, tau):
    """
    Fit exponentially modified gaussian function to simulate the
    shape of chromatographic peaks.

    Parameters
    ----------
    x: np.ndarray
        The peak for fitting
    h: float
        The amplitude of Gaussian
    sigma: float
        Gaussian variance
    mu: float
        Gaussian mean
    tau: float
        Exponent relaxation time

    Returns
    -------
    y: np.ndarray
        Fitted peak

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
    [2] Kalambet Y, et al. J Chemometrics. 2011; 25, 352–356.
    [3] Jeansonne MS, et al. J Chromatogr Sci. 1991, 29, 258-266.
    """
    # set up parameters
    m = (x - mu) / sigma
    c = sigma / tau
    z = (c - m) / np.sqrt(2)
    y = np.empty(z.shape)
    ix = np.zeros(z.shape, dtype=bool)
    # using log transform to calculate the values
    if (z < 0).any():
        ix[z < 0] = True
        # gaussian
        g = c * c / 2 - (m[ix] * c)
        # exponential
        e = np.log(h * c * np.sqrt(np.pi / 2)) + np.log(special.erfc(z[ix]))
        y[ix] = g + e

    if (z > 6.71e7).any():
        ix[z > 6.71e7] = True
        m2 = m[z > 6.71e7]
        y[z > 6.71e7] = np.log(h / (1 - m2 / c)) - m2 * m2 / 2

    ix = np.logical_not(ix)
    if ix.any() > 0:
        m2 = m[ix]
        g = -m2 * m2 / 2
        e = np.log(h * c * np.sqrt(np.pi / 2)) + np.log(special.erfcx(z[ix]))
        y[ix] = g + e

    return np.exp(np.clip(y, -200, None))


def _gaussian(x, h, mu, sigma):
    """ Fit gaussian curve. """
    y = ((x - mu) / sigma) ** 2 / 2
    # avoid underflow error
    return h * np.exp(-np.clip(y, None, 200))


def fit_curve(rt, intensity, shape='emg'):
    """
    Fit the curve using scipy's curve_fit optimization.

    Parameters
    ----------
    rt: np.ndarray
        Retention time
    intensity: np.ndarray
        Peak intensities
    shape: str
        Function for fitting the peak. Accepted functions are:
        "emg": exponentially modified gaussian function
        "gau": gaussian function
        Default is "emg".

    Returns
    -------
    param: np.ndarray
        Parameters for fitting the peak

    """
    # initialization
    j = np.argmax(intensity)
    h0, m0 = intensity[j], rt[j]
    s0 = np.std(rt)
    t0 = 1
    # fit curve using optimization
    if shape == 'emg':
        param, _ = optimize.curve_fit(_emg, rt, intensity, p0=(h0, s0, m0, t0))
    elif shape == 'gaussian':
        param, _ = optimize.curve_fit(_gaussian, rt, intensity, p0=(h0, m0, s0))
    else:
        param = None
    return param
