"""GE peripheral nerve stimulation (PNS) model.

Port of ``legacy/libs/toppe-main/+toppe/pns.m`` (Rolf Schulte, GE, 2014; modified by J.-F. Nielsen for TOPPE), model 1
("convolution in time domain"). Reference: T. Linz, "Introduction to dB/dt" (2006); IEC/DIN EN 60601-2-33.

The original evaluates the convolution sum point by point (O(n^2)). Here it is the same causal convolution done with
an FFT, so whole sequences can be checked. GE's scanner shows the result as a percentage of the first-level
threshold: below 80 % is "normal mode", and 80-100 % needs first controlled mode.
"""

import numpy as np
from scipy.signal import fftconvolve

# Coil parameters from the pns.m header (values from the scanner's Scandbdt.cfg / GRSubsystemHWO.xml):
#   coil: (chronaxie [s], rheobase [T/s], alpha = effective coil length [m])
COILS = {
    "xrmw": (360e-6, 20.0, 0.324),      # MR750w
    "xrm": (334e-6, 23.4, 0.333),       # MR750
    "whole": (370e-6, 23.7, 0.344),     # HDx TRM whole
    "zoom": (354e-6, 29.1, 0.309),      # HDx TRM zoom
    "hrmb": (359e-6, 26.5, 0.370),      # UHP
    "hrmw": (642.4e-6, 17.9, 0.310),    # Premier
    "magnus": (611e-6, 55.2, 0.324),    # MAGNUS
}


def pns(grad, coil, gdt=4e-6):
    """PNS threshold (%) over time for a gradient waveform.

    Parameters
    ----------
    grad : (n_axes, n_points) array in T/m, sampled every ``gdt`` seconds.
    coil : one of ``COILS``.

    Returns
    -------
    p_thresh : (n_points,) combined PNS threshold in percent, sqrt of the sum of squares over axes (pns.m Eq. 31).
    p_axes : (n_axes, n_points) per-axis threshold in percent.
    """
    chronaxie, rheobase, alpha = COILS[coil.lower()]
    grad = np.atleast_2d(np.asarray(grad, dtype=float))
    sr_min = rheobase / alpha                                    # min slew rate causing stimulation [T/m/s]
    slew = np.concatenate([np.zeros((grad.shape[0], 1)), np.diff(grad, axis=1)], axis=1) / gdt
    n = grad.shape[1]
    kernel = 1.0 / (chronaxie + (np.arange(n) + 0.5) * gdt) ** 2   # pns.m: (chronaxie + ((l2:-1:1) - 0.5)*gdt).^2
    p_axes = 100 * gdt * chronaxie / sr_min * np.stack([fftconvolve(s, kernel)[:n] for s in slew])
    return np.sqrt((p_axes ** 2).sum(axis=0)), p_axes
