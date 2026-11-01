"""Global B0 drift between dynamics, from the three navigator lines of each imaging TR.

Port of ``legacy/libs/libs_EPTI_Recon/library_EPTI_recon/DriftB0_esti_nav.m`` (the ``B0_update = 0`` path of
recon_EPTI_P1.m L183-185). The phase accrued between navigators 1 and 3 (2 echo spacings apart) is compared between the
current dynamic and the first one; the difference, averaged over x and coils, is the field change in Hz. It has the
same sign convention as the calibration B0 map (phase = 2 pi B0 t), so it is simply added to it.
"""

import numpy as np


def estimate(nav_first, nav_now, esp):
    """``nav_* [3, coil, kx]`` ghost-corrected navigator k-space; returns the drift in Hz."""
    def hyb(k):
        return np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(k, axes=-1), axis=-1), axes=-1)
    n1, n2 = hyb(nav_first), hyb(nav_now)
    t = n2 * np.conj(n1) / (np.abs(n2) * np.abs(n1) + np.finfo(float).eps)
    t = t[2] * np.conj(t[0]) * np.abs(n1[0])
    return float(np.angle(t.mean()) / (2 * esp) / (2 * np.pi))
