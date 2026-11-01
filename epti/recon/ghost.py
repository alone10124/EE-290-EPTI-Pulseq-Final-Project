"""EPI Nyquist-ghost correction from the three navigator lines.

Replaces the encrypted P-code ``mrir_artifact_ghost_compute`` / ``mrir_artifact_ghost_correct`` (MGH mrir toolbox),
as called from ``legacy/libs/libs_EPTI_Recon/library_EPTI_recon/EPTI_SMS_Estimate_Linear_Coeff_WJG.m`` L16 and
``EPTI_SMS_Deghost_WJG.m`` L11-13.

Odd and even EPI lines are read with opposite gradient polarity. Gradient delays and eddy currents shift them
relative to each other, which appears in hybrid (x, ky) space as a phase difference linear in x. The navigators are
three lines at ky = 0 read +, -, +. Comparing the reversed line with the forward ones gives that phase difference, which
is fit as ``phi(x) = phi0 + phi1 * x`` and removed from every reversed line. The fit is the standard magnitude-weighted
estimate (Ahn & Cho, IEEE TMI 1987): the slope comes from adjacent-sample correlations, and the offset from the
de-sloped sum.

Conventions matched to the P-code with tests/ref/fixtures.mat (black-box input/output pairs):

* hybrid transform: ``mrir_iDFT_freqencode(x) = N * fftshift(ifft(ifftshift(x)))`` along the readout, and
  ``mrir_fDFT_freqencode`` is its inverse, so x = 0 is sample N/2;
* coefficients ``(a, b)`` are the *correction*: reversed lines are multiplied by ``exp(i (a + b (n - N/2)))``.
"""

import numpy as np


def to_hybrid(k):
    n = k.shape[-1]
    return n * np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(k, axes=-1), axis=-1), axes=-1)


def to_kspace(h):
    n = h.shape[-1]
    return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(h, axes=-1), axis=-1), axes=-1) / n


def estimate(nav_fwd, nav_rev):
    """Correction coefficients from navigators.

    Parameters
    ----------
    nav_fwd : ``[..., n]`` k-space of the forward navigators (nav1 and nav3 averaged), kx increasing.
    nav_rev : ``[..., n]`` k-space of the reversed navigator (nav2), already in kx-increasing order.

    Returns ``(a, b)`` with shape ``nav_fwd.shape[:-1]``: constant (rad) and slope (rad per readout sample).
    """
    d = to_hybrid(nav_rev) * np.conj(to_hybrid(nav_fwd))          # reversed vs forward phase, weighted by |signal|^2
    n = d.shape[-1]
    slope = np.angle(np.sum(d[..., 1:] * np.conj(d[..., :-1]), axis=-1))
    xs = np.arange(n) - n / 2
    const = np.angle(np.sum(d * np.exp(-1j * slope[..., None] * xs), axis=-1))
    return -const, -slope


def correct(k_rev, a, b):
    """Apply the correction to reversed lines ``[..., n]``; ``a`` and ``b`` broadcast against ``k_rev.shape[:-1]``."""
    n = k_rev.shape[-1]
    xs = np.arange(n) - n / 2
    return to_kspace(to_hybrid(k_rev) * np.exp(1j * (np.asarray(a)[..., None] + np.asarray(b)[..., None] * xs)))


def calib_coefficients(nav, group_of_slice):
    """Per-coil coefficients for each imaging slice group, from the calibration navigators.

    Parameters
    ----------
    nav : ``[pe_line, slice, 3, coil, n]`` regridded calibration navigators (nav1 +, nav2 -, nav3 +).
    group_of_slice : ``[slice]`` index of the imaging slice group each calibration slice belongs to.

    As in ``EPTI_SMS_Estimate_Linear_Coeff_WJG.m``, the navigators of the slices in one SMS group are averaged before
    fitting, and the per-shot coefficients are reduced with a median over PE lines.
    Returns ``(a, b)`` with shape ``[group, coil]``.
    """
    groups = np.unique(group_of_slice)
    a = np.zeros((len(groups), nav.shape[3]))
    b = np.zeros_like(a)
    for gi, g in enumerate(groups):
        nv = nav[:, group_of_slice == g].mean(axis=1)              # [pe_line, 3, coil, n]
        aa, bb = estimate(0.5 * (nv[:, 0] + nv[:, 2]), nv[:, 1])   # [pe_line, coil]
        a[gi], b[gi] = np.median(aa, axis=0), np.median(bb, axis=0)
    return a, b
