"""Ramp-sampling regridding of EPI readouts onto a uniform k-space grid.

Port of ``legacy/libs/libs_EPTI_Recon/library_EPTI_recon/EPI_trapezoid_regridding_WJG.m`` (MGH, Jian Wu; itself
modeled on Gadgetron's EPI reconstruction). The ADC samples part of each trapezoid ramp, so the samples aren't evenly
spaced in k. With the trapezoid's k(t), the uniform k-space values are found by least squares through a sinc
interpolation matrix: ``M = F @ pinv(Q)``, with ``Q[p, q] = sinc(k_p - k_q)``. ``F`` maps the uniform grid to image x,
and an FFT along the readout returns to k-space.

The off-center correction in the original multiplies by 1 (position and read direction are hard-coded to 0), so it is
left out.

Deliberate deviation (bug fix): the original evaluates the image on ``x = linspace(-0.5, (N-1)/(2N), N)``, inherited
from Gadgetron's EPI code. That spacing is (2N-1)/(2N(N-1)), not 1/N, so every regridded readout is stretched by
~0.7% and misinterpolated by up to ~8% of the peak (tests/test_recon.py, band-limited test object). We use the
uniform grid ``x = (n - N/2) / N``, which matches the centered FFT used everywhere else and reproduces the test
object exactly.
"""

import numpy as np


def _trajectory(n_samples, rise, flat, fall, dwell, nx):
    """k index (units of the uniform grid, centered on 0) of each ADC sample, for positive and negative lobes."""
    tot_time = rise + flat + fall
    read_time = dwell * n_samples
    acq_delay = 0.5 * (tot_time - read_time)         # "balanced": ADC centered on the lobe
    tot_area = 0.5 * rise + flat + 0.5 * fall
    read_area = tot_area
    if rise > 0:
        read_area -= 0.5 * acq_delay ** 2 / rise
    if fall > 0:
        read_area -= 0.5 * (tot_time - (acq_delay + read_time)) ** 2 / fall
    scale = nx / read_area
    t = np.arange(1, n_samples + 1) * dwell + acq_delay      # end of each dwell period, as in the original
    k = np.where(t <= rise, 0.5 / rise * t ** 2,
                 np.where(t <= rise + flat, 0.5 * rise + (t - rise),
                          tot_area - 0.5 / fall * (tot_time - t) ** 2))
    pre = 0.5 * tot_area
    return scale * (k - pre), scale * (-k + tot_area - pre)


def regrid_operators(n_samples, rise, flat, fall, dwell, nx):
    """``(M_pos, M_neg)``: ``[nx, n_samples]`` matrices mapping raw samples of a positive/negative lobe (in acquisition
    order) to image-space x on an ``nx``-point grid."""
    traj_pos, traj_neg = _trajectory(n_samples, rise, flat, fall, dwell, nx)
    km = nx // 2
    k_even = np.arange(-km, km + 1)                          # Ne = 2*Km + 1 uniform points
    x = (np.arange(nx) - nx // 2) / nx                      # original: linspace(-0.5, (nx-1)/(2nx), nx), see above
    F = np.exp(1j * 2 * np.pi * np.outer(x, k_even)) / np.sqrt(len(k_even))
    q_pos = np.sinc(traj_pos[:, None] - k_even[None, :])
    q_neg = np.sinc(traj_neg[:, None] - k_even[None, :])
    return F @ np.linalg.pinv(q_pos), F @ np.linalg.pinv(q_neg)


def regrid(lines, polarity, p):
    """Regrid readout lines to uniform k-space.

    Parameters
    ----------
    lines : ``[..., n_lobe, n_coil, n_samples]`` raw data in acquisition order (negative lobes NOT flipped).
    polarity : ``[n_lobe]`` +1/-1 readout sign of each lobe.
    p : ``EptiParams`` with the readout timing.

    Returns ``[..., n_lobe, n_coil, nx]`` k-space with kx increasing along the last axis for every lobe.
    """
    m_pos, m_neg = regrid_operators(lines.shape[-1], p.rise, p.flat, p.fall, p.dwell, p.nx)
    x = np.empty(lines.shape[:-1] + (p.nx,), dtype=np.complex64)
    pos = np.asarray(polarity) > 0
    x[..., pos, :, :] = lines[..., pos, :, :] @ m_pos.T
    x[..., ~pos, :, :] = lines[..., ~pos, :, :] @ m_neg.T
    # MATLAB: fftshift(fft(ifftshift(x, 1), [], 1), 1) along the readout
    return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(x, axes=-1), axis=-1), axes=-1).astype(np.complex64)
