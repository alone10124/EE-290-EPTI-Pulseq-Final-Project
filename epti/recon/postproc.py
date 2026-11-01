"""Post-processing of reconstructed EPTI echo images: apodization, T2* fit and optimal echo combination.

* Apodization replaces the P-code ``mrir_filter_raw_apodize_1d(.., 0.1)`` (SubspaceRecon_Process_v5.m L98). The
  MATLAB fixture shows it is exactly ``scipy.signal.windows.tukey(N, 0.1)`` applied to centered k-space.
* ``t2star_fit`` / ``combine_echoes`` port ``legacy/libs/libs_EPTI_Recon/library_EPTI_preprocessing/
  optimalCombineEchoes.m`` and ``T2s_GE_EPTI_adaptiveMask.m`` (Fuyixue Wang, MGH): a log-linear mono-exponential fit
  that skips the first 3 echoes and, per voxel, any echo below 10% of the first fitted echo; then echoes are combined
  with weights ``TE * exp(-TE / T2*)`` (the optimal BOLD weighting of Posse et al., MRM 1999).
"""

import numpy as np
import sigpy as sp
from scipy.ndimage import binary_fill_holes
from scipy.signal import windows


def apodize(img, factor=0.1):
    """Tukey-window the k-space of ``img [..., y, x]`` along both in-plane axes."""
    ny, nx = img.shape[-2:]
    w = windows.tukey(ny, factor)[:, None] * windows.tukey(nx, factor)[None, :]
    return sp.ifft(sp.fft(img, axes=(-2, -1)) * w, axes=(-2, -1))


def t2star_fit(mag, tes, mask, first=3, threshold=0.1):
    """Log-linear T2* fit. ``mag [echo, y, x]``; returns ``(S0, T2* in s, n_echoes_used)``."""
    m = mag[first:]
    te = np.asarray(tes)[first:]
    use = (m > threshold * m[:1]) & mask[None] & (m > 0)
    logm = np.log(np.where(use, m, 1.0))
    w = use.astype(float)
    # weighted normal equations for log S = log S0 - TE / T2*, per voxel
    s0, s1, s2 = w.sum(0), (w * te[:, None, None]).sum(0), (w * te[:, None, None] ** 2).sum(0)
    b0, b1 = (w * logm).sum(0), (w * logm * te[:, None, None]).sum(0)
    det = s0 * s2 - s1 ** 2
    ok = mask & (det > 0) & (s0 >= 2)
    slope = np.where(ok, (s0 * b1 - s1 * b0) / np.where(ok, det, 1), 0)
    inter = np.where(ok, (s2 * b0 - s1 * b1) / np.where(ok, det, 1), 0)
    with np.errstate(divide="ignore"):
        t2s = np.where(ok & (slope < 0), -1 / np.where(slope < 0, slope, -1), 0.0)
    t2s = np.clip(np.where(ok & (t2s <= 0), 4 * te[-1], t2s), 0, 4 * te[-1]) * ok    # as optimalCombineEchoes.m
    return np.exp(inter) * ok, t2s, s0.astype(int)


def fit_mask(mag_mean_over_echo):
    """Fitting mask: > 20% of the mean, holes filled (optimalCombineEchoes.m)."""
    m = mag_mean_over_echo > 0.2 * mag_mean_over_echo.mean()
    return binary_fill_holes(m)


def combine_echoes(mag, tes, t2s):
    """Optimally combined image from ``mag [echo, y, x]`` with weights TE*exp(-TE/T2*)."""
    te = np.asarray(tes)[:, None, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        w = te * np.exp(-te / np.where(t2s > 0, t2s, np.inf)[None])
        w = np.nan_to_num(w / w.sum(0, keepdims=True))
    return (mag * w).sum(0)
