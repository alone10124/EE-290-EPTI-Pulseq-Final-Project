"""Calibration-scan processing: coil sensitivities, B0 map and odd/even background phase.

Follows ``legacy/libs/libs_EPTI_Recon/library_EPTI_recon/EPTI_Preprocessing_B0_sens_BARTv5.m`` (MGH, Zijing Dong,
Fuyixue Wang), minus coil compression (not needed for the 8-16 simulated coils; see the plan), with library pieces:

* GRAPPA: ``pygrappa.mdgrappa`` (Tikhonov ``lamda=0.01``, kernel from echo 1's ACS lines reused for every echo, as the
  original does) in place of the P-code ``mrir_array_GRAPPA_2d_kernel_improved_dev`` / ``mrir_array_GRAPPA_conv_kernel_v2``.
* ESPIRiT: ``sigpy.mri.app.EspiritCalib`` per slice (calibration width 24, eigenvalue crop 0.1, one map) in place of
  ``bart ecalib -r 24 -c 0.1 -m1``. The original runs a 3D ESPIRiT across all slices, which isn't meaningful for a few
  slices 63 mm apart.
* B0: temporal phase unwrapping plus a least-squares line fit of phase against TE (``dB_fitting_JumpCorrect.m``).

The B0 map is in the data's own phase convention: ``phase(TE) = phase0 + 2*pi*b0_hz*TE``. KomaMRI precesses an
off-resonance of +df as ``exp(-i 2 pi df t)``, so for simulated data ``b0_hz = -df`` (tests/test_recon_sim.py checks this).
"""

from dataclasses import dataclass

import numpy as np
import pygrappa
import sigpy as sp
import sigpy.mri as mr

from epti.seq.patterns import cal_pe_index


def ifft2c(k):
    return sp.ifft(k, axes=(-2, -1))       # centered, orthonormal


def fft2c(x):
    return sp.fft(x, axes=(-2, -1))


def zpad_ky(k, n):
    """Zero-pad the ky axis (second to last) to ``n`` lines, centered as MATLAB's zpad.m."""
    m = k.shape[-2]
    start = n // 2 + int(np.ceil(-m / 2))
    out = np.zeros(k.shape[:-2] + (n, k.shape[-1]), dtype=k.dtype)
    out[..., start:start + m, :] = k
    return out


@dataclass
class CalibResult:
    positions: np.ndarray      # [slice] m
    sens: np.ndarray           # [slice, coil, y, x]
    b0_hz: np.ndarray          # [slice, y, x]
    phase0_odd: np.ndarray     # [slice, y, x] rad, background phase of positive-readout echoes
    phase0_even: np.ndarray    # [slice, y, x] rad, background phase of negative-readout echoes
    images: np.ndarray         # [slice, echo, y, x] coil-combined echo images
    mask: np.ndarray           # [slice, y, x]
    tes: np.ndarray            # [echo] s


def echo_times(p, n_echo):
    """TE of each calibration echo: TE + (i - N/2) * ESP for i = 1..N (BARTv5 L8-10)."""
    return p.te + (np.arange(1, n_echo + 1) - n_echo / 2) * p.esp


def grappa_fill(k, acs_lines, lamda=0.01, kernel=(5, 5)):
    """Fill missing ky lines of ``k`` ``[echo, coil, ky, kx]`` with GRAPPA; weights fit on echo 0's ACS region."""
    calib = np.moveaxis(k[0][:, acs_lines, :], 0, -1)                     # [ky_acs, kx, coil]... pygrappa: (ky, kx, coil)
    out = np.empty_like(k)
    weights = None
    for e in range(k.shape[0]):
        ks = np.moveaxis(k[e], 0, -1)                                        # [ky, kx, coil]
        if weights is None:
            filled, weights = pygrappa.mdgrappa(ks, calib, kernel_size=kernel, coil_axis=-1, lamda=lamda,
                                                ret_weights=True)
        else:
            filled = pygrappa.mdgrappa(ks, calib, kernel_size=kernel, coil_axis=-1, lamda=lamda, weights=weights)
        out[e] = np.moveaxis(filled, -1, 0)
    return out


def process(kcal, p, ny, acs_npe=18, espirit_width=24, espirit_crop=0.1):
    """Coil maps, B0 and odd/even phase from the calibration scan.

    Parameters
    ----------
    kcal : ``[pe_line, slice, echo, coil, kx]`` regridded, ghost-corrected k-space of the calibration echoes
        (navigators removed).
    p : calibration ``EptiParams``.
    ny : PE matrix of the imaging scan (the calibration's ``cnseg`` lines are zero-padded to it).
    """
    n_pe, n_slc, n_echo, n_coil, nx = kcal.shape
    pe_idx = cal_pe_index(acs_npe, p.cnseg, p.pat_factor).astype(int) + (p.cnseg - 1) // 2
    acs = np.arange(acs_npe) + (p.cnseg - acs_npe) // 2 + ((p.cnseg - acs_npe) % 2)
    acs = np.intersect1d(acs, pe_idx)
    tes = echo_times(p, n_echo)

    sens, b0, ph_odd, ph_even, images, masks = [], [], [], [], [], []
    for s in range(n_slc):
        k = np.zeros((n_echo, n_coil, p.cnseg, nx), dtype=np.complex64)
        k[:, :, pe_idx, :] = np.moveaxis(kcal[:, s], 0, -2)                # [echo, coil, ky, kx]
        if p.pat_factor > 1:
            k = grappa_fill(k, acs)
        k = zpad_ky(k, ny)
        im = ifft2c(k)                                                       # [echo, coil, y, x]
        smap = mr.app.EspiritCalib(k[0], calib_width=espirit_width, crop=espirit_crop, show_pbar=False).run()
        comb = np.sum(np.conj(smap) * im, axis=1) / (np.sum(np.abs(smap) ** 2, axis=0) + np.finfo(np.float32).eps)
        mag0 = np.abs(comb[0])
        mask = mag0 > 0.2 * mag0.mean()
        b0_s = fit_b0(np.angle(comb), tes, mask)
        demod = comb * np.exp(-1j * 2 * np.pi * b0_s[None] * tes[:, None, None])
        # odd/even = 1-based echo parity = positive/negative readout lobes. Averaging as recon_EPTI_P1.m L204-205.
        ph_odd.append(np.angle(np.mean(np.exp(1j * np.angle(demod[0::2])), axis=0)) * mask)
        ph_even.append(np.angle(np.mean(np.exp(1j * np.angle(demod[1::2])), axis=0)) * mask)
        sens.append(smap), b0.append(b0_s), images.append(comb), masks.append(mask)
    return CalibResult(positions=np.asarray(p.slice_positions), sens=np.stack(sens), b0_hz=np.stack(b0),
                       phase0_odd=np.stack(ph_odd), phase0_even=np.stack(ph_even), images=np.stack(images),
                       mask=np.stack(masks), tes=tes)


def fit_b0(phase, tes, mask):
    """Unwrap ``phase [echo, y, x]`` along echoes and fit ``phase / 2pi = c + b0 * TE`` in ``mask`` (b0 in Hz)."""
    ph = np.unwrap(phase, axis=0) / (2 * np.pi)
    A = np.stack([np.ones_like(tes), tes], axis=1)
    sol, *_ = np.linalg.lstsq(A, ph[:, mask], rcond=None)
    out = np.zeros(mask.shape)
    out[mask] = sol[1]
    return out
