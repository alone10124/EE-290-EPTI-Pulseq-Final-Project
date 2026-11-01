"""EPTI forward model for one SMS slice group, as a SigPy linear operator.

Replaces the custom BART ``picsS`` + ``linopScript`` forward model of
``legacy/libs/libs_EPTI_Recon/library_EPTI_recon/SubspaceRecon_Process_v5.m`` (Zijing Dong, Fuyixue Wang, MGH):

    y[t, c](ky_t, kx) = F_xy{ sum_s S[s, c] * sum_k Phi'[t, k, s] * coeff[k, s] }  sampled on ky line ky_t,
    Phi'[t, k, s] = Phi[t, k] * exp(i 2 pi B0_s TE_t) * Phase0_s(t) * E_s(t)

``coeff [K, slice, y, x]`` are the subspace coefficient maps. ``Phi`` is the temporal basis, B0 and Phase0 come from
the calibration scan (Phase0 = odd/even background phase, chosen by each lobe's readout polarity), ``S`` holds the coil
sensitivities, and ``E_s(t)`` is the blipped-CAIPI phase of slice ``s`` at lobe ``t``.

``E_s(t)`` is derived from the sequence physics rather than copied from ``Generate_SMS_encode_Matrix2.m``. The kz blips
put lobe ``t`` at ``kz(t) = deltakz * cumsum(kz_step)`` (0 for the navigators), giving phase ``exp(-i 2 pi kz z)``
(KomaMRI: ``Mxy ~ exp(-i gamma int Bz dt)``). The sequence sets an ADC phase ``phi = 2 pi kz z_ref`` with
``z_ref = group_center + slice_sep (MB-1)/2`` (seq_EPTI_imaging_P1.m L208-210), and Pulseq/KomaMRI demodulate with
``exp(-i phi)`` (KomaMRIBase ADC.jl ``get_adc_phase_compensation``). So the residual is
``E_s(t) = exp(-i 2 pi kz(t) (z_s + sign * z_ref))`` with ``sign = +1``. For MB = 2 at the center group, the upper
slice flips sign on alternate kz states and the lower one doesn't. The opposite convention (``sign = -1``) is
available for scanner data; confirm it on the first scanner data set. Getting it wrong swaps the SMS slices, which
the simulation test caught.
"""

import numpy as np
import sigpy as sp

from epti.seq.patterns import calc_pattern_sms_blip


def ky_indices(p, ny):
    """0-based ky grid line of each lobe. ``pattern_GE`` is 1-based with the k-space center at (Rseg+PFshift)/2;
    the centered FFT puts ky = 0 at ``ny // 2``."""
    return (np.round(p.pattern_ge).astype(int) - (p.rseg + p.pf_shift) // 2 + ny // 2)


def caipi_phase(p, slice_z, group_center, n_nav=3, adc_phase_sign=1):
    """``E [t, s]``: blipped-CAIPI phase of each slice at each lobe (see module docstring)."""
    n_echo = p.n_lobes - n_nav
    steps, dkz = calc_pattern_sms_blip(p.mb, n_echo, p.slice_sep)
    kz = dkz * np.concatenate([np.zeros(n_nav), np.cumsum(steps)])
    z_ref = group_center + p.slice_sep * (p.mb - 1) / 2
    return np.exp(-1j * 2 * np.pi * kz[:, None] * (np.asarray(slice_z)[None, :] + adc_phase_sign * z_ref))


def time_phase(b0_hz, phase0_odd, phase0_even, tes, polarity):
    """``Phase_T [t, s, y, x] = exp(i 2 pi B0 TE_t) * exp(i Phase0_{odd/even})`` (SubspaceRecon_Process_v5.m L35-38)."""
    p0 = np.where(np.asarray(polarity)[:, None, None, None] > 0, np.exp(1j * phase0_odd)[None],
                  np.exp(1j * phase0_even)[None])
    return np.exp(1j * 2 * np.pi * b0_hz[None] * np.asarray(tes)[:, None, None, None]) * p0


def build(phi, phase_t, caipi, sens, ky_idx):
    """Forward operator ``A: coeff [K, s, y, x] -> k-space [t, coil, ky, kx]`` (zero except on each lobe's ky line).

    phi [t, K]; phase_t [t, s, y, x]; caipi [t, s]; sens [s, coil, y, x]; ky_idx [t].
    Returns ``(A, mask)`` with ``mask [t, 1, ky, 1]``.
    """
    nt, k = phi.shape
    ns, nc, ny, nx = sens.shape
    phi_prime = (phi[:, :, None, None, None] * (phase_t * caipi[:, :, None, None])[:, None]).astype(np.complex64)
    mask = np.zeros((nt, 1, ny, 1), dtype=np.complex64)
    mask[np.arange(nt), 0, ky_idx, 0] = 1
    ops = [
        sp.linop.Reshape([1, k, ns, ny, nx], [k, ns, ny, nx]),
        sp.linop.Multiply([1, k, ns, ny, nx], phi_prime),                  # -> [t, K, s, y, x]
        sp.linop.Sum([nt, k, ns, ny, nx], axes=(1,)),                       # -> [t, s, y, x]
        sp.linop.Reshape([nt, ns, 1, ny, nx], [nt, ns, ny, nx]),
        sp.linop.Multiply([nt, ns, 1, ny, nx], sens[None].astype(np.complex64)),   # -> [t, s, coil, y, x]
        sp.linop.Sum([nt, ns, nc, ny, nx], axes=(1,)),                      # SMS: slices add up -> [t, coil, y, x]
        sp.linop.FFT([nt, nc, ny, nx], axes=(-2, -1)),                      # centered, orthonormal
        sp.linop.Multiply([nt, nc, ny, nx], mask),
    ]
    return sp.linop.Compose(ops[::-1]), mask


def images(coeff, phi, phase_t):
    """Echo images ``[t, s, y, x] = sum_k Phi[t,k] Phase_T[t,s] coeff[k,s]`` (SubspaceRecon_Process_v5.m L79)."""
    return np.einsum("tk,ksyx->tsyx", phi, coeff) * phase_t
