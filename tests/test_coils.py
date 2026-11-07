"""Coil maps (mrtools.recon.coils) on a synthetic multi-coil phantom."""

import numpy as np
import sigpy as sp
import sigpy.mri as mr

from mrtools.recon.coils import coil_combine, espirit_maps


def test_espirit_coil_combination_recovers_object():
    """Fully sampled 8-coil data -> ESPIRiT maps -> SENSE-1 combination reproduces the object's magnitude."""
    n, nc = 72, 8
    obj = sp.shepp_logan((n, n)).astype(np.complex64)
    sens = mr.birdcage_maps((nc, n, n))
    kspace = sp.fft(sens * obj, axes=(-2, -1))
    maps = espirit_maps(kspace)
    combined = coil_combine(sp.ifft(kspace, axes=(-2, -1)), maps)
    mask = np.abs(obj) > 0.1
    err = np.abs(np.abs(combined[mask]) - np.abs(obj[mask])) / np.abs(obj).max()
    assert np.median(err) < 0.02
