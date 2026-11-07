"""SMS building blocks (mrtools.seq.sms)."""

import numpy as np
import pytest
import sigpy.mri.rf as sigpy_rf

from mrtools.hardware import GE_MR750W
from mrtools.seq.sms import caipi_delta_kz, caipi_kz_index, multiband_sinc_pulse


@pytest.mark.parametrize("mb, sep", [(2, 0.063), (3, 0.042)])
def test_multiband_slice_profile(mb, sep):
    """Bloch-simulate the pulse on its plateau: mb slices of the right width and flip angle at the right positions."""
    sys = GE_MR750W.opts(slow=True)
    thickness, flip = 3e-3, 60
    rf, gz, _ = multiband_sinc_pulse(flip, 4e-3, thickness, sep, mb, sys)
    dt, n = sys.rf_raster_time, len(rf.signal)
    z = np.linspace(-0.1, 0.1, 8001)
    a, b = sigpy_rf.sim.abrm(rf.signal * dt * 2 * np.pi, gz.amplitude * z * n * dt, True)
    mxy = np.abs(2 * np.conj(a) * b)
    above = z[mxy > 0.5 * mxy.max()]
    bands = np.split(above, np.nonzero(np.diff(above) > 1e-3)[0] + 1)
    centers = (np.arange(mb) - (mb - 1) / 2) * sep
    assert len(bands) == mb
    for band, c in zip(bands, centers):
        assert abs(band.mean() - c) < 0.1e-3
        assert abs((band[-1] - band[0]) - thickness) < 0.3e-3
        peak = mxy[(z >= band[0]) & (z <= band[-1])].max()
        assert np.degrees(np.arcsin(min(peak, 1))) == pytest.approx(flip, abs=3)
    outside = np.min(np.abs(z[:, None] - centers[None]), axis=1) > 3 * thickness
    assert mxy[outside].max() < 0.02


def test_caipi_pattern_gives_fov_over_mb_shift():
    """Slice s at s*sep from slice 0 gains 2*pi*s/mb of phase per line: a FOV/mb shift in the image."""
    mb, sep = 3, 0.042
    k = caipi_kz_index(mb, 9) * caipi_delta_kz(mb, sep)
    for s in range(mb):
        phase = 2 * np.pi * k * s * sep
        assert np.allclose(np.exp(1j * np.diff(phase[:mb])), np.exp(1j * 2 * np.pi * s / mb))

