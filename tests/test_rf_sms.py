"""SMS excitation pulse (epti.seq.rf_sms)."""

import numpy as np
import pytest
import sigpy.mri.rf as sigpy_rf
from scipy import signal

from epti.hardware import PROFILES
from epti.seq.rf_sms import matlab_resample, matlab_smooth, slice_select_trapezoid, slr_pulse, sms_pulse


def _build(profile):
    hw = PROFILES[profile]
    sys, sys2 = hw.opts(), hw.opts(slow=True)
    rf, gz_wave, _, _ = sms_pulse(60, 3e-3, 4, 4e-3, 2, 0.063, sys2,
                                  0.5 * sys.max_grad / sys.gamma * 100, 0.4 * sys.max_slew / sys.gamma / 10,
                                  matlab_compat=hw.matlab_compat)
    return rf, slice_select_trapezoid(gz_wave, sys), sys2


# ---------------------------------------------------------------- Tier 1: physics, no MATLAB needed
@pytest.mark.parametrize("profile", ["siemens_prisma_ref", "ge_mr750w"])
def test_sms_slice_profile(profile):
    """Bloch-simulate the Pulseq pulse: two 3 mm slices at +-31.5 mm with a 60 deg flip, nothing in between."""
    rf, gz, sys2 = _build(profile)
    dt, n = sys2.rf_raster_time, len(rf.signal)
    z = np.linspace(-0.08, 0.08, 6401)
    a, b = sigpy_rf.sim.abrm(rf.signal * dt * 2 * np.pi, gz.amplitude * z * n * dt, True)
    mxy = np.abs(2 * np.conj(a) * b)
    above = z[mxy > 0.5 * mxy.max()]
    bands = np.split(above, np.nonzero(np.diff(above) > 1e-3)[0] + 1)
    assert len(bands) == 2
    for band, center in zip(bands, (-0.0315, 0.0315)):
        assert abs(band.mean() - center) < 0.1e-3
        assert abs((band[-1] - band[0]) - 3e-3) < 0.15e-3
    assert np.degrees(np.arcsin(min(mxy.max(), 1))) == pytest.approx(60, abs=2)
    outside = np.abs(np.abs(z) - 0.0315) > 6e-3
    assert mxy[outside].max() < 0.01


def test_matlab_resample_preserves_band_limited_signal():
    t = np.arange(200)
    x = np.sin(2 * np.pi * 3 * t / 200) * np.hanning(200)
    y = matlab_resample(x, 1000, 200)
    assert len(y) == 1000
    assert np.allclose(y[::5], x, atol=2e-3)


def test_matlab_smooth_ends_and_middle():
    x = np.array([1.0, 2, 3, 10, 5, 6, 7])
    y = matlab_smooth(x)
    assert y[0] == 1 and y[-1] == 7
    assert y[1] == pytest.approx(2) and y[3] == pytest.approx((2 + 3 + 10 + 5 + 6) / 5)


def test_slr_base_pulse_flip_and_gradient():
    rf, gz, _ = slr_pulse(60, 0.3, 4, 4, 3.2, 8)
    gamma_khz_per_g = 4.2576
    flip = 2 * np.pi * gamma_khz_per_g * 4e-3 * rf.sum()       # rad (rf in Gauss on a 4 us raster)
    assert np.degrees(flip) == pytest.approx(60, rel=1e-6)
    plateau = (4 / 4) / (gamma_khz_per_g * 0.3)                  # TBW/dur / (gamma * thickness), G/cm
    assert gz.max() == pytest.approx(plateau, rel=1e-6)
    assert gz.min() < 0                                          # rephaser lobe present


# ---------------------------------------------------------------- Tier 2: MATLAB fixtures
def test_sigpy_dzrf_matches_pauly_dzrf(fixtures):
    """SigPy designs the 200-tap filter via a shifted 201-tap firls; MATLAB solves the even-length problem directly."""
    ours = np.real(sigpy_rf.dzrf(200, 4, "st", "ls"))
    ref = np.ravel(fixtures["dzrf_200_4_st_ls"])
    assert np.abs(ours - ref).max() / np.abs(ref).max() < 1e-2


def test_base_pulse_matches_toppe_makeslr(fixtures):
    sys_ge = fixtures["sysGE"]
    rf, gz, center = slr_pulse(60, 0.3, 4, 4, float(sys_ge.maxGrad), float(sys_ge.maxSlew))
    rf_ref, gz_ref = np.ravel(fixtures["makeslr_rf_G"]), np.ravel(fixtures["makeslr_gz_Gcm"])
    # MATLAB's exactly symmetric filter has two equal central maxima; their midpoint is the center getsmspulse uses
    peaks = np.nonzero(np.abs(rf_ref) == np.abs(rf_ref).max())[0]
    assert len(peaks) == 2 and peaks.mean() == center
    assert rf.shape == rf_ref.shape and gz.shape == gz_ref.shape
    assert np.abs(rf - rf_ref).max() / np.abs(rf_ref).max() < 1e-2
    assert np.abs(gz - gz_ref).max() / np.abs(gz_ref).max() < 1e-6


def test_sms_pulse_matches_getsmspulse(fixtures):
    rf, _, _ = _build("siemens_prisma_ref")
    ref = np.ravel(fixtures["sms_rf_signal"])
    assert rf.delay == pytest.approx(float(fixtures["sms_rf_delay"]), abs=1e-9)
    assert len(rf.signal) == len(ref)
    assert np.abs(rf.signal - ref).max() / np.abs(ref).max() < 1e-2
    assert rf.center == pytest.approx(float(fixtures["sms_rf_center"]), abs=2e-6)
