"""GE PNS model (epti.seq.ge_pns)."""

import numpy as np
import pytest

from epti.seq.ge_pns import COILS, pns


def _pns_direct(grad, coil, gdt=4e-6):
    """Literal transcription of the time-domain loop in legacy/libs/toppe-main/+toppe/pns.m (model 1)."""
    chronaxie, rheobase, alpha = COILS[coil]
    sr_min = rheobase / alpha
    slew = np.concatenate([np.zeros((grad.shape[0], 1)), np.diff(grad, axis=1)], axis=1) / gdt
    n = grad.shape[1]
    pt = np.zeros_like(grad)
    for l2 in range(1, n + 1):
        denom = (chronaxie + (np.arange(l2, 0, -1) - 0.5) * gdt) ** 2
        pt[:, l2 - 1] = 100 * gdt * chronaxie / sr_min * (slew[:, :l2] / denom).sum(axis=1)
    return np.sqrt((pt ** 2).sum(axis=0)), pt


def test_fft_convolution_equals_direct_loop():
    rng = np.random.default_rng(0)
    grad = np.cumsum(rng.standard_normal((3, 400)), axis=1) * 1e-4      # T/m, smooth-ish random walk
    ours, ours_axes = pns(grad, "xrmw")
    ref, ref_axes = _pns_direct(grad, "xrmw")
    assert np.allclose(ours_axes, ref_axes, rtol=1e-9, atol=1e-9)
    assert np.allclose(ours, ref, rtol=1e-9, atol=1e-9)


def test_long_ramp_approaches_slew_over_rheobase():
    """A long constant-slew ramp stimulates at 100 * (slew / SRmin) * tau / (chronaxie + tau)."""
    chronaxie, rheobase, alpha = COILS["xrmw"]
    gdt, slew, n = 4e-6, 50.0, 5000
    grad = np.zeros((1, n))
    grad[0, 1:] = slew * gdt * np.arange(1, n)
    p, _ = pns(grad, "xrmw", gdt)
    tau = (n - 1) * gdt
    expected = 100 * slew / (rheobase / alpha) * tau / (chronaxie + tau)
    assert p[-1] == pytest.approx(expected, rel=2e-3)
