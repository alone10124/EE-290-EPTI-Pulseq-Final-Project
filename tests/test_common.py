"""Basic building blocks (mrtools.seq.common) and the GE checks on a tiny EPI train (mrtools.seq.checks)."""

import numpy as np
import pypulseq as pp
import pytest

from mrtools.hardware import GE_MR750W
from mrtools.seq.checks import check_sequence
from mrtools.seq.common import interleaved_slice_order, make_epi_readout, make_fatsat, slice_positions, trid


def test_readout_samples_only_flat_top_and_covers_kspace():
    sys = GE_MR750W.opts()
    fov, nx = 0.216, 72
    gx, adc = make_epi_readout(fov, nx, sys, dwell=4e-6)
    assert adc.num_samples == nx
    assert adc.delay >= gx.rise_time - 1e-12                                     # starts on the flat top
    assert adc.delay + adc.num_samples * adc.dwell <= gx.rise_time + gx.flat_time + 1e-12
    k_covered = gx.amplitude * adc.num_samples * adc.dwell                       # 1/m
    assert k_covered == pytest.approx(nx / fov, rel=1e-6)                         # one k-space width per lobe


def test_fatsat_frequency():
    sys = GE_MR750W.opts(slow=True)
    rf, spoiler = make_fatsat(sys)
    assert rf.freq_offset == pytest.approx(-3.4e-6 * 3.0 * 42.576e6, rel=1e-3)   # about -434 Hz at 3 T
    assert spoiler.channel == "z"


def _tiny_epi(hw):
    sys = hw.opts()
    gx, adc = make_epi_readout(0.216, 64, sys, dwell=4e-6)
    gy = pp.make_trapezoid("y", system=sys, area=1 / 0.216)
    seq = pp.Sequence(system=sys)
    rf = pp.make_block_pulse(np.deg2rad(30), duration=0.5e-3, system=sys, delay=sys.rf_dead_time, use="excitation")
    seq.add_block(rf, trid(1))
    for i in range(16):
        seq.add_block(pp.scale_grad(gx, (-1) ** i), adc)
        seq.add_block(gy)
    seq.add_block(pp.make_delay(10e-3))
    return seq, pp.calc_duration(gx) + pp.calc_duration(gy)


def test_ge_checks_pass_on_tiny_epi():
    seq, _ = _tiny_epi(GE_MR750W)
    report = check_sequence(seq, GE_MR750W, "tiny EPI")
    assert report.ok, str(report)
    assert [c.name for c in report.checks][-1] == "TRID labels present"


def test_slice_positions_and_order():
    assert np.allclose(slice_positions(12e-3, 3e-3), [-4.5e-3, -1.5e-3, 1.5e-3, 4.5e-3])
    assert interleaved_slice_order(5).tolist() == [0, 2, 4, 1, 3]
    assert sorted(interleaved_slice_order(42)) == list(range(42))
