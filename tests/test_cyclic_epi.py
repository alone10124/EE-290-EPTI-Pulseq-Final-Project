"""The cyclic EPI sequence script, its recon script and the run config, without KomaMRI."""

import numpy as np
import pytest
from pypulseq.supported_labels_rf_use import get_supported_labels

from conftest import REPO, load_script
from mrtools import config
from mrtools.seq.checks import check_sequence

seq_script = load_script("sequences/scripts/cyclic_epi.py", "cyclic_epi_seq")
recon_script = load_script("recon/cyclic_epi.py", "cyclic_epi_recon")
FOV = (0.22, 0.22)


@pytest.fixture(scope="module")
def cyc_cfg():
    return config.load(REPO / "configs" / "cyclic_epi.toml")


def adc_lines(seq, n):
    """k-space and time since excitation of every ADC sample, as [shot, echo, sample] arrays (one slice, no dummy)."""
    k, _, t_exc, _, t_adc = seq.calculate_kspace()
    t_adc, t_exc = np.asarray(t_adc), np.asarray(t_exc)
    t_rel = t_adc - t_exc[np.searchsorted(t_exc, t_adc, side="right") - 1]
    return k[0].reshape(n, n, n), k[1].reshape(n, n, n), t_rel.reshape(n, n, n)


def test_config(cyc_cfg):
    assert cyc_cfg.general["sequence"] == cyc_cfg.general["recon"] == "cyclic_epi"
    assert cyc_cfg.hardware.ncoils == 32 and cyc_cfg.simulation["ncoils"] == 8
    assert cyc_cfg.sim_seq_path.name == "cyclic_epi_sim.seq"


@pytest.mark.parametrize("wrap_timing", ["ignore", "pad"])
def test_all_slices_pass_checks_with_one_trid_per_shot(cyc_cfg, wrap_timing):
    protocol = {k: v for k, v in cyc_cfg.sequence.items() if k != "sim_slices"}
    protocol |= dict(matrix=(16, 16), wrap_timing=wrap_timing)          # full matrix: see the scanner .seq instead
    seq = seq_script.build_sequence(cyc_cfg.hardware, **protocol)
    report = check_sequence(seq, cyc_cfg.hardware, "cyclic_epi")
    assert report.ok, str(report)
    assert len(seq.definitions["SlicePositions"]) == 40
    trid_id = get_supported_labels().index("TRID") + 1
    trids = {v[0] for v in seq.label_set_library.data.values() if v[1] == trid_id}
    assert trids == set(range(1, 18))                                   # shots 0..15 -> 1..16, dummy -> 17


@pytest.mark.parametrize("wrap_timing", ["ignore", "pad"])
def test_rows_cycle_and_wrap_delay(hw, wrap_timing):
    n = 16
    seq = seq_script.build_sequence(hw, fov=FOV, matrix=(n, n), slices=[20], wrap_timing=wrap_timing)
    d = seq.definitions
    kx, ky, t = adc_lines(seq, n)
    rows = (np.arange(n)[:, None] + np.arange(n)[None, :]) % n             # [shot, echo]
    assert np.allclose(ky.mean(-1) * FOV[1], rows - n / 2, atol=1e-6)      # shot s, echo j reads row (s + j) mod n
    assert np.allclose(np.ptp(ky, axis=-1), 0, atol=1e-6)                 # no ky motion during a line
    assert np.all(np.diff(kx[:, 0], axis=-1) > 0) and np.all(np.diff(kx[:, 1], axis=-1) < 0)

    te = t.mean(-1)                                                        # [shot, echo], center of each line
    expected = d["TE"] + np.arange(n) * d["EchoSpacing"] + (np.arange(n)[None, :] > n - 1 - np.arange(n)[:, None]) * d["WrapDelay"]
    assert np.allclose(te, expected, atol=1e-7)
    if wrap_timing == "pad":
        assert d["WrapDelay"] == 0 and np.allclose(np.ptp(te, axis=0), 0, atol=1e-7)   # one TE per echo
    else:
        assert d["WrapDelay"] > 1e-4


def test_slice_subset_keeps_shot_period(hw):
    full = seq_script.build_sequence(hw, matrix=(16, 16), slab=12e-3)
    sub = seq_script.build_sequence(hw, matrix=(16, 16), slab=12e-3, slices=[2])
    assert sub.definitions["ShotPeriod"] == pytest.approx(full.definitions["ShotPeriod"])
    assert sub.duration()[0] == pytest.approx(full.duration()[0])
    assert full.duration()[0] == pytest.approx(17 * full.definitions["ShotPeriod"])  # 16 shots + 1 dummy


@pytest.mark.parametrize("wrap_timing, t2s_tol, b0_tol", [("pad", 1e-3, 0.05), ("ignore", 0.05, 2.0)])
def test_recon_recovers_t2s_and_b0(hw, wrap_timing, t2s_tol, b0_tol):
    """Signal from the sequence's own trajectory and sample times with a known T2* and B0 (KomaMRI's sign:
    +df Hz -> exp(-i 2 pi df t)), 4 coils -> sort, ESPIRiT combination, fit. With "pad" the fit is exact; with
    "ignore" the late rows after the wrap bias it slightly."""
    n = 32
    seq = seq_script.build_sequence(hw, fov=FOV, matrix=(n, n), fat_sat=False, slices=[20], wrap_timing=wrap_timing)
    kx, ky, t = (a.ravel() for a in adc_lines(seq, n))
    grid = (np.arange(n) - n // 2) * FOV[0] / n
    yy, xx = np.meshgrid(grid, grid, indexing="ij")
    disk = xx ** 2 + (yy / 0.8) ** 2 < 0.08 ** 2
    t2s = np.where(xx > 0, 0.03, 0.06)
    b0 = 20 + 300 * yy                                                     # Hz
    angles = np.arange(4) * np.pi / 2
    sens = np.stack([np.exp(-((xx - 0.1 * np.cos(a)) ** 2 + (yy - 0.1 * np.sin(a)) ** 2) / 0.02) * np.exp(1j * a)
                     for a in angles])
    p = np.flatnonzero(disk)
    decay = np.exp(np.outer(t, -1 / t2s.ravel()[p] - 2j * np.pi * b0.ravel()[p]))
    phase = np.exp(-2j * np.pi * (np.outer(kx, xx.ravel()[p]) + np.outer(ky, yy.ravel()[p])))
    signal = (decay * phase) @ sens.reshape(4, -1)[:, p].T                 # [samples, coils]

    d = seq.definitions
    te = d["TE"] + np.arange(n) * d["EchoSpacing"]
    images = recon_script.combine(recon_script.sort_kspace(signal, 1, n, n), "espirit")[0]
    fit = recon_script.fit_maps(images, te)
    m = recon_script.interior(disk) & fit["support"]
    assert m.sum() > 200
    assert np.median(np.abs(fit["t2s"][m] - t2s[m]) / t2s[m]) < t2s_tol
    assert np.median(np.abs(fit["b0_hz"][m] - b0[m])) < b0_tol
