"""The simple EPI sequence script, its recon script and the run config, without KomaMRI."""

import numpy as np
import pytest

from conftest import load_script
from mrtools import config
from mrtools.hardware import HardwareProfile
from mrtools.seq.checks import check_sequence

seq_script = load_script("sequences/scripts/simple_epi.py", "simple_epi_seq")
recon_script = load_script("recon/simple_epi.py", "simple_epi_recon")


def test_config_paths_and_hardware(cfg):
    assert cfg.run == "simple_epi" and cfg.general["recon"] == "simple_epi"
    assert isinstance(cfg.hardware, HardwareProfile) and cfg.hardware.ncoils == 32
    assert cfg.seq_path.name == "simple_epi.seq" and cfg.sim_seq_path.name == "simple_epi_sim.seq"
    assert cfg.sim_input == cfg.sim_seq_path == cfg.recon_seq                   # recon reads what was simulated
    assert cfg.recon_raw == config.REPO / "out" / "simple_epi_sim" / "raw.npz"


def test_unknown_hardware_key_is_an_error(cfg):
    d = {f: getattr(cfg.hardware, f) for f in cfg.hardware.__dataclass_fields__}
    with pytest.raises(ValueError, match="max_slw"):
        HardwareProfile.from_dict(d | {"max_slw": 1.0})


def test_full_sequence_from_config_passes_checks(cfg):
    protocol = {k: v for k, v in cfg.sequence.items() if k != "sim_slices"}
    seq = seq_script.build_sequence(cfg.hardware, **protocol)
    report = check_sequence(seq, cfg.hardware, "simple_epi")
    assert report.ok, str(report)
    d = seq.definitions
    positions = np.asarray(d["SlicePositions"])
    assert len(positions) == 40                                     # 120 mm / 3 mm
    assert np.allclose(np.diff(np.sort(positions)), 3e-3)
    assert positions[0] == pytest.approx(-58.5e-3)                  # interleaved: slice 0 first, then 2, ...
    assert positions[1] == pytest.approx(-52.5e-3)


def test_slice_subset_keeps_only_those_slices(hw):
    d = seq_script.build_sequence(hw, slices=[3, 20]).definitions
    assert np.allclose(d["SlicePositions"], [1.5e-3, -49.5e-3])     # acquisition order: even first
    assert np.atleast_1d(d["SliceIndices"]).astype(int).tolist() == [20, 3]


def test_recon_recovers_object_from_sequence_trajectory(hw):
    """Signal computed directly from the sequence's own k-space trajectory (no simulator) -> recon == object.

    Checks the sample sorting, the reversal of odd lines and the FFT conventions together.
    """
    n, fov = 32, (0.22, 0.22)
    seq = seq_script.build_sequence(hw, fov=fov, matrix=(n, n), fat_sat=False, slices=[0])
    k_adc = seq.calculate_kspace()[0]                               # [3, n_samples], 1/m
    grid = (np.arange(n) - n // 2) * fov[0] / n                     # pixel positions of sigpy's centered FFT
    yy, xx = np.meshgrid(grid, grid, indexing="ij")
    obj = np.zeros((n, n))
    obj[8:20, 10:26] = 1.0
    obj[22:26, 5:9] = 0.5
    phase = np.exp(-2j * np.pi * (np.outer(k_adc[0], xx.ravel()) + np.outer(k_adc[1], yy.ravel())))
    signal = (phase @ obj.ravel())[:, None]                         # one coil
    k = recon_script.sort_kspace(signal, 1, n, n)
    img = np.abs(recon_script.reconstruct(k, "rss")[0])
    img /= img.max()
    assert np.abs(img - obj).max() < 1e-3
