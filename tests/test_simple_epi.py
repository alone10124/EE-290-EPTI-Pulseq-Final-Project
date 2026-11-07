"""The simple EPI sequence script and its recon script, without KomaMRI."""

import numpy as np
import pytest

from conftest import load_script
from mrtools.seq.checks import check_sequence

seq_script = load_script("sequences/scripts/simple_epi.py", "simple_epi_seq")
recon_script = load_script("recon/simple_epi.py", "simple_epi_recon")


def test_full_sequence_passes_checks_and_records_definitions():
    seq = seq_script.build_sequence()
    report = check_sequence(seq, seq_script.HARDWARE, "simple_epi")
    assert report.ok, str(report)
    d = seq.definitions
    positions = np.asarray(d["SlicePositions"])
    assert len(positions) == 40                                     # 120 mm / 3 mm
    assert np.allclose(np.diff(np.sort(positions)), 3e-3)
    assert positions[0] == pytest.approx(-58.5e-3)                  # interleaved: slice 0 first, then 2, ...
    assert positions[1] == pytest.approx(-52.5e-3)


def test_slice_subset_keeps_only_those_slices():
    seq = seq_script.build_sequence(slices=[3, 20])
    assert np.allclose(seq.definitions["SlicePositions"], [1.5e-3, -49.5e-3])     # acquisition order: even first


def test_recon_recovers_object_from_sequence_trajectory():
    """Signal computed directly from the sequence's own k-space trajectory (no simulator) -> recon == object.

    Checks the sample sorting, the reversal of odd lines and the FFT conventions together.
    """
    n = 32
    seq = seq_script.build_sequence(matrix=(n, n), fat_sat=False, slices=[0])
    k_adc = seq.calculate_kspace()[0]                               # [3, n_samples], 1/m
    fov = seq_script.FOV
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
