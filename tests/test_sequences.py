"""Full EPTI calibration + imaging sequences (epti.seq.build)."""

import numpy as np
import pypulseq as pp
import pytest

from epti.config import CalibProtocol, ImagingProtocol, load_config
from epti.hardware import PROFILES
from epti.seq.build import build_all
from epti.seq.calib import build_calib
from epti.seq.compare import compare_sequences
from epti.seq.imaging import build_imaging, excited_slice_positions

from conftest import REPO


@pytest.fixture(scope="module")
def ge_full():
    return build_all(load_config(REPO / "config/epti_3mm_ge.yaml"))


@pytest.fixture(scope="module")
def sim_small():
    return build_all(load_config(REPO / "config/sim_small.yaml"))


# ---------------------------------------------------------------- Tier 1
def test_ge_protocol_passes_all_hardware_checks(ge_full):
    for rep in ge_full.reports:
        assert rep.ok, str(rep)


def test_ge_protocol_timing(ge_full):
    # One TR per volume. The per-group delay is rounded *up* to the raster (as in seq_EPTI_imaging_P1.m L167), so the
    # volume can run up to one raster per slice group long, plus the 2-raster preamble delay.
    slack = (21 + 2) * 4e-6
    assert 2.0 <= ge_full.imaging.duration()[0] <= 2.0 + slack
    assert 29 * 1.2 <= ge_full.calib.duration()[0] <= 29 * 1.2 + (29 * 42 + 2) * 4e-6   # 29 PE lines x TR
    assert ge_full.imaging_info.esp == pytest.approx(600e-6)


def test_sim_subset_sizes(sim_small):
    """Two slices at +-31.5 mm in both scans; ADC counts match the protocol exactly."""
    n_ro = sim_small.imaging_info.adc.num_samples
    assert np.allclose(np.sort(sim_small.calib_info.positions), [-0.0315, 0.0315])
    assert np.allclose(excited_slice_positions(ImagingProtocol(), [10]), [-0.0315, 0.0315])
    _, _, _, t_adc_img, _ = sim_small.imaging.waveforms_and_times()
    _, _, _, t_adc_cal, _ = sim_small.calib.waveforms_and_times()
    assert len(t_adc_img) == 2 * 1 * 128 * n_ro                # dynamics x groups x lobes x samples
    assert len(t_adc_cal) == 29 * 2 * 17 * n_ro                # PE lines x slices x (3 nav + 14 echoes) x samples
    assert 4.0 <= sim_small.imaging.duration()[0] <= 4.0 + 4 * 4e-6          # TR preserved per dynamic


def test_definitions_needed_by_recon(ge_full):
    img, cal = ge_full.imaging.definitions, ge_full.calib.definitions
    for key in ("TrapezoidGriddingParameters", "TargetGriddedSamples", "Rseg", "Nseg", "Rpe", "PFshift", "SMS",
                "pattern_GE", "GRE_Lines", "FOV", "Resolution", "SlicePositions", "actual_SlicePositions", "TE",
                "AdcDwell", "SliceSep", "NumDynamics", "PlayedSliceGroups"):
        assert key in img, key
    for key in ("ETL", "CNseg", "pe_lines", "PAT_factor", "SlicePositions", "actual_SlicePositions", "AdcDwell"):
        assert key in cal, key
    assert len(np.atleast_1d(img["pattern_GE"])) == 128


def test_write_read_roundtrip_and_comparator(tmp_path):
    """Writing and re-reading a .seq file changes nothing the scanner would play (also validates compare_sequences)."""
    seq, _ = build_imaging(ImagingProtocol(slice_groups=(10,)), PROFILES["siemens_prisma_ref"])
    path = tmp_path / "img.seq"
    seq.write(str(path))
    back = pp.Sequence()
    back.read(str(path))
    d = compare_sequences(back, seq)
    assert d["duration"] < 1e-9
    assert d["n_adc_samples"][0] == d["n_adc_samples"][1]
    assert d["adc_time_maxdiff"] < 1e-9
    assert d["grad_rel_maxdiff"] < 1e-5      # .seq files store amplitudes to ~6 significant digits
    assert d["rf_time_maxdiff"] < 1e-9
    assert d["labels"] == []


def test_counter_labels_at_each_adc():
    """LIN counts 0..127 within each slice group and SLC counts groups (Siemens-style labels)."""
    seq, _ = build_imaging(ImagingProtocol(slice_groups=(9, 10, 11)), PROFILES["siemens_prisma_ref"])
    labels = seq.evaluate_labels(evolution="adc")
    lin = np.asarray(labels["LIN"]).reshape(3, 128)
    assert np.array_equal(lin, np.tile(np.arange(128), (3, 1)))
    assert np.array_equal(np.asarray(labels["SLC"]).reshape(3, 128)[:, 0], [0, 1, 2])


# ---------------------------------------------------------------- Tier 2: parity with the MATLAB scripts
@pytest.mark.parametrize("builder, protocol, ref_name, window", [
    (build_imaging, ImagingProtocol(), "P1_EPTI_sms2.seq", (0, 0.2)),
    (build_calib, CalibProtocol(), "P1_Calib.seq", (0, 0.1)),
])
def test_parity_with_matlab(ref_seq, tmp_path, builder, protocol, ref_name, window):
    ref = ref_seq(ref_name)
    built, _ = builder(protocol, PROFILES["siemens_prisma_ref"])
    built.write(str(tmp_path / ref_name))     # parse ours the same way as the MATLAB file (adds raster definitions)
    ours = pp.Sequence()
    ours.read(str(tmp_path / ref_name))
    d = compare_sequences(ref, ours, time_range=window)
    assert d["duration"] < 1e-6
    assert d["n_adc_samples"][0] == d["n_adc_samples"][1]
    assert d["adc_time_maxdiff"] < 1e-7
    assert d["grad_rel_maxdiff"] < 1e-3
    assert d["rf_time_maxdiff"] < 2e-6
    assert d["definitions"] == []
    full = compare_sequences(ref, ours)
    assert full["labels"] == [], full["labels"]
