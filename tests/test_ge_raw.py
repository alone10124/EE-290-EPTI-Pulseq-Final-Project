"""GE raw data path: tools/ge/ge_raw_export.m -> epti.io.ge_raw -> recon reader.

The export tests use the KomaMRI simulation in out/sim/ (run sim/run_sim.py config/sim_small.yaml first) and are
skipped without it. The P-file round trip also needs MATLAB and the cloned converters in tools/ge/.
"""

import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pytest

from epti.io.ge_raw import is_ge_export, read_ge_export
from epti.io.ismrmrd_io import read_acquisitions, read_calib
from epti.io.seqdefs import read_params
from epti.recon.pipeline import coil_compression_matrix, compress

REPO = Path(__file__).resolve().parents[1]
SIM = REPO / "out" / "sim"
GE = REPO / "tools" / "ge"


def write_export(path, data):
    """Write ``[n_adc, n_coil, n_samples]`` the way ge_raw_export.m does (MATLAB [ns nc nadc] = row-major reversed)."""
    with h5py.File(path, "w") as f:
        f["adc_real"] = data.real.astype(np.float32)
        f["adc_imag"] = data.imag.astype(np.float32)


@pytest.fixture(scope="module")
def sim_calib():
    path = SIM / "sim_small_calib.h5"
    if not path.exists():
        pytest.skip("run sim/run_sim.py config/sim_small.yaml first")
    return read_acquisitions(path)


def test_export_reads_like_ismrmrd(tmp_path, sim_calib):
    exp = tmp_path / "calib_export.h5"
    write_export(exp, sim_calib)
    assert is_ge_export(exp) and not is_ge_export(SIM / "sim_small_calib.h5")
    np.testing.assert_array_equal(read_acquisitions(exp), sim_calib)
    np.testing.assert_array_equal(read_acquisitions(exp, conjugate=True), np.conj(sim_calib))
    p = read_params(REPO / "sequences" / "sim_small_calib.seq")
    assert read_calib(exp, p).shape[:3] == (p.pe_lines, len(p.slice_positions), p.etl)


def test_wrong_adc_count_is_a_clear_error(tmp_path, sim_calib):
    exp = tmp_path / "short.h5"
    write_export(exp, sim_calib[:-1])
    p = read_params(REPO / "sequences" / "sim_small_calib.seq")
    with pytest.raises(ValueError, match="ADC events, but the .seq has"):
        read_calib(exp, p)


def test_coil_compression_keeps_energy_when_not_compressing(sim_calib):
    raw = sim_calib.reshape(-1, 17, *sim_calib.shape[1:])
    nc = raw.shape[-2]
    cc = coil_compression_matrix(raw, nc)                 # full rank: a unitary rotation of the coils
    np.testing.assert_allclose(np.linalg.norm(compress(raw, cc)), np.linalg.norm(raw), rtol=1e-4)
    cc4 = coil_compression_matrix(raw, 4)
    assert compress(raw, cc4).shape[-2] == 4


@pytest.mark.skipif(shutil.which("matlab") is None or not (GE / "PulCeq_tv6").exists(),
                    reason="needs MATLAB and the converters cloned in tools/ge (see tools/ge/README.md)")
def test_pfile_round_trip(tmp_path, sim_calib):
    """Fake P-file (ADCs in the slots PulCeq tv6 assigns, wrapping across 'slices' every 100 views) -> export ->
    same data in the same ADC order. Takes about 2 minutes (MATLAB start-up and the tv6 conversion)."""
    seq = tmp_path / "sim_small_calib.seq"
    shutil.copy(REPO / "sequences" / "sim_small_calib.seq", seq)
    adc = tmp_path / "adc.h5"
    out = tmp_path / "export.h5"
    write_export(adc, sim_calib)
    tar = tmp_path / "sim_small_calib_tv6.tar"
    script = (f"cd('{GE.as_posix()}'); ge_dryrun_tv6({{'{seq.as_posix()}'}}, 100); "
              f"test_raw_export('{adc.as_posix()}', '{tar.as_posix()}', '{out.as_posix()}', 100)")
    res = subprocess.run(["matlab", "-batch", script], capture_output=True, text=True, timeout=900, shell=True)
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    np.testing.assert_array_equal(read_ge_export(out), sim_calib)
