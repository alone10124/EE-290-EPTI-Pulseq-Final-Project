"""Read EPTI raw data: ISMRMRD files (KomaMRI simulations) or GE scanner data exported by tools/ge/ge_raw_export.m.

Replaces ``legacy/libs/utils/ReadRaw_EPTI_calib_pulseq_ISMRMRD.m`` and ``ReadRaw_EPTI_pulseq_ISMRMRD.m``. Like the
originals, acquisitions are ordered purely by arrival and the labels/idx fields are ignored. Unlike them:

* sizes come from the data, not the XML header (KomaMRI's header has no ``receiverChannels`` or ``patientPosition``,
  and the ``ismrmrd`` package's strict parser rejects it);
* lines are NOT flipped here. Each lobe's readout polarity is returned instead, and regridding uses the matching
  trajectory (the original flips negative lines on read, then flips them back inside the regridding).
"""

import ismrmrd
import numpy as np

from epti.io.ge_raw import is_ge_export, read_ge_export

NOISE_FLAG = 1 << (19 - 1)      # ISMRMRD_ACQ_IS_NOISE_MEASUREMENT


def read_acquisitions(path, noise_std=0.0, seed=0, conjugate=False):
    """All non-noise acquisitions as ``[n_acq, n_coil, n_samples]`` complex64.

    ``path`` is an ISMRMRD file or a GE export (``epti.io.ge_raw``); the format is detected from the file.
    ``noise_std`` adds complex Gaussian noise (standard deviation per real/imaginary part) for simulation studies.
    ``conjugate`` conjugates the data, for scanners whose receiver uses the opposite phase convention.
    """
    if is_ge_export(path):
        data = read_ge_export(path)
    else:
        ds = ismrmrd.Dataset(str(path), "dataset", create_if_needed=False)
        acqs = [ds.read_acquisition(i) for i in range(ds.number_of_acquisitions())]
        data = np.stack([a.data for a in acqs if not (int(a.flags) & NOISE_FLAG)]).astype(np.complex64)
    if conjugate:
        data = np.conj(data)
    if noise_std > 0:
        rng = np.random.default_rng(seed)
        data += (noise_std * (rng.standard_normal(data.shape) + 1j * rng.standard_normal(data.shape))).astype(np.complex64)
    return data


def calib_polarity(nav_etl, etl):
    """Readout sign of the calibration scan's lobes: navigators +,-,+ ...; echoes restart at + after the extra
    phase-encoding lobe (seq_EPTI_calib_P1.m L172-186)."""
    nav = [(-1) ** ((i + 1) % 2) for i in range(1, nav_etl + 1)]
    echoes = [(-1) ** ((i + 1) % 2) for i in range(1, etl + 1)]
    return np.array(nav + echoes)


def imaging_polarity(n_lobes):
    """Readout sign of the imaging scan's lobes: +,-,+,- ... continuously (seq_EPTI_imaging_P1.m L193-224)."""
    return np.array([1 if t % 2 == 0 else -1 for t in range(n_lobes)])


def read_calib(path, p, **kw):
    """Calibration k-space ``[pe_line, slice, lobe, coil, sample]`` (acquisition order: PE line, slice, lobe)."""
    data = read_acquisitions(path, **kw)
    n_slc = len(p.slice_positions)
    _check_count(data, p.pe_lines * n_slc * p.etl, path, f"{p.pe_lines} PE lines x {n_slc} slices x {p.etl} lobes")
    return data.reshape(p.pe_lines, n_slc, p.etl, data.shape[1], data.shape[2])


def read_imaging(path, p, **kw):
    """Imaging k-space ``[dynamic, slice_group, lobe, coil, sample]`` (acquisition order: dynamic, group, lobe)."""
    data = read_acquisitions(path, **kw)
    n_groups = len(p.actual_slice_positions)
    per_dyn = n_groups * p.n_lobes
    if data.shape[0] % per_dyn:
        raise ValueError(f"{path}: {data.shape[0]} ADC events is not a whole number of dynamics "
                         f"({n_groups} slice groups x {p.n_lobes} lobes = {per_dyn} per dynamic). "
                         "Does the data match the imaging .seq?")
    n_dyn = data.shape[0] // per_dyn
    return data.reshape(n_dyn, n_groups, p.n_lobes, data.shape[1], data.shape[2])


def _check_count(data, expected, path, what):
    if data.shape[0] != expected:
        raise ValueError(f"{path}: {data.shape[0]} ADC events, but the .seq has {expected} ({what}). "
                         "Does the data match the .seq file (and, for tv6, the .tar used for the export)?")
