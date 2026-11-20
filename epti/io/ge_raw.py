"""GE scanner raw data, after ``tools/ge/ge_raw_export.m`` has put it in ADC order.

The export is an HDF5 file with ``/adc_real`` and ``/adc_imag``, shape ``[n_adc, n_coil, n_samples]``, one row per
ADC event in the order the events appear in the ``.seq`` file. That is the same layout
``epti.io.ismrmrd_io.read_acquisitions`` returns for KomaMRI data, so the recon treats both the same way.

Reading GE formats directly in Python isn't practical: P-file headers change with each GE software release (the
export uses TOPPE's maintained ``loadpfile``), and ScanArchive files need GE's Orchestra SDK.
"""

import h5py
import numpy as np


def is_ge_export(path) -> bool:
    """True if ``path`` is an HDF5 file written by ``ge_raw_export.m`` (not an ISMRMRD dataset)."""
    try:
        with h5py.File(path, "r") as f:
            return "adc_real" in f and "adc_imag" in f
    except OSError:
        return False


def read_ge_export(path) -> np.ndarray:
    """``[n_adc, n_coil, n_samples]`` complex64, in ADC order."""
    with h5py.File(path, "r") as f:
        return (f["adc_real"][()] + 1j * f["adc_imag"][()]).astype(np.complex64)
