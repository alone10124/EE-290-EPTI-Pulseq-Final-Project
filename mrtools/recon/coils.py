"""Coil sensitivity maps from calibration data (ESPIRiT), and coil combination.

The calibration data are a fully sampled k-space of each slice, e.g. single-band EPI with the same readout as the
imaging, so the maps share its distortion. ESPIRiT (Uecker et al., MRM 2014) estimates each coil's sensitivity from the
central ``calib_width`` lines; ``sigpy.mri.app.EspiritCalib`` does the work. Because the data are fully sampled, no
GRAPPA step is needed (the EPTI calibration needed one only because it skipped outer k-space lines; see
docs/LESSONS_LEARNED.md).

Conventions: k-space and images are ``[coil, y, x]``, using SigPy's centered, orthonormal FFT.
"""

import numpy as np
import sigpy as sp
import sigpy.mri as mr


def espirit_maps(kspace, calib_width=24, thresh=0.02, crop=0.95, device=-1):
    """Sensitivity maps ``[coil, y, x]`` from one slice's fully sampled multi-coil k-space ``[coil, ky, kx]``.

    ``thresh`` and ``crop`` are SigPy's defaults: the kernel singular-value threshold, and the eigenvalue below which a
    voxel's maps are set to zero (outside the object).
    """
    return mr.app.EspiritCalib(np.asarray(kspace), calib_width=calib_width, thresh=thresh, crop=crop,
                               device=sp.Device(device), show_pbar=False).run()


def coil_combine(images, maps):
    """SENSE-1 coil combination of ``images [coil, y, x]``: ``sum(conj(S) * I) / sum(|S|^2)``."""
    num = np.sum(np.conj(maps) * images, axis=0)
    den = np.sum(np.abs(maps) ** 2, axis=0)
    return np.divide(num, den, out=np.zeros_like(num), where=den > 1e-12)
