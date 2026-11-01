"""NIfTI output with nibabel.

Replaces ``legacy/libs/utils/saveniftidata_EPTI_pulseq.m`` (which needs the Siemens protocol header for its geometry,
so it can't work for GE or simulated data). The affine is built from the ``.seq`` definitions instead: in-plane voxel
size from FOV/matrix, slice positions from ``SlicePositions``, with the image centered on the scanner isocenter
(the sequence has no in-plane offsets). Axes follow the recon arrays (x = readout, y = phase encoding, z = slice), in
scanner coordinates. The original also wrote TR in mm units by mistake (``pixdim(5)``); here the time unit is seconds.
"""

import nibabel as nib
import numpy as np


def save(path, volume, voxel_size, slice_positions, tr=None):
    """Write ``volume [slice, y, x]`` or ``[slice, y, x, t]`` (t = echoes or dynamics) as NIfTI.

    ``voxel_size`` = (dx, dy) in m; ``slice_positions`` ascending z in m. Unevenly spaced slices (e.g. a simulated
    subset) are written with the median spacing, and a warning is printed, since NIfTI needs one spacing.
    """
    vol = np.moveaxis(np.asarray(volume), 0, 2)                     # -> [y, x, z, ...]
    vol = np.swapaxes(vol, 0, 1)                                    # -> [x, y, z, ...]
    z = np.asarray(slice_positions)
    dz = float(np.median(np.diff(z))) if len(z) > 1 else 3e-3
    if len(z) > 2 and np.ptp(np.diff(z)) > 1e-6:
        print(f"warning: {path}: uneven slice spacing, writing median spacing {dz * 1e3:.2f} mm")
    nx, ny = vol.shape[:2]
    dx, dy = voxel_size
    aff = np.diag([dx * 1e3, dy * 1e3, dz * 1e3, 1.0])
    aff[:3, 3] = [-(nx / 2) * dx * 1e3, -(ny / 2) * dy * 1e3, z[0] * 1e3]
    img = nib.Nifti1Image(vol.astype(np.float32), aff)
    img.header.set_xyzt_units("mm", "sec")
    if tr is not None and vol.ndim == 4:
        img.header["pixdim"][4] = tr
    nib.save(img, str(path))
