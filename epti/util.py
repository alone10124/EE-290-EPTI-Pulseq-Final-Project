"""Small numerical helpers shared across the package."""

import numpy as np


def mround(x):
    """Round half away from zero, like MATLAB's ``round`` (Python/NumPy round half to even).

    The MGH scripts use ``round`` when snapping times to rasters, so ports must use this to reproduce them exactly.
    """
    x = np.asarray(x, dtype=float)
    out = np.sign(x) * np.floor(np.abs(x) + 0.5)
    return out.item() if out.ndim == 0 else out


def ceil_to_raster(t, raster, tol=1e-9):
    """Round a duration up to a multiple of ``raster``, ignoring floating-point noise below ``tol`` raster units."""
    return np.ceil(t / raster - tol) * raster


def round_to_raster(t, raster):
    """Round a duration to the nearest multiple of ``raster`` (MATLAB-style half away from zero)."""
    return mround(t / raster) * raster
