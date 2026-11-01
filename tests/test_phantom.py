"""Ground-truth maps of the simulation phantom (sim/phantom.py). Doesn't start Julia."""

import sys

import numpy as np

from conftest import REPO

sys.path.insert(0, str(REPO / "sim"))
from phantom import SpinPhantom, ground_truth_maps  # noqa: E402


def _uniform_slab(spacing, df_hz=7.0, t2s=0.05):
    g = np.arange(-0.05, 0.05 + 1e-9, spacing)
    x, y = (a.ravel() for a in np.meshgrid(g, g))
    n = len(x)
    return SpinPhantom(x=x, y=y, z=np.zeros(n), rho=np.ones(n), t1=np.ones(n), t2=np.full(n, 0.08),
                       t2s=np.full(n, t2s), df_hz=np.full(n, df_hz), slice_id=np.zeros(n, int))


def test_uniform_slab_gives_flat_maps():
    """2 mm spins in 3 mm voxels: density must be flat (no checkerboard), and T2*/B0 must equal the spin values."""
    truth = ground_truth_maps(_uniform_slab(0.002), (0.216, 0.216), (72, 72), [0.0])
    inside = truth["rho"][0] > 0
    assert np.allclose(truth["rho"][0][inside], 1.0)
    assert np.allclose(truth["t2s"][0][inside], 0.05)
    assert np.allclose(truth["b0_hz"][0][inside], 7.0)
    # the slab spans +-5 cm, i.e. about 100/3 voxels across
    assert abs(inside.sum(axis=1).max() - 100 / 3) <= 2
