"""Ground-truth maps of the simulation phantom (mrtools.sim.phantom). Doesn't start Julia."""

import numpy as np
import pytest

from mrtools.sim.phantom import (SpinPhantom, bandlimited_image, ground_truth_maps, load_truth, reference,
                                 save_truth)


def _uniform_slab(spacing, df_hz=7.0, t2s=0.05):
    g = np.arange(-0.05, 0.05 + 1e-9, spacing)
    x, y = (a.ravel() for a in np.meshgrid(g, g))
    n = len(x)
    return SpinPhantom(x=x, y=y, z=np.zeros(n), rho=np.ones(n), t1=np.ones(n), t2=np.full(n, 0.08),
                       t2s=np.full(n, t2s), df_hz=np.full(n, df_hz), slice_id=np.zeros(n, int))


def test_uniform_slab_gives_flat_maps():
    """2 mm spins in 3 mm voxels: density must be flat (no checkerboard), and T2*/B0 must equal the spin values."""
    truth = ground_truth_maps(_uniform_slab(0.002), (0.216, 0.216), (72, 72), 1)
    inside = truth["rho"][0] > 0
    assert np.allclose(truth["rho"][0][inside], 1.0)
    assert np.allclose(truth["t2s"][0][inside], 0.05)
    assert np.allclose(truth["b0_hz"][0][inside], 7.0)
    # the slab spans +-5 cm, i.e. about 100/3 voxels across
    assert abs(inside.sum(axis=1).max() - 100 / 3) <= 2


def test_truth_grid_matches_fft_pixel_grid():
    """A spin at (x, y) = (0, 0) must land in pixel n // 2, where the centered FFT (sigpy.ifft) puts x = 0, and a spin
    one voxel to the right in pixel n // 2 + 1. (The old grid was off by half a voxel.)"""
    n, fov = 64, 0.2
    x = np.array([0.0, fov / n])
    one = np.ones(2)
    ph = SpinPhantom(x=x, y=np.zeros(2), z=np.zeros(2), rho=np.array([1.0, 2.0]), t1=one, t2=one, t2s=one,
                     df_hz=np.zeros(2), slice_id=np.zeros(2, int))
    rho = ground_truth_maps(ph, (fov, fov), (n, n), 1)["rho"][0]
    assert rho[n // 2, n // 2] == 1.0 and rho[n // 2, n // 2 + 1] == 2.0


def test_bandlimited_image_of_uniform_square():
    """Spins on a fine grid filling a square: the band-limited image integrates to the total spin weight, is ~the spin
    weight inside (up to Gibbs ripple) and ~0 far outside, and is centered where the FFT pixel grid puts x = y = 0."""
    n, fov, spacing, side, w = 32, 0.064, 0.5e-3, 0.04, 0.7
    g = np.arange(-side / 2, side / 2, spacing) + spacing / 2      # square of side 4 cm centered at 0
    x, y = (a.ravel() for a in np.meshgrid(g, g))
    img = bandlimited_image(x, y, np.full(len(x), w), (fov, fov), (n, n)).real
    assert img.sum() * (fov / n) ** 2 == pytest.approx(w * side ** 2, rel=1e-6)
    c = n // 2
    assert np.allclose(img[c - 6:c + 7, c - 6:c + 7], w, rtol=0.05)  # inside: weight +- Gibbs ripple
    assert abs(img[0, 0]) < 0.05 * w
    assert np.allclose(img, img[::-1, ::-1][np.r_[-1, 0:n - 1]][:, np.r_[-1, 0:n - 1]], atol=1e-9)  # symmetric about c


def test_truth_file_round_trip_and_reference(tmp_path):
    """run_sim's truth.npz holds spins; reference() on a recon grid gives the same maps and band-limited image as
    computing them directly, with the recon's own weighting. An old-format (maps) file gives None."""
    spins = _uniform_slab(0.002)
    spins.is_fat = spins.x > 0.03
    rss = 1 + spins.y
    save_truth(tmp_path / "truth.npz", spins, [1.5e-3], rss)
    loaded, rss2, positions = load_truth(tmp_path / "truth.npz")
    assert np.array_equal(loaded.x, spins.x) and np.array_equal(loaded.is_fat, spins.is_fat)
    assert np.array_equal(rss2, rss) and positions.tolist() == [1.5e-3]

    fov, n = (0.216, 0.216), (72, 72)
    ref = reference(tmp_path / "truth.npz", fov, n, lambda s, c: s.rho * c * ~s.is_fat)
    direct = ground_truth_maps(spins, fov, n, 1)
    assert all(np.array_equal(ref[k], direct[k]) for k in ("rho", "t2s", "b0_hz"))
    w = spins.rho * rss * ~spins.is_fat
    assert np.allclose(ref["expected"][0], np.abs(bandlimited_image(spins.x, spins.y, w, fov, n)))

    np.savez(tmp_path / "old.npz", rho=direct["rho"])
    assert load_truth(tmp_path / "old.npz") is None and load_truth(tmp_path / "missing.npz") is None
