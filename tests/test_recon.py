"""Reconstruction building blocks (epti.recon) and the end-to-end recon on simulated data."""

from types import SimpleNamespace

import numpy as np
import pygrappa  # noqa: F401  (import check)
import pytest
import sigpy as sp
import sigpy.mri as mr
from scipy.signal import windows

from epti.recon import calib, ghost, llr, model, postproc, regrid, subspace

from conftest import REPO


# ---------------------------------------------------------------- Tier 2: MGH P-code / MATLAB fixtures
def test_hybrid_transform_matches_mrir_dft(fixtures):
    x = np.moveaxis(fixtures["dft_in"], 0, -1)
    assert np.allclose(ghost.to_hybrid(x), np.moveaxis(fixtures["idft_freqencode_out"], 0, -1), atol=1e-10)
    assert np.allclose(ghost.to_kspace(x), np.moveaxis(fixtures["fdft_freqencode_out"], 0, -1), atol=1e-12)


def test_ghost_estimate_and_correct_match_pcode(fixtures):
    nav = fixtures["ghost_nav_in"]                                     # [Nx, 3, coil, pe, SEG]
    fwd = 0.5 * (nav[:, 0, :, :, 0] + nav[:, 2, :, :, 0])
    a, b = ghost.estimate(np.moveaxis(fwd, 0, -1), np.moveaxis(nav[:, 1, :, :, 1], 0, -1))
    ref = fixtures["ghost_linear_fit_coeff"]                           # [2, coil]
    assert np.allclose(np.median(a, axis=1), ref[0], atol=1e-6)
    assert np.allclose(np.median(b, axis=1), ref[1], atol=1e-8)
    k = np.moveaxis(fixtures["ghost_data_in"], 0, -1)                  # [echo, coil, pe, SEG, Nx]
    out = k[..., 0, :].copy()
    rev = np.arange(k.shape[0]) % 2 == 1
    out[rev] = ghost.correct(k[rev][..., 1, :], ref[0][None, :, None], ref[1][None, :, None])
    ref_out = fixtures["ghost_deghost_out"]
    assert np.abs(np.moveaxis(out, -1, 0) - ref_out).max() < 1e-10 * np.abs(ref_out).max()


def test_apodization_window_matches_pcode(fixtures):
    assert np.allclose(fixtures["apod_ones72_dim1"][:, 0], windows.tukey(72, 0.1))
    assert np.allclose(fixtures["apod_ones72_dim2"][0, :], windows.tukey(72, 0.1))


def test_subspace_basis_matches_matlab(fixtures):
    """Same span as MATLAB's basis (SVD signs are arbitrary, so compare projectors) and same scaling."""
    tes = np.ravel(fixtures["TEs_all"])
    dt = tes[1] - tes[0]
    phi, s = subspace.basis(dt * np.arange(1, len(tes) + 1), 4)      # relative times, as gen_GE_basis_T2B0 (T0 = 0)
    ref = fixtures["Phi"]
    assert np.allclose(s[:10], np.ravel(fixtures["S"])[:10], rtol=1e-6)
    p_ours = phi @ np.linalg.pinv(phi)
    p_ref = ref @ np.linalg.pinv(ref)
    assert np.abs(p_ours - p_ref).max() < 1e-8
    # the fixture is EPTI_gen_Subspace_Basis_GE's unscaled U; recon_EPTI_P1.m L99 scales it by S^(1/3) afterwards
    assert np.allclose(np.linalg.norm(phi, axis=0), s[:4] ** (1 / 3), rtol=1e-9)


# ---------------------------------------------------------------- Tier 1: known answers
def test_regrid_recovers_uniform_kspace_from_ramp_samples():
    """Sample a known band-limited object along the trapezoid k(t); regridding must return its uniform k-space.

    The object is a few Gaussian blobs (2-pixel width), so its k-space has decayed before the grid edge. With the
    original x grid this fails at ~8% (see the regrid.py docstring); with the uniform grid it is exact.
    """
    nx, n_samp = 72, 216
    p = SimpleNamespace(rise=268e-6, flat=64e-6, fall=268e-6, dwell=2e-6, nx=nx)
    traj_pos, traj_neg = regrid._trajectory(n_samp, p.rise, p.flat, p.fall, p.dwell, nx)
    x_src, amp, sigma = np.array([-0.21, 0.03, 0.17]), np.array([1.0, 0.6, -0.8]), 2 / nx

    def kspace(k):
        return (amp[None] * np.exp(-2 * np.pi ** 2 * sigma ** 2 * k[:, None] ** 2)
                * np.exp(-1j * 2 * np.pi * k[:, None] * x_src[None])).sum(1)

    lines = np.stack([kspace(traj_pos), kspace(traj_neg)])[:, None, :]     # [lobe, coil, sample]; lobe 1 negative
    out = regrid.regrid(lines, np.array([1, -1]), p)[:, 0]
    expected = kspace(np.arange(nx) - nx // 2)

    def image(k):                       # centered inverse FFT, normalized to peak 1 (regrid has its own scaling)
        img = np.abs(np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(k))))
        return img / img.max()

    for lobe in range(2):
        assert np.abs(image(out[lobe]) - image(expected)).max() < 1e-3


def test_grappa_fills_undersampled_calibration():
    n, nc = 72, 8
    img = sp.shepp_logan((n, n)).astype(np.complex64)
    maps = mr.birdcage_maps((nc, n, n))
    k_full = sp.fft(maps * img, axes=(-2, -1))[None]                 # [echo, coil, ky, kx]
    keep = np.zeros(n, bool)
    keep[::3] = True
    keep[n // 2 - 9:n // 2 + 9] = True
    k_us = k_full * keep[None, None, :, None]
    filled = calib.grappa_fill(k_us, np.arange(n // 2 - 9, n // 2 + 9))
    err = np.linalg.norm(filled - k_full) / np.linalg.norm(k_full)
    zero = np.linalg.norm(k_us - k_full) / np.linalg.norm(k_full)
    assert err < 0.05 and err < zero / 5


def test_model_adjoint():
    rng = np.random.default_rng(0)
    nt, k, ns, nc, n = 12, 3, 2, 4, 12
    phi = rng.standard_normal((nt, k))
    phase_t = np.exp(1j * rng.uniform(-3, 3, (nt, ns, n, n)))
    caipi = np.exp(1j * np.pi * rng.integers(0, 2, (nt, ns)))
    sens = rng.standard_normal((ns, nc, n, n)) + 1j * rng.standard_normal((ns, nc, n, n))
    A, _ = model.build(phi, phase_t, caipi, sens, rng.integers(0, n, nt))
    x = (rng.standard_normal(A.ishape) + 1j * rng.standard_normal(A.ishape)).astype(np.complex64)
    y = (rng.standard_normal(A.oshape) + 1j * rng.standard_normal(A.oshape)).astype(np.complex64)
    lhs, rhs = np.vdot(A(x), y), np.vdot(x, A.H(y))
    assert abs(lhs - rhs) < 1e-4 * abs(lhs)


def test_llr_prox():
    rng = np.random.default_rng(0)
    shape = (4, 2, 12, 12)
    x = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    assert np.allclose(llr.LLR(shape, 0.0)(1.0, x), x)                          # lambda = 0: identity
    # equals a direct per-block singular-value soft threshold (6x6 blocks, never across slices)
    out = llr.LLR(shape, 0.7, randshift=False)(2.0, x)                         # threshold = 0.7 * 2.0
    for s in range(2):
        for by in range(2):
            for bx in range(2):
                blk = x[:, s, 6 * by:6 * by + 6, 6 * bx:6 * bx + 6].reshape(4, 36).T
                u, sv, vh = np.linalg.svd(blk, full_matrices=False)
                ref = (u * np.maximum(sv - 1.4, 0)) @ vh
                assert np.allclose(out[:, s, 6 * by:6 * by + 6, 6 * bx:6 * bx + 6].reshape(4, 36).T, ref)
    # rank-1 blocks plus noise below the threshold come out exactly rank 1
    low = np.einsum("k,syx->ksyx", rng.standard_normal(4), rng.standard_normal((2, 12, 12)))
    out = llr.LLR(shape, 0.5, randshift=False)(1.0, low + 0.01 * rng.standard_normal(shape))
    sv = np.linalg.svd(out[:, 0, :6, :6].reshape(4, 36).T, compute_uv=False)
    assert sv[1] < 1e-9 * sv[0]


def test_t2star_fit_recovers_known_decay():
    tes = np.linspace(5e-3, 80e-3, 40)
    t2s_true = np.full((6, 6), 0.045)
    t2s_true[:3] = 0.07
    mag = 100 * np.exp(-tes[:, None, None] / t2s_true[None])
    s0, t2s, _ = postproc.t2star_fit(mag, tes, np.ones((6, 6), bool))
    assert np.allclose(t2s, t2s_true, rtol=1e-6) and np.allclose(s0, 100, rtol=1e-6)


# ---------------------------------------------------------------- end to end on the KomaMRI simulation
def test_recon_on_simulation_matches_truth():
    """Needs out/sim/sim_small_* (python sim/run_sim.py config/sim_small.yaml). ~1-2 min with 10 iterations."""
    paths = [REPO / "out/sim/sim_small_calib.h5", REPO / "out/sim/sim_small_imaging.h5",
             REPO / "sequences/sim_small_calib.seq", REPO / "sequences/sim_small_imaging.seq",
             REPO / "out/sim/sim_small_truth.npz"]
    if not all(p.exists() for p in paths):
        pytest.skip("simulation outputs missing: run make_epti.py and sim/run_sim.py on config/sim_small.yaml")
    import sys
    sys.path.insert(0, str(REPO / "recon"))
    from run_recon import compare_with_truth
    from epti.recon.pipeline import ReconOptions, reconstruct
    from epti.recon.solve import SolverOptions

    res = reconstruct(*paths[:4], ReconOptions(solver=SolverOptions(max_iter=10), verbose=False))
    for row in compare_with_truth(res, np.load(paths[4])):
        assert row["b0_median_abs_err_hz"] < 1.0                     # calib B0 (sign convention b0 = -df)
        assert row["t2s_median_rel_err"] < 0.10
        assert row["corr_own"] > row["corr_other"] + 0.2             # SMS slices not swapped or mixed
