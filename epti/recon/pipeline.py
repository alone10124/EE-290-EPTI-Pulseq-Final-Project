"""End-to-end EPTI reconstruction (core path): calibration -> per-dynamic subspace recon -> T2* and combined images.

Follows ``legacy/recon/recon_EPTI_P1.m`` with ``B0_update = 0`` and ``moco_recon = 0`` (the high-order B0 update and
motion-aware recon are Phase 4). Line numbers below refer to that script.
"""

import time
from dataclasses import dataclass, field

import numpy as np

from epti.io.ismrmrd_io import calib_polarity, imaging_polarity, read_calib, read_imaging
from epti.io.seqdefs import read_params
from epti.recon import calib, drift, ghost, model, postproc, regrid, subspace
from epti.recon.solve import SolverOptions, solve


@dataclass
class ReconOptions:
    k: int = 4                       # subspace size (L96)
    echo_cut: tuple = (2, 0)         # lobes dropped from the start/end of the output (L106-107)
    apodize: float = 0.1             # Tukey factor (L56-57); 0 = off
    drift_correction: bool = True    # navigator-based global B0 drift (L183-185)
    noise_std: float = 0.0           # add complex Gaussian noise to the raw data (simulation studies)
    adc_phase_sign: int = 1          # ADC phase convention for the CAIPI model (+1: Pulseq/KomaMRI); see model.py
    solver: SolverOptions = field(default_factory=SolverOptions)
    verbose: bool = True


@dataclass
class ReconResult:
    slice_positions: np.ndarray      # [slice] m, ascending
    tes: np.ndarray                  # [echo] s, after the echo cut
    images: np.ndarray               # [dynamic, slice, echo, y, x] complex echo images
    t2s: np.ndarray                  # [slice, y, x] s
    s0: np.ndarray                   # [slice, y, x]
    combined: np.ndarray             # [dynamic, slice, y, x] optimally combined
    mask: np.ndarray                 # [slice, y, x]
    calib: calib.CalibResult
    drift_hz: np.ndarray             # [dynamic, group]
    timings: dict


def group_slices(p_img):
    """Positions of the slices excited by each played slice group: ``[group, mb]`` (getsmspulse.m band offsets)."""
    centers = np.asarray(p_img.actual_slice_positions)
    bands = (-p_img.mb / 2 + 0.5 + np.arange(p_img.mb)) * p_img.slice_sep
    return centers[:, None] + bands[None, :]


def _match(positions, z, tol=1e-5):
    i = int(np.argmin(np.abs(np.asarray(positions) - z)))
    if abs(positions[i] - z) > tol:
        raise ValueError(f"no calibration slice at z = {z * 1e3:.2f} mm (have {np.round(np.asarray(positions) * 1e3, 2)})")
    return i


def reconstruct(calib_h5, imaging_h5, calib_seq, imaging_seq, opts: ReconOptions = None) -> ReconResult:
    opts = opts or ReconOptions()
    log = print if opts.verbose else (lambda *a, **k: None)
    tim = {}
    t0 = time.time()
    pc, pi = read_params(calib_seq), read_params(imaging_seq)
    if pi.slice_sep is None:
        pi.slice_sep = (pi.thickness + float(pi.extra.get("SliceGap", 0))) * len(pi.slice_positions) / pi.mb
    ny = pi.rseg * pi.nseg
    zs = group_slices(pi)                                                # [group, mb]
    group_of_calib_slice = np.array([int(np.argmin(np.abs(zs - z).min(axis=1))) for z in pc.slice_positions])

    # --- calibration (L61-88)
    raw = read_calib(calib_h5, pc, noise_std=opts.noise_std, seed=1)
    pol_c = calib_polarity(3, pc.etl - 3)
    kc = regrid.regrid(raw, pol_c, pc)
    ga, gb = ghost.calib_coefficients(kc[:, :, :3], group_of_calib_slice)        # [group, coil]
    ech = kc[:, :, 3:].copy()
    rev = pol_c[3:] < 0
    for s, g in enumerate(group_of_calib_slice):
        ech[:, s][:, rev] = ghost.correct(ech[:, s][:, rev], ga[g], gb[g])
    cal = calib.process(ech, pc, ny)
    tim["calibration"] = time.time() - t0
    log(f"calibration: {len(pc.slice_positions)} slices, ghost slope {np.median(gb)*pi.nx:.2f} rad/FOV, "
        f"{tim['calibration']:.1f} s")

    # --- imaging setup
    tes = subspace.imaging_echo_times(pi)
    phi, _ = subspace.basis(subspace.basis_times(pi), opts.k)
    pol_i = imaging_polarity(pi.n_lobes)
    ky_idx = model.ky_indices(pi, ny)
    raw = read_imaging(imaging_h5, pi, noise_std=opts.noise_std, seed=2)
    n_dyn, n_grp = raw.shape[:2]
    ki = regrid.regrid(raw, pol_i, pi)                                           # [dyn, group, t, coil, kx]
    for g in range(n_grp):
        ki[:, g][:, pol_i < 0] = ghost.correct(ki[:, g][:, pol_i < 0], ga[g], gb[g])

    cut = slice(opts.echo_cut[0], pi.n_lobes - opts.echo_cut[1])
    all_z = np.sort(zs.ravel())
    images = np.zeros((n_dyn, len(all_z), cut.stop - cut.start, ny, pi.nx), dtype=np.complex64)
    drifts = np.zeros((n_dyn, n_grp))
    t_solve = 0.0
    for g in range(n_grp):
        cs = [_match(cal.positions, z) for z in zs[g]]
        sens, b0 = cal.sens[cs], cal.b0_hz[cs]
        caipi = model.caipi_phase(pi, zs[g], pi.actual_slice_positions[g], adc_phase_sign=opts.adc_phase_sign)
        for d in range(n_dyn):
            if opts.drift_correction and d > 0:
                drifts[d, g] = drift.estimate(ki[0, g, :3], ki[d, g, :3], pi.esp)
            phase_t = model.time_phase(b0 + drifts[d, g], cal.phase0_odd[cs], cal.phase0_even[cs], tes, pol_i)
            A, mask = model.build(phi, phase_t, caipi, sens, ky_idx)
            y = np.zeros(A.oshape, dtype=np.complex64)
            y[np.arange(pi.n_lobes), :, ky_idx, :] = ki[d, g]                    # each lobe's line at its ky
            t1 = time.time()
            coeff = solve(A, y, A.ishape, opts.solver)
            t_solve += time.time() - t1
            img = model.images(coeff, phi, phase_t)[cut]                         # [t, s, y, x]
            if opts.apodize:
                img = postproc.apodize(img, opts.apodize)
            for s, z in enumerate(zs[g]):
                images[d, int(np.searchsorted(all_z, z))] = img[:, s]
            log(f"group {g} (z = {np.round(zs[g] * 1e3, 1)} mm) dynamic {d}: drift {drifts[d, g]:+.2f} Hz, "
                f"solve {time.time() - t1:.1f} s")
    tim["subspace solves"] = t_solve

    # --- T2* and optimal echo combination on the dynamics' mean magnitude (optimalCombineEchoes.m)
    tes_out = tes[cut]
    mag_mean = np.abs(images).mean(axis=0)                                       # [slice, echo, y, x]
    t2s, s0, masks, combined = [], [], [], np.zeros(images.shape[:2] + images.shape[3:])
    for s in range(len(all_z)):
        msk = postproc.fit_mask(np.sqrt((mag_mean[s] ** 2).sum(0)))
        a, b, _ = postproc.t2star_fit(mag_mean[s], tes_out, msk)
        s0.append(a), t2s.append(b), masks.append(msk)
        for d in range(n_dyn):
            combined[d, s] = postproc.combine_echoes(np.abs(images[d, s]), tes_out, b)
    tim["total"] = time.time() - t0
    return ReconResult(slice_positions=all_z, tes=tes_out, images=images, t2s=np.stack(t2s), s0=np.stack(s0),
                       combined=combined, mask=np.stack(masks), calib=cal, drift_hz=drifts, timings=tim)
