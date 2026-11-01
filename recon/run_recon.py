"""Reconstruct EPTI data (core path) and write results, NIfTI files and a summary figure.

    python recon/run_recon.py config/sim_small.yaml                 # uses sequences/<config>_*.seq, out/sim/<config>_*.h5
    python recon/run_recon.py config/sim_small.yaml --iters 20      # faster, slightly less converged
    python recon/run_recon.py config/sim_small.yaml --noise 200     # add complex noise to the raw data first

Writes out/recon/<config>_recon.npz, out/recon/<config>_*.nii.gz and out/recon/<config>_summary.png. If simulation
ground truth exists (out/sim/<config>_truth.npz), it also prints and plots the comparison with it.
"""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from epti.io import nifti  # noqa: E402
from epti.io.seqdefs import read_params  # noqa: E402
from epti.recon.pipeline import ReconOptions, reconstruct  # noqa: E402
from epti.recon.solve import SolverOptions  # noqa: E402

REPO = Path(__file__).resolve().parents[1]


def compare_with_truth(res, truth):
    """Per-slice metrics against the phantom ground truth (simulation only)."""
    rows = []
    for s, z in enumerate(res.slice_positions):
        ts = int(np.argmin(np.abs(truth["slice_positions"] - z)))
        others = [i for i in range(len(truth["slice_positions"])) if i != ts]
        m = res.mask[s] & (truth["rho"][ts] > 0.5)
        rel = np.abs(res.t2s[s][m] - truth["t2s"][ts][m]) / truth["t2s"][ts][m]
        cs = int(np.argmin(np.abs(res.calib.positions - z)))
        mb = res.calib.mask[cs] & (truth["rho"][ts] > 0)
        b0_err = np.abs(res.calib.b0_hz[cs][mb] + truth["b0_hz"][ts][mb])            # b0 = -df (KomaMRI)
        comb = res.combined[0, s].ravel()
        rows.append(dict(
            z_mm=z * 1e3, t2s_median_rel_err=float(np.median(rel)), t2s_iqr=tuple(np.percentile(rel, [25, 75])),
            b0_median_abs_err_hz=float(np.median(b0_err)),
            corr_own=float(np.corrcoef(comb, truth["rho"][ts].ravel())[0, 1]),
            corr_other=max([float(np.corrcoef(comb, truth["rho"][o].ravel())[0, 1]) for o in others], default=np.nan)))
    return rows


def figure(res, truth, path, title):
    n = len(res.slice_positions)
    cols = 6 if truth is not None else 4
    fig, ax = plt.subplots(n, cols, figsize=(3.1 * cols, 3.0 * n), squeeze=False)
    e_early, e_late = 10, len(res.tes) - 10
    for s, z in enumerate(res.slice_positions):
        cs = int(np.argmin(np.abs(res.calib.positions - z)))
        panels = [
            (np.abs(res.images[0, s, e_early]), f"echo TE={res.tes[e_early]*1e3:.1f} ms", "gray", None),
            (np.abs(res.images[0, s, e_late]), f"echo TE={res.tes[e_late]*1e3:.1f} ms", "gray", None),
            (res.combined[0, s], "optimal echo combination", "gray", None),
            (res.t2s[s] * 1e3, "recon T2* (ms)", "magma", (0, 120)),
        ]
        if truth is not None:
            ts = int(np.argmin(np.abs(truth["slice_positions"] - z)))
            panels += [(truth["t2s"][ts] * 1e3 * (truth["rho"][ts] > 0), "truth T2* (ms)", "magma", (0, 120)),
                       (res.calib.b0_hz[cs], "calib B0 (Hz)", "RdBu_r", (-50, 50))]
        for a, (img, t, cmap, lim) in zip(ax[s], panels):
            im = a.imshow(img, cmap=cmap, origin="lower", vmin=None if lim is None else lim[0],
                          vmax=None if lim is None else lim[1])
            a.set_title(f"z={z*1e3:+.1f} mm: {t}", fontsize=8)
            a.axis("off")
            if lim is not None:
                fig.colorbar(im, ax=a, fraction=0.046)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=110)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path, help="config used to make the sequences (only its name is used)")
    parser.add_argument("--seqdir", type=Path, default=REPO / "sequences")
    parser.add_argument("--datadir", type=Path, default=REPO / "out" / "sim")
    parser.add_argument("--outdir", type=Path, default=REPO / "out" / "recon")
    parser.add_argument("--iters", type=int, default=70, help="ADMM iterations (original: 70)")
    parser.add_argument("--lamda", type=float, default=1e-3, help="LLR weight (0 = no regularization)")
    parser.add_argument("--noise", type=float, default=0.0, help="std of complex noise added to the raw data")
    parser.add_argument("--device", type=int, default=-1, help="-1 = CPU, 0 = first CUDA GPU (needs cupy)")
    args = parser.parse_args(argv)

    stem = args.config.stem
    seq_c, seq_i = args.seqdir / f"{stem}_calib.seq", args.seqdir / f"{stem}_imaging.seq"
    h5_c, h5_i = args.datadir / f"{stem}_calib.h5", args.datadir / f"{stem}_imaging.h5"
    opts = ReconOptions(noise_std=args.noise,
                        solver=SolverOptions(max_iter=args.iters, lamda=args.lamda, device=args.device))
    res = reconstruct(h5_c, h5_i, seq_c, seq_i, opts)
    print("timings (s):", {k: round(v, 1) for k, v in res.timings.items()})

    args.outdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.outdir / f"{stem}_recon.npz", slice_positions=res.slice_positions, tes=res.tes,
                        images=res.images, t2s=res.t2s, s0=res.s0, combined=res.combined, mask=res.mask,
                        b0_hz=res.calib.b0_hz, sens=res.calib.sens, drift_hz=res.drift_hz)
    p = read_params(seq_i)
    vox = (p.fov[0] / p.nx, p.fov[1] / (p.rseg * p.nseg))
    nifti.save(args.outdir / f"{stem}_t2s.nii.gz", res.t2s * 1e3, vox, res.slice_positions)
    nifti.save(args.outdir / f"{stem}_combined.nii.gz", np.moveaxis(res.combined, 0, -1), vox, res.slice_positions)
    nifti.save(args.outdir / f"{stem}_echoes_dyn0.nii.gz", np.moveaxis(np.abs(res.images[0]), 1, -1), vox,
               res.slice_positions)

    truth_path = args.datadir / f"{stem}_truth.npz"
    truth = np.load(truth_path) if truth_path.exists() else None
    if truth is not None:
        for r in compare_with_truth(res, truth):
            print(f"z={r['z_mm']:+.1f} mm: T2* median error {r['t2s_median_rel_err']*100:.1f}% "
                  f"(IQR {r['t2s_iqr'][0]*100:.1f}-{r['t2s_iqr'][1]*100:.1f}%), calib B0 median error "
                  f"{r['b0_median_abs_err_hz']:.2f} Hz, correlation with own slice {r['corr_own']:.2f} / other {r['corr_other']:.2f}")
    figure(res, truth, args.outdir / f"{stem}_summary.png",
           f"{stem}: EPTI subspace recon ({args.iters} ADMM iterations, lambda={args.lamda:g}, noise={args.noise:g})")
    print(f"wrote {args.outdir}/{stem}_*")
    return 0


if __name__ == "__main__":
    sys.exit(main())
