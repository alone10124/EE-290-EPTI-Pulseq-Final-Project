"""Reconstruct data from the simple EPI sequence (sequences/scripts/simple_epi.py).

    python recon/simple_epi.py sequences/simple_epi_sim.seq       # reads out/simple_epi_sim/raw.npz

Steps: sort the ADC samples into k-space, reverse the lines read with a negative readout gradient, inverse 2D FFT per
coil, combine the coils. No ghost correction, ramp-sampling regridding or B0 correction yet
(see docs/LESSONS_LEARNED.md). Matrix size, slice positions and TE come from the .seq file's [DEFINITIONS].

Writes ``out/<seq name>/recon.npz`` (complex images [slice, y, x], slices in increasing z) and ``recon.png``. If the
data came from the simulator (``truth.npz`` exists), the figure also shows the image the phantom should give at this
TE, rho * exp(-TE / T2*), and the difference, with its normalized RMS error (NRMSE) inside the head.
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pypulseq as pp
import sigpy as sp

from mrtools.recon.coils import coil_combine, espirit_maps

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]


def sort_kspace(signal, n_slices, nx, ny):
    """ADC samples ``[n_slices * ny * nx, n_coils]`` in acquisition order -> k-space ``[slice, coil, ky, kx]``.

    Slices stay in acquisition order. Odd lines (1, 3, ...) were read with a negative readout gradient, so their
    samples are reversed to put kx in increasing order.
    """
    n_coils = signal.shape[1]
    k = np.asarray(signal).reshape(n_slices, ny, nx, n_coils).transpose(0, 3, 1, 2).copy()
    k[:, :, 1::2] = k[:, :, 1::2, ::-1]
    return k


def reconstruct(kspace, coil_combine_method="espirit"):
    """k-space ``[slice, coil, ky, kx]`` -> images ``[slice, y, x]``.

    ``"espirit"`` estimates coil maps from the center of each slice's own k-space and combines with them (keeps the
    image phase); ``"rss"`` is root-sum-of-squares (magnitude only).
    """
    coil_images = sp.ifft(kspace, axes=(-2, -1))
    if coil_combine_method == "rss":
        return np.sqrt((np.abs(coil_images) ** 2).sum(axis=1))
    if coil_combine_method == "espirit":
        return np.stack([coil_combine(img, espirit_maps(k)) for img, k in zip(coil_images, kspace)])
    raise ValueError(f"unknown coil_combine_method {coil_combine_method!r}")


def plot(images, positions, te, truth, path, label):
    """Recon (and, with ``truth``, the expected image and the difference) for each slice; returns NRMSE per slice."""
    n, cols = len(images), 3 if truth is not None else 1
    fig, axes = plt.subplots(n, cols, figsize=(4 * cols, 4 * n), squeeze=False)
    errors = []
    for i, (img, z) in enumerate(zip(images, positions)):
        mag = np.abs(img)
        axes[i, 0].imshow(mag, cmap="gray", origin="lower")
        axes[i, 0].set_title(f"recon ({label}), z = {z*1e3:.1f} mm")
        if truth is not None:
            j = int(np.argmin(np.abs(truth["slice_positions"] - z)))
            expected = truth["rho"][j] * np.exp(-te / np.maximum(truth["t2s"][j], 1e-6))
            head = truth["rho"][j] > 0
            scale = (mag[head] @ expected[head]) / (mag[head] @ mag[head])   # least-squares intensity match
            err = scale * mag - expected
            errors.append(np.sqrt(np.mean(err[head] ** 2) / np.mean(expected[head] ** 2)))
            vmax = expected.max()
            axes[i, 0].images[0].set_clim(0, vmax / scale)
            axes[i, 1].imshow(expected, cmap="gray", origin="lower", vmin=0, vmax=vmax)
            axes[i, 1].set_title(f"expected: rho * exp(-TE/T2*), TE = {te*1e3:.0f} ms")
            axes[i, 2].imshow(np.abs(err), cmap="magma", origin="lower", vmin=0, vmax=0.5 * vmax)
            axes[i, 2].set_title(f"|difference|, NRMSE {errors[-1]:.1%}")
    for a in axes.ravel():
        a.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=100, bbox_inches="tight")
    plt.close(fig)
    return errors


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq", type=Path, help="the .seq file the data was acquired (or simulated) with")
    parser.add_argument("--raw", type=Path, help="raw data (default: out/<seq name>/raw.npz)")
    parser.add_argument("--coil-combine", choices=["espirit", "rss"], default="espirit")
    args = parser.parse_args(argv)

    out = REPO / "out" / args.seq.stem
    seq = pp.Sequence()
    seq.read(str(args.seq))
    d = seq.definitions
    nx, ny = (int(v) for v in d["MatrixSize"])
    positions = np.atleast_1d(d["SlicePositions"])          # acquisition order
    te = float(d["TE"])

    signal = np.load(args.raw or out / "raw.npz")["signal"]
    images = reconstruct(sort_kspace(signal, len(positions), nx, ny), args.coil_combine)
    order = np.argsort(positions)                           # acquisition order -> increasing z
    images, positions = images[order], positions[order]
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / "recon.npz", images=images, slice_positions=positions)

    truth = np.load(out / "truth.npz") if (out / "truth.npz").exists() else None
    errors = plot(images, positions, te, truth, out / "recon.png", args.coil_combine)
    for z, e in zip(positions, errors):
        print(f"z = {z*1e3:.1f} mm: NRMSE vs expected image {e:.1%}")
    print(f"wrote {out.relative_to(REPO)}/recon.npz and recon.png")


if __name__ == "__main__":
    main()
