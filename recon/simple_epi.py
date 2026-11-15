"""Reconstruct data from the simple EPI sequence (sequences/scripts/simple_epi.py).

    python recon/simple_epi.py configs/simple_epi.toml     # reads out/simple_epi_sim/raw.npz

Which data: ``reconstruction.seq`` / ``reconstruction.raw`` in the config, by default what the simulator wrote.
Option: ``reconstruction.coil_combine`` = "espirit" or "rss".

Steps: sort the ADC samples into k-space, reverse the lines read with a negative readout gradient, inverse 2D FFT per
coil, combine the coils. No ghost correction, ramp-sampling regridding, B0 correction or coil-intensity normalization
yet. Matrix size, slice positions and TE come from the .seq file's [DEFINITIONS].

Writes ``recon.npz`` (complex images [slice, y, x], slices in increasing z) and ``recon.png`` next to the raw data.
If the data came from the simulator (``truth.npz`` with ``expected``), the figure also shows the image an ideal
Cartesian scan would give (``rho * exp(-TE/T2*)``, band-limited to the same k-space and shaded by the same receive
coils; see sim/run_sim.py) and the difference, with its normalized RMS error (NRMSE) inside the head. What's left in
that error is EPI physics (B0 distortion, T2* decay during the readout) plus any recon error.
"""

import argparse

import matplotlib
import numpy as np
import pypulseq as pp
import sigpy as sp

from mrtools import config
from mrtools.recon.coils import coil_combine, espirit_maps

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


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
            expected = truth["expected"][j]
            head = truth["rho"][j] > 0
            scale = (mag[head] @ expected[head]) / (mag[head] @ mag[head])   # least-squares intensity match
            err = scale * mag - expected
            errors.append(np.sqrt(np.mean(err[head] ** 2) / np.mean(expected[head] ** 2)))
            vmax = expected.max()
            axes[i, 0].images[0].set_clim(0, vmax / scale)
            axes[i, 1].imshow(expected, cmap="gray", origin="lower", vmin=0, vmax=vmax)
            axes[i, 1].set_title(f"expected (ideal scan), TE = {te*1e3:.0f} ms")
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
    parser.add_argument("config", help="run config, e.g. configs/simple_epi.toml")
    cfg = config.load(parser.parse_args(argv).config)
    method = cfg.reconstruction.get("coil_combine", "espirit")

    out = cfg.recon_raw.parent
    seq = pp.Sequence()
    seq.read(str(cfg.recon_seq))
    d = seq.definitions
    nx, ny = (int(v) for v in d["MatrixSize"])
    positions = np.atleast_1d(d["SlicePositions"])          # acquisition order
    te = float(d["TE"])

    signal = np.load(cfg.recon_raw)["signal"]
    images = reconstruct(sort_kspace(signal, len(positions), nx, ny), method)
    order = np.argsort(positions)                           # acquisition order -> increasing z
    images, positions = images[order], positions[order]
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / "recon.npz", images=images, slice_positions=positions)

    truth = np.load(out / "truth.npz") if (out / "truth.npz").exists() else None
    if truth is not None and "expected" not in truth:
        truth = None                                        # older truth.npz without the expected image
    errors = plot(images, positions, te, truth, out / "recon.png", method)
    for z, e in zip(positions, errors):
        print(f"z = {z*1e3:.1f} mm: NRMSE vs expected image {e:.1%}")
    print(f"wrote {out.relative_to(config.REPO)}/recon.npz and recon.png")


if __name__ == "__main__":
    main()
