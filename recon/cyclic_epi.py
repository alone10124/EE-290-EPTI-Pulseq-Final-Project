"""Reconstruct data from the cyclic EPI sequence (sequences/scripts/cyclic_epi.py): one image per echo, then T2* and B0.

    python recon/cyclic_epi.py configs/cyclic_epi.toml     # reads out/cyclic_epi_sim/raw.npz

Which data: ``reconstruction.seq`` / ``reconstruction.raw`` in the config, by default what the simulator wrote.
Option: ``reconstruction.coil_combine`` = "espirit" or "rss" (magnitude only, so no B0 map).

Steps, per slice:
1. Sort the ADC samples into [shot, echo, coil, kx] and reverse kx on odd echoes (negative readout lobes).
2. Regroup by echo: echo j of shot s is ky row (s + j) mod ny, so echo j from all shots is a full k-space.
3. Inverse FFT per coil and combine the coils (ESPIRiT maps from echo 0, which has the most signal and no late rows).
4. Fit each pixel: |I_j| = |S0| exp(-TE_j / T2*) (log-linear, magnitude-weighted) and B0 from the phase change between
   consecutive echoes, with TE_j = TE + j * EchoSpacing from the .seq [DEFINITIONS].

The rows that come after the wrap blip (rows 0 ... j-1 of echo j, ``WrapDelay`` late; 0 with ``wrap_timing = "pad"``)
are treated as on time; the draft's delay-aware solve (cyclic_epi_draft/recon_cyclic_epi.py) is not used.

B0 sign: KomaMRI precesses an off-resonance of +df Hz as exp(-i 2 pi df t) (see mrtools.sim.phantom), so the phase
*decreases* with TE where B0 is positive: B0 = -(phase change per echo) / (2 pi EchoSpacing). The same sign is assumed
for scanner data; check it on the first scan. B0 is unambiguous for |B0| < 1 / (2 EchoSpacing) (~600 Hz).

Writes ``recon.npz`` (``echo_images`` [slice, echo, y, x], ``te``, ``s0``, ``t2s`` (s), ``b0_hz``, ``support``; slices in
increasing z), ``recon.png`` and ``decay.png`` next to the raw data. If the data came from the simulator (``truth.npz``),
the figure compares echo 0 with the image an ideal Cartesian scan at that TE would give (as recon/simple_epi.py) and
the fitted maps with the phantom's T2* and B0, and prints the errors.

T2* and B0 fitting adapted from the team's draft recon_cyclic_epi.py (``fit_maps``, 2026-11); changes: B0 sign as
above, B0 in Hz, slices and coils.
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
    """ADC samples ``[n_shots * n_slices * ny * nx, n_coils]`` (acquisition order) -> k-space
    ``[slice, echo, coil, ky, kx]``, slices in acquisition order.

    The sequence loops shot -> slice -> echo -> sample; shot s, echo j is row (s + j) mod ny. Odd echoes were read with
    a negative readout gradient, so their samples are reversed to put kx in increasing order.
    """
    n_coils = signal.shape[1]
    d = np.asarray(signal, dtype=np.complex64).reshape(ny, n_slices, ny, nx, n_coils)   # [shot, slice, echo, kx, coil]
    d = d.transpose(1, 0, 2, 4, 3).copy()                                               # [slice, shot, echo, coil, kx]
    d[:, :, 1::2] = d[:, :, 1::2, :, ::-1]
    k = np.empty((n_slices, ny, n_coils, ny, nx), np.complex64)
    shots = np.arange(ny)
    for j in range(ny):
        k[:, j, :, (shots + j) % ny] = d[:, :, j].transpose(1, 0, 2, 3)                 # [row, slice, coil, kx]
    return k


def combine(kspace, method="espirit"):
    """k-space ``[slice, echo, coil, ky, kx]`` -> images ``[slice, echo, y, x]``.

    ``"espirit"``: SENSE-1 combination with ESPIRiT maps estimated from each slice's echo 0 (complex, keeps the phase).
    ``"rss"``: root-sum-of-squares (magnitude only).
    """
    out = np.empty(kspace.shape[:2] + kspace.shape[3:], np.complex64)
    for z, k in enumerate(kspace):
        coil_images = sp.ifft(k, axes=(-2, -1))                                          # [echo, coil, y, x]
        if method == "rss":
            out[z] = np.sqrt((np.abs(coil_images) ** 2).sum(axis=1))
        elif method == "espirit":
            maps = espirit_maps(k[0])
            out[z] = np.stack([coil_combine(img, maps) for img in coil_images])
        else:
            raise ValueError(f"unknown coil_combine method {method!r}")
    return out


def fit_maps(images, te, threshold=0.1):
    """Per-pixel fit of echo images ``[echo, y, x]`` at echo times ``te`` (s).

    Returns ``s0`` (complex, at TE = 0), ``t2s`` (s; 0 outside the support or where the fit gives no decay), ``b0_hz``
    (``None`` for magnitude images) and ``support`` (mean magnitude > ``threshold`` x its maximum).

    T2*: weighted least squares of log|I| = log|S0| - TE / T2*, weights |I|^2 (late, noisy echoes count less).
    B0: angle of sum_j I_{j+1} conj(I_j) over all echo pairs, divided by -2 pi * echo spacing (KomaMRI sign, see the
    module doc).
    """
    mag = np.abs(images).astype(float) + 1e-12
    support = mag.mean(0) > threshold * mag.mean(0).max()
    w = mag ** 2
    t = te[:, None, None]
    log = np.log(mag)
    sw, st, stt = w.sum(0), (w * t).sum(0), (w * t * t).sum(0)
    sl, stl = (w * log).sum(0), (w * t * log).sum(0)
    r2s = -(sw * stl - st * sl) / (sw * stt - st ** 2 + 1e-30)
    log_s0 = (sl + r2s * st) / sw
    t2s = np.where(support & (r2s > 0), 1 / np.maximum(r2s, 1e-9), 0.0)
    if np.iscomplexobj(images):
        dt = float(np.mean(np.diff(te)))
        rate = np.angle((images[1:] * np.conj(images[:-1])).sum(0)) / dt              # rad/s of phase change
        b0_hz = np.where(support, -rate / (2 * np.pi), 0.0)
        s0 = np.exp(log_s0) * np.exp(1j * (np.angle(images[0]) - rate * te[0]))
    else:
        b0_hz, s0 = None, np.exp(log_s0)
    return dict(s0=np.where(support, s0, 0), t2s=t2s, b0_hz=b0_hz, support=support)


def interior(mask, n=2):
    """``mask`` eroded by ``n`` pixels (drops partial-volume edges)."""
    m = mask.copy()
    for _ in range(n):
        m &= np.roll(m, 1, 0) & np.roll(m, -1, 0) & np.roll(m, 1, 1) & np.roll(m, -1, 1)
    return m


def compare_with_truth(images, fits, truth, j):
    """Errors against the simulator's ground truth for slice ``j`` of ``truth``: NRMSE of echo 0 vs the expected
    image, and median |error| of T2* (ms) and B0 (Hz) inside the head (eroded by 2 pixels)."""
    expected = truth["expected"][j]
    head = truth["rho"][j] > 0
    mag = np.abs(images[0])
    scale = (mag[head] @ expected[head]) / (mag[head] @ mag[head])                   # least-squares intensity match
    nrmse = np.sqrt(np.mean((scale * mag - expected)[head] ** 2) / np.mean(expected[head] ** 2))
    m = interior(head) & fits["support"] & (fits["t2s"] > 0)
    err = dict(nrmse=nrmse, scale=scale, mask=m,
               t2s_ms=np.median(np.abs(fits["t2s"] - truth["t2s"][j])[m]) * 1e3,
               t2s_rel=np.median(np.abs(fits["t2s"] - truth["t2s"][j])[m] / truth["t2s"][j][m]))
    if fits["b0_hz"] is not None:
        err["b0_hz"] = np.median(np.abs(fits["b0_hz"] - truth["b0_hz"][j])[m])
    return err


def plot(images, fits, te, positions, truth, path, label):
    """Two rows per slice. Top: four echo images, then (with ``truth``) echo 0's expected image and the difference.
    Bottom: fitted T2* and B0 next to the phantom's maps and the fit errors. Returns the errors per slice."""
    n_echo = images.shape[1]
    show = [0, n_echo // 4, n_echo // 2, n_echo - 1]
    cols = 6
    fig, axes = plt.subplots(2 * len(images), cols, figsize=(3.4 * cols, 3.6 * 2 * len(images)), squeeze=False)
    errors = []

    def panel(ax, img, title, cmap="gray", lim=(None, None), bar=False):
        h = ax.imshow(img, cmap=cmap, origin="lower", vmin=lim[0], vmax=lim[1])
        ax.set_title(title, fontsize=10)
        if bar:
            fig.colorbar(h, ax=ax, fraction=0.046)

    for i, (img, fit, z) in enumerate(zip(images, fits, positions)):
        top, bottom = axes[2 * i], axes[2 * i + 1]
        vmax = np.abs(img[0]).max()
        for c, e in enumerate(show):
            panel(top[c], np.abs(img[e]), f"echo {e}, TE {te[e]*1e3:.1f} ms" + (f"\nz = {z*1e3:.1f} mm, {label}" if c == 0 else ""),
                  lim=(0, vmax))
        t2s_ms = np.where(fit["support"], fit["t2s"] * 1e3, np.nan)
        b0 = None if fit["b0_hz"] is None else np.where(fit["support"], fit["b0_hz"], np.nan)
        err = None
        if truth is not None:
            j = int(np.argmin(np.abs(truth["slice_positions"] - z)))
            err = compare_with_truth(img, fit, truth, j)
            errors.append(err)
            ev = truth["expected"][j].max()
            panel(top[4], truth["expected"][j], f"expected (ideal scan), TE {te[0]*1e3:.1f} ms", lim=(0, ev))
            panel(top[5], np.abs(err["scale"] * np.abs(img[0]) - truth["expected"][j]),
                  f"|echo 0 - expected|, NRMSE {err['nrmse']:.1%}", cmap="magma", lim=(0, 0.5 * ev))
            head = truth["rho"][j] > 0
            true_t2s = np.where(head, truth["t2s"][j] * 1e3, np.nan)
            true_b0 = np.where(head, truth["b0_hz"][j], np.nan)
        t2_title = "T2* fit (ms)" + (f"\nmedian error {err['t2s_ms']:.1f} ms ({err['t2s_rel']:.1%})" if err else "")
        panel(bottom[0], t2s_ms, t2_title, "viridis", (0, 120), bar=True)
        if b0 is not None:
            b0_title = "B0 fit (Hz)" + (f"\nmedian error {err['b0_hz']:.2f} Hz" if err else "")
            panel(bottom[3], b0, b0_title, "RdBu_r", (-60, 60), bar=True)
        if err is not None:
            panel(bottom[1], true_t2s, "T2* truth (ms)", "viridis", (0, 120), bar=True)
            panel(bottom[2], np.where(err["mask"], t2s_ms - true_t2s, np.nan), "T2* fit - truth (ms), interior",
                  "RdBu_r", (-10, 10), bar=True)
            panel(bottom[4], true_b0, "B0 truth (Hz)", "RdBu_r", (-60, 60), bar=True)
            if b0 is not None:
                panel(bottom[5], np.where(err["mask"], b0 - true_b0, np.nan), "B0 fit - truth (Hz), interior",
                      "RdBu_r", (-2, 2), bar=True)
    for a in axes.ravel():
        a.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=90, bbox_inches="tight")
    plt.close(fig)
    return errors


def plot_decay(images, fits, te, truth, positions, path):
    """Signal vs TE at three pixels spanning the fitted T2* range (first slice), with the fitted curve."""
    img, fit = images[0], fits[0]
    t2s = fit["t2s"]
    ok = fit["support"] & (t2s > 0)
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    for q in (10, 50, 90):
        y, x = np.unravel_index(np.argmin(np.where(ok, np.abs(t2s - np.percentile(t2s[ok], q)), np.inf)), t2s.shape)
        pts, = ax.semilogy(te * 1e3, np.abs(img[:, y, x]), "o", ms=3)
        label = f"({y},{x}): T2* fit {t2s[y, x]*1e3:.0f} ms"
        if truth is not None:
            j = int(np.argmin(np.abs(truth["slice_positions"] - positions[0])))
            label += f", truth {truth['t2s'][j][y, x]*1e3:.0f} ms"
        ax.semilogy(te * 1e3, np.abs(fit["s0"][y, x]) * np.exp(-te / t2s[y, x]), "-", color=pts.get_color(), label=label)
    ax.set_xlabel("TE (ms)")
    ax.set_ylabel("|signal|")
    ax.set_title("T2* decay (points: echo images, lines: fit)", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="run config, e.g. configs/cyclic_epi.toml")
    cfg = config.load(parser.parse_args(argv).config)
    method = cfg.reconstruction.get("coil_combine", "espirit")

    out = cfg.recon_raw.parent
    seq = pp.Sequence()
    seq.read(str(cfg.recon_seq))
    d = seq.definitions
    nx, ny = (int(v) for v in d["MatrixSize"])
    positions = np.atleast_1d(d["SlicePositions"])          # acquisition order
    te = float(d["TE"]) + np.arange(ny) * float(d["EchoSpacing"])
    print(f"{ny} echoes, TE {te[0]*1e3:.2f} to {te[-1]*1e3:.1f} ms; wrap timing {d.get('WrapTiming', '?')}, "
          f"rows after the wrap {float(d.get('WrapDelay', 0))*1e6:.0f} us late (treated as on time)")

    signal = np.load(cfg.recon_raw)["signal"]
    images = combine(sort_kspace(signal, len(positions), nx, ny), method)
    order = np.argsort(positions)                           # acquisition order -> increasing z
    images, positions = images[order], positions[order]
    fits = [fit_maps(img, te) for img in images]
    out.mkdir(parents=True, exist_ok=True)
    stack = lambda key: np.stack([f[key] for f in fits]) if fits[0][key] is not None else np.array([])  # noqa: E731
    np.savez(out / "recon.npz", echo_images=images, te=te, slice_positions=positions, s0=stack("s0"),
             t2s=stack("t2s"), b0_hz=stack("b0_hz"), support=stack("support"))

    truth = np.load(out / "truth.npz") if (out / "truth.npz").exists() else None
    if truth is not None and "expected" not in truth:
        truth = None
    errors = plot(images, fits, te, positions, truth, out / "recon.png", method)
    plot_decay(images, fits, te, truth, positions, out / "decay.png")
    for z, e in zip(positions, errors):
        b0 = f", B0 median error {e['b0_hz']:.2f} Hz" if "b0_hz" in e else ""
        print(f"z = {z*1e3:.1f} mm: echo 0 NRMSE vs expected image {e['nrmse']:.1%}; "
              f"T2* median error {e['t2s_ms']:.2f} ms ({e['t2s_rel']:.1%}){b0}")
    print(f"wrote {out.relative_to(config.REPO)}/recon.npz, recon.png and decay.png")


if __name__ == "__main__":
    main()
