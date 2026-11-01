"""Quick look at simulated calibration data: low-res multi-echo images by plain FFT, next to the phantom ground truth.

    python sim/preview_calib.py config/sim_small.yaml      # after sim/run_sim.py; writes out/sim/<config>_calib_preview.png

This is a sanity check of the simulation, NOT the EPTI reconstruction:
* the calibration scan reads one ky line per TR (17 echoes of it), so each echo's lines form a 49-line k-space;
* the outer lines are 3x undersampled (no GRAPPA here) and the ramp-sampled readout isn't regridded, so expect some
  aliasing and blur;
* only positive-polarity echoes are used, so there is no odd/even ghost correction to do.
The phase change between two same-polarity echoes divided by 2*pi*dTE estimates the B0 map.
"""

import argparse
import sys
from pathlib import Path

import ismrmrd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pypulseq as pp  # noqa: E402

from epti.config import load_config  # noqa: E402

REPO = Path(__file__).resolve().parents[1]


def load_acquisitions(path):
    ds = ismrmrd.Dataset(str(path), "dataset", create_if_needed=False)
    return np.stack([ds.read_acquisition(i).data for i in range(ds.number_of_acquisitions())])  # [acq, coil, sample]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path)
    args = parser.parse_args(argv)
    cfg = load_config(args.config)
    stem = args.config.stem
    sim = REPO / "out" / "sim"

    seq = pp.Sequence()
    seq.read(str(REPO / "sequences" / f"{stem}_calib.seq"))
    pe_lines = int(seq.definitions["pe_lines"])
    n_echo = int(seq.definitions["ETL"])                   # 3 navigators + echoes
    positions = np.atleast_1d(seq.definitions["SlicePositions"])
    n_slc = len(positions)
    te0 = float(seq.definitions["TE"])                      # TE of the center of the echo train
    esp = sum(float(v) for v in np.atleast_1d(seq.definitions["TrapezoidGriddingParameters"])[:3])

    data = load_acquisitions(sim / f"{stem}_calib.h5")   # acquisition order: PE line, slice, echo
    n_acq, n_coil, n_ro = data.shape
    assert n_acq == pe_lines * n_slc * n_echo, (n_acq, pe_lines, n_slc, n_echo)
    data = data.reshape(pe_lines, n_slc, n_echo, n_coil, n_ro)

    from epti.seq.patterns import cal_pe_index
    pe = cal_pe_index(cfg.calib.acs_npe, cfg.calib.nseg, cfg.calib.pat_factor).astype(int) + cfg.calib.nseg // 2
    # positive-polarity imaging echoes: echo i (1-based, after the 3 navigators) has polarity (-1)^((i+1)%2)
    pos_echoes = [cfg.calib.nav_etl + i - 1 for i in range(1, cfg.calib.etl + 1) if (i + 1) % 2 == 0]
    etl_pre = int(np.floor(cfg.calib.etl / 2 - 1 + 0.5))
    te = {e: te0 + (e - cfg.calib.nav_etl - etl_pre) * esp for e in pos_echoes}

    def image(slc, echo):
        k = np.zeros((cfg.calib.nseg, n_coil, n_ro), complex)
        k[pe] = data[:, slc, echo]
        img = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(k, axes=(0, 2)), axes=(0, 2)), axes=(0, 2))
        crop = slice(n_ro // 4, n_ro - n_ro // 4)                # remove the 2x readout oversampling
        return img[:, :, crop]                                   # [y, coil, x]

    truth = np.load(sim / f"{stem}_truth.npz")
    e1, e2 = pos_echoes[0], pos_echoes[1]
    fig, ax = plt.subplots(n_slc, 5, figsize=(16, 3.4 * n_slc), squeeze=False)
    for s in range(n_slc):
        i1, i2 = image(s, e1), image(s, e2)
        sos1, sos3 = np.sqrt((np.abs(i1) ** 2).sum(1)), np.sqrt((np.abs(image(s, pos_echoes[-1])) ** 2).sum(1))
        # coil-combined phase difference (sum over coils of i2 * conj(i1)), as a B0 estimate in Hz
        dphi = np.angle((i2 * np.conj(i1)).sum(1))
        b0 = dphi / (2 * np.pi * (te[e2] - te[e1]))
        mask = sos1 > 0.15 * sos1.max()
        panels = [
            (truth["rho"][s], "truth: proton density", "gray", None),
            (sos1, f"sim echo TE={te[e1]*1e3:.1f} ms", "gray", None),
            (sos3, f"sim echo TE={te[pos_echoes[-1]]*1e3:.1f} ms", "gray", None),
            (truth["b0_hz"][s] * (truth["rho"][s] > 0), "truth: B0 (Hz)", "RdBu_r", (-60, 60)),
            (np.where(mask, b0, np.nan), "sim: B0 from echo phase (Hz)", "RdBu_r", (-60, 60)),
        ]
        for a, (img, title, cmap, lim) in zip(ax[s], panels):
            im = a.imshow(img, cmap=cmap, origin="lower", aspect="auto",
                          vmin=None if lim is None else lim[0], vmax=None if lim is None else lim[1])
            a.set_title(f"z={positions[s]*1e3:+.1f} mm  {title}", fontsize=9)
            a.axis("off")
            if lim is not None:
                fig.colorbar(im, ax=a, fraction=0.046)
    fig.suptitle(f"{stem}: KomaMRI calibration scan, plain-FFT preview ({n_coil} coils, sum of squares)")
    fig.tight_layout()
    out = sim / f"{stem}_calib_preview.png"
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
