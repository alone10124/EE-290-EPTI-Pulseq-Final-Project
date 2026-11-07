"""Simulate any .seq file with KomaMRI on the brain phantom (the same for every sequence).

    python sim/run_sim.py sequences/simple_epi_sim.seq            # -> out/simple_epi_sim/raw.npz and truth.npz

The phantom gets one brain slice at each z position the sequence excites. They are read from the .seq file's
[DEFINITIONS] (``SlicePositions`` in m, plus ``FOV`` and ``MatrixSize`` for the ground-truth maps); command-line
options override them or fill them in for .seq files that don't define them.

Simulating the whole sequence costs time proportional to its duration times the number of spins, so simulate a
short version: e.g. write the sequence with only one or two slices (``--slices`` in the sequence scripts).

Outputs:
* ``raw.npz``: ``signal`` [n_adc_samples, n_coils], every ADC sample in acquisition order (complex).
* ``truth.npz``: ``rho``, ``t2s`` (s), ``b0_hz`` maps [slice, y, x] on the recon grid, and ``slice_positions``
  (increasing z). See ``mrtools.sim.phantom``.
"""

import argparse
from pathlib import Path

import numpy as np
import pypulseq as pp

from mrtools.hardware import GE_MR750W
from mrtools.sim import phantom as ph
from mrtools.sim.koma import simulate

REPO = Path(__file__).resolve().parents[1]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seq", type=Path, help=".seq file to simulate")
    parser.add_argument("--slices-mm", type=float, nargs="+", help="phantom slice z positions (mm)")
    parser.add_argument("--fov-mm", type=float, nargs=2, help="ground-truth map FOV x y (mm)")
    parser.add_argument("--matrix", type=int, nargs=2, help="ground-truth map matrix nx ny")
    parser.add_argument("--ncoils", type=int, default=8)
    parser.add_argument("--spin-subsampling", type=int, default=1,
                        help="phantom spin spacing = 0.5 mm * this; keep it well below the voxel size")
    parser.add_argument("--fat", choices=["water", "3T"], default="water",
                        help='"water": no fat chemical shift; "3T": fat at -3.4 ppm (so fat-sat has something to do)')
    parser.add_argument("--b0-scale", type=float, default=1.0, help="scale of the synthetic B0 map (0 = off)")
    parser.add_argument("--gpu", action="store_true", help="simulate on an NVIDIA GPU (CUDA)")
    args = parser.parse_args(argv)

    seq = pp.Sequence()
    seq.read(str(args.seq))
    d = seq.definitions
    positions = (np.asarray(args.slices_mm) * 1e-3 if args.slices_mm is not None
                 else np.atleast_1d(d["SlicePositions"]))
    positions = np.unique(positions)                       # increasing z, each slice once
    fov = np.asarray(args.fov_mm) * 1e-3 if args.fov_mm else np.atleast_1d(d["FOV"])[:2]
    matrix = args.matrix or [int(v) for v in np.atleast_1d(d["MatrixSize"])[:2]]

    hw = GE_MR750W
    spins = ph.make_phantom(positions, ss=args.spin_subsampling, fat=args.fat, b0_scale=args.b0_scale, B0=hw.B0)
    scanner = ph.make_scanner(hw, ncoils=args.ncoils)
    print(f"phantom: {len(spins)} spins at z = {np.round(positions * 1e3, 1).tolist()} mm; {args.ncoils} coils")
    signal, secs = simulate(args.seq, ph.to_koma(spins), scanner, gpu=args.gpu)

    out = REPO / "out" / args.seq.stem
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / "raw.npz", signal=signal)
    truth = ph.ground_truth_maps(spins, fov, matrix, positions)
    np.savez(out / "truth.npz", slice_positions=positions, **truth)
    print(f"simulated in {secs:.0f} s: {signal.shape[0]} samples x {signal.shape[1]} coils "
          f"-> {out.relative_to(REPO)}/raw.npz, truth.npz")


if __name__ == "__main__":
    main()
