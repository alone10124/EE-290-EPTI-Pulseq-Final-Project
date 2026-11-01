"""Simulate the EPTI calibration + imaging scans with KomaMRI and write ISMRMRD raw data.

    python sequences/scripts/make_epti.py config/sim_small.yaml       # 1. make the .seq files
    python sim/run_sim.py config/sim_small.yaml                       # 2. simulate them

Writes, under out/sim/:
    <config>_calib.h5, <config>_imaging.h5   ISMRMRD raw data (same input format the recon will take from the scanner)
    <config>_truth.npz                        ground-truth proton density, T2* and B0 maps on the recon grid

Both scans use the same phantom and receive coils, as on a real scanner. Use --gpu on a machine with an NVIDIA GPU
(CUDA). On the Radeon 780M laptop, run on the CPU threads (the default).
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import phantom as ph  # noqa: E402  (sets the juliacall environment variables before komamripy is imported)

from epti.config import load_config  # noqa: E402
from epti.seq.imaging import excited_slice_positions  # noqa: E402

REPO = Path(__file__).resolve().parents[1]


def simulate(obj, seq_path, scanner, gpu):
    import komamripy as km
    from juliacall import Main as jl

    seq = km.read_seq(str(seq_path))
    params = jl.seval("(gpu) -> Dict{String,Any}(\"return_type\" => \"raw\", \"gpu\" => gpu)")(gpu)
    t0 = time.time()
    raw = km.simulate(obj, seq, scanner, sim_params=params)
    return raw, time.time() - t0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path)
    parser.add_argument("--seqdir", type=Path, default=REPO / "sequences", help="where make_epti.py wrote the .seq files")
    parser.add_argument("--outdir", type=Path, default=REPO / "out" / "sim")
    parser.add_argument("--ss", type=int, default=2, help="phantom subsampling of the 0.5 mm grid (2 -> 1 mm spins)")
    parser.add_argument("--ncoils", type=int, default=8)
    parser.add_argument("--fat", choices=["water", "3T"], default="water")
    parser.add_argument("--b0-scale", type=float, default=1.0, help="scale of the synthetic B0 map (0 = off)")
    parser.add_argument("--gpu", action="store_true", help="simulate on an NVIDIA GPU (CUDA)")
    parser.add_argument("--only", choices=["calib", "imaging"], help="simulate just one of the two scans")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    stem = args.config.stem
    if cfg.imaging.slice_groups is None:
        print("warning: the config plays every slice group; the simulation will be slow. See config/sim_small.yaml.")
    slices = excited_slice_positions(cfg.imaging, cfg.imaging.slice_groups)
    args.outdir.mkdir(parents=True, exist_ok=True)

    import komamripy as km
    if args.gpu:
        km.load_cuda()

    spins = ph.make_phantom(slices, ss=args.ss, fat=args.fat, b0_scale=args.b0_scale, B0=cfg.hardware.B0)
    obj = ph.to_koma(spins)
    scanner = ph.make_scanner(cfg.hardware, ncoils=args.ncoils)
    print(f"phantom: {len(spins)} spins in {len(slices)} slices at z = {np.round(np.asarray(slices) * 1e3, 1)} mm; "
          f"{args.ncoils} receive coils")

    truth = ph.ground_truth_maps(spins, cfg.imaging.fov, (cfg.imaging.nx, cfg.imaging.ny), slices)
    np.savez(args.outdir / f"{stem}_truth.npz", slice_positions=slices, **truth)

    for scan in ("calib", "imaging"):
        if args.only and scan != args.only:
            continue
        seq_path = args.seqdir / f"{stem}_{scan}.seq"
        if not seq_path.exists():
            raise SystemExit(f"{seq_path} not found: run sequences/scripts/make_epti.py {args.config} first")
        raw, secs = simulate(obj, seq_path, scanner, args.gpu)
        out = args.outdir / f"{stem}_{scan}.h5"
        out.unlink(missing_ok=True)
        km.save(km.ISMRMRDFile(str(out)), raw)
        print(f"{scan}: simulated in {secs:.0f} s -> {out} ({len(raw.profiles)} acquisitions)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
