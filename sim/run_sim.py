"""Simulate a .seq file with KomaMRI on the brain phantom (the same for every sequence).

    python sim/run_sim.py configs/simple_epi.toml     # -> out/simple_epi_sim/raw.npz and truth.npz

Which .seq file: ``simulation.seq`` in the config, else the sim sequence (``sequences/<run>_sim.seq``), else the full
one. The phantom gets one brain slice at each z position the sequence excites. The geometry comes from the .seq file's
[DEFINITIONS], since only the .seq knows which slices were actually written: ``SlicePositions`` (m), ``SliceIndices``
(optional), ``FOV`` and ``MatrixSize`` (for the ground-truth maps), ``TE`` (optional, for the expected image).
Phantom options come from ``[simulation]``, the scanner (incl. number of coils) from ``[hardware]``; an optional
``simulation.ncoils`` simulates fewer receive channels than the hardware has (raw data grows with samples x coils).

Simulating the whole sequence costs time proportional to its duration times the number of spins, so simulate a
short version (``sequence.sim_slices``).

Outputs (in the folder of ``simulation.raw``, default ``out/<seq name>/``):
* ``raw.npz``: ``signal`` [n_adc_samples, n_coils], every ADC sample in acquisition order (complex).
* ``truth.npz`` (slices in increasing z, maps [slice, y, x] on the recon grid):
  ``rho``, ``t2s`` (s), ``b0_hz`` voxel maps, ``slice_positions``, and, if TE is defined, ``expected``: the image an
  ideal, fully sampled Cartesian scan at that TE would give after RSS-type coil combination,
  i.e. ``rho * exp(-TE/T2*) * coil_rss`` band-limited to the k-space grid (see ``mrtools.sim.phantom``). Comparing a
  recon with it leaves only EPI physics (B0 distortion, T2* decay during the readout) and recon errors.
"""

import argparse
import dataclasses

import numpy as np
import pypulseq as pp

from mrtools import config
from mrtools.sim import phantom as ph
from mrtools.sim.koma import julia_threads, simulate


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="run config, e.g. configs/simple_epi.toml")
    cfg = config.load(parser.parse_args(argv).config)
    opts, hw = cfg.simulation, cfg.hardware
    if "ncoils" in opts:
        hw = dataclasses.replace(hw, ncoils=int(opts["ncoils"]))

    seq_path = cfg.sim_input
    seq = pp.Sequence()
    seq.read(str(seq_path))
    d = seq.definitions
    positions = np.atleast_1d(d["SlicePositions"]).astype(float)
    indices = np.atleast_1d(d.get("SliceIndices", np.arange(len(positions)))).astype(int)
    order = np.argsort(positions)                          # increasing z
    positions, indices = positions[order], indices[order]
    fov = np.atleast_1d(d["FOV"])[:2]
    matrix = [int(v) for v in np.atleast_1d(d["MatrixSize"])[:2]]

    spins = ph.make_phantom(positions, indices, ss=opts.get("spin_subsampling", 1), fat=opts.get("fat", "water"),
                            b0_scale=opts.get("b0_scale", 1.0), B0=hw.B0, relaxation=opts.get("relaxation", True))
    scanner = hw.make_scanner()
    print(f"phantom: {len(spins)} spins, slices {indices.tolist()} at z = {np.round(positions * 1e3, 1).tolist()} mm; "
          f"{hw.ncoils} coils; {julia_threads()} Julia threads")
    signal, secs = simulate(seq_path, ph.to_koma(spins), scanner)

    out = cfg.sim_raw.parent
    out.mkdir(parents=True, exist_ok=True)
    np.savez(cfg.sim_raw, signal=signal)
    truth = ph.ground_truth_maps(spins, fov, matrix, len(positions))
    if "TE" in d:
        w = spins.rho * np.exp(-float(d["TE"]) / spins.t2s) * ph.coil_rss(scanner, spins)
        truth["expected"] = np.stack([np.abs(ph.bandlimited_image(spins.x[m], spins.y[m], w[m], fov, matrix))
                                      for m in (spins.slice_id == s for s in range(len(positions)))])
    np.savez(out / "truth.npz", slice_positions=positions, **truth)
    print(f"simulated {seq_path.name} in {secs:.1f} s: {signal.shape[0]} samples x {signal.shape[1]} coils "
          f"-> {cfg.sim_raw.relative_to(config.REPO)}, truth.npz")


if __name__ == "__main__":
    main()
