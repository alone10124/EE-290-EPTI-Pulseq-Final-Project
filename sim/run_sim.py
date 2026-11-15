"""Simulate a .seq file with KomaMRI on the brain phantom (the same for every sequence).

    python sim/run_sim.py configs/simple_epi.toml     # -> out/simple_epi_sim/raw.npz and truth.npz

Which .seq file: ``simulation.seq`` in the config, else the sim sequence (``sequences/<run>_sim.seq``), else the full
one. The phantom gets one brain slice at each z position the sequence excites. The geometry comes from the .seq file's
[DEFINITIONS], since only the .seq knows which slices were actually written: ``SlicePositions`` (m) and
``SliceIndices`` (optional). Nothing else about the sequence is needed, so this script is the same for every sequence.
Phantom options come from ``[simulation]``, the scanner (incl. number of coils) from ``[hardware]``; an optional
``simulation.ncoils`` simulates fewer receive channels than the hardware has (raw data grows with samples x coils),
and an optional ``simulation.threads`` caps KomaMRI's threads (default: all; each thread holds a full copy of the
signal, see ``mrtools.sim.koma.simulate``).

Simulating the whole sequence costs time proportional to its duration times the number of spins, so simulate a
short version (``sequence.sim_slices``).

Outputs (in the folder of ``simulation.raw``, default ``out/<seq name>/``):
* ``raw.npz``: ``signal`` [n_adc_samples, n_coils], every ADC sample in acquisition order (complex).
* ``truth.npz``: the ground truth as spins, not maps (``mrtools.sim.phantom.save_truth``): every spin's position,
  rho, T1, T2, T2*, off-resonance (Hz), slice id and fat flag, its RSS receive sensitivity ``coil_rss`` (KomaMRI's coil
  model, computed here so the recon doesn't need Julia) and ``slice_positions`` (increasing z). A recon script turns it
  into maps and an "expected" image on its own grid and with its own contrast: ``mrtools.sim.phantom.reference``.
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

    spins = ph.make_phantom(positions, indices, ss=opts.get("spin_subsampling", 1), fat=opts.get("fat", "water"),
                            b0_scale=opts.get("b0_scale", 1.0), B0=hw.B0, relaxation=opts.get("relaxation", True))
    scanner = hw.make_scanner()
    print(f"phantom: {len(spins)} spins, slices {indices.tolist()} at z = {np.round(positions * 1e3, 1).tolist()} mm; "
          f"{hw.ncoils} coils; {opts.get('threads', julia_threads())} of {julia_threads()} Julia threads")
    signal, secs = simulate(seq_path, ph.to_koma(spins), scanner, threads=opts.get("threads"))

    out = cfg.sim_raw.parent
    out.mkdir(parents=True, exist_ok=True)
    np.savez(cfg.sim_raw, signal=signal)
    ph.save_truth(out / "truth.npz", spins, positions, ph.coil_rss(scanner, spins))
    print(f"simulated {seq_path.name} in {secs:.1f} s: {signal.shape[0]} samples x {signal.shape[1]} coils "
          f"-> {config.display_path(cfg.sim_raw)}, truth.npz")


if __name__ == "__main__":
    main()
