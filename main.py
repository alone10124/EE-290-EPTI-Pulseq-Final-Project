"""Run the whole workflow for one config: PyPulseq -> KomaMRI simulation -> reconstructed image.

    python main.py configs/simple_epi.toml             # or just: python main.py simple_epi
    python main.py configs/simple_epi.toml --skip-sim  # rebuild the sequence and redo the recon on existing data

Everything else (protocol, slices to simulate, phantom options, coil combination, hardware) is in the config; see
mrtools/config.py. The config's [general] section names the scripts; this file runs the same commands you would type:

    python sequences/scripts/<sequence>.py <config>   # -> sequences/<run>.seq (+ <run>_sim.seq if sim_slices is set)
    python sim/run_sim.py <config>                    # -> out/<run>_sim/raw.npz, truth.npz
    python recon/<recon>.py <config>                  # -> out/<run>_sim/recon.npz, recon.png
"""

import argparse
import subprocess
import sys
from pathlib import Path

from mrtools import config


def run(*args):
    print(f"\n$ python {' '.join(args)}", flush=True)
    subprocess.run([sys.executable, *args], cwd=config.REPO, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="config file, or a name in configs/ (e.g. simple_epi)")
    parser.add_argument("--skip-sim", action="store_true", help="reuse the existing raw data instead of simulating")
    args = parser.parse_args(argv)

    path = Path(args.config)
    if not path.suffix:
        path = config.REPO / "configs" / f"{args.config}.toml"
    path = str(path.resolve())
    cfg = config.load(path)

    run(f"sequences/scripts/{cfg.general['sequence']}.py", path)
    if not args.skip_sim:
        run("sim/run_sim.py", path)
    run(f"recon/{cfg.general['recon']}.py", path)


if __name__ == "__main__":
    main()
