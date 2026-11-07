"""Run the whole workflow for one sequence: PyPulseq -> KomaMRI simulation -> reconstructed image.

    python main.py simple_epi                      # uses sequences/scripts/simple_epi.py and recon/simple_epi.py
    python main.py simple_epi --sim-slices 10 20   # simulate two slices instead of one
    python main.py simple_epi --skip-sim           # rebuild the sequence and redo the recon on existing data

Each sequence <name> has a script that writes it (sequences/scripts/<name>.py, accepting --name and --slices) and a
recon script (recon/<name>.py, taking the .seq file). The simulator (sim/run_sim.py) is shared. This file runs the
same commands you would type by hand:

    python sequences/scripts/<name>.py                                  # full sequence, for the scanner
    python sequences/scripts/<name>.py --name <name>_sim --slices 20    # short version to simulate
    python sim/run_sim.py sequences/<name>_sim.seq
    python recon/<name>.py sequences/<name>_sim.seq                     # -> out/<name>_sim/recon.png
"""

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent


def run(*args):
    print(f"\n$ python {' '.join(args)}", flush=True)
    subprocess.run([sys.executable, *args], cwd=REPO, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sequence", help="sequence name, e.g. simple_epi")
    parser.add_argument("--sim-slices", nargs="+", default=["20"], help="slice indices to simulate (default: 20)")
    parser.add_argument("--skip-sim", action="store_true", help="reuse out/<name>_sim/raw.npz instead of simulating")
    args = parser.parse_args(argv)

    name, sim_name = args.sequence, f"{args.sequence}_sim"
    seq_script = f"sequences/scripts/{name}.py"
    run(seq_script)
    run(seq_script, "--name", sim_name, "--slices", *args.sim_slices)
    if not args.skip_sim:
        run("sim/run_sim.py", f"sequences/{sim_name}.seq")
    run(f"recon/{name}.py", f"sequences/{sim_name}.seq")


if __name__ == "__main__":
    main()
