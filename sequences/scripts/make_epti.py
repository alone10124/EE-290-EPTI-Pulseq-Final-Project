"""Generate the EPTI calibration + imaging .seq files from a config, and check them against the hardware profile.

    python sequences/scripts/make_epti.py config/epti_3mm_ge.yaml          # full protocol for the GE scanner
    python sequences/scripts/make_epti.py config/sim_small.yaml            # small subset for KomaMRI simulation

Writes sequences/<config name>_calib.seq and sequences/<config name>_imaging.seq. The exit code is 1 if any hardware
check fails; the files are still written so you can inspect them.
"""

import argparse
import sys
import time
from pathlib import Path

from epti.config import load_config
from epti.seq.build import build_all, write_all

REPO = Path(__file__).resolve().parents[2]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path, help="YAML config (see config/)")
    parser.add_argument("--outdir", type=Path, default=REPO / "sequences", help="output folder (default: sequences/)")
    parser.add_argument("--no-check", action="store_true", help="skip the hardware checks")
    parser.add_argument("--plot", action="store_true", help="plot the first slice-group TR of the imaging scan")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    t0 = time.time()
    built = build_all(cfg, check=not args.no_check)
    paths = write_all(built, args.outdir, args.config.stem)

    img, cal = built.imaging_info, built.calib_info
    print(f"built in {time.time() - t0:.1f} s with hardware profile '{cfg.hardware.name}'")
    print(f"  imaging: {paths.imaging}  ({built.imaging.duration()[0]:.2f} s, TE1 {img.te1*1e3:.2f} ms, "
          f"ESP {img.esp*1e6:.0f} us, {img.n_groups_played} slice group(s), {cfg.imaging.n_dynamics} dynamic(s))")
    print(f"  calib:   {paths.calib}  ({built.calib.duration()[0]:.2f} s, TE {cal.te*1e3:.2f} ms, "
          f"{len(cal.positions)} slice(s), {len(cal.pe_index)} PE lines)")
    for rep in built.reports:
        print(rep)
    if args.plot:
        built.imaging.plot(time_range=(0, img.dur_group + img.delay_te2))
    return 0 if all(r.ok for r in built.reports) else 1


if __name__ == "__main__":
    sys.exit(main())
