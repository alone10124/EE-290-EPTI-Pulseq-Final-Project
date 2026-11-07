"""Simple EPI: multi-slice, single-shot 2D gradient-echo EPI, the baseline for building and debugging the workflow.

    python sequences/scripts/simple_epi.py                                  # -> sequences/simple_epi.seq (all slices)
    python sequences/scripts/simple_epi.py --name simple_epi_sim --slices 20   # only slice 20, for simulation

Per slice (one TRID segment):

    [fat-sat + spoiler] -> slice-selective excitation -> prephasers (x, y) + slice rephaser
    -> ny readout lobes of alternating polarity, with a ky blip between lobes -> spoiler [-> TR padding]

Slices are acquired in interleaved order (even indices first), so neighbouring slices aren't excited back to back.
Simplifying assumptions (see docs/LESSONS_LEARNED.md): the ADC samples only on the readout flat top, the ky blips sit in
their own short blocks between lobes, and there is no acceleration, partial Fourier or simultaneous multi-slice.

The [DEFINITIONS] section of the .seq file records what the simulation and recon need: FOV, matrix, the z position of
each slice in acquisition order, TE and echo spacing.
"""

import argparse
from pathlib import Path

import numpy as np
import pypulseq as pp

from mrtools.hardware import GE_MR750W
from mrtools.seq.checks import check_sequence
from mrtools.seq.common import interleaved_slice_order, make_epi_readout, make_fatsat, slice_positions, trid

# ---- Protocol ------------------------------------------------------------------------------------------------------
FOV = (220e-3, 220e-3)          # m, x (readout) and y (phase encode)
MATRIX = (128, 128)             # nx, ny
SLAB = 120e-3                   # m, z coverage; contiguous slices
SLICE_THICKNESS = 3e-3          # m  -> 40 slices
FLIP_ANGLE = 90                 # deg
RF_DURATION = 3e-3              # s
TIME_BW_PRODUCT = 4
DWELL = 4e-6                    # s, ADC dwell; multiple of GE's 2 us ADC raster
FAT_SAT = True
TR = 0.0                        # s, volume TR; 0 = as short as possible
HARDWARE = GE_MR750W

REPO = Path(__file__).resolve().parents[2]


def build_sequence(fov=FOV, matrix=MATRIX, slab=SLAB, slice_thickness=SLICE_THICKNESS, flip_angle=FLIP_ANGLE,
                   rf_duration=RF_DURATION, time_bw_product=TIME_BW_PRODUCT, dwell=DWELL, fat_sat=FAT_SAT, tr=TR,
                   hw=HARDWARE, slices=None):
    """Build the sequence and return it. ``slices`` keeps only those slice indices (for simulation); the timing,
    including the per-slice TR padding, stays as in the full sequence."""
    sys_fast = hw.opts()                    # EPI readout and blips
    sys = hw.opts(slow=True)                # everything else
    nx, ny = matrix
    positions = slice_positions(slab, slice_thickness)
    order = interleaved_slice_order(len(positions))
    keep = [int(s) for s in order if slices is None or s in slices]

    rf, gz, gz_reph = pp.make_sinc_pulse(
        flip_angle=np.deg2rad(flip_angle), duration=rf_duration, slice_thickness=slice_thickness, apodization=0.5,
        time_bw_product=time_bw_product, system=sys, return_gz=True, delay=sys.rf_dead_time, use="excitation")
    gx, adc = make_epi_readout(fov[0], nx, sys_fast, dwell=dwell)
    blip = pp.make_trapezoid("y", system=sys_fast, area=1 / fov[1])
    gx_pre = pp.make_trapezoid("x", system=sys, area=-gx.area / 2)
    gy_pre = pp.make_trapezoid("y", system=sys, area=-ny / 2 / fov[1])
    spoiler = pp.make_trapezoid("z", system=sys, area=4 / slice_thickness)
    fatsat = make_fatsat(sys) if fat_sat else None

    # Timing
    t_ro, t_blip = pp.calc_duration(gx), pp.calc_duration(blip)
    t_exc, t_pre = pp.calc_duration(rf, gz), pp.calc_duration(gx_pre, gy_pre, gz_reph)
    t_sat = pp.calc_duration(fatsat[0]) + pp.calc_duration(fatsat[1]) if fatsat else 0.0
    slice_min = t_sat + t_exc + t_pre + ny * t_ro + (ny - 1) * t_blip + pp.calc_duration(spoiler)
    pad = max(0.0, tr / len(positions) - slice_min)
    pad = np.ceil(pad / hw.block_duration_raster - 1e-9) * hw.block_duration_raster
    rf_center = rf.delay + pp.calc_rf_center(rf)[0]
    te = (t_exc - rf_center) + t_pre + (ny // 2) * (t_ro + t_blip) + t_ro / 2   # RF center -> ky = 0, mid-readout

    seq = pp.Sequence(system=sys)
    for s in keep:
        # Shift the excitation frequency to the slice position; the phase offset keeps the phase at the RF center
        # at zero, so every slice starts with the same phase.
        rf.freq_offset = gz.amplitude * positions[s]
        rf.phase_offset = -2 * np.pi * rf.freq_offset * pp.calc_rf_center(rf)[0]
        if fatsat:
            seq.add_block(fatsat[0], trid(1))
            seq.add_block(fatsat[1])
            seq.add_block(rf, gz)
        else:
            seq.add_block(rf, gz, trid(1))
        seq.add_block(gx_pre, gy_pre, gz_reph)
        for line in range(ny):
            seq.add_block(pp.scale_grad(gx, (-1) ** line), adc)
            if line < ny - 1:
                seq.add_block(blip)
        seq.add_block(spoiler)
        if pad > 0:
            seq.add_block(pp.make_delay(pad))

    seq.set_definition("Name", "simple_epi")
    seq.set_definition("FOV", [fov[0], fov[1], slice_thickness])
    seq.set_definition("MatrixSize", [nx, ny])
    seq.set_definition("SlicePositions", [float(positions[s]) for s in keep])   # m, in acquisition order
    seq.set_definition("TE", float(te))
    seq.set_definition("EchoSpacing", float(t_ro + t_blip))
    return seq


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--name", default="simple_epi", help="output file: sequences/<name>.seq")
    parser.add_argument("--slices", type=int, nargs="+", help="keep only these slice indices (0 = most inferior)")
    args = parser.parse_args(argv)

    seq = build_sequence(slices=args.slices)
    print(check_sequence(seq, HARDWARE, f"{args.name}.seq on {HARDWARE.name}"))
    d = seq.definitions
    n_all = len(slice_positions(SLAB, SLICE_THICKNESS))
    print(f"{len(d['SlicePositions'])} of {n_all} slices; TE {d['TE']*1e3:.1f} ms, echo spacing "
          f"{d['EchoSpacing']*1e6:.0f} us, duration {seq.duration()[0]:.2f} s")
    out = REPO / "sequences" / f"{args.name}.seq"
    seq.write(str(out))
    print(f"wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
