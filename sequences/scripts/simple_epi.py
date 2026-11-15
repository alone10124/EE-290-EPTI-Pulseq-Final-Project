"""Simple EPI: multi-slice, single-shot 2D gradient-echo EPI, the baseline for building and debugging the workflow.

    python sequences/scripts/simple_epi.py configs/simple_epi.toml
    # -> sequences/simple_epi.seq (sequence.slices, for the scanner) and sequences/simple_epi_sim.seq (sim_slices)

Per slice (one TRID segment):

    [fat-sat + spoiler] -> slice-selective excitation -> prephasers (x, y) + slice rephaser
    -> ny readout lobes of alternating polarity, with a ky blip between lobes -> spoiler [-> TR padding]

Slices are acquired in interleaved order (even indices first), so neighbouring slices aren't excited back to back.
Simplifying assumptions: the ADC samples only on the readout flat top (no ramp sampling), the ky blips sit in their
own short blocks between lobes, and there is no acceleration, partial Fourier or simultaneous multi-slice.

The protocol defaults are the keyword defaults of ``build_sequence``; the config's ``[sequence]`` table overrides
them. The [DEFINITIONS] section of the .seq file records what the simulation and recon need: FOV, matrix, the index
and z position of each slice in acquisition order, TE and echo spacing.
"""

import argparse

import numpy as np
import pypulseq as pp

from mrtools import config
from mrtools.seq.common import (interleaved_slice_order, make_epi_readout, make_fatsat, slice_positions, trid,
                                write_seq)


def build_sequence(hw, fov=(220e-3, 220e-3), matrix=(128, 128), slab=120e-3, slice_thickness=3e-3, flip_angle=90,
                   rf_duration=3e-3, time_bw_product=4, dwell=4e-6, fat_sat=True, tr=0.0, slices="all", name="simple_epi"):
    """Build the sequence for hardware profile ``hw`` and return it.

    Units are SI (m, s), ``flip_angle`` in degrees, ``tr`` the volume TR (0 = shortest). ``slices`` is ``"all"`` or
    a list of slice indices to keep (0 = most inferior), e.g. for simulation; the timing, including the per-slice TR
    padding, stays as in the full sequence.
    """
    sys = hw.opts()
    nx, ny = matrix
    positions = slice_positions(slab, slice_thickness)
    order = interleaved_slice_order(len(positions))
    keep = [int(s) for s in order if isinstance(slices, str) or s in slices]
    if not keep:
        raise ValueError(f"no slices left: {slices} vs {len(positions)} slices")

    rf, gz, gz_reph = pp.make_sinc_pulse(
        flip_angle=np.deg2rad(flip_angle), duration=rf_duration, slice_thickness=slice_thickness, apodization=0.5,
        time_bw_product=time_bw_product, system=sys, return_gz=True, delay=sys.rf_dead_time, use="excitation")
    gx, adc = make_epi_readout(fov[0], nx, sys, dwell=dwell)
    blip = pp.make_trapezoid("y", system=sys, area=1 / fov[1])
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
        # Every slice is the same segment (TRID 1); only the RF frequency/phase change, which GE pge2 allows.
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

    seq.set_definition("Name", name)
    seq.set_definition("FOV", [fov[0], fov[1], slice_thickness])
    seq.set_definition("MatrixSize", [nx, ny])
    seq.set_definition("SliceIndices", keep)                                     # acquisition order
    seq.set_definition("SlicePositions", [float(positions[s]) for s in keep])   # m, acquisition order
    seq.set_definition("TE", float(te))
    seq.set_definition("EchoSpacing", float(t_ro + t_blip))
    return seq


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="run config, e.g. configs/simple_epi.toml")
    cfg = config.load(parser.parse_args(argv).config)

    protocol = dict(cfg.sequence)
    sim_slices = protocol.pop("sim_slices", [])
    write_seq(build_sequence(cfg.hardware, name=cfg.run, **protocol), cfg.seq_path, cfg.hardware)
    if sim_slices:
        protocol["slices"] = sim_slices
        write_seq(build_sequence(cfg.hardware, name=f"{cfg.run}_sim", **protocol), cfg.sim_seq_path, cfg.hardware)


if __name__ == "__main__":
    main()
