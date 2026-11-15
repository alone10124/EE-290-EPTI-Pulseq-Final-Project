"""Cyclic EPI: multi-shot, multi-slice 2D gradient-echo EPI in which every echo index gives a full k-space.

    python sequences/scripts/cyclic_epi.py configs/cyclic_epi.toml
    # -> sequences/cyclic_epi.seq (sequence.slices, for the scanner) and sequences/cyclic_epi_sim.seq (sim_slices)

For an ny-row k-space there are ny shots, each reading ny lines. Shot ``s`` starts at ky row ``s`` and moves down one
row per echo; after the last row a long "wrap" blip (-(ny - 1) rows) takes it back to row 0:

    row(s, j) = (s + j) mod ny          (echo j = 0 ... ny - 1)

So echo ``j`` taken from all shots fills every row once: one full k-space per echo, i.e. ny images along the T2* decay.
Within an echo, every row is acquired at the same time after its excitation, so the images have no phase-encode
distortion or T2* blur, and since the readout polarity depends only on the echo index there is no even/odd ghost.

Per shot and slice (one TRID segment):

    [fat-sat + spoiler] -> slice-selective excitation -> prephasers (x, y to row s) + slice rephaser
    -> ny readout lobes (polarity (-1)^j), a +1-row blip or the wrap blip between lobes -> spoiler -> padding

The wrap blip is longer than a normal blip (with the MR750w limits: 532 vs 48 us), so the echoes after the wrap come
late. ``wrap_timing`` chooses what to do about it:

* ``"ignore"``: normal blips; in echo j, rows 0 ... j-1 (the shots that already wrapped) are ``WrapDelay`` late.
  The recon treats them as on time.
* ``"pad"``: every blip lasts as long as the wrap blip, so every row of echo j is at the same TE (``WrapDelay`` = 0),
  at the cost of a longer echo spacing.

Shots are the outer loop and slices the inner one (interleaved order), so each slice sees one excitation per volume TR.
``dummy_shots`` extra shots without ADC come first, to reach the T1 steady state. TRID = s + 1 for shot s and ny + 1 for
the dummy shots: the wrap sits at a different echo in every shot, so every shot is its own segment. A slice subset
(``slices``, e.g. for simulation) keeps the shot period of the full sequence (a delay stands in for the omitted slices),
so the simulated T1 history matches the scanner's; ky-dependent intensity changes would otherwise show up as artifacts.

Simplifying assumptions as in simple_epi.py: ADC only on the readout flat top, blips in their own blocks between lobes,
no acceleration.

Adapted from the team's draft ``make_cyclic_epi_ge.py`` (2026-11; Stage-2 copy in ``cyclic_epi_draft/``). Changes: built
from ``mrtools.seq.common`` (readout, fat-sat, TRID, slice order) and the config's hardware profile, fat-sat added,
``wrap_timing`` option, z-only spoiler as in simple_epi.py, slice subsets keep the shot period, and the timing the recon
needs goes into the .seq [DEFINITIONS] instead of a separate meta file.
"""

import argparse

import numpy as np
import pypulseq as pp

from mrtools import config
from mrtools.seq.common import (interleaved_slice_order, make_epi_readout, make_fatsat, slice_positions, trid,
                                write_seq)


def build_sequence(hw, fov=(220e-3, 220e-3), matrix=(128, 128), slab=120e-3, slice_thickness=3e-3, flip_angle=90,
                   rf_duration=3e-3, time_bw_product=4, dwell=4e-6, fat_sat=True, tr=0.0, dummy_shots=1,
                   wrap_timing="ignore", slices="all", name="cyclic_epi"):
    """Build the sequence for hardware profile ``hw`` and return it.

    Units are SI (m, s), ``flip_angle`` in degrees, ``tr`` the volume TR (0 = shortest). ``matrix`` is (nx, ny); there
    are ny shots of ny echoes. ``slices`` is ``"all"`` or a list of slice indices to keep (0 = most inferior); the shot
    period stays that of the full sequence. ``wrap_timing`` is ``"ignore"`` or ``"pad"`` (see the module doc).
    """
    if wrap_timing not in ("ignore", "pad"):
        raise ValueError(f"wrap_timing must be 'ignore' or 'pad', not {wrap_timing!r}")
    sys = hw.opts()
    nx, ny = matrix
    dk = 1 / fov[1]
    positions = slice_positions(slab, slice_thickness)
    order = interleaved_slice_order(len(positions))
    keep = [int(s) for s in order if isinstance(slices, str) or s in slices]
    if not keep:
        raise ValueError(f"no slices left: {slices} vs {len(positions)} slices")

    rf, gz, gz_reph = pp.make_sinc_pulse(
        flip_angle=np.deg2rad(flip_angle), duration=rf_duration, slice_thickness=slice_thickness, apodization=0.5,
        time_bw_product=time_bw_product, system=sys, return_gz=True, delay=sys.rf_dead_time, use="excitation")
    gx, adc = make_epi_readout(fov[0], nx, sys, dwell=dwell)
    gx_pre = pp.make_trapezoid("x", system=sys, area=-gx.area / 2)
    wrap = pp.make_trapezoid("y", system=sys, area=-(ny - 1) * dk)
    t_wrap = pp.calc_duration(wrap)
    if wrap_timing == "pad":    # same area, stretched to the wrap blip's duration
        blip = pp.make_trapezoid("y", system=sys, area=dk, duration=t_wrap)
    else:
        blip = pp.make_trapezoid("y", system=sys, area=dk)
    spoiler = pp.make_trapezoid("z", system=sys, area=4 / slice_thickness)
    fatsat = make_fatsat(sys) if fat_sat else None

    # The y prephaser goes to row s: area (s - ny/2) dk. All prephaser blocks get the duration of the largest one, so
    # every shot has the same timing up to the wrap.
    t_pre = max(pp.calc_duration(g) for g in
                (gx_pre, gz_reph, pp.make_trapezoid("y", system=sys, area=ny / 2 * dk)))
    t_pre = np.ceil(t_pre / hw.grad_raster_time - 1e-9) * hw.grad_raster_time
    gx_pre = pp.make_trapezoid("x", system=sys, area=-gx.area / 2, duration=t_pre)
    gz_reph = pp.make_trapezoid("z", system=sys, area=gz_reph.area, duration=t_pre)

    def gy_pre(s):
        area = (s - ny / 2) * dk
        return [pp.make_trapezoid("y", system=sys, area=area, duration=t_pre)] if area else []

    # Timing. Shot 0 never wraps; every other shot wraps once, so it is the longest.
    t_ro, t_blip = pp.calc_duration(gx, adc), pp.calc_duration(blip)
    t_exc = pp.calc_duration(rf, gz)
    t_sat = pp.calc_duration(fatsat[0]) + pp.calc_duration(fatsat[1]) if fatsat else 0.0
    t_spoil = pp.calc_duration(spoiler)

    def segment(s):                                                  # one shot, one slice, before padding
        wraps = 1 if s % ny else 0
        return t_sat + t_exc + t_pre + ny * t_ro + (ny - 1 - wraps) * t_blip + wraps * t_wrap + t_spoil

    seg_min = segment(1)
    slot = max(seg_min, tr / len(positions))                         # time per slice and shot
    slot = np.ceil(slot / hw.block_duration_raster - 1e-9) * hw.block_duration_raster
    rf_center = rf.delay + pp.calc_rf_center(rf)[0]
    te = (t_exc - rf_center) + t_pre + t_ro / 2                     # RF center -> center of echo 0 (kx = 0)
    esp = t_ro + t_blip

    seq = pp.Sequence(system=sys)
    shots = [(0, False)] * dummy_shots + [(s, True) for s in range(ny)]
    for s, acquire in shots:
        segment_id = s + 1 if acquire else ny + 1
        for z in keep:
            # Excitation at the slice position; the phase offset keeps the phase at the RF center at zero.
            rf.freq_offset = gz.amplitude * positions[z]
            rf.phase_offset = -2 * np.pi * rf.freq_offset * pp.calc_rf_center(rf)[0]
            if fatsat:
                seq.add_block(fatsat[0], trid(segment_id))
                seq.add_block(fatsat[1])
                seq.add_block(rf, gz)
            else:
                seq.add_block(rf, gz, trid(segment_id))
            seq.add_block(gx_pre, gz_reph, *gy_pre(s))
            for j in range(ny):
                lobe = pp.scale_grad(gx, (-1) ** j)
                if acquire:
                    seq.add_block(lobe, adc)
                else:
                    seq.add_block(lobe, pp.make_delay(t_ro))
                if j < ny - 1:
                    seq.add_block(wrap if (s + j) % ny == ny - 1 else blip)
            seq.add_block(spoiler)
            pad = slot - segment(s)
            if pad > 1e-9:
                seq.add_block(pp.make_delay(round(pad / hw.block_duration_raster) * hw.block_duration_raster))
        if len(keep) < len(positions):     # the omitted slices' time, so the shot period matches the full sequence
            seq.add_block(pp.make_delay(round((len(positions) - len(keep)) * slot / hw.block_duration_raster)
                                        * hw.block_duration_raster))

    seq.set_definition("Name", name)
    seq.set_definition("FOV", [fov[0], fov[1], slice_thickness])
    seq.set_definition("MatrixSize", [nx, ny])
    seq.set_definition("SliceIndices", keep)                                     # acquisition order
    seq.set_definition("SlicePositions", [float(positions[z]) for z in keep])   # m, acquisition order
    seq.set_definition("TE", float(te))                                          # echo 0; echo j: TE + j * EchoSpacing
    seq.set_definition("EchoSpacing", float(esp))
    seq.set_definition("WrapDelay", round(float(t_wrap - t_blip), 9))                    # how late rows 0..j-1 of echo j are
    seq.set_definition("WrapTiming", wrap_timing)
    seq.set_definition("DummyShots", int(dummy_shots))
    seq.set_definition("ShotPeriod", float(slot * len(positions)))               # = volume TR
    return seq


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="run config, e.g. configs/cyclic_epi.toml")
    cfg = config.load(parser.parse_args(argv).config)

    protocol = dict(cfg.sequence)
    sim_slices = protocol.pop("sim_slices", [])
    write_seq(build_sequence(cfg.hardware, name=cfg.run, **protocol), cfg.seq_path, cfg.hardware)
    if sim_slices:
        protocol["slices"] = sim_slices
        write_seq(build_sequence(cfg.hardware, name=f"{cfg.run}_sim", **protocol), cfg.sim_seq_path, cfg.hardware)


if __name__ == "__main__":
    main()
