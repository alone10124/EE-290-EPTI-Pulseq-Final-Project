"""EPTI calibration scan: multi-echo gradient-echo EPI used for coil sensitivity, B0 and odd/even phase maps.

Port of ``legacy/sequences/seq_EPTI_calib_P1.m`` (MGH, Jian Wu, Zijing Dong, Fuyixue Wang; version 09/2025).
Line numbers in comments refer to that script.

Each TR excites one slice with a sinc pulse and reads a *single* ky line ``etl`` times with bipolar lobes (no blips),
so every echo of the line is a different TE. The ky line changes between repetitions (``cal_pe_index``: fully
sampled center, every ``pat_factor``-th line outside). Three navigator lobes before the phase-encoding step measure the
odd/even phase.
"""

from types import SimpleNamespace

import numpy as np
import pypulseq as pp

from epti.config import CalibProtocol
from epti.hardware import HardwareProfile
from epti.seq.common import BlockWriter, make_fatsat, make_readout, stretch, trapezoid
from epti.seq.patterns import cal_pe_index, interleaved_slice_order
from epti.util import ceil_to_raster, mround

TRID_SLICE = 1
TRID_PREAMBLE = 2


def slice_positions(p: CalibProtocol):
    """All slice positions (m) in acquisition order (L140-145)."""
    order = interleaved_slice_order(p.n_slices)
    return -((p.thickness + p.slice_gap) * (np.arange(p.n_slices) - (p.n_slices - 1) / 2))[order]


def build_calib(p: CalibProtocol, hw: HardwareProfile):
    """Build the calibration sequence. Returns ``(seq, info)``."""
    sys, sys2 = hw.opts(), hw.opts(slow=True)
    if abs(round(p.blip_dur / sys.grad_raster_time) * sys.grad_raster_time - p.blip_dur) > 1e-9:
        raise ValueError("blip_dur must be on the gradient raster (L70)")

    rf_fs, gz_fs = make_fatsat(p.fatsat, sys2, hw)                                   # L71-78
    rf, gz, gz_reph = pp.make_sinc_pulse(np.deg2rad(p.rf_angle_deg), system=sys2, duration=p.rf_dur,  # L80-81
                                         slice_thickness=p.thickness, apodization=p.rf_apodization,
                                         time_bw_product=p.rf_tbw, return_gz=True, use="excitation",
                                         delay=sys2.rf_dead_time)
    deltak = 1 / p.fov
    min_blip = np.ceil(2 * np.sqrt(deltak / sys.max_slew) / sys.grad_raster_time / 2) * sys.grad_raster_time * 2
    if min_blip > p.blip_dur + 1e-12:                                                 # L86-87
        raise ValueError("blip_dur too short")
    gx, adc = make_readout(p.fov, p.nx, p.esp, p.blip_dur, p.ro_os, sys, hw)        # L89-100

    # prephasers (L102-108); gy_pre is the largest PE step and is scaled per ky line
    gx_pre = pp.make_trapezoid("x", system=sys2, area=-gx.area / 2)
    pe_index = cal_pe_index(p.acs_npe, p.nseg, p.pat_factor)
    gy_pre = trapezoid("y", np.abs(pe_index).max() * deltak, p.esp, sys2)
    dur_pre = pp.calc_duration(gx_pre, gy_pre, gz_reph) * 2
    gx_pre = stretch(gx_pre, dur_pre, sys2)
    gz_reph = stretch(gz_reph, dur_pre, sys2)

    # spoilers, tripled in duration to lower PNS (L110-117)
    gx_sp = pp.make_trapezoid("x", system=sys2, area=-2 * p.nx * deltak)
    gy_sp = pp.make_trapezoid("y", system=sys2, area=-2 * p.ny * deltak)
    gz_sp = pp.make_trapezoid("z", system=sys2, area=-3 / p.thickness)
    dur_sp = pp.calc_duration(gx_sp, gy_sp, gz_sp) * 3
    gx_sp, gy_sp, gz_sp, gz_fs = (stretch(g, dur_sp, sys2) for g in (gx_sp, gy_sp, gz_sp, gz_fs))

    # TE / TR (L119-138). The TE formulas are kept exactly as written, including two quirks: the min-TE value leaves
    # out delay_te1, and the explicit-TE branch subtracts the prephaser duration twice.
    esp = pp.calc_duration(gx)
    etl_pre = mround(p.etl / 2 - 1)
    dur_to_center = (etl_pre + 0.5) * esp
    dur_epi = p.etl * esp
    dur_nav = p.nav_etl * esp
    rf_center_incl_delay = rf.delay + pp.calc_rf_center(rf)[0]
    raster = sys.grad_raster_time
    if p.min_te:
        delay_te1 = sys.block_duration_raster * 2
        te = ceil_to_raster(pp.calc_duration(rf, gz) - rf_center_incl_delay + dur_pre + dur_nav + p.esp
                            + dur_to_center, raster)
    else:
        delay_te1 = ceil_to_raster(p.te - pp.calc_duration(rf, gz) + rf_center_incl_delay - dur_nav - dur_to_center
                                   - dur_pre - dur_pre, raster)
        te = ceil_to_raster(pp.calc_duration(rf, gz) - rf_center_incl_delay + dur_pre + dur_nav + dur_to_center
                            + p.esp, raster) + delay_te1
    if delay_te1 < 0:
        raise ValueError("TE too short (L131)")
    dur_slice = (pp.calc_duration(rf_fs) + dur_sp + pp.calc_duration(rf, gz) + dur_pre + delay_te1 + dur_nav
                 + p.esp + dur_epi + dur_sp)

    positions_all = slice_positions(p)
    positions = positions_all
    if p.slice_positions is not None:
        wanted = np.asarray(p.slice_positions)
        keep = [np.any(np.abs(wanted - z) < 1e-6) for z in positions_all]
        positions = positions_all[keep]
        if len(positions) != len(wanted):
            raise ValueError(f"calib slice_positions {wanted} are not all on the slice grid {np.sort(positions_all)}")
    n_played = len(positions)
    if p.min_tr:
        delay_te2 = sys.block_duration_raster * 2
    else:
        delay_te2 = ceil_to_raster((p.tr - n_played * dur_slice) / n_played, raster)
    if delay_te2 < 0:
        raise ValueError(f"TR too short: need >= {n_played * dur_slice * 1e3:.1f} ms (L138)")

    # assemble (L147-192)
    if hw.matlab_compat:
        sys.adc_raster_time = adc.dwell      # L147: the recon reads the ADC dwell back as AdcRasterTime
    seq = pp.Sequence(system=sys)
    w = BlockWriter(seq, hw)
    w.add(pp.make_delay(sys.block_duration_raster * 2), trid=TRID_PREAMBLE)
    w.counter("SET", "REP", 0)
    w.counter("SET", "AVG", len(pe_index))   # the original hard-codes 29 (= pe_lines for P1); nothing reads it
    rf_phase, rf_inc = 0.0, 0.0
    for pe in pe_index:
        w.counter("SET", "SLC", 0)
        scale = pe * deltak / gy_pre.area
        for pos in positions:
            w.add(rf_fs, trid=TRID_SLICE)
            w.add(gz_fs, gx_sp, gy_sp)
            rf.freq_offset = gz.amplitude * pos
            rf.phase_offset = np.deg2rad(rf_phase)          # RF spoiling (increment 0 in P1)
            adc.phase_offset = np.deg2rad(rf_phase)
            rf_inc = (rf_inc + p.rf_spoil_inc_deg) % 360.0
            rf_phase = (rf_phase + rf_inc) % 360.0
            w.add(rf, gz)
            w.add(gx_pre, gz_reph)
            w.add(pp.make_delay(delay_te1))
            w.counter("SET", "LIN", 0)
            for i in range(1, p.nav_etl + 1):                       # navigators: +, -, +
                w.add(pp.scale_grad(gx, (-1) ** ((i + 1) % 2)), adc)
                w.counter("INC", "LIN", 1)
            # phase-encoding step, played during one extra (negative) readout lobe without ADC (L179)
            w.add(pp.scale_grad(gy_pre, scale), pp.scale_grad(gx, (-1) ** ((p.nav_etl + 2) % 2)))
            for i in range(1, p.etl + 1):                           # echoes, all on the same ky line
                w.add(pp.scale_grad(gx, (-1) ** ((i + 1) % 2)), adc)
                w.counter("INC", "LIN", 1)
            w.counter("INC", "SLC", 1)
            w.add(pp.scale_grad(gx_sp, (-1) ** ((p.etl + 2) % 2)), gy_sp, gz_sp)
            w.add(pp.make_delay(delay_te2))
        w.counter("INC", "REP", 1)

    # definitions read back by the recon (L203-218), plus our additions
    seq.set_definition("FOV", [p.fov, p.fov, positions_all.max() - positions_all.min() + p.thickness])
    seq.set_definition("Resolution", [p.fov / p.nx, p.fov / p.ny, p.thickness])
    seq.set_definition("Name", p.name)
    seq.set_definition("TargetGriddedSamples", p.nx)
    seq.set_definition("TE", te)
    seq.set_definition("ETL", p.etl + p.nav_etl)
    seq.set_definition("CNseg", p.nseg)
    seq.set_definition("pe_lines", len(pe_index))
    seq.set_definition("PAT_factor", p.pat_factor)
    seq.set_definition("ReadoutOversamplingFactor", p.ro_os)
    seq.set_definition("ReceiverGainHigh", 1)
    seq.set_definition("TrapezoidGriddingParameters", [gx.rise_time, gx.flat_time, gx.fall_time, adc.delay - gx.delay])
    seq.set_definition("SlicePositions", positions)
    seq.set_definition("actual_SlicePositions", positions)
    seq.set_definition("SliceThickness", p.thickness)
    seq.set_definition("SliceGap", p.slice_gap)
    seq.set_definition("AdcDwell", adc.dwell)

    info = SimpleNamespace(te=te, tr=p.tr, esp=esp, adc=adc, gx=gx, rf=rf, gz=gz, pe_index=pe_index,
                           delay_te1=delay_te1, delay_te2=delay_te2, dur_slice=dur_slice, positions=positions)
    return seq, info
