"""EPTI imaging scan: 2D single-shot SMS EPTI with gradient-echo contrast.

Port of ``legacy/sequences/seq_EPTI_imaging_P1.m`` (MGH, Jian Wu, Zijing Dong, Fuyixue Wang; version 09/2025).
Line numbers in comments refer to that script.

Each TR excites one slice group (``mb`` slices, ``slice_sep`` apart) with an SMS pulse and plays a 128-lobe EPI train:
3 navigator lobes at the k-space center (for odd/even ghost correction) followed by ``etl`` lobes that zig-zag through
ky (Gy blips, ``calc_pattern_wjg3``) while alternating kz (Gz blipped-CAIPI blips).

Additions to the original, all off by default:

* ``protocol.slice_groups`` plays only some slice groups (for simulation), keeping the full-volume geometry and TR.
* ``protocol.n_dynamics`` repeats the volume inside the file.
* TRID labels and dropping LIN/SLC labels or the trigger, controlled by the hardware profile (GE pge2).
"""

from types import SimpleNamespace

import numpy as np
import pypulseq as pp

from epti.config import ImagingProtocol
from epti.hardware import HardwareProfile
from epti.seq.blips import epti_blip
from epti.seq.common import BlockWriter, make_fatsat, make_readout, stretch, trapezoid
from epti.seq.patterns import calc_pattern_sms_blip, calc_pattern_wjg3, interleaved_slice_order
from epti.seq.rf_sms import sms_pulse, slice_select_trapezoid
from epti.util import ceil_to_raster

TRID_SLICE_GROUP = 1
TRID_PREAMBLE = 2
TRID_TRIGGER = 3


def group_positions(p: ImagingProtocol):
    """Excitation offset (m) of every slice group, by spatial group index k (L174-175)."""
    k = np.arange(p.n_groups)
    return -(p.thickness + p.slice_gap) * (k - (p.n_groups - 1) / 2)


def excited_slice_positions(p: ImagingProtocol, groups=None):
    """Physical slice positions (m) excited by the given spatial slice groups (all groups if None).

    The SMS pulse for group k is centered at ``group_positions(p)[k]`` with bands at
    ``(-mb/2 + 0.5 + sl) * slice_sep`` around it (getsmspulse.m L58). Use this to give the calibration scan the same
    slices as a simulated subset of the imaging scan.
    """
    groups = range(p.n_groups) if groups is None else groups
    centers = group_positions(p)[list(groups)]
    bands = (-p.mb / 2 + 0.5 + np.arange(p.mb)) * p.slice_sep
    return np.sort((centers[:, None] + bands[None, :]).ravel())


def build_imaging(p: ImagingProtocol, hw: HardwareProfile):
    """Build the EPTI imaging sequence. Returns ``(seq, info)``; ``info`` holds derived timing for reports."""
    sys, sys2 = hw.opts(), hw.opts(slow=True)
    if p.mb > 4:
        raise ValueError("only MB <= 4 is supported (L33)")
    if p.n_slices % p.mb:
        raise ValueError("n_slices must be a multiple of mb (L70)")
    if abs(round(p.blip_dur / sys.grad_raster_time) * sys.grad_raster_time - p.blip_dur) > 1e-9:
        raise ValueError("blip_dur must be on the gradient raster (L69)")
    slice_sep = p.slice_sep

    # --- fat saturation (L72-78)
    rf_fs, gz_fs = make_fatsat(p.fatsat, sys2, hw)

    # --- SMS excitation (L80-102)
    rf, gz_wave, _, _ = sms_pulse(
        p.rf_angle_deg, p.thickness, p.rf_tbw, p.rf_dur, p.mb, slice_sep, sys2,
        ge_max_grad_gcm=0.5 * sys.max_grad / sys.gamma * 100,      # sysGE limits as in L80-82
        ge_max_slew_gcm_ms=0.4 * sys.max_slew / sys.gamma / 10,
        matlab_compat=hw.matlab_compat)
    gz = slice_select_trapezoid(gz_wave, sys)
    gz_reph = trapezoid("z", -gz.area / 2, p.prephase_dur, sys2)

    # --- EPTI k-space pattern and blip durations (L105-118)
    deltakx, deltaky = 1 / p.fov[0], 1 / p.fov[1]
    ky_acq, _ = calc_pattern_wjg3(p.nseg, p.rseg, p.etl + p.nav_etl, p.rpe, 2, p.pf_shift)
    ky_step = np.diff(ky_acq[0])[2:]          # one step per imaging lobe (the 3 navigators sit at the center)
    kz_step, deltakz = calc_pattern_sms_blip(p.mb, p.etl, slice_sep)

    def min_blip(area):
        return np.ceil(2 * np.sqrt(area / sys.max_slew) / sys.grad_raster_time / 2) * sys.grad_raster_time * 2

    needed = max(min_blip(deltaky * np.abs(ky_step).max()), min_blip(deltakz * np.abs(kz_step).max()))
    if needed > p.blip_dur + 1e-12:
        raise ValueError(f"blip_dur {p.blip_dur*1e6:.0f} us too short; blips need {needed*1e6:.0f} us (L118)")

    # --- readout lobe + ADC (L119-131)
    gx, adc = make_readout(p.fov[0], p.nx, p.esp, p.blip_dur, p.ro_os, sys, hw)

    # --- prephasers (L133-139); the EPTI pattern starts at the ky center, so gyPre has zero area
    gx_pre = pp.make_trapezoid("x", system=sys2, area=-gx.area / 2)
    dur_pre = pp.calc_duration(gz_reph)
    gz_reph = stretch(gz_reph, dur_pre, sys2)
    gx_pre = stretch(gx_pre, dur_pre, sys2)
    gy_pre = trapezoid("y", 0.0, dur_pre, sys2)

    # --- spoilers, doubled in duration to lower PNS (L141-148)
    gx_sp = pp.make_trapezoid("x", system=sys2, area=-2 * p.nx * deltakx)
    gy_sp = pp.make_trapezoid("y", system=sys2, area=-2 * p.ny * deltaky)
    gz_sp = pp.make_trapezoid("z", system=sys2, area=-4 / p.thickness)
    dur_sp = pp.calc_duration(gx_sp, gy_sp, gz_sp) * 2
    gx_sp, gy_sp, gz_sp, gz_fs = (stretch(g, dur_sp, sys2) for g in (gx_sp, gy_sp, gz_sp, gz_fs))

    # --- TE and TR delays (L150-169)
    n_lobes = p.etl + p.nav_etl
    dur_epi = n_lobes * pp.calc_duration(gx)
    rf_center_incl_delay = rf.delay + pp.calc_rf_center(rf)[0]
    if p.min_te:
        delay_te1 = sys.block_duration_raster * 2
        te1 = (pp.calc_duration(rf, gz) - rf_center_incl_delay + pp.calc_duration(gx_pre, gy_pre, gz_reph)
               + sys.block_duration_raster * 2 + (n_lobes / 2 - 0.5) * p.esp)
    else:
        te1 = p.te1
        delay_te1 = ceil_to_raster(te1 - pp.calc_duration(rf, gz) + rf_center_incl_delay
                                   - (n_lobes / 2 - 0.5) * p.esp - dur_pre, sys.grad_raster_time)
    if delay_te1 < 0:
        raise ValueError("TE too short for this protocol (L160)")
    trig = pp.make_digital_output_pulse("ext1", duration=100e-6, delay=500e-6) if hw.trigger_out else None
    dur_trig = pp.calc_duration(trig) if trig is not None else 0.0
    dur_group = (pp.calc_duration(rf_fs) + dur_sp + pp.calc_duration(rf, gz) + dur_pre + delay_te1
                 + dur_epi + dur_sp)

    # --- slice-group order (L171-176) and the subset to play
    acq_order = interleaved_slice_order(p.n_groups)
    if p.slice_groups is not None:
        bad = set(p.slice_groups) - set(range(p.n_groups))
        if bad:
            raise ValueError(f"slice_groups out of range 0..{p.n_groups - 1}: {sorted(bad)}")
        acq_order = np.array([k for k in acq_order if k in set(p.slice_groups)])
    positions = group_positions(p)[acq_order]
    n_played = len(acq_order)
    if p.min_tr:
        delay_te2 = sys.block_duration_raster * 2
    else:
        # per-group delay so that the volume still takes TR when only some groups are played
        delay_te2 = ceil_to_raster((p.tr - dur_trig - n_played * dur_group) / n_played, sys.grad_raster_time)
    if delay_te2 < 0:
        raise ValueError(f"TR too short: need >= {(dur_trig + n_played * dur_group)*1e3:.1f} ms (L169)")

    # --- assemble (L178-228)
    if hw.matlab_compat:
        sys.adc_raster_time = adc.dwell      # L178: the recon reads the ADC dwell back as AdcRasterTime
    seq = pp.Sequence(system=sys)
    w = BlockWriter(seq, hw)
    w.add(pp.make_delay(sys.block_duration_raster * 2), trid=TRID_PREAMBLE)
    for dyn in range(p.n_dynamics):
        w.counter("SET", "SLC", 0)
        if p.n_dynamics > 1:
            w.counter("SET", "REP", dyn)
        if trig is not None:
            w.add(trig, trid=TRID_TRIGGER)
        for pos in positions:
            w.add(rf_fs, trid=TRID_SLICE_GROUP)
            w.add(gz_fs, gx_sp, gy_sp)
            rf.freq_offset = gz.amplitude * pos
            w.add(rf, gz)
            w.add(gx_pre, gy_pre, gz_reph)
            w.add(pp.make_delay(delay_te1))
            w.counter("SET", "LIN", 0)
            adc.phase_offset = 0.0
            for i in range(1, p.nav_etl + 1):                     # navigators: +, -, + lobes at ky = center
                lobe = pp.scale_grad(gx, (-1) ** ((i + 1) % 2))
                if i == p.nav_etl:                                # rising half of the first blips (state 1)
                    gyb = epti_blip("y", ky_step[0], ky_step[0], deltaky, p.blip_dur, gx, sys, 1)
                    gzb = epti_blip("z", kz_step[0], kz_step[0], deltakz, p.blip_dur, gx, sys, 1)
                    w.add(lobe, gyb, gzb, adc)
                else:
                    w.add(lobe, adc)
                w.counter("INC", "LIN", 1)
            phi = 0.0
            for i in range(1, p.etl + 1):                         # imaging lobes
                # ADC phase that undoes the kz-blip phase at this slice-group offset (L208-210)
                phi += 2 * np.pi * deltakz * kz_step[i - 1] * (pos + slice_sep * (p.mb - 1) / 2)
                adc.phase_offset = phi
                lobe = pp.scale_grad(gx, (-1) ** (i % 2))
                if i == p.etl:
                    gyb = epti_blip("y", ky_step[i - 1], ky_step[i - 1], deltaky, p.blip_dur, gx, sys, 3)
                    gzb = epti_blip("z", kz_step[i - 1], kz_step[i - 1], deltakz, p.blip_dur, gx, sys, 3)
                else:
                    gyb = epti_blip("y", ky_step[i - 1], ky_step[i], deltaky, p.blip_dur, gx, sys, 2)
                    gzb = epti_blip("z", kz_step[i - 1], kz_step[i], deltakz, p.blip_dur, gx, sys, 2)
                w.add(lobe, gyb, gzb, adc)
                w.counter("INC", "LIN", 1)
            w.counter("INC", "SLC", 1)
            w.add(gx_sp, gy_sp, gz_sp)
            w.add(pp.make_delay(delay_te2))

    # --- definitions read back by the recon (L240-264), plus our additions at the end
    order_all = interleaved_slice_order(p.n_slices)
    slice_positions = -((p.thickness + p.slice_gap) * (np.arange(p.n_slices) - (p.n_slices - 1) / 2))[order_all]
    seq.set_definition("FOV", [p.fov[0], p.fov[1], slice_positions.max() - slice_positions.min() + p.thickness])
    seq.set_definition("Resolution", [p.fov[0] / p.nx, p.fov[1] / p.ny, p.thickness])
    seq.set_definition("Name", p.name)
    seq.set_definition("TargetGriddedSamples", p.nx)
    seq.set_definition("TE", te1)
    seq.set_definition("GRE_Lines", n_lobes)
    seq.set_definition("Nseg", p.nseg)
    seq.set_definition("Rpe", p.rpe)
    seq.set_definition("Rseg", p.rseg)
    seq.set_definition("PFshift", p.pf_shift)
    seq.set_definition("SMS", p.mb)
    seq.set_definition("pattern_GE", ky_acq[0])
    seq.set_definition("ReadoutOversamplingFactor", p.ro_os)
    seq.set_definition("ReceiverGainHigh", 1)
    seq.set_definition("TrapezoidGriddingParameters", [gx.rise_time, gx.flat_time, gx.fall_time, adc.delay - gx.delay])
    seq.set_definition("SlicePositions", slice_positions)
    seq.set_definition("actual_SlicePositions", positions)
    seq.set_definition("SliceThickness", p.thickness)
    seq.set_definition("SliceGap", p.slice_gap)
    # additions (not in the original): explicit values instead of the AdcRasterTime hack, and the sim subset
    seq.set_definition("AdcDwell", adc.dwell)
    seq.set_definition("SliceSep", slice_sep)
    seq.set_definition("NumDynamics", p.n_dynamics)
    seq.set_definition("PlayedSliceGroups", acq_order)

    info = SimpleNamespace(
        te1=te1, tr=p.tr, esp=pp.calc_duration(gx), adc=adc, gx=gx, rf=rf, gz=gz,
        delay_te1=delay_te1, delay_te2=delay_te2, dur_group=dur_group, n_groups_played=n_played,
        ky_acq=ky_acq[0], kz_step=kz_step, deltakz=deltakz, positions=positions, acq_order=acq_order,
    )
    return seq, info
