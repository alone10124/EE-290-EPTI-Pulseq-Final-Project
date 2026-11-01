"""Building blocks shared by the EPTI calibration and imaging sequences.

Both MGH scripts (``legacy/sequences/seq_EPTI_calib_P1.m`` and ``seq_EPTI_imaging_P1.m``) construct these events with
identical code; line numbers below refer to the imaging script.
"""

from types import SimpleNamespace

import numpy as np
import pypulseq as pp

from epti.util import mround


def make_fatsat(fatsat, sys_slow, hw):
    """Gaussian fat-saturation pulse and its spoiler (L72-78).

    The spoiler is returned with a 2 ms duration; callers stretch it to the common spoiler duration later, as the
    scripts do (L148).
    """
    sat_freq = fatsat.ppm * 1e-6 * sys_slow.B0 * sys_slow.gamma
    rf = pp.make_gauss_pulse(np.deg2rad(fatsat.flip_deg), system=sys_slow, duration=fatsat.duration,
                             delay=sys_slow.rf_dead_time,
                             dwell=hw.fatsat_dwell or hw.rf_raster_time, bandwidth=abs(sat_freq),
                             freq_offset=sat_freq, use="saturation")
    rf.phase_offset = -2 * np.pi * rf.freq_offset * pp.calc_rf_center(rf)[0]   # undo the offset-induced phase
    gz = trapezoid("z", fatsat.spoil_area, 2e-3, sys_slow)
    return rf, gz


def make_readout(fov_x, nx, esp, blip_dur, ro_os, sys, hw):
    """EPI readout lobe with ramp sampling and its ADC (L119-131).

    The lobe lasts exactly one echo spacing. Its area is enlarged so that the area *between* the blips (the part not
    overlapped by the triangular blips at either end) equals the k-space width. The ADC covers ``esp - blip_dur`` and
    is centered on the lobe.
    """
    deltak = 1 / fov_x
    k_width = nx * deltak
    extra_area = blip_dur / 2 * blip_dur / 2 * sys.max_slew
    gx = trapezoid("x", k_width + extra_area, esp, sys)
    actual_area = (gx.area
                   - gx.amplitude / gx.rise_time * blip_dur / 2 * blip_dur / 2 / 2
                   - gx.amplitude / gx.fall_time * blip_dur / 2 * blip_dur / 2 / 2)
    gx = pp.scale_grad(gx, k_width / actual_area)

    readout_time = esp - blip_dur
    dwell_nyquist = deltak / gx.amplitude / ro_os
    dwell = mround(dwell_nyquist / hw.adc_raster_time) * hw.adc_raster_time
    n_div = hw.adc_samples_divisor
    n_samples = int(mround(readout_time / dwell / n_div) * n_div)
    adc = pp.make_adc(n_samples, dwell=dwell, delay=blip_dur / 2, system=sys)
    time_to_center = dwell * ((n_samples - 1) / 2 + 0.5)
    adc.delay = mround((gx.rise_time + gx.flat_time / 2 - time_to_center) / sys.rf_raster_time) * sys.rf_raster_time
    return gx, adc


def trapezoid(channel, area, duration, system):
    """Trapezoid with a given area and total duration, built the way MATLAB Pulseq's ``mr.makeTrapezoid`` does.

    MATLAB (``legacy/libs/pulseq-master151/matlab/+mr/makeTrapezoid.m`` L164-197) solves for the *lowest* amplitude
    that fits the duration, gives it the fastest ramp the slew limit allows (rounded up to the raster), and puts the
    rest in the plateau. PyPulseq 1.5.0's ``make_trapezoid(area=..., duration=...)`` instead reuses the ramp time of
    the *shortest possible* gradient for that area, giving long ramps and a lower plateau. Both have the right area
    and duration, but the shapes differ, so the readout lobe, prephasers and spoilers wouldn't match the MGH
    sequence. Found by the MATLAB parity test (tests/test_sequences.py).
    """
    max_slew, max_grad, raster = system.max_slew, system.max_grad, system.grad_raster_time
    dc = 1 / abs(2 * max_slew) + 1 / abs(2 * max_slew)
    if duration ** 2 <= 4 * abs(area) * dc:
        raise ValueError(f"area {area:.4g} 1/m does not fit in {duration * 1e6:.0f} us at this slew rate")
    amplitude = (duration - np.sqrt(duration ** 2 - 4 * abs(area) * dc)) / (2 * dc)
    rise = np.ceil(abs(amplitude) / max_slew / raster) * raster
    if rise == 0:
        rise = raster
    flat = duration - 2 * rise
    amplitude = area / (rise / 2 + rise / 2 + flat)
    if abs(amplitude) > max_grad * (1 + 1e-9):
        raise ValueError(f"trapezoid amplitude {abs(amplitude) / max_grad:.0%} of max_grad")
    return pp.make_trapezoid(channel, system=system, amplitude=amplitude, rise_time=rise, fall_time=rise,
                             flat_time=flat)


def stretch(grad, duration, system):
    """Same-area trapezoid with a fixed duration (used to lengthen spoilers/prephasers and lower PNS)."""
    return trapezoid(grad.channel, grad.area, duration, system)


def label(name, kind, value):
    """``pp.make_label`` with the MATLAB argument order ``(kind, name, value)`` swapped and plain Python types."""
    if isinstance(value, (bool, np.bool_)):
        value = bool(value)
    elif isinstance(value, (int, np.integer)):
        value = int(value)
    return pp.make_label(name, kind, value)


class BlockWriter:
    """Adds blocks to a sequence, applying the hardware profile's labeling conventions.

    * Counter labels (LIN/SLC/REP/AVG) are written only when ``hw.write_counter_labels`` is set. The recon doesn't use
      them, and GE's interpreter doesn't need them. The MATLAB scripts put each label in its own zero-duration block.
      PyPulseq 1.5 can't do that here: its gradient-continuity check compares a block only with the immediately
      preceding one, while MATLAB's skips zero-duration blocks, and the EPTI blips cross lobe boundaries at full
      amplitude. So a label is held and attached to the *next* block instead. Pulseq applies labels at the start of
      the block that carries them, before its ADC, so every ADC still sees the same counter values.
    * TRID labels (``hw.write_trid``) are attached to the first block of each repeating segment, as GE pge2 requires.
    """

    def __init__(self, seq, hw):
        self.seq = seq
        self.hw = hw
        self.pending = []

    def add(self, *events, trid=None):
        events = [e for e in events if e is not None] + self.pending
        self.pending = []
        if trid is not None and self.hw.write_trid:
            events.append(label("TRID", "SET", trid))
        self.seq.add_block(*events)

    def counter(self, kind, name, value):
        if self.hw.write_counter_labels:
            self.pending.append(label(name, kind, value))


def summarize(seq):
    """Small timing summary used in reports."""
    return SimpleNamespace(duration=seq.duration()[0], n_blocks=len(seq.block_events))
