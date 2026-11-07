"""Basic sequence building blocks in plain PyPulseq: fat saturation, an EPI readout lobe, TRID labels, and slice
ordering.

Current simplifying assumption (see docs/LESSONS_LEARNED.md): the ADC samples only on the readout's flat top, so the
k-space samples are evenly spaced and no ramp-sampling regridding is needed. Ramp sampling, which shortens the echo
spacing, can be added later; the ``epti`` branch has a tested regridding implementation.
"""

import numpy as np
import pypulseq as pp


def make_fatsat(system, ppm=-3.4, flip_deg=110.0, duration=5.12e-3, spoil_area=-1000.0):
    """Gaussian fat-saturation pulse (at ``ppm`` from water) and the gradient spoiler that follows it.

    The fat-sat phase is set so that the frequency offset doesn't leave a phase at the pulse center.
    ``spoil_area`` is in 1/m (-1000 1/m dephases by one cycle over 1 mm).
    """
    sat_freq = ppm * 1e-6 * system.B0 * system.gamma
    rf = pp.make_gauss_pulse(np.deg2rad(flip_deg), system=system, duration=duration, delay=system.rf_dead_time,
                             bandwidth=abs(sat_freq), freq_offset=sat_freq, use="saturation")
    rf.phase_offset = -2 * np.pi * rf.freq_offset * pp.calc_rf_center(rf)[0]
    spoiler = pp.make_trapezoid("z", system=system, area=spoil_area)
    return rf, spoiler


def make_epi_readout(fov, nx, system, dwell=4e-6):
    """One EPI readout lobe and its ADC, sampling only on the flat top.

    Parameters
    ----------
    fov : readout field of view (m).
    nx : readout matrix size.
    dwell : ADC dwell (s); rounded to the ADC raster.

    Returns ``(gx, adc)``. The lobe's flat-top area is ``nx / fov`` (the k-space width), and the ADC is centered on the
    flat top with ``nx`` samples. Negative lobes are ``pp.scale_grad(gx, -1)`` with the same ADC;
    their samples come out in reverse kx order.
    """
    dwell = round(dwell / system.adc_raster_time) * system.adc_raster_time
    adc_duration = nx * dwell
    flat_time = np.ceil(adc_duration / system.grad_raster_time - 1e-9) * system.grad_raster_time
    gx = pp.make_trapezoid("x", system=system, flat_area=nx / fov * flat_time / adc_duration, flat_time=flat_time)
    delay = gx.rise_time + (flat_time - adc_duration) / 2
    delay = round(delay / system.adc_raster_time) * system.adc_raster_time
    adc = pp.make_adc(nx, dwell=dwell, delay=delay, system=system)
    return gx, adc


def trid(value):
    """TRID label: marks the first block of each repeating segment, which GE's pge2 interpreter requires."""
    return pp.make_label("TRID", "SET", int(value))


def slice_positions(slab, thickness):
    """z position (m) of each contiguous slice center covering ``slab``, inferior to superior, centered on z = 0."""
    n = int(round(slab / thickness))
    return (np.arange(n) - (n - 1) / 2) * thickness


def interleaved_slice_order(n_slices):
    """Slice indices in acquisition order: even indices first, then odd (0, 2, 4, ..., 1, 3, 5, ...).

    Neighbouring slices are then excited half a TR apart, so the imperfect edges of one slice profile don't
    saturate the next slice just before it is excited.
    """
    return np.concatenate([np.arange(0, n_slices, 2), np.arange(1, n_slices, 2)])
