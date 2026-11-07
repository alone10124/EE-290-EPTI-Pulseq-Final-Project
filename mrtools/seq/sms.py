"""Simultaneous multi-slice (SMS / multiband) building blocks in plain PyPulseq + NumPy.

* ``multiband_sinc_pulse``: one RF pulse that excites ``mb`` slices at once. It's a single-band sinc pulse
  multiplied by a sum of complex exponentials, one per slice, each shifting the excitation by that slice's offset
  (``frequency = slice-select gradient * distance``). Each band keeps the single-band flip angle; the peak RF grows
  with ``mb``.
* ``caipi_kz_index``: the blipped-CAIPI pattern. Small z-gradient blips between EPI lines give slice ``s`` an extra
  phase of ``2 pi s / mb`` per line, which shifts it by ``FOV/mb`` in the phase-encode direction. The collapsed slices
  then overlap less, which makes unfolding them with coil maps better conditioned (Setsompop et al., MRM 2012).

The default sub-pulse phases (Wong, ISMRM 2012, p. 2209, Table 1) lower the peak RF power. The MGH EPTI code used the
same table (``getsmsphase.m`` on the ``epti`` branch).
"""

import numpy as np
import pypulseq as pp

WONG_PHASES = {1: [0.0], 2: [0.0, np.pi], 3: [0.0, 0.730, 4.602], 4: [0.0, 3.875, 5.940, 6.197]}


def multiband_sinc_pulse(flip_deg, duration, slice_thickness, slice_sep, mb, system, time_bw_product=4.0,
                         apodization=0.5, phases=None):
    """Multiband sinc excitation centered on z = 0 (set ``rf.freq_offset = gz.amplitude * z_center`` to move it).

    Parameters
    ----------
    flip_deg : flip angle of each band (degrees).
    duration, slice_thickness, time_bw_product, apodization : as ``pp.make_sinc_pulse``.
    slice_sep : center-to-center distance between simultaneously excited slices (m).
    mb : multiband factor (number of slices).
    phases : per-band phases (rad); defaults to the Wong table for ``mb <= 4``, else zeros.

    Returns ``(rf, gz, gz_reph)``, like ``pp.make_sinc_pulse(..., return_gz=True)``.
    """
    rf, gz, gz_reph = pp.make_sinc_pulse(np.deg2rad(flip_deg), system=system, duration=duration,
                                         slice_thickness=slice_thickness, time_bw_product=time_bw_product,
                                         apodization=apodization, delay=system.rf_dead_time, return_gz=True,
                                         use="excitation")
    if phases is None:
        phases = WONG_PHASES.get(mb, [0.0] * mb)
    offsets = (np.arange(mb) - (mb - 1) / 2) * slice_sep             # band positions relative to the group center
    t = rf.t - rf.center                                             # time relative to the pulse center
    modulation = sum(np.exp(1j * (2 * np.pi * gz.amplitude * z * t + ph)) for z, ph in zip(offsets, phases))
    rf_mb = pp.make_arbitrary_rf(rf.signal * modulation, np.deg2rad(flip_deg), no_signal_scaling=True,
                                 delay=rf.delay, system=system, use="excitation", center=rf.center)
    return rf_mb, gz, gz_reph


def caipi_kz_index(mb, n_lines):
    """Blipped-CAIPI kz state of each EPI line: 0, 1, ..., mb-1, 0, 1, ... (in units of ``caipi_delta_kz``).

    The z blip before line ``n`` has area ``(index[n] - index[n-1]) * delta_kz``: a small +1 step, or a rewind of
    ``-(mb-1)`` every ``mb`` lines.
    """
    return np.arange(n_lines) % mb


def caipi_delta_kz(mb, slice_sep):
    """kz step (1/m) that gives neighboring simultaneous slices a phase difference of ``2 pi / mb`` per step."""
    return 1.0 / (mb * slice_sep)

