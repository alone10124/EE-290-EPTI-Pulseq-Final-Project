"""Simultaneous multi-slice (SMS / multiband) excitation pulse for the EPTI imaging scan.

Ports the chain ``getsmspulse.m`` -> ``toppe.utils.rf.makeslr`` -> John Pauly's SLR tools that the MGH imaging script
uses (seq_EPTI_imaging_P1.m L80-102). Sources:

* ``legacy/libs/utils/GE/getsmspulse.m``, ``getsmsphase.m``, ``rf2pulseq.m``, ``g2pulseq.m`` (MIT, (c) 2020 harmonizedMRI)
* ``legacy/libs/toppe-main/+toppe/+utils/+rf/makeslr.m`` (``sub_myslrrf``), ``+toppe/+utils/trapwave2.m``,
  ``+toppe/makeGElength.m`` (TOPPE, J.-F. Nielsen)
* The SLR/least-squares filter itself comes from ``sigpy.mri.rf.dzrf`` instead of TOPPE's copy of Pauly's ``dzrf``.
  SigPy builds the 200-tap filter from a 201-tap ``firls`` design plus a half-sample shift. MATLAB's ``firls(199)``
  solves the even-length problem directly, so the two are close but not bit-identical; tests/test_rf_sms.py measures
  the difference against ``fixtures.mat``.

The base pulse is designed on GE's 4 us raster (TOPPE's native raster), modulated to MB bands with the phase schedule
of Wong (ISMRM 2012, p. 2209), then interpolated to the Pulseq RF/gradient rasters.
"""

from math import gcd

import numpy as np
import pypulseq as pp
from scipy import signal
from sigpy.mri import rf as sigpy_rf

from epti.util import ceil_to_raster, mround

GE_RASTER = 4e-6            # s, TOPPE waveform raster
GAMMA_HZ_PER_T = 42.576e6   # toppe.systemspecs sys.gamma
GAMMA_HZ_PER_G = 4257.6     # rf2pulseq.m / g2pulseq.m

# Phases (rad) of the SMS sub-pulses, from Table 1 of Wong, ISMRM 2012 p. 2209 (getsmsphase.m)
SMS_PHASES = {
    1: [0.0],
    2: [0.0, np.pi],
    3: [0.0, 0.730, 4.602],
    4: [0.0, 3.875, 5.940, 6.197],
    5: [0.0, 3.778, 5.335, 0.872, 0.471],
    6: [0.0, 2.005, 1.674, 5.012, 5.736, 4.123],
    7: [0.0, 3.002, 5.998, 5.909, 2.624, 2.528, 2.440],
    8: [0.0, 1.036, 3.414, 3.778, 3.215, 1.756, 4.555, 2.467],
}


# ----------------------------------------------------------------------------------------------- MATLAB helpers
def matlab_resample(x, p, q, n=10, beta=5.0):
    """MATLAB ``resample(x, p, q)`` with its default anti-aliasing filter.

    MATLAB designs ``firls(L-1, [0 2fc 2fc 1], [1 1 0 0]) .* kaiser(L, 5)`` with ``L = 2*n*max(p,q) + 1`` and
    ``fc = 1/(2*max(p,q))``, normalized to gain ``p``. ``scipy.signal.resample_poly`` multiplies any filter it is given
    by ``p`` itself, so we pass the unit-gain version. Both then apply it zero-phase with output length ceil(len*p/q).
    """
    g = gcd(int(p), int(q))
    p, q = int(p) // g, int(q) // g
    if p == q == 1:
        return np.array(x, dtype=float)
    pq = max(p, q)
    fc = 1 / (2 * pq)
    length = 2 * n * pq + 1
    h = signal.firls(length, [0, 2 * fc, 2 * fc, 1], [1, 1, 0, 0]) * signal.windows.kaiser(length, beta)
    h = h / h.sum()
    return signal.resample_poly(np.asarray(x, dtype=float), p, q, window=h)


def matlab_smooth(x, span=5):
    """MATLAB ``smooth(x)``: moving average with the span shrinking symmetrically at the ends."""
    x = np.asarray(x, dtype=float)
    half = span // 2
    out = np.empty_like(x)
    for i in range(len(x)):
        k = min(i, len(x) - 1 - i, half)
        out[i] = x[i - k:i + k + 1].mean()
    return out


def _make_ge_length(w):
    """Zero-pad to a multiple of 4 samples (toppe.makeGElength)."""
    pad = (-len(w)) % 4
    return np.concatenate([w, np.zeros(pad)])


def _trapwave2(area, mxg, mxs, dt):
    """Shortest trapezoid/triangle with the given area (toppe.utils.trapwave2).

    area in G/cm*s, mxg in G/cm, mxs in G/cm/ms, dt in ms. Returns the waveform in G/cm on the dt raster.
    """
    area = area * 1e3                 # G/cm*ms
    invert = area < 0
    area = abs(area)
    mxg, mxs = 0.995 * mxg, 0.995 * mxs
    tr = mxg / mxs                    # ramp time to full scale
    a_crit = mxs * tr ** 2            # area of the largest triangle
    dg = mxs * dt
    if area < a_crit:
        n = int(np.ceil(np.sqrt(area / mxs) / dt))
        ramp = np.arange(n + 1) * dg
        wave = np.concatenate([ramp, ramp[::-1]])
    else:
        nr = int(np.ceil(tr / dt))
        ramp = np.arange(nr) * mxs * dt
        n_plat = int(np.ceil((area - 2 * ramp.sum() * dt) / mxg / dt))
        wave = np.concatenate([ramp, np.full(n_plat, mxg), ramp[::-1]])
    wave = wave / (wave.sum() * dt) * area
    return -wave if invert else wave


def _interp_linear(xq, x, y, extrapolate):
    out = np.interp(xq, x, y)
    if extrapolate:
        lo, hi = xq < x[0], xq > x[-1]
        out[lo] = y[0] + (xq[lo] - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
        out[hi] = y[-1] + (xq[hi] - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
    return out


def _colon_count(start, step, stop):
    """Number of elements in MATLAB's ``start:step:stop``."""
    return int(np.floor((stop - start) / step + 1e-9)) + 1


def rf_to_pulseq(rf_gauss, raster_in, raster_out):
    """Gauss -> Hz and linear interpolation (with extrapolation) to the Pulseq RF raster (rf2pulseq.m)."""
    rf = np.asarray(rf_gauss) * GAMMA_HZ_PER_G
    t_in = np.arange(1, len(rf) + 1) * raster_in - raster_in / 2
    t_out = raster_out / 2 + raster_out * np.arange(_colon_count(raster_out / 2, raster_out, len(rf) * raster_in))
    return (_interp_linear(t_out, t_in, rf.real, True) + 1j * _interp_linear(t_out, t_in, rf.imag, True))


def grad_to_pulseq(g_gcm, raster_in, raster_out):
    """G/cm -> Hz/m and linear interpolation to the Pulseq gradient raster (g2pulseq.m)."""
    if raster_out < raster_in:
        raise ValueError("raster_out must be >= raster_in")
    g = np.asarray(g_gcm, dtype=float) * GAMMA_HZ_PER_G * 100
    t_in = np.arange(1, len(g) + 1) * raster_in - raster_in / 2
    t_out = raster_out / 2 + raster_out * np.arange(_colon_count(raster_out / 2, raster_out, t_in[-1]))
    return np.interp(t_out, t_in, g)


# ----------------------------------------------------------------------------------------------- pulse design
def slr_pulse(flip_deg, slthick_cm, tbw, dur_ms, max_grad_gcm, max_slew_gcm_ms):
    """Small-tip ('st') least-squares ('ls') SLR pulse and slice-select gradient on the 4 us GE raster.

    Port of ``toppe.utils.rf.makeslr(flip, slthick, tbw, dur, 1e-6, sysGE, 'type','st','ftype','ls')`` as called by
    getsmspulse.m. ``ncycles=1e-6`` means no prephaser and no spoiler bridge, so only ``sub_myslrrf`` matters:
    makeslr.m L45-55 and the ``sub_myslrrf`` body.

    Returns ``(rf, gz, center)``: RF in Gauss and gradient in G/cm (same length, with a rephaser lobe at the end of
    ``gz``), and ``center``, the 0-based (fractional) sample index of the pulse's center of symmetry.

    Note on peaks: MATLAB's even-length ``firls`` filter is exactly symmetric, so after 5x resampling its two central
    taps are *exactly* equal maxima, and the MATLAB code (``find(rf == max(rf))``) picks up both. SigPy's filter is
    symmetric only to ~1e-3, so a literal ``argmax`` would pick one of them and shift every downstream time
    reference by half a tap spacing (10 us here). We therefore use the two central taps by construction.
    """
    dt = GE_RASTER * 1e3                   # ms
    res = int(mround(dur_ms / dt))
    dur = res * dt
    mxs = 0.99 * max_slew_gcm_ms           # makeslr passes 0.99*maxSlew into sub_myslrrf
    nrf = 200
    if round(tbw) != tbw:
        raise ValueError("tbw must be an integer for small-tip design")

    rf = np.real(sigpy_rf.dzrf(nrf, tbw, "st", "ls"))
    rf = matlab_resample(rf, res, nrf)
    rf = np.deg2rad(flip_deg) * rf / rf.sum()          # flip angle (rad) = sum(rf)
    gamma = GAMMA_HZ_PER_T * 1e-7                      # kHz/Gauss
    rf = rf / gamma / dt / 2 / np.pi                   # Gauss
    if len(rf) % 2:
        rf = np.append(rf, 0.0)
    npix = len(rf)
    up = res / nrf
    center_rf = up * (nrf - 1) / 2                     # symmetry center of the resampled filter (0-based)

    # reference point (peak) for the rephaser area: MATLAB uses the first of the two central maxima
    i_peak = int(mround(center_rf - up / 2))           # 0-based; MATLAB I = i_peak + 1
    s = matlab_smooth(rf[i_peak - 20:i_peak + 20])     # MATLAB rf((I-20):(I+19))
    # MATLAB: iref = find(s==max(s)) + I - 20 (1-based; one sample past the smoothed peak, kept as in the original)
    iref = int(np.argmax(s)) + i_peak - 19             # 0-based index of that same sample

    gplateau = (tbw / dur) / (gamma * slthick_cm)      # G/cm
    if gplateau > max_grad_gcm:
        raise ValueError("slice-select gradient exceeds max_grad")
    step = mxs * dt * 0.995
    ramp = step * np.arange(1, _colon_count(step, step, gplateau) + 1)
    if ramp.size == 0:
        ramp = np.array([0.0])
    if gplateau - ramp[-1] > step:
        ramp = np.append(ramp, (gplateau + ramp[-1]) / 2)
    trap = np.concatenate([ramp, np.full(npix, gplateau), ramp[::-1]])
    iref += len(ramp)
    area_rep = trap[iref:].sum() * dt * 1e-3            # G/cm*s
    gzrep = -_trapwave2(area_rep, max_grad_gcm, mxs, dt)

    gz = _make_ge_length(np.concatenate([trap, gzrep]))
    rf = _make_ge_length(np.concatenate([np.zeros(len(ramp)), rf, np.zeros(len(ramp) + len(gzrep))]))
    assert gz.sum() > 0, "makeslr would insert a spoiler bridge here; not ported"
    rf = np.concatenate([[0.0], rf, [0.0]])
    gz = np.concatenate([[0.0], gz, [0.0]])
    n = max(len(rf), len(gz))
    center = 1 + len(ramp) + center_rf                 # leading [0] and the ramp-length zero pad come first
    return np.pad(rf, (0, n - len(rf))), np.pad(gz, (0, n - len(gz))), center


def sms_pulse(flip_deg, thickness, tbw, dur, mb, slice_sep, system, ge_max_grad_gcm, ge_max_slew_gcm_ms,
              matlab_compat=False):
    """Multiband RF event plus the raw slice-select waveform (port of getsmspulse.m).

    Parameters
    ----------
    flip_deg, thickness (m), tbw, dur (s), mb, slice_sep (m): pulse spec.
    system : pp.Opts used for the Pulseq conversion (``sys2`` in the imaging script).
    ge_max_grad_gcm, ge_max_slew_gcm_ms : limits for the TOPPE design (``sysGE`` in the script).
    matlab_compat : reproduce two quirks of the original, for the parity test only:

        * the zero-trimming test ``find(abs(diff(rfp) > 0))`` compares *real parts* of complex values in MATLAB (the
          MB=2 pulse is purely imaginary, so the trim point is set by round-off). Otherwise we use ``abs(diff)``.
        * the RF center is left to Pulseq's peak finder, which may lock onto one lobe of the two-lobe MB envelope.
          Otherwise we set it explicitly to the symmetric center of the base pulse.

    Returns
    -------
    rf : PyPulseq RF event (signal in Hz; ``use='excitation'``).
    gz_waveform : slice-select gradient (Hz/m) on ``system.grad_raster_time``, including the TOPPE rephaser.
    gz_delay : delay (s) the original applies to ``gz`` when the RF delay is raised to the dead time.
    freq : frequency offset (Hz) corresponding to ``slice_sep`` (unused by the scripts).
    """
    rf1, gz, center = slr_pulse(flip_deg, thickness * 100, tbw, dur * 1e3, ge_max_grad_gcm, ge_max_slew_gcm_ms)

    # modulate: one band per slice, referenced to the base-pulse peak (getsmspulse.m L50-64).
    # MATLAB: I = find(abs(rf1) == max(abs(rf1))) finds the two central taps, whose mean is the symmetry center;
    # iPeak = mean(I) + 2 in 1-based indexing. See slr_pulse for why we don't search for the maximum.
    i_peak = (center + 1) + 2
    gplateau = gz.max()
    t = GE_RASTER * np.arange(1, len(rf1) + 1) - GE_RASTER * i_peak
    rf = np.zeros(len(rf1), dtype=complex)
    for sl, phs in enumerate(SMS_PHASES[mb]):
        offset_cm = (-mb / 2 + 0.5 + sl) * slice_sep * 100
        f = GAMMA_HZ_PER_T * gplateau * 1e-4 * offset_cm
        rf += rf1 * np.exp(1j * 2 * np.pi * f * t) * np.exp(1j * phs)
    freq = GAMMA_HZ_PER_T * gplateau * 1e-4 * slice_sep * 100

    # convert to Pulseq units and rasters (getsmspulse.m L100-133)
    rfp = rf_to_pulseq(rf, GE_RASTER, system.rf_raster_time)
    gzp = grad_to_pulseq(gz, GE_RASTER, system.grad_raster_time)

    def first_change(w):
        d = np.diff(w)
        changed = np.real(d) > 0 if matlab_compat else np.abs(d) > 0
        return int(np.nonzero(changed)[0][0]) + 1              # 1-based, as MATLAB find()

    ratio = int(mround(system.grad_raster_time / system.rf_raster_time))
    i1 = first_change(rfp)
    n_delay = i1 - (i1 % ratio)
    rfp = rfp[n_delay:]
    delay = n_delay * system.rf_raster_time
    j1 = first_change(rfp[::-1])
    rfp = rfp[:len(rfp) - j1 + 1]
    wavdur = len(rfp) * system.rf_raster_time
    target = ceil_to_raster(wavdur, system.grad_raster_time)
    rfp = np.concatenate([rfp, np.zeros(int(mround((target - wavdur) / system.rf_raster_time)))])
    gz_delay = 0.0
    if delay < system.rf_dead_time:
        gz_delay = system.rf_dead_time - delay
        delay = system.rf_dead_time

    rf_center = None
    if not matlab_compat:
        # symmetry center of the base pulse (sample k of rf1 sits at (k + 0.5) * raster), relative to the trimmed start
        rf_center = (center + 0.5) * GE_RASTER - n_delay * system.rf_raster_time

    rf_event = pp.make_arbitrary_rf(rfp, flip_angle=np.deg2rad(flip_deg), no_signal_scaling=True, delay=delay,
                                    system=system, use="excitation", center=rf_center)
    return rf_event, gzp, gz_delay, freq


def slice_select_trapezoid(gz_waveform, system):
    """Rebuild the TOPPE slice-select waveform as a plain trapezoid (seq_EPTI_imaging_P1.m L86-99).

    The original keeps only the plateau amplitude and ramp/flat durations found on the waveform; the TOPPE rephaser
    is discarded and a new one is made by the caller.
    """
    w = np.asarray(gz_waveform)
    ind1 = int(np.argmax(w))                  # start of plateau (0-based)
    ind2 = int(np.argmin(w))                  # peak of the (negative) rephaser
    ind3 = int(np.argmin(np.abs(w[ind1:ind2 + 1]))) + 1   # 1-based, relative to ind1 (zero crossing)
    raster = system.grad_raster_time
    flat_time = (ind3 - (ind1 + 1)) * raster  # MATLAB: (ind3 - ind1) * raster with 1-based ind1
    rise = ind1 * raster                      # MATLAB: (ind1 - 1) * raster
    return pp.make_trapezoid("z", system=system, amplitude=w[ind1], flat_time=flat_time,
                             rise_time=rise, fall_time=rise)
