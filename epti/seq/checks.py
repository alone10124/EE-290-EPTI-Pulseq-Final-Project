"""Hardware checks for a built sequence: timing/rasters, gradient limits, forbidden echo spacings, GE PNS.

These replace the Siemens-only checks at the end of the MGH scripts (``gradSpectrum_bb_v1`` against Prisma acoustic
resonances, and ``check_pns`` with the Siemens SAFE model). GE gives forbidden *echo spacings* directly (the scanner's
``epiesp*.dat``), so the ESP is checked against those ranges instead of computing a gradient spectrum.
"""

from dataclasses import dataclass, field

import numpy as np
from pypulseq.supported_labels_rf_use import get_supported_labels

from epti.hardware import HardwareProfile
from epti.seq.ge_pns import pns as ge_pns


@dataclass
class Check:
    name: str
    ok: bool
    detail: str


@dataclass
class Report:
    title: str
    checks: list = field(default_factory=list)

    @property
    def ok(self):
        return all(c.ok for c in self.checks)

    def add(self, name, ok, detail):
        self.checks.append(Check(name, bool(ok), detail))

    def __str__(self):
        lines = [f"== {self.title}: {'PASS' if self.ok else 'FAIL'}"]
        lines += [f"  [{'ok' if c.ok else '!!'}] {c.name}: {c.detail}" for c in self.checks]
        return "\n".join(lines)


def gradient_waveforms(seq, dt, time_range=None):
    """Gradients resampled on a uniform ``dt`` grid. Returns ``(t, g)`` with ``g`` of shape (3, n) in T/m."""
    wave_data = seq.waveforms_and_times(time_range=time_range)[0]
    t_end = time_range[1] if time_range else seq.duration()[0]
    t0 = time_range[0] if time_range else 0.0
    t = np.arange(t0, t_end, dt)
    g = np.zeros((3, len(t)))
    for axis, wd in enumerate(wave_data[:3]):
        if wd.shape[1]:
            g[axis] = np.interp(t, wd[0], wd[1], left=0.0, right=0.0)
    return t, g / seq.system.gamma


def check_sequence(seq, hw: HardwareProfile, esp, title, pns_window=2.0):
    """Run all checks. ``esp`` is the readout echo spacing (s); PNS is evaluated over the first ``pns_window`` s."""
    rep = Report(title)

    ok, errors = seq.check_timing()
    detail = "no errors" if ok else f"{len(errors)} errors, first: {errors[0]}"
    rep.add("timing and rasters (PyPulseq check_timing)", ok, detail)

    adc_dwells = {round(seq.adc_library.data[k][1], 12) for k in seq.adc_library.data}
    bad = [d for d in adc_dwells if abs(d / hw.adc_raster_time - round(d / hw.adc_raster_time)) > 1e-6]
    rep.add(f"ADC dwell on {hw.adc_raster_time*1e6:g} us raster", not bad,
            f"dwells {sorted(d*1e6 for d in adc_dwells)} us")

    window = (0.0, min(pns_window, seq.duration()[0]))
    t, g = gradient_waveforms(seq, hw.grad_raster_time, window)
    gmax = np.abs(g).max(axis=1) * 1e3
    smax = np.abs(np.diff(g, axis=1)).max(axis=1) / hw.grad_raster_time
    rep.add(f"max gradient <= {hw.max_grad:g} mT/m", np.all(gmax <= hw.max_grad * 1.001),
            "x/y/z = " + "/".join(f"{v:.1f}" for v in gmax) + " mT/m")
    rep.add(f"max slew <= {hw.max_slew:g} T/m/s", np.all(smax <= hw.max_slew * 1.01),
            "x/y/z = " + "/".join(f"{v:.0f}" for v in smax) + " T/m/s")

    hits = [r for r in hw.forbidden_esp if r[0] <= esp <= r[1]]
    rep.add("echo spacing outside forbidden ranges", not hits,
            f"ESP {esp*1e6:.1f} us; forbidden {[(lo*1e6, hi*1e6) for lo, hi in hw.forbidden_esp] or 'none listed'}")

    if hw.pns_coil:
        p_thresh, _ = ge_pns(g, hw.pns_coil, gdt=hw.grad_raster_time)
        peak = p_thresh.max()
        rep.add(f"GE PNS ({hw.pns_coil}) < 80% (normal mode)", peak < 80,
                f"peak {peak:.1f}% over the first {window[1]:.1f} s")

    if hw.write_trid:
        # PyPulseq stores a label as (value, 1 + index into get_supported_labels())
        trid_id = get_supported_labels().index("TRID") + 1
        n_trid = sum(1 for v in seq.label_set_library.data.values() if v[1] == trid_id)
        rep.add("TRID labels present", n_trid > 0, f"{n_trid} distinct TRID values (GE pge2 uses them to find segments)")
    return rep
