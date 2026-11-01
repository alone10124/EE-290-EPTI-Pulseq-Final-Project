"""Scanner hardware profiles: gradient/RF limits, raster times, and interpreter-specific options.

Two built-in profiles:

* ``siemens_prisma_ref`` -- exactly the settings in ``legacy/sequences/seq_EPTI_*_P1.m``. Used only to check that the
  Python port reproduces the MATLAB-generated ``.seq`` files (parity tests). Not for the GE scanner.
* ``ge_mr750w`` -- GE rasters and dead times for the pge2 Pulseq interpreter. Gradient limits are the MR750w/XRMW
  values from TOPPE's PNS table; confirm them with the scanner facility (see docs/EPTI_PORT_PLAN.md, Phase 5).
"""

from dataclasses import dataclass, replace

import pypulseq as pp


@dataclass(frozen=True)
class HardwareProfile:
    name: str
    max_grad: float                  # mT/m
    max_slew: float                  # T/m/s, "fast" limit for readout lobes and blips (`sys` in the MATLAB scripts)
    max_slew_slow: float             # T/m/s, used for RF, spoilers and prephasers to lower PNS (`sys2` in the scripts)
    rf_dead_time: float              # s
    rf_ringdown_time: float          # s
    adc_dead_time: float             # s
    grad_raster_time: float          # s
    rf_raster_time: float            # s
    adc_raster_time: float           # s
    block_duration_raster: float     # s
    B0: float                        # T
    adc_samples_divisor: int = 4     # ADC sample count is rounded to a multiple of this
    fatsat_dwell: float | None = None   # RF dwell for the Gaussian fat-sat pulse; None -> rf_raster_time
    forbidden_esp: tuple = ()        # ((lo, hi), ...) forbidden echo-spacing ranges in s
    pns_coil: str | None = None      # GE gradient coil name for the PNS model (epti.seq.ge_pns)
    write_counter_labels: bool = True   # LIN/SLC/REP/AVG labels in zero-duration blocks, as the MATLAB scripts do
    write_trid: bool = False         # TRID labels marking segment starts (required by GE pge2)
    trigger_out: bool = True         # 'ext1' digital output pulse at the start of each volume (Siemens)
    matlab_compat: bool = False      # reproduce MATLAB quirks bit-for-bit (see epti.seq.rf_sms and the builders)

    def opts(self, slow: bool = False) -> pp.Opts:
        """PyPulseq system limits. ``slow=True`` gives the reduced-slew variant (`sys2` in the MATLAB scripts)."""
        return pp.Opts(
            max_grad=self.max_grad, grad_unit="mT/m",
            max_slew=self.max_slew_slow if slow else self.max_slew, slew_unit="T/m/s",
            rf_ringdown_time=self.rf_ringdown_time, rf_dead_time=self.rf_dead_time,
            adc_dead_time=self.adc_dead_time, adc_raster_time=self.adc_raster_time,
            grad_raster_time=self.grad_raster_time, block_duration_raster=self.block_duration_raster,
            rf_raster_time=self.rf_raster_time, B0=self.B0,
        )

    def with_overrides(self, **kwargs) -> "HardwareProfile":
        if "forbidden_esp" in kwargs:
            kwargs["forbidden_esp"] = tuple(tuple(r) for r in kwargs["forbidden_esp"])
        return replace(self, **kwargs)


PROFILES = {
    # legacy/sequences/seq_EPTI_imaging_P1.m L52, L57-68 (identical in the calib script)
    "siemens_prisma_ref": HardwareProfile(
        name="siemens_prisma_ref",
        max_grad=64.0, max_slew=200.0, max_slew_slow=160.0,
        rf_dead_time=100e-6, rf_ringdown_time=60e-6, adc_dead_time=20e-6,
        grad_raster_time=10e-6, rf_raster_time=1e-6, adc_raster_time=2.5e-6, block_duration_raster=10e-6,
        B0=2.89, adc_samples_divisor=4, fatsat_dwell=10e-6,
        write_counter_labels=True, write_trid=False, trigger_out=True, matlab_compat=True,
    ),
    # GE rasters: grad 4 us, RF 2 us, ADC 2 us (HarmonizedMRI SequenceExamples-GE/pge2 README).
    # Dead/ringdown times: the original values, which satisfy TOPPE's documented GE minimums of rfDead >= 72 us,
    # ringdown >= 54 us, adcDead >= 40 us (legacy/libs/toppe-main/+toppe/systemspecs.m L31-33).
    # Forbidden ESP band: TOPPE's default (systemspecs.m L24); the real list is in the scanner's epiesp*.dat.
    # max_grad / max_slew: the MR750w XRMW coil limits listed in TOPPE's PNS model table
    # (legacy/libs/toppe-main/+toppe/pns.m header: gmax 33 mT/m, smax 120 T/m/s). The "slow" slew keeps the original
    # 160/200 ratio. Still confirm with the scanner facility (derating, gradient mode).
    "ge_mr750w": HardwareProfile(
        name="ge_mr750w",
        max_grad=33.0, max_slew=120.0, max_slew_slow=96.0,
        rf_dead_time=100e-6, rf_ringdown_time=60e-6, adc_dead_time=40e-6,
        grad_raster_time=4e-6, rf_raster_time=2e-6, adc_raster_time=2e-6, block_duration_raster=4e-6,
        B0=3.0, adc_samples_divisor=4, fatsat_dwell=None,
        forbidden_esp=((410e-6, 510e-6),), pns_coil="xrmw",
        write_counter_labels=False, write_trid=True, trigger_out=False, matlab_compat=False,
    ),
}


def get_profile(spec) -> HardwareProfile:
    """Look up a profile by name, or by a mapping ``{"profile": name, <field>: <override>, ...}`` (from YAML)."""
    if isinstance(spec, HardwareProfile):
        return spec
    if isinstance(spec, str):
        return PROFILES[spec]
    spec = dict(spec)
    base = PROFILES[spec.pop("profile")]
    return base.with_overrides(**spec)
