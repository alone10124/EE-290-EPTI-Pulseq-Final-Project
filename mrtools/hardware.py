"""Scanner hardware profiles: gradient/RF limits, raster times, and Pulseq-interpreter options.

The ``ge_mr750w`` profile holds the GE rasters and dead times for the pge2 Pulseq interpreter, and MR750w gradient
limits. Where the values come from:

* Gradient limits: GE's published MR750w specs, 44 mT/m and 200 T/m/s per axis. ``max_slew_slow`` (120 T/m/s) is used
  for everything except the EPI readout and blips. Confirm with the scanner facility which gradient mode is used.
* Rasters: grad 4 us, RF 2 us, ADC 2 us (HarmonizedMRI SequenceExamples-GE/pge2 README).
* Dead/ringdown times: rfDead >= 72 us, ringdown >= 54 us, adcDead >= 40 us are GE's minimums as documented by TOPPE
  (``legacy/libs/toppe-main/+toppe/systemspecs.m`` L31-33 on the ``epti`` branch); we use slightly larger values.
* Not modeled for now (see docs/LESSONS_LEARNED.md): forbidden echo spacings and peripheral nerve stimulation (PNS).
"""

from dataclasses import dataclass

import pypulseq as pp


@dataclass(frozen=True)
class HardwareProfile:
    name: str
    max_grad: float                  # mT/m
    max_slew: float                  # T/m/s, used for readout lobes and blips
    max_slew_slow: float             # T/m/s, a gentler limit for RF, spoilers and prephasers (lowers nerve stimulation)
    rf_dead_time: float              # s
    rf_ringdown_time: float          # s
    adc_dead_time: float             # s
    grad_raster_time: float          # s
    rf_raster_time: float            # s
    adc_raster_time: float           # s
    block_duration_raster: float     # s
    B0: float                        # T
    write_trid: bool = True          # TRID labels marking segment starts (required by GE pge2)

    def opts(self, slow: bool = False) -> pp.Opts:
        """PyPulseq system limits; ``slow=True`` uses ``max_slew_slow``."""
        return pp.Opts(
            max_grad=self.max_grad, grad_unit="mT/m",
            max_slew=self.max_slew_slow if slow else self.max_slew, slew_unit="T/m/s",
            rf_ringdown_time=self.rf_ringdown_time, rf_dead_time=self.rf_dead_time,
            adc_dead_time=self.adc_dead_time, adc_raster_time=self.adc_raster_time,
            grad_raster_time=self.grad_raster_time, block_duration_raster=self.block_duration_raster,
            rf_raster_time=self.rf_raster_time, B0=self.B0,
        )


GE_MR750W = HardwareProfile(
    name="ge_mr750w",
    max_grad=44.0, max_slew=200.0, max_slew_slow=120.0,
    rf_dead_time=100e-6, rf_ringdown_time=60e-6, adc_dead_time=40e-6,
    grad_raster_time=4e-6, rf_raster_time=2e-6, adc_raster_time=2e-6, block_duration_raster=4e-6,
    B0=3.0, write_trid=True,
)
