"""Scanner hardware profile: gradient/RF limits, raster times, receive channels, and the matching KomaMRI scanner.

The values for a particular scanner live in the ``[hardware]`` table of a run's config file (``configs/<run>.toml``),
not in this package, so the same code can target different scanners. Build a profile with
``HardwareProfile.from_dict(config["hardware"])`` (``mrtools.config.load`` does this for you).

Not modeled for now: GE's forbidden echo spacings and peripheral nerve stimulation (PNS). The ``epti`` branch has
both (``epti/seq/checks.py``, ``epti/seq/ge_pns.py``).
"""

from dataclasses import dataclass, fields

import pypulseq as pp


@dataclass(frozen=True)
class HardwareProfile:
    name: str
    max_grad: float                  # mT/m
    max_slew: float                  # T/m/s
    rf_dead_time: float              # s
    rf_ringdown_time: float          # s
    adc_dead_time: float             # s
    grad_raster_time: float          # s
    rf_raster_time: float            # s
    adc_raster_time: float           # s
    block_duration_raster: float     # s
    B0: float                        # T
    ncoils: int                      # receive channels (the simulator's coil array uses this many)
    write_trid: bool = True          # TRID labels marking segment starts (required by GE pge2)

    @classmethod
    def from_dict(cls, d):
        """Profile from a config's ``[hardware]`` table. Missing or unknown keys raise an error, so typos don't pass
        silently."""
        names = {f.name for f in fields(cls)}
        unknown = set(d) - names
        if unknown:
            raise ValueError(f"unknown [hardware] keys: {sorted(unknown)}; allowed: {sorted(names)}")
        return cls(**d)

    def opts(self) -> pp.Opts:
        """PyPulseq system limits."""
        return pp.Opts(
            max_grad=self.max_grad, grad_unit="mT/m",
            max_slew=self.max_slew, slew_unit="T/m/s",
            rf_ringdown_time=self.rf_ringdown_time, rf_dead_time=self.rf_dead_time,
            adc_dead_time=self.adc_dead_time, adc_raster_time=self.adc_raster_time,
            grad_raster_time=self.grad_raster_time, block_duration_raster=self.block_duration_raster,
            rf_raster_time=self.rf_raster_time, B0=self.B0,
        )

    def make_scanner(self, coil_radius=0.15, coil_half_length=0.12):
        """KomaMRI ``Scanner`` with this profile's limits and an ``ncoils``-channel receive array.

        KomaMRI's ``BirdcageCoilSens`` gives one receive channel per rung (so it's an array, despite the name): ideal
        wires on a cylinder of ``coil_radius`` (m) and half-length ``coil_half_length`` (m). It's a stand-in for the
        real head coil, so its geometry is a simulation choice rather than part of the hardware description. Its
        sensitivity is symmetric in z, so two slices at +-z see the same coil maps. Starts Julia on first use.
        """
        import mrtools.sim  # noqa: F401  (sets the juliacall environment variables before Julia starts)
        from juliacall import Main as jl

        jl.seval("using KomaMRI")
        limits = jl.seval("(B0, Gmax, Smax, adc, dur, gr, rf, rd, rfd, adcd) -> HardwareLimits(B0=B0, Gmax=Gmax, "
                          "Smax=Smax, ADC_Δt=adc, DUR_Δt=dur, GR_Δt=gr, RF_Δt=rf, RF_ring_down_time=rd, "
                          "RF_dead_time=rfd, ADC_dead_time=adcd)")(
            self.B0, self.max_grad * 1e-3, self.max_slew, self.adc_raster_time, self.block_duration_raster,
            self.grad_raster_time, self.rf_raster_time, self.rf_ringdown_time, self.rf_dead_time, self.adc_dead_time)
        receiver = jl.seval("(n, r, L) -> BirdcageCoilSens(ncoils=n, radius=r, L=L)")(
            self.ncoils, coil_radius, coil_half_length)
        return jl.seval("(lim, rx) -> Scanner(limits=lim, receiver=rx)")(limits, receiver)
