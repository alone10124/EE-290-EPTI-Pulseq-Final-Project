# Lessons learned from the EPTI port

We first ported MGH's EPTI (echo planar time-resolved imaging) MATLAB code to Python. Then we switched to a simpler
EPI-based method of our own. We start from a **simple single-shot, multi-slice EPI** to build and
debug the workflow. The full EPTI work (PyPulseq sequences matched to MATLAB, KomaMRI
simulation, a subspace reconstruction that recovered T2* to about 4.5%, and the MGH sources in `legacy/`) is on the
**`epti` branch**. This note keeps what carries over.

## GE MR750w hardware (in `mrtools/hardware.py`, `mrtools/seq/checks.py`)
- **Gradient limits:** GE publishes 44 mT/m and 200 T/m/s per axis for the MR750w. The EPTI port had used 33 mT/m
  and 120 T/m/s from a table in TOPPE's PNS code (`pns.m`), which doesn't match GE's specs. We use 200 T/m/s for the
  EPI readout and blips and 120 T/m/s for everything else.
- **Rasters for GE's pge2 Pulseq interpreter:** gradients 4 µs, RF 2 µs, ADC 2 µs, blocks 4 µs. Every echo spacing
  (ESP), RF duration and delay has to land on these. MGH's Siemens values (10 µs gradients, 2.5 µs ADC) didn't.
- **TRID labels** must mark the first block of each repeating segment (`mrtools.seq.common.trid`).
- **Forbidden echo spacings:** GE scanners have them (readout frequencies that hit mechanical resonances of the
  gradient coil), but they aren't in the datasheets; they're on the scanner in `/usr/g/bin/epiesp.dat`. TOPPE guesses
  410–510 µs. We don't check them for now (see the table below); with a 128 matrix and flat-top-only sampling the
  ESP is about 830 µs, well above that band anyway.
- **Confirm with the scanner facility:** the gradient mode and limits, which interpreter (pge2 or tv6) and software
  version, `epiesp.dat`, the receive channel count, the raw-data format and how to get it into Python, and whether
  trigger output works.

## PyPulseq pitfalls
- `pp.make_trapezoid(area=..., duration=...)` picks longer ramps (a lower plateau) than MATLAB's `mr.makeTrapezoid`.
  That only matters when matching MATLAB; it's fine for new designs.
- A zero-duration block (e.g. one holding only a label) between two gradients that are non-zero at the boundary
  fails PyPulseq's continuity check; MATLAB skips such blocks. Put labels on a real block.
- `pp.make_label(label, type, value)` has the opposite argument order from MATLAB's `mr.makeLabel(type, label, value)`.

## KomaMRI simulation (in `sim/`)
- Set `PYTHON_JULIACALL_THREADS=auto` **and** `PYTHON_JULIACALL_HANDLE_SIGNALS=yes` before importing komamripy; with
  several threads and no signal handling it core-dumps. `sim/phantom.py` sets both.
- **The simulator ignores the phantom's `T2s` field** (it uses T2 only). For gradient-echo decay, set `T2 := T2*`.
  `sim/phantom.py` does this.
- **Sign conventions:**
  - Off-resonance +Δf precesses as `exp(−i2πΔf t)`, so a B0 map fitted from echo phase comes out as `−Δf`.
  - Gradient phase is `exp(−i2π k·r)`.
  - The ADC phase offset is applied as `exp(−iφ)`.
- **SMS:** with blipped-CAIPI, getting the ADC-phase sign wrong **swaps the simultaneous slices** in the recon. This
  happened in EPTI and was caught by comparing with the phantom. Keep a phantom with distinct slices (`sim/phantom.py`
  rotates every other slice) and check against its ground truth. On the scanner, check the sign on the first data set.
- The `BirdcageCoilSens` receive array is symmetric in z, so slices at ±z see identical coil maps. Unfolding then
  relies on the CAIPI shift. Use `ArbitraryCoilSens` for more realistic coils.
- KomaMRI's ISMRMRD header has no `receiverChannels` or `patientPosition`, and the strict `ismrmrd` parser rejects it.
  Our simulation stage now saves plain NumPy arrays instead.

## Calibration and coil maps
- GRAPPA only fills k-space lines that were *skipped*. A calibration scan that's **fully sampled at low or imaging
  resolution** needs only ESPIRiT for coil maps (`mrtools.recon.coils`).
- EPTI's calibration was larger and undersampled (needing GRAPPA) because it also measured a high-resolution B0 map
  and odd/even phase maps for its forward model. We don't model B0 voxel by voxel, so we don't need that.
- Plan: put the calibration **in the same .seq file** as single-band EPI with the **same readout** as the imaging. The
  coil maps then share the EPI distortion, which helps any coil-based unfolding (SENSE, SMS). Tag the calibration acquisitions (e.g. with the
  `REF` label) or count them, and expect the first dynamics to start out of steady state.

## Simplifying assumptions for now, and where to get each piece back
Everything below exists, tested, on the `epti` branch.

| Assumption now | To add back | On the `epti` branch |
|---|---|---|
| ADC samples only on the readout flat top | Ramp sampling + regridding (shorter echo spacing) | `epti/recon/regrid.py`, readout in `epti/seq/common.py` |
| Odd/even (Nyquist) ghosts ignored | Navigator-based linear phase correction | `epti/recon/ghost.py` (matches MGH's P-code exactly) |
| No B0 drift between dynamics | Navigator-based drift estimate | `epti/recon/drift.py` |
| Calibration fully sampled | GRAPPA for undersampled calibration | `epti/recon/calib.py` (`grappa_fill`) |
| All coils used | Geometric coil compression for 32+ channels | MGH/Lustig code in `legacy/` |
| No distortion or B0 modeling | B0 map from multi-echo calibration | `epti/recon/calib.py` (`fit_b0`) |
| Simulation output as NumPy | ISMRMRD reading, NIfTI writing | `epti/io/` |
| No noise in simulation | Complex Gaussian noise | `read_acquisitions(noise_std=...)` in `epti/io/ismrmrd_io.py` |
| No forbidden echo-spacing check | Check the ESP (most common gap between ADCs within an EPI train) against the scanner's `epiesp.dat` | `epti/seq/checks.py` (`forbidden_esp` in `epti/hardware.py`) |
| No multi-echo post-processing | T2* fit, BOLD-optimal echo combination, k-space apodization | `epti/recon/postproc.py` |
| No peripheral nerve stimulation (PNS) check | GE's PNS model; see the note below | `epti/seq/ge_pns.py` (TOPPE `pns.m` port), `tests/test_ge_pns.py` |
| Full ky sampling, single shot | In-plane acceleration (SENSE/GRAPPA) and/or partial Fourier to shorten the readout | GRAPPA in `epti/recon/calib.py`; SENSE via SigPy |

**PNS, when it comes back:** GE's scanner refuses (or needs "first controlled mode" for) sequences above 80% of its
PNS threshold. With GE's model for the MR750w coil (`xrmw`), EPI readouts at 120 T/m/s already reached 80–92%;
derating the readout slew to about 90–100 T/m/s brought them under 80%. At 200 T/m/s it will be worse, so expect to
derate the readout slew before scanning.

**Readout length:** with a 128 × 128 matrix, 220 mm FOV and flat-top-only sampling, the EPI train is about 106 ms and
TE about 56 ms, longer than gray/white-matter T2* (about 60–70 ms). That blurs the image along the phase-encode
direction and lengthens TE beyond the usual ~30 ms for BOLD at 3 T. Ramp sampling, in-plane acceleration (R = 2–3)
and partial Fourier are the usual ways to shorten it.

**One bug to remember if regridding comes back:** MGH's regridding (inherited from Gadgetron) evaluates the image on
`linspace(-0.5, (N-1)/(2N), N)`, which stretches readouts by 0.7% and misinterpolates by up to 8%. The `epti` branch
version uses the uniform grid.
