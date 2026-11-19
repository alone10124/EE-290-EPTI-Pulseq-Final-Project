# GE 3T session: Thu 2026-11-19, 5-7 PM (with Hannah)

The goal of the session is to load our `.seq` files on the scanner, get them to run, and bring back raw data to
check what it looks like. This page lists what we checked beforehand, what to bring, and what to do if something fails.

## Already checked locally (2026-11-19)

Both GE conversion paths were run on our files with the real converters, using MATLAB R2023b. The scripts are in
[`tools/ge/`](../tools/ge/README.md).

| Check | pge2 (`.pge`) | tv6 TOPPE (`.tar`) |
|---|---|---|
| `seq2ceq` conversion (TRID segments, same block durations in every segment instance) | pass | pass |
| GE block timing (4 µs gradient, 2 µs RF/ADC raster, dead times) | pass (`pge2.check`) | pass (`ceq2ge` / `writemod`) |
| Peak B1 / gradient / slew | 0.076 G, 2.32 G/cm, 10.0 G/cm/ms | pass |
| GE PNS per segment, under 80% | pass for XRMW (MR750w), XRM (MR750), HRMW (Premier), HRMB (UHP) | not computed by tv6 |
| Converted waveforms match the `.seq` | pass for every segment (`pge2.validate`) | n/a |
| SAR / gradient heating (`toppe.preflightcheck`) | checked on the scanner | pass. Peak 10 s SAR: 1.0 W/kg (calib), 0.1-0.2 W/kg (imaging) |

**Bug found and fixed:** both sequences started with an 8 µs delay block that had its own TRID. That made a segment
with no waveforms, and `pge2.check` crashed on it ("g must be size 3xn"). The scanner would likely have rejected it too.
The block is now left out for GE (`epti/seq/calib.py`, `epti/seq/imaging.py`). All 37 tests pass.

## Files to bring (USB stick + laptop), all in `sequences/`

Bring the laptop too, so we can regenerate or reconvert during the session (MATLAB + `tools/ge/` + `.venv`).

| Scan | `.seq` (source) | pge2 | tv6 | Duration |
|---|---|---|---|---|
| Calibration (42 slices, 29 PE lines, no blips) | `epti_3mm_ge_calib.seq` | `epti_3mm_ge_calib.pge` | `epti_3mm_ge_calib_tv6.tar` | 34.8 s |
| EPTI imaging, 1 dynamic | `epti_3mm_ge_imaging.seq` | `epti_3mm_ge_imaging.pge` | `epti_3mm_ge_imaging_tv6.tar` | 2.0 s |
| EPTI imaging, 3 dynamics (preferred) | `epti_3mm_ge_3dyn_imaging.seq` | `epti_3mm_ge_3dyn_imaging.pge` | `epti_3mm_ge_3dyn_imaging_tv6.tar` | 6.0 s |
| Fallback: calibration, 2 slices | `sim_small_calib.seq` | `sim_small_calib.pge` | `sim_small_calib_tv6.tar` | 34.8 s |
| Fallback: imaging, 1 slice group (2 slices), 2 dynamics | `sim_small_imaging.seq` | `sim_small_imaging.pge` | `sim_small_imaging_tv6.tar` | 4.0 s |

The 1-dynamic imaging scan is all first-TR data, before steady state. Use the 3-dynamic one and discard dynamic 1.
The calibration and imaging scans must use the same readout. `make_epti.py` checks this, so don't change ESP or
blip duration in only one of them.

## Ask Hannah (ideally before Thursday, by email)

1. **pge2 or tv6?** Her email says "TOPPE package", which may mean tv6 (`.tar`). We have both. Which PulCeq/TOPPE
   versions does she use? Current TOPPE `main` writes the v7 entry file, which doesn't match PulCeq's tv6 branch
   (see `tools/ge/README.md`).
2. **`psd_rf_wait` and `psd_grd_wait`** for this scanner. They only change the conversion, not our `.seq`. We used 100 µs.
3. **Scanner model and gradient coil** (MR750w, MR750, Premier or UHP?). PNS passes on all four. If it isn't a 750w,
   we could relax `max_slew` later.
4. **Forbidden echo spacings** (`epiesp*.dat`). We use 600 µs, outside TOPPE's default 410-510 µs band.
5. **Head coil and channel count**, and the **raw data format** (P-file or ScanArchive) and how we get it off the
   scanner. Our recon reads ISMRMRD, so we will need a converter afterwards.

## Order of the session

1. Use a **phantom** (a human scan needs IRB approval). Use the spherical or head-shaped phantom with some structure.
2. GE localizer. Center the landmark on the middle of the phantom. Our slices cover ±63 mm in z around isocenter,
   in logical axial orientation.
3. **Shim well** (GE auto-shim, or higher order if available). EPTI's B0 maps make shimming errors visible, but bad
   shimming still distorts the images and loses signal.
4. Load and run **`epti_3mm_ge_calib` first**: no blips, simple Cartesian GRE-EPI, longest TR. If it runs, the
   conversion setup is right.
   Then run `epti_3mm_ge_3dyn_imaging` (or `epti_3mm_ge_imaging`).
5. Auto prescan (center frequency, transmit gain, receive gain) before each scan. The pge2 receive gain is set from the
   first ADC. That is the navigator of the first slice group, which has the strongest signal, so it shouldn't clip.
6. Copy the raw data for every scan that ran. Also save the scan.xml / Pulse View output if Hannah can produce it:
   `pge2.validate(ceq, sysGE, seq, xmlPath)` compares the scanner's waveforms with our `.seq`.

## If a scan won't run

Isolate the problem one step at a time:
- **Calibration fails too:** the problem is the setup (interpreter version, file path, entry file), not EPTI. Try one
  of Hannah's known-good example sequences, then the `sim_small_*` files.
- **Calibration runs, imaging doesn't:** try `sim_small_imaging` (1 slice group, smaller loop). If that runs, the
  problem is the scan length or loop size, not the waveforms.
- **Timing or PNS error on the console:** rerun `tools/ge/ge_dryrun.m` with the scanner's real limits
  (`'g_max', 'slew_max', 'coil'`). To lower PNS, reduce `max_slew` in `config/epti_3mm_ge.yaml`, then
  `python sequences/scripts/make_epti.py config/epti_3mm_ge.yaml` and reconvert. A full rebuild and conversion takes
  about 5 minutes.
- **The scan runs but slices are swapped, or the images ghost badly:** don't debug that at the scanner. Save the data:
  this is the `ReconOptions.adc_phase_sign` item in `docs/EPTI_PORT_PLAN.md`.
- If the tv6 conversion is done on Windows, it needs the line-ending patch in `tools/ge/tv6_patch/` (CRLF breaks the
  entry file).

## After the session

- Raw data loader: P-file / ScanArchive → ISMRMRD, or a direct reader. ADC order is the same as in the `.seq`
  (dynamic, slice group, lobe for imaging).
- First recon: check `adc_phase_sign` and the ghost-correction sign.
- Write the confirmed hardware values back into `config/epti_3mm_ge.yaml` and the `ge_mr750w` profile.
