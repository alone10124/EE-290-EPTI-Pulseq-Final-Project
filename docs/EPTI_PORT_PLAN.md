# Plan: Port EPTI (MATLAB) into the Python PyPulseq → KomaMRI → recon workflow

## Status (2026-10-31)
| Phase | State |
|---|---|
| 0. Cleanup & environment | **Done.** `legacy/` pruned to the files actually used (386 kept, 545 removed; see `legacy/README.md`). New dependencies installed and verified on Python 3.13 / numpy 2.5. `pyproject.toml`, `.gitignore` and the README are updated. |
| 1. Sequences in PyPulseq | **Done.** `epti/seq/` builds both scans, and `sequences/scripts/make_epti.py` writes and checks them. **26/26 tests pass, including MATLAB parity**: with the `siemens_prisma_ref` profile, every gradient sample, ADC time, RF event, definition and per-ADC label value matches the MGH-generated `tests/ref/*.seq`. |
| 2. Simulation | **Done.** `sim/run_sim.py` simulates both scans (52k spins, 8 coils, two slices, synthetic B0) in about 1 min on CPU, and writes ISMRMRD + ground truth to `out/sim/`. `sim/preview_calib.py` makes a plain-FFT sanity figure. Noise isn't added yet; plan to add it when the recon reads the data. |
| 3. Core recon | **Done.** `epti/recon/` + `recon/run_recon.py`. On the simulation: T2* median error 4.4–4.6%, calibration B0 0.16–0.18 Hz, SMS slices correctly separated; about 2.5 min per dynamic on CPU (70 ADMM iterations). 37/37 tests pass. |

## Remaining work (as of 2026-11-01)
Suggested order: 1 → 2 → 3, and send the Phase 5 questions to the scanner facility now (the answers take time).

1. **Leftovers from Phases 2–3**
   - [ ] Add realistic noise to the simulation and pick an SNR (`recon/run_recon.py --noise`, or add it in `sim/run_sim.py`).
   - [ ] Tune the recon on noisy data (LLR `--lamda`, `--iters`). The T2* map currently shows horizontal striping.
   - [ ] GPU check: run the recon on an NVIDIA machine (`pip install cupy-cuda12x`, `--device 0`), compare it with the CPU result, and time both.
   - [ ] Commit, including `legacy/` and the MATLAB references in `tests/ref/`.
2. **Phase 4: extras**
   - [ ] KomaMRI motion experiments (`Translate` / `Rotate` across dynamics) in `sim/`.
   - [ ] High-order B0 update (`B0_update = 1`: `recon_lowresI2`, `DriftB0_esti_nav_WJG3_4`), with AFNI registration replaced by Python.
   - [ ] Motion-aware recon (`moco_recon`), which builds on the high-order B0 update.
   - [ ] GCC coil compression (needed for 32+ channel head coils, not for 8 simulated coils).
   - [ ] Optional: build the EPTI BART fork to cross-check `epti/recon/solve.py`.
3. **Phase 5: GE scanner**
   - [ ] Questions for the scanner facility:
     - Which interpreter (pge2 or tv6), and which software version?
     - What are the real Gmax/Smax and gradient coil? We currently use 33 mT/m and 120 T/m/s from TOPPE's table.
     - Can we get `epiesp*.dat` (the forbidden echo spacings)?
     - How many receive channels?
     - Is trigger output supported?
     - What is the raw format (ScanArchive or P-file), and how is it converted to ISMRMRD?
   - [ ] Update `config/epti_3mm_ge.yaml` and the `ge_mr750w` profile with the answers, then rerun `make_epti.py`'s checks.
   - [ ] Add a raw-data loader or converter if the scanner doesn't give ISMRMRD.
   - [ ] On the first scanner data: check `ReconOptions.adc_phase_sign`. If the two SMS slices come out swapped, flip it.

Findings from Phase 1 that change details below:
- **Two porting pitfalls found by the parity tests (both fixed):**
  1. PyPulseq 1.5.0's `make_trapezoid(area=..., duration=...)` uses different ramp times from MATLAB's `makeTrapezoid`.
     MATLAB uses the lowest amplitude with the fastest ramp; PyPulseq reuses the ramps of the shortest gradient for that area.
     Our `epti.seq.common.trapezoid` reproduces MATLAB's choice. This affected the readout lobe itself (ramps of 210 µs vs 120 µs).
  2. MATLAB's even-length `firls` filter is exactly symmetric, so `getsmspulse` finds two tied peaks and uses their midpoint as its
     time reference. SigPy's filter is symmetric only to about 1e-3, so a literal `argmax` shifted the multiband modulation by 10 µs.
     We now take the symmetry center from the filter layout.
- **GE protocol (`config/epti_3mm_ge.yaml`)**: ESP 600 µs, readout/blip slew derated to 100 T/m/s, blip 168 µs, TR 2.0 s, calibration RF 2.048 ms.
  TE₁ becomes 41.5 ms (the original is 34.4 ms).
  - The binding constraint is GE PNS, not timing: 120 T/m/s gives 80.3–80.5%, just over the 80% normal-mode limit; 100 T/m/s gives 78.3–78.5%.
  - The calibration and imaging scans must share the readout (ESP, blip_dur), because the ghost correction is measured on one and applied
    to the other. `make_epti.py` checks this.
- **GE limits**: 33 mT/m and 120 T/m/s come from the MR750w/XRMW row of TOPPE's PNS table (`+toppe/pns.m`). Still confirm on site.
- **Counter labels** ride on the next real block instead of a zero-duration block, because PyPulseq 1.5 can't reproduce MATLAB's
  label-only blocks between blips. The per-ADC label values are unchanged, and the parity test checks this.
- **Phase 2 findings, for the recon:**
  - KomaMRI's simulator ignores the phantom's `T2s` field, so the phantom sets T2 := T2* to get gradient-echo decay.
  - Koma's ISMRMRD header has no `patientPosition` or `receiverChannels`, so the strict `ismrmrd` XML parser fails.
    The reader must take sizes from the acquisitions, which the original MATLAB reader effectively does too.
  - The B0 map estimated from the calibration echo phase has the right shape and size but the **opposite sign** from the
    phantom's `Δw`. The Phase 3 B0 code must fix the sign convention against `out/sim/*_truth.npz` (add a test).
  - The default `BirdcageCoilSens` is symmetric in z, so both SMS slices see the same coil maps; separating them relies on the
    CAIPI FOV/2 shift. Consider `ArbitraryCoilSens` later.
- **Phase 3 findings:**
  - Ghost correction matches the encrypted P-code exactly. The P-code uses `mrir_iDFT_freqencode = N * fftshift(ifft(ifftshift(x)))`, and the
    correction is `exp(i(a + b(n - N/2)))`. Even ideal simulated data shows a 2.7 rad/FOV odd/even slope, from ADC sample-timing conventions,
    and the ghost correction removes it.
  - **CAIPI sign convention:** Pulseq/KomaMRI demodulate the ADC phase with `exp(-i phi)`, so the residual blip phase is
    `exp(-i 2 pi kz (z_s + z_ref))`. The upper slice flips sign, not the lower one. The wrong sign swapped the two SMS slices.
    `ReconOptions.adc_phase_sign` exists for scanner data: **verify it on the first GE data set.**
  - The subspace basis is built on times relative to the start of the echo train (`gen_GE_basis_T2B0(..., T0 = 0)`), not absolute TEs
    (found by the MATLAB fixture test).
  - **Bug fixed relative to the original:** the regridding's image grid `linspace(-0.5, (N-1)/(2N), N)` (from Gadgetron) stretches readouts
    by 0.7% and causes up to 8% interpolation error. We use the uniform grid.
  - Not ported yet: coil compression (not needed for 8 simulated coils), the high-order B0 update and the motion-aware recon (Phase 4).
    The sim currently adds no noise (`--noise`).
- **`tdr_sliceorder`** (encrypted): our Siemens-style interleave guess matches all 18 MATLAB cases.
- **KomaMRI**: with more than one Julia thread, `PYTHON_JULIACALL_HANDLE_SIGNALS=yes` is needed as well as `PYTHON_JULIACALL_THREADS=auto`,
  or it core-dumps.

## Context
The team has a working Python sim loop (PyPulseq `.seq` → komamripy → numpy FFT) for a toy EPI. It has just received
the official MGH EPTI code: MATLAB sequence scripts and recon, a patched BART, and bundled copies of Pulseq, TOPPE and ISMRMRD.

Goal: run EPTI end to end in Python, in simulation now and on a GE scanner later. Use the GPU where practical, drop
unneeded dependencies, and keep the `.seq` files GE-compatible. **Prefer existing libraries (SigPy, pygrappa, scipy) over
translating code**, especially the obfuscated P-code.

User decisions:
- **MATLAB**: used only for small runs that need no data (the two seq scripts plus one fixture script). We don't depend on MATLAB recon outputs.
- **GE conversion**: turning `.seq` into the GE interpreter's format (pge2) happens outside this repo. Our `.seq` files must still follow GE hardware rules.
- **Recon scope**: the core path first. The high-order B0 update and the motion-aware recon come later.
- **Legacy code**: the MATLAB files actually used go into `legacy/` as a reference. The rest is deleted.
- **GE raw-data I/O**: a later phase, because the GE interpreter and raw format are still unknown.

## Glossary (also goes in README)
- **EPI / EPTI**: EPI reads many k-space lines after one excitation by alternating the readout gradient. EPTI zig-zags
  through ky over time, so every echo time gets its own image, and a model-based recon fills the gaps.
- **PE (phase encoding)**: the ky direction. PE lines are the ky rows (npe = 72).
- **ESP (echo spacing)**: the time between EPI readout lines. GE forbids ESPs whose readout frequency hits mechanical resonances of the gradient coil.
- **SMS / MB (simultaneous multi-slice / multiband factor)**: one RF pulse excites MB slices at once. With MB=2 they are 63 mm apart.
  The slices are separated with coil sensitivities plus blipped-CAIPI kz blips, which shift one slice by FOV/2.
- **K**: the number of subspace basis vectors (4). The recon solves for K coefficient maps instead of ~125 echo images.
  Each echo image is Φ(t)·coefficients, where Φ comes from an SVD of T2 decay curves.
- **GRAPPA**: k-space parallel imaging (Griswold, MRM 2002). Missing ky lines are filled from neighbouring acquired lines
  across coils, with weights fit on the fully sampled center (ACS) lines. Here R=3 with 18 ACS lines.
- **Ghost (Nyquist) correction**: odd and even lines are read with opposite gradient polarity, which gives them a phase mismatch and an N/2 ghost.
  Three navigator lines (+,−,+) measure the mismatch as a function of x. It is fit as constant + linear and then removed.
- **Hybrid / x-space**: data after an inverse FFT along the readout only, i.e. (x, ky).
- **Even-line flip**: even lines are read right-to-left, so their samples are reversed.

## Findings that shape the design
- **The raw readers ignore labels and ISMRMRD idx fields.** `libs/utils/ReadRaw_EPTI_*_ISMRMRD.m` reshape purely by arrival order.
  - KomaMRI writes ISMRMRD in ADC order, so its output plugs straight in.
  - The NAV/REV/SEG `mr.makeLabel` calls are no-ops.
- **The ISMRMRD path reads data correctly, but its NIfTI step only works for Siemens.**
  - `nifti_header` is assigned only in `ReadRaw_EPTI_calib_pulseq_Siemens.m:57`.
  - `recon_EPTI_P1.m:75` and `saveniftidata_EPTI_pulseq.m:7` read it, so `ISMRMRD_data=1` should fail with "Unrecognized field name".
    This comes from reading the code, not running it. The default Siemens path is complete.
- **Recon parameters come from the `.seq` `[DEFINITIONS]`** (`EPTI_SMS_Ext_parameters_pulseq.m`).
  We will write an explicit `AdcDwell` instead of the `adcRasterTime` hack.
- **BART:**
  - The EPTI fork adds `picsS`, `linopScript` and `opScript` (`bart-0.8.00_EPTI/Makefile:207-209`).
  - It hard-codes LLR to skip SMS slices (`src/lowrank/lrthresh.c:516`) and defaults to CUDA with an MGH-specific CUDA path.
  - **Stock BART cannot do the EPTI solve.** `pics -B` asserts that the basis has no FFT/coil dims, but EPTI's basis
    Φ'(x,y,t) = Φ(t)·exp(i2πΔB0(x,y)·TE_t)·Phase0 is voxel-dependent.
  - BART's Python wrapper only shells out to a compiled `bart` binary through .cfl files.
  - No conda-forge package was found. The official installs are apt (`bart`, `bart-cuda`), MacPorts/Homebrew, or building from source.
  - Decision: **no BART**. SigPy covers `ecalib` (ESPIRiT), and `picsS` becomes SigPy linops + ADMM.
- **TOPPE is used only for the SMS RF design** (`getsmspulse.m` → Pauly `dzrf`). We replace it with `sigpy.mri.rf`.
- **P-code**: `libs/libs_EPTI_Recon/library_EPTI_rawprocessing/library_EPI/mrir_toolbox/*.p` (MGH mrir toolbox).
  It runs in MATLAB but can't be read. Each core-path P-code function is replaced by a library call where one exists
  (see the table). Ghost correction has no library equivalent and gets about 30 lines of our own code, checked against a MATLAB fixture.
- **Hidden dependencies**: AFNI (B0-update/moco path only) and three MATLAB toolboxes. Siemens readers are not needed.
- **Sequence gotchas**:
  - ESP = 490 µs is off GE's 4 µs raster and sits inside TOPPE's default GE forbidden band [410, 510] µs.
  - The raster times are Siemens values. There are no TRID labels, B0 is 2.89 T, and the PNS and forbidden-frequency checks are Prisma-only.
  - The `../libs` path is broken.
  - `getsmspulse.m` L114/L118 compares complex values.
- **PyPulseq API differences**: `make_label(label, type, value)`, `calc_rf_center()[0]`, `align(right=..., left=[...])`,
  and `mround` (MATLAB-style rounding).
- **KomaMRI v0.14** already supports labels/TRID, multi-coil receive, Δw, T2*, motion, ISMRMRD output and CPU threads/CUDA.
  The sim currently runs on one thread.
- **GPU**: SigPy uses CuPy, which runs on NVIDIA GPUs (e.g. the team's RTX 3050 laptop). The Radeon 780M laptop (gfx1103) runs SigPy on CPU,
  because neither CuPy-ROCm nor ROCm officially supports it. The problem is small (72×72, ~125 echoes, K=4, MB=2), so CPU is adequate.
  - **Why not PyTorch:** SigPy already has ESPIRiT, the RF tools, the linear operators and an ADMM solver. That gives one array
    stack (numpy/CuPy) and no hand-written solver. PyTorch's main advantage would be an unofficial ROCm path to the 780M.
  - **Fallback:** if CPU runtime becomes a bottleneck, rewrite only `model.py` + `solve.py` in PyTorch. The other stages stay the same.
- **License**: the EPTI recon library is BSD-style (© 2023 Wang & Dong). Lustig's coil-compression code has its own LICENSE.
  Ports keep those notices. The headers also require the MGH access request, so check its terms before pushing `legacy/` to a public repo.

## Library-first substitutions (what replaces each original piece)
| Original | Kind | Replacement |
|---|---|---|
| `bart picsS` + `linopScript` (fork) | C, CUDA | SigPy linops (`Multiply`, `Sum`, `FFT`, `Reshape`) + `sp.app.LinearLeastSquares(solver='ADMM', rho=1e-2, max_iter=70, max_cg_iter=8)` + **custom LLR `sp.prox.Prox`** (~30 lines; SigPy has none) |
| `bart ecalib -r24 -c0.1 -m1` | stock BART | `sigpy.mri.app.EspiritCalib` |
| `mrir_iDFT_/fDFT_freqencode` | P-code | `sigpy.ifft/fft(axes=(readout,))` (centered FFT) |
| `mrir_array_GRAPPA_2d_kernel_improved_dev` + `_conv_kernel_v2` | P-code | `pygrappa` (Tikhonov `lamda`, 3×3 kernel). Fallback: ESPIRiT on ACS + `sigpy.mri.app.SenseRecon` |
| `mrir_filter_raw_apodize_1d` | P-code | `scipy.signal.windows.tukey` (shape matched to fixture) |
| `mrir_artifact_ghost_compute/_correct` | P-code | **ours** (~30 lines, standard linear fit), checked against fixture |
| `tdr_sliceorder` | P-code | 2-line interleave, checked against fixture |
| `getsmspulse` / TOPPE `makeslr` | MATLAB | `sigpy.mri.rf.dzrf` + MB modulation |
| GCC coil compression (Lustig) | readable MATLAB | deferred: not needed for 8–16 sim coils; later port or simple SVD compression |
| `EPI_trapezoid_regridding_WJG` | readable MATLAB | port (sinc matrix + `np.linalg.pinv`) |
| Mask (`imbinarize`, `bwareaopen`, `imfill`) | MATLAB toolbox | `skimage.filters.threshold_otsu`, `remove_small_objects`, `scipy.ndimage.binary_fill_holes` |
| NIfTI (Jimmy Shen toolbox + Siemens vox2ras) | MATLAB | `nibabel`, affine built from `.seq` FOV/resolution/slice positions |
| `ismrmrd.Dataset` (MATLAB) | MATLAB | `ismrmrd` (pip) |
| `mr.Sequence.read` (recon side) | MATLAB Pulseq | `pp.Sequence().read()` |

## Conventions
- **Provenance**: every ported or replaced function's docstring cites `legacy/<path>:<lines>`, the original authors, and the paper
  or algorithm it implements. For library substitutions it also names the library call and how equivalence was checked.
  The original LICENSE files are copied into `legacy/`.
- **Array shapes**: named, low-dimensional arrays such as `[slc, echo, coil, ky, kx]`, never the 10-D column-major mrir layout.
- **Config classes**: `Protocol` (scan parameters from the top of the MATLAB scripts) and `HardwareProfile` (the `mr.opts` limits and
  rasters, plus the GE forbidden-ESP band and PNS coil). Both are new dataclasses in `epti/config.py`, loaded from YAML in `config/`.

## Target layout
```
pyproject.toml            # pip install -e .
epti/
  config.py hardware.py   # Protocol + HardwareProfile; profiles: ge_mr750w, siemens_prisma_ref
  seq/   rf_sms.py patterns.py blips.py calib.py imaging.py checks.py ge_pns.py
  io/    seqdefs.py ismrmrd_io.py nifti.py
  recon/ device.py regrid.py ghost.py calib.py subspace.py sampling.py drift.py
         model.py (SigPy linop builder) llr.py (LLR Prox) solve.py postproc.py pipeline.py
sequences/scripts/make_epti.py    # config → calib.seq + imaging.seq
sim/phantom.py sim/run_sim.py sim/sim_test.ipynb   # Koma → out/sim/{calib,imaging}.h5 (ISMRMRD)
recon/run_recon.py                # (2× .h5, 2× .seq, config) → out/recon/*.npz, *.nii.gz, figures
config/epti_3mm.yaml config/sim_small.yaml
tests/ (ref/*.mat, ref/*.seq from MATLAB)
legacy/ (+ make_reference_fixtures.m, written by us, run once by the team)
```

## MATLAB runs the team does once (no data, BART or AFNI needed; Signal Processing + Curve Fitting toolboxes)
1. The two seq scripts (path fixed) → `tests/ref/P1_Calib.seq`, `tests/ref/P1_EPTI_sms2.seq`.
2. `legacy/make_reference_fixtures.m` → `tests/ref/fixtures.mat` (`-v7`, readable with `scipy.io.loadmat`). It saves:
   - Seq helper outputs:
     - `calc_Pattern_WJG3(1,72,128,6,2,0)` and `calc_Pattern_WJG4(...)`
     - `calc_Pattern_SMS_blip(2,125,0.063)` and `cal_PEindex(18,49,3)`
     - `tdr_sliceorder(21,2)` and `tdr_sliceorder(42,2)`
     - `getsmspulse(...)`: rf signal/t/delay and the gz waveform
   - Subspace basis: Φ, TEs and S from `EPTI_gen_Subspace_Basis_GE`, built from the reference `.seq`.
   - P-code black boxes:
     - ghost compute/correct on a synthetic navigator with a known linear odd/even phase
     - `mrir_filter_raw_apodize_1d(ones(72,72),1,0.1)` and the same along dim 2, which gives the window directly
   - Expected runtime: about a minute.

## Phase 0: Cleanup and environment
1. `git mv` the used files into `legacy/`:
   - the seq scripts (path fixed) and `recon_EPTI_P1.m`
   - the seq helpers, `libs/utils/GE/*`, `idea/asc/*`, `calc_Pattern_WJG3.m`, `tdr_sliceorder.p`, the ~13 TOPPE RF files and Pulseq `+mr`
   - the recon-closure `.m` files and the core-path P-code (for the fixtures)
   - the `bart-0.8.00_EPTI/src` overlay (reference only) and the LICENSE files
2. Delete:
   - `libs/ismrmrd-master` and `libs/utils/ISMRMRD`
   - the Siemens readers, `Nifti_Analyze`, `mrir_save_mgz_measdat` and the unused P-code
   - the Pulseq demos/tests/docs/src and the unused TOPPE files
   - dead code: the extra `calc_Pattern*` variants, `SubspaceRecon_v6-8` and `*.asv`
3. Requirements: `sigpy`, `pygrappa`, `ismrmrd`, `h5py`, `nibabel`, `pyyaml`, `scikit-image`, `pytest`, plus a `pyproject.toml`.
   - On NVIDIA machines, optionally add `cupy-cuda12x`.
   - Check sigpy and pygrappa on Python 3.13 / numpy 2.5 right after installing. If sigpy fails, vendor the needed functions.
     If pygrappa fails, use the ESPIRiT+SENSE fallback.
4. Set `PYTHON_JULIACALL_THREADS=auto` before importing komamripy. Update the README (folder structure, glossary).

## Phase 1: Sequences in PyPulseq (`epti/seq/`)
1. Port the helpers 1:1 and add pytest checks against `fixtures.mat`.
2. `rf_sms.py`:
   - `sigpy.mri.rf.dzrf(200, tb, 'st', 'ls')`, then resample and scale the flip angle the way `makeslr` does.
   - Modulate with the `getsmsphase` table, then `pp.make_arbitrary_rf(no_signal_scaling=True)` with an explicit RF center.
   - Build the gz plateau trapezoid as in `seq_EPTI_imaging_P1.m:87-102`.
   - Check against the fixture pulse and a Bloch slice profile (`sigpy.mri.rf.sim.abrm`).
3. `calib.py` / `imaging.py`: ports driven by `Protocol` + `HardwareProfile`.
   - Write all the recon definitions plus `AdcDwell`, `SliceSep` and `NumDynamics`.
   - Put TRID labels at segment starts.
   - **Sim-subset option**: play only the selected slice group(s) and N dynamics, keeping the full geometry and the per-slice TR.
4. Parity: with the `siemens_prisma_ref` profile, Python output must match `tests/ref/*.seq` block by block
   (`waveforms_and_times`, `calculate_kspace`, ADC times, definitions).
5. `ge_mr750w` profile:
   - Rasters: gradients 4 µs, RF 2 µs, ADC 2 µs, blocks 4 µs. B0 = 3.0 T.
   - Dead and ringdown times use TOPPE's minimums from `systemspecs.m:31-33` (72/54/40 µs, read from source) until pge2's `sysGE` values are known.
   - Gmax/Smax are placeholders until the scanner facility confirms them.
   - ESP must be on the raster and outside the forbidden band. The ADC dwell must be a multiple of 2 µs.
6. `checks.py`: timing and raster alignment, Gmax/Smax, gradient spectrum against the forbidden ESPs,
   GE PNS (port of `+toppe/pns.m`, coil `xrmw`), and a TE/TR report.

## Phase 2: Simulation (`sim/`)
1. `phantom.py` builds on `km.brain_phantom2D`:
   - two different brain slices placed at ±sliceSep/2
   - a smooth synthetic Δw map and noise
   - `ss` as a knob for spin density
   - analytic ground truth (echo images, T2*, B0)
2. Scanner: GE-profile limits, with `receiver=BirdcageCoilSens(ncoils=8–16)`. This is a multi-channel array; EPTI needs coil diversity.
3. `run_sim.py`: simulate calib and imaging with the same phantom and coils, using `return_type="raw"`, and save ISMRMRD.
   Use CUDA on NVIDIA machines and CPU threads elsewhere (e.g. the 780M laptop).

## Phase 3: Core recon (`epti/recon/`, numpy/SigPy; CuPy when available), in pipeline order
1. `io/seqdefs.py`: `pp.Sequence().read()` → `EptiParams`.
   `io/ismrmrd_io.py`: reshape by arrival order, take Nx and the coil count from the data, flip even lines, split off the 3 navigator echoes.
2. `regrid.py`: port of `EPI_trapezoid_regridding_WJG.m`.
3. `ghost.py` (ours; replaces `mrir_artifact_ghost_*` at `EPTI_SMS_Estimate_Linear_Coeff_WJG.m:16`, `EPTI_SMS_Deghost_WJG.m:11-13`):
   - FFT along the readout into hybrid space, then a magnitude-weighted constant + linear phase fit between the reversed line and the mean of the forward navigator lines, per coil.
   - Take the median over calibration shots, apply the fit to the reversed lines, and FFT back.
4. `calib.py` (from `EPTI_Preprocessing_B0_sens_BARTv5.m`):
   - Fill the PE grid and run `pygrappa` (R=3, 18 ACS lines, 3×3 kernel, λ=0.01), then zero-pad to npe.
   - `EspiritCalib` gives the sensitivity maps, followed by a SENSE-1 combine.
   - B0 map: `np.unwrap` along echoes plus a least-squares fit. Then compute Phase0 for odd and even echoes, and the mask.
5. `subspace.py`: Φ from 91 T2 values (10–300 ms), SVD, K=4, scaled by S^(1/3); checked against the fixture.
   `sampling.py`: `pattern_GE` placement (`putRawtoPattern_SMS_zigzag_PF_selfnav_WJ4.m`).
   `drift.py`: B0 drift from the 3-line navigators (`DriftB0_esti_nav.m`).
6. `model.py`: SigPy linop for y = M·F·Σ_s E_s·S_{s,c}·Σ_k Φ'_{s,t,k}·c_{s,k}.
   - Build it from `sp.linop.Multiply`, `Sum`, `FFT` and `Reshape`, one slice group at a time.
   - Add an adjoint dot-product test.
7. `llr.py`: a `sp.prox.Prox` subclass that soft-thresholds singular values on 6×6×K blocks with random shifts, never across slices
   (BART's `-N` fully overlapping blocks can be added later as an option).
   `solve.py`: `LinearLeastSquares(A, y, proxg=LLR, solver='ADMM', rho=1e-2, max_iter=70, max_cg_iter=8)`, with data scaling as with `-w`.
   The linear option is CG with no regularization.
8. `postproc.py`: Σ_k c·Φ·Phase_T, flips and echo cut, Tukey apodization, `optimalCombineEchoes` and the T2* fit.
   `io/nifti.py`: write with nibabel (fixing the `pixdim[4]` bug).
9. `device.py`: `sp.Device(0)` if CuPy is available, else CPU, with a config override.
   `pipeline.py` and `recon/run_recon.py` tie it together.

## Phase 4 (later): extras
- High-order B0 update: replace AFNI with Python rigid registration.
- `moco_recon`.
- Koma motion experiments.
- GCC coil compression (needed for real 32+ channel data).
- Optionally, build the EPTI BART fork as an independent cross-check of `solve.py`.

## Phase 5 (later): GE scanner — questions for the scanner facility
- Which interpreter (pge2 or tv6), and which software version?
- What are Gmax, Smax and the gradient coil?
- What does `epiesp*.dat` list as forbidden ESPs?
- How many receive channels?
- Is trigger output supported?
- What is the raw format (ScanArchive or P-file), and how is it converted to ISMRMRD?

## Verification
- **How the P-code replacements are checked.** Library replacements are different implementations, so we don't expect them to match the P-code exactly.
  - **Tier 1: synthetic ground truth, no MATLAB.** These tests guard against silent convention mismatches (FFT centering and scaling, axis order, phase sign):
    - FFT helpers: centered round trip, plus a check of the scaling convention.
    - `pygrappa`: fully sampled multi-coil phantom → undersample to R=3 with 18 ACS lines → fill → NRMSE against the truth.
      Also confirms pygrappa's coil-axis and `lamda` conventions.
    - `EspiritCalib`: SENSE-combined image vs. ground truth, using sim data with known birdcage sensitivities.
    - LLR prox: λ=0 is the identity, and it lowers the rank of synthetic low-rank-plus-noise blocks.
    - Ghost correction: recovers a known injected linear phase.
  - **Tier 2: MATLAB fixtures.** Only where we wrote the code ourselves or the exact shape matters: ghost compute/correct,
    `tdr_sliceorder`, the apodization window, the subspace basis, and the seq helpers.
- `pytest tests/`:
  - helpers and ghost correction against `fixtures.mat`
  - `.seq` parity against `tests/ref`
  - GE profile: raster, timing, forbidden-ESP and PNS checks pass
  - adjoint test with relative error < 1e-5
  - regridding round-trip
  - pygrappa fill vs. fully sampled k-space on a synthetic multi-coil phantom
- End to end: `make_epti.py --config config/sim_small.yaml` → `sim/run_sim.py` → `recon/run_recon.py`, then compare to the phantom ground truth:
  - echo NRMSE, B0 and T2* error inside the mask
  - SMS leakage between the two slabs below a few %
- CPU and CuPy (NVIDIA) runs agree within float tolerance; time both.
