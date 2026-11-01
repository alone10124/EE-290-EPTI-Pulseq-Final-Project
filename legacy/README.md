# legacy/ — original MGH EPTI code, used as a porting reference

This folder holds the official EPTI release (version 09/2025, by Jian Wu, Zijing Dong and Fuyixue Wang,
[Advanced MRI Acquisition Lab, Martinos Center](https://mriacq-lab.martinos.org)). It has been pruned down to the files
our pipeline actually depends on. **Don't edit it.** The Python port lives in `epti/`, and every ported function's
docstring points back to the file and lines here that it came from.

> The MGH headers ask every user to complete the access request at https://martinos.org/~fw089/sequences.html
> (it's used for BRAIN Initiative grant tracking). Check those terms before publishing this folder in a public repo.

## What's here

| Path | What it is | Status |
|---|---|---|
| `sequences/seq_EPTI_calib_P1.m`, `seq_EPTI_imaging_P1.m` | MATLAB Pulseq scripts for the calibration and imaging scans | Run once in MATLAB to make reference `.seq` files (see below) |
| `recon/recon_EPTI_P1.m` | Top-level MATLAB recon | Reference only. It needs the compiled BART fork, plus AFNI when `B0_update=1`, and Siemens raw data by default |
| `libs/pulseq-master151/matlab/+mr` | MATLAB Pulseq 1.5.x (writes file format 1.5.1) | Needed by the scripts above |
| `libs/toppe-main/+toppe` | TOPPE. Only its SLR RF tools are used, via `libs/utils/GE/getsmspulse.m` | Kept whole because of runtime `import toppe.utils.*` calls. `pns.m` is the GE PNS model we port |
| `libs/utils/` | EPTI sequence helpers (blips, patterns, PE index), Siemens PNS/acoustic checks + `idea/asc/*.asc`, ISMRMRD raw readers, `.seq` parameter reader | Reference / port source |
| `libs/libs_EPTI_Recon/library_EPTI_recon/` | Core recon functions (regridding, ghost correction drivers, calibration, subspace recon) | Port source |
| `libs/libs_EPTI_Recon/library_EPTI_preprocessing/` | Echo combination and T2* fitting | Port source |
| `libs/libs_EPTI_Recon/library_EPTI_rawprocessing/library_EPI/mrir_toolbox/` | MGH "mrir" toolbox. **P-code (`.p`): runs in MATLAB but the source can't be read** | Black-box fixtures only |
| `bart-0.8.00_EPTI/` | Overlay files that turn stock BART 0.8.00 into the EPTI fork (`picsS`, `linopScript`, LLR change) | Reference for the solver's semantics; never compiled here |

## What was removed (still in git history: commit `d9f8b19`)
- `libs/ismrmrd-master`: the full ISMRMRD C++ library. Replaced by `pip install ismrmrd`.
- `libs/utils/ISMRMRD/`: a `siemens_to_ismrmrd` parameter map that no code referenced.
- Siemens raw readers: `mapVBVD`, `read_meas_*`, the XA30/XA60 readers, `*_Siemens.m`, and the Siemens-I/O P-code.
- NIfTI/geometry code that relies on Siemens headers (`Nifti_Analyze`, `mrir_save_mgz_measdat`, `MyNIFTIWrite`). Replaced by nibabel.
- Unused versions: other `calc_Pattern*` variants, `SubspaceRecon_Process_v6-8`, spin-echo bases, `*.asv` autosaves,
  and unused helpers under `library_EPTI_recon/utils` and `coilCompression_code`.
- Pulseq demos, tests, docs and C++ source; TOPPE images, uncompiled MEX sources, and a duplicate `jpauly/matlab` folder.

P-code dependencies can't be traced, so every non-I/O `.p` file in `mrir_toolbox` (including `private/`) was kept.

## Running the reference MATLAB scripts
You need MATLAB with the Signal Processing and Curve Fitting toolboxes (used by `getsmspulse`'s `firls`, `resample` and `smooth`).
You don't need scanner data, BART or AFNI.

One command does everything:
```matlab
cd legacy
make_reference_fixtures        % about a minute; close the plot windows it opens
```
It runs the two unchanged sequence scripts. They call `addpath(genpath('../libs'))`, which now resolves to `legacy/libs`.
It copies their output to `tests/ref/P1_Calib.seq` and `tests/ref/P1_EPTI_sms2.seq`, then saves `tests/ref/fixtures.mat`
with the helper outputs and the black-box P-code input/output pairs the Python tests compare against. Commit all three files.

These `.seq` files use Siemens Prisma settings. They're the parity reference for our Python port's `siemens_prisma_ref`
profile, **not** files to run on the GE scanner.

## Licenses (they must be kept with any ported code)
- `libs/libs_EPTI_Recon/*/license.txt`: BSD-style, © 2023 Fuyixue Wang, Zijing Dong.
- `libs/utils/LICENSE`: MIT. `libs/utils/GE/LICENSE.m`: MIT, © 2020 harmonizedMRI.
- `.../coilCompression_code/LICENSE.txt`: © 2013 UC Regents (Lustig). Educational, research and not-for-profit use.
- `libs/pulseq-master151/LICENSE`: MIT (Pulseq authors).
- `bart-0.8.00_EPTI/LICENSE`: BART's BSD-style license.
- TOPPE: this copy ships no license file; see the upstream repo (toppeMRI/toppe).
