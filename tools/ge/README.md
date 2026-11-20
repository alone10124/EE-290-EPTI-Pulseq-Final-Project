# GE conversion dry run

These scripts convert our `.seq` files to GE scanner format and run the converters' own checks locally, before a
scanner session. They need MATLAB. See [docs/SCANNER_SESSION.md](../../docs/SCANNER_SESSION.md) for the session plan.

```powershell
cd tools/ge
matlab -batch "ge_dryrun"        # pge2: seq2ceq -> pge2.check -> pge2.validate -> sequences/<name>.pge
matlab -batch "ge_dryrun_tv6"    # TOPPE v6: seq2ge -> sequences/<name>_tv6.tar (includes preflightcheck: SAR, heating)
```

Both take an optional list of `.seq` files: `ge_dryrun({'../../sequences/sim_small_imaging.seq'})`. `ge_dryrun` also
takes hardware options, for example `ge_dryrun([], 'coil', 'xrm', 'g_max', 5, 'slew_max', 20, 'psd_rf_wait', 100e-6)`.

## Third-party code (not ours, cloned 2026-11-19, git-ignored)

| Folder | Source | Version |
|---|---|---|
| `PulCeq/` | github.com/HarmonizedMRI/PulCeq | `main` 434f871 (pge2) |
| `PulCeq_tv6/` | git worktree of `PulCeq`, branch `tv6_v1.10` | for tv6 |
| `toppe/` | github.com/toppeMRI/toppe | `main` f3ae4f3 |

Pulseq itself comes from `legacy/libs/pulseq-master151`.

PulCeq says it is superseded by PulSeg + pge2 (github.com/HarmonizedMRI), but Hannah uses PulCeq.
If her versions differ from these, rerun with hers.

## `tv6_patch/`: workarounds for the tv6 path

These replace `toppe.writeentryfile` / `toppe.writecoresfile`, because they come first on the path.
- **Entry file version:** current TOPPE writes the v7 `toppeN.entry` format by default. PulCeq's tv6 `ceq2ge` doesn't
  ask for v6, so `preflightcheck` fails ("Path in toppeN.entry must end with a slash"). The patch defaults to v6.
- **Line endings:** `fopen(..., 'wt')` writes CRLF on Windows. The tv6 interpreter runs on Linux, and the entry
  file path check fails. The patch writes LF.
- TOPPE branch `dev/v6` doesn't work either: its `writemod` rejects any RF block whose gradients outlast the RF.

## Raw data back from the scanner

`ge_raw_export.m` turns a P-file (tv6, needs the `_tv6.tar` that was run) or a ScanArchive (pge2, needs GE's Orchestra
SDK) into a plain HDF5 file in ADC order, which `recon/run_recon.py` reads with `--calib-data` / `--imaging-data`.
See "Reconstructing the scanner data" in [docs/SCANNER_SESSION.md](../../docs/SCANNER_SESSION.md).
`test_raw_export.m` exists only for `tests/test_ge_raw.py`.
