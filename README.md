# EE 290 Final Project - EPI Pulse Sequence

We're designing our own EPI-based method. To build and debug the workflow first, this repo runs a **simple EPI**:
single-shot, multi-slice 2D gradient-echo EPI with no acceleration or simultaneous multi-slice tricks. The workflow is
Python only:

1. **Sequence:** a PyPulseq script builds the `.seq` file from the building blocks in `mrtools.seq` and checks it
   against the scanner's limits.
2. **Simulation:** KomaMRI (through `komamripy`) simulates it on a multi-coil brain phantom with known ground truth.
3. **Reconstruction:** a recon script turns the raw signal into images (NumPy/SigPy, with helpers from `mrtools.recon`).

Later, the same `.seq` file goes to the scanner. It gets converted to GE format (pge2) outside this repo.

We first ported MGH's EPTI method; that work (incl. SMS pulses, ramp-sampling regridding, PNS and forbidden echo
spacing checks) is on the **`epti` branch**.

## Quick start

```bash
python main.py configs/simple_epi.toml     # sequence -> simulation -> image, about 50 s on a laptop CPU
```

Each run is described by one config file, `configs/<run>.toml`, with sections `[general]` (which scripts),
`[sequence]` (protocol, SI units), `[simulation]`, `[reconstruction]` and `[hardware]` (the scanner, incl. number of
receive coils). Output files are named after the config: `sequences/<run>.seq` (full, for the scanner),
`sequences/<run>_sim.seq` (the `sim_slices` subset, for the simulator) and `out/<run>_sim/`. To try a different FOV or
scanner with the same scripts, copy the config. See `mrtools/config.py` for the details and optional path overrides.

`main.py` runs these three commands, which you can also run one at a time:

```bash
python sequences/scripts/simple_epi.py configs/simple_epi.toml   # -> sequences/simple_epi.seq and simple_epi_sim.seq
python sim/run_sim.py configs/simple_epi.toml                    # -> out/simple_epi_sim/raw.npz, truth.npz
python recon/simple_epi.py configs/simple_epi.toml               # -> out/simple_epi_sim/recon.npz, recon.png
```

`python main.py configs/simple_epi.toml --skip-sim` reuses the existing raw data.

`recon.png` shows each reconstructed slice next to the image an ideal, fully sampled Cartesian scan would give
(`rho * exp(-TE/T2*)`, band-limited to the same k-space and shaded by the same receive coils), and their difference.
What remains in that difference is EPI physics (B0 distortion, T2* decay during the readout) plus any recon error. For
debugging, `simulation.b0_scale = 0` and `simulation.relaxation = false` turn those off (the error then drops to ~0.5%).

## Adding a new sequence

Shared code lives in `mrtools` and should change rarely. A new sequence `<name>` needs two scripts and a config:

- `sequences/scripts/<name>.py`: takes the config path; builds `sequences/<run>.seq` (and `<run>_sim.seq` for
  `sim_slices`) from the `[sequence]` table, whose keys are its `build_sequence` arguments. It should set the
  `MatrixSize`, `FOV`, `SlicePositions` and `SliceIndices` definitions (plus whatever its recon needs).
- `recon/<name>.py`: takes the config path and reads `config.load(...).recon_raw` / `.recon_seq`.
- `configs/<name>.toml`: `[general] sequence = "<name>"`.

The simulator (`sim/run_sim.py`) is shared. Then `python main.py configs/<name>.toml` runs the whole workflow.

## Folder Structure

```
./
├── main.py                # Runs the whole workflow for one config: python main.py configs/<run>.toml
├── configs/               # One TOML file per run (protocol, simulation, recon, hardware)
├── mrtools/               # Shared Python package (pip install -e .)
│   ├── config.py          # Reads configs, file-naming rules
│   ├── hardware.py        # HardwareProfile (limits, rasters, coils; builds the KomaMRI scanner)
│   ├── seq/               # PyPulseq building blocks (EPI readout, fat-sat, slice order, TRID), hardware checks
│   ├── sim/               # KomaMRI: brain phantom with ground truth, simulator wrapper
│   └── recon/             # ESPIRiT coil maps, coil combination
├── sequences/             # Generated .seq files (not tracked by git; rebuilt by the scripts)
│   └── scripts/           # One script per sequence; builds its .seq file
├── sim/                   # Simulator script, shared by all sequences (+ notebook)
├── recon/                 # One recon script per sequence
├── out/                   # Simulated data, images, figures (not tracked by git)
├── tests/                 # pytest tests
├── pyproject.toml         # Package definition for mrtools
├── requirements.txt       # Pinned environment (pip freeze)
├── LICENSE                # Project license details
└── README.md              # Project overview (this file)
```

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .                        # makes `import mrtools` work everywhere
```

KomaMRI runs on a single Julia thread unless you tell it otherwise. These two settings must be in place **before**
`komamripy` is imported (the second is required with more than one thread, or Julia can crash with a core dump).
Importing `mrtools.sim` sets them for you (`sim/run_sim.py` prints the thread count); in a notebook, import
`mrtools.sim` before `komamripy`, or set them yourself:

```bash
export PYTHON_JULIACALL_THREADS=auto          # use all CPU cores
export PYTHON_JULIACALL_HANDLE_SIGNALS=yes
```

Run the tests with `pytest`.

## Glossary

- **EPI**: reads many k-space lines after a single excitation by alternating the readout gradient.
- **PE (phase encoding)**: the ky direction. PE lines are the ky rows of k-space.
- **ESP (echo spacing)**: the time between EPI readout lines. GE forbids ESPs whose readout frequency hits mechanical
  resonances of the gradient coil.
- **Coil sensitivity maps / ESPIRiT**: how strongly each receive coil sees each voxel; ESPIRiT estimates them from
  fully sampled central k-space.
- **GRAPPA**: fills in *skipped* k-space lines from neighbouring acquired lines across coils. Not needed while the
  calibration is fully sampled.
- **T2\***: the gradient-echo decay constant, S(TE) = S0·exp(−TE/T2\*). Fitting it from several TEs gives T2\* maps and
  the BOLD-optimal echo combination.
- **PNS (peripheral nerve stimulation)**: fast gradient switching can stimulate nerves. GE limits it; not modeled
  yet (the `epti` branch has a GE PNS model).
- **TRID**: a Pulseq label marking the start of a repeating segment. GE's interpreter requires it.
