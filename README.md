# EE 290 Final Project - EPI Pulse Sequence

We're designing our own EPI-based method. To build and debug the workflow first, this repo runs a **simple EPI**:
single-shot, multi-slice 2D gradient-echo EPI with no acceleration or simultaneous multi-slice tricks. The workflow is
Python only:

1. **Sequence:** a PyPulseq script builds the `.seq` file from the building blocks in `mrtools.seq` and checks it
   against the GE MR750w's limits.
2. **Simulation:** KomaMRI (through `komamripy`) simulates it on a multi-coil brain phantom with known ground truth.
3. **Reconstruction:** a recon script turns the raw signal into images (NumPy/SigPy, with helpers from `mrtools.recon`).

Later, the same `.seq` file goes to the scanner. It gets converted to GE format (pge2) outside this repo.

We first ported MGH's EPTI method; that work is on the **`epti` branch**. [docs/LESSONS_LEARNED.md](docs/LESSONS_LEARNED.md)
lists what carried over, the simplifying assumptions we're making now, and where to get each removed piece back.

## Quick start

```bash
python main.py simple_epi          # sequence -> simulation -> image, about 25 s on a laptop CPU
```

`main.py` runs these four commands, which you can also run one at a time:

```bash
python sequences/scripts/simple_epi.py                                     # -> sequences/simple_epi.seq (all 40 slices, for the scanner)
python sequences/scripts/simple_epi.py --name simple_epi_sim --slices 20   # -> sequences/simple_epi_sim.seq (one slice, to simulate)
python sim/run_sim.py sequences/simple_epi_sim.seq                         # -> out/simple_epi_sim/raw.npz, truth.npz
python recon/simple_epi.py sequences/simple_epi_sim.seq                    # -> out/simple_epi_sim/recon.npz, recon.png
```

`recon.png` shows each reconstructed slice next to the image the phantom should give at that TE, and their difference.
Options: `python main.py simple_epi --sim-slices 19 20` simulates two slices; `--skip-sim` reuses the existing data.

The protocol (FOV, matrix, slice thickness, ...) is set by the constants at the top of
`sequences/scripts/simple_epi.py`. The script records what the simulator and recon need (FOV, matrix, slice positions,
TE) in the `.seq` file's `[DEFINITIONS]`, so there is no separate config file.

## Adding a new sequence

Shared code lives in `mrtools` and should change rarely. A new sequence `<name>` needs two scripts:

- `sequences/scripts/<name>.py`: builds the `.seq` file. To work with `main.py` and the simulator it should accept
  `--name` and `--slices`, and set the `MatrixSize`, `FOV` and `SlicePositions` definitions (plus whatever its recon
  needs).
- `recon/<name>.py`: takes the `.seq` file and reads `out/<seq name>/raw.npz`.

The simulator (`sim/run_sim.py`) is shared. Then `python main.py <name>` runs the whole workflow.

## Folder Structure

```
./
├── main.py                # Runs the whole workflow for one sequence: python main.py <name>
├── docs/                  # Design notes (LESSONS_LEARNED.md)
├── mrtools/               # Shared Python package (pip install -e .)
│   ├── hardware.py        # GE MR750w hardware profile
│   ├── seq/               # PyPulseq building blocks (EPI readout, fat-sat, slice order, SMS pulses, CAIPI), GE checks
│   ├── sim/               # KomaMRI: brain phantom with ground truth, simulator wrapper
│   └── recon/             # ESPIRiT coil maps, coil combination
├── sequences/             # .seq files (for the scanner or KomaMRI)
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
# NVIDIA GPU only (optional): pip install cupy-cuda12x, then add --gpu to sim/run_sim.py
```

KomaMRI runs on a single Julia thread unless you tell it otherwise. These two settings must be in place **before**
`komamripy` is imported (the second is required with more than one thread, or Julia can crash with a core dump).
`mrtools.sim` sets them for you; set them yourself if you import `komamripy` directly, e.g. in a notebook.

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
  yet (see docs/LESSONS_LEARNED.md).
- **TRID**: a Pulseq label marking the start of a repeating segment. GE's interpreter requires it.
