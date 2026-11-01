# EE 290 Final Project - EPTI Pulse Sequence Improvements

We're modifying and improving the EPTI (echo planar time-resolved imaging) pulse sequence. The workflow is Python only:

1. Generate the `.seq` files with PyPulseq: a calibration scan plus the EPTI imaging scan.
2. Simulate them with KomaMRI (through `komamripy`), which writes raw data in ISMRMRD format.
3. Reconstruct with the `epti` package (NumPy/SigPy), and save figures and NIfTI images.

Later, the same `.seq` files go to a GE scanner. They get converted to GE format (pge2) outside this repo, and the raw data
comes back into the same recon. The porting plan and status are in [docs/EPTI_PORT_PLAN.md](docs/EPTI_PORT_PLAN.md).
The original MGH MATLAB code lives in [legacy/](legacy/README.md), kept for reference only.

## Folder Structure

```
./
├── config/                # .yaml configuration files (protocol + hardware profile) that drive the flow
├── docs/                  # Porting plan and design notes
├── epti/                  # Python package (pip install -e .), shared by sequences, sim and recon
│   ├── seq/               # PyPulseq sequence builders and hardware checks
│   ├── io/                # .seq definitions, ISMRMRD raw data, NIfTI output
│   └── recon/             # EPTI reconstruction
├── legacy/                # Original MGH EPTI MATLAB/C code (reference only, do not edit)
├── out/                   # Plots, simulated raw data, recon results (not tracked by git)
├── recon/                 # Recon entry-point scripts
├── sequences/             # .seq files. Can be imported to GE scanner or KomaMRI sim scripts
│   └── scripts/           # PyPulseq scripts that generate .seq files
├── sim/                   # Simulation scripts, using KomaMRI Python wrapper
├── tests/                 # pytest tests; tests/ref/ holds MATLAB-generated reference outputs
├── pyproject.toml         # Package definition for epti
├── requirements.txt       # Pinned environment (pip freeze)
├── LICENSE                # Project license details
└── README.md              # Project overview (this file)
```

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .                        # makes `import epti` work everywhere
# NVIDIA GPU only (optional): pip install cupy-cuda12x
```

KomaMRI runs on a single Julia thread unless you tell it otherwise. Set both of these **before** importing `komamripy`
(the second is required with more than one thread, or Julia can crash with a core dump):

```bash
export PYTHON_JULIACALL_THREADS=auto          # use all CPU cores
export PYTHON_JULIACALL_HANDLE_SIGNALS=yes
```
In a notebook, set them with `os.environ[...] = ...` in the first cell, before `import komamripy`.

### End-to-end run (simulation)

```bash
python sequences/scripts/make_epti.py config/sim_small.yaml   # 1. sequences (+ hardware checks)
python sim/run_sim.py config/sim_small.yaml                   # 2. KomaMRI simulation -> out/sim/   (~1 min)
python recon/run_recon.py config/sim_small.yaml               # 3. EPTI recon -> out/recon/         (~5 min; --iters 20 for ~2)
```
Step 3 writes echo images, T2* and combined images (`.npz` and NIfTI) and `out/recon/sim_small_summary.png`. It also compares the
result with the phantom's ground truth.

Run the tests with `pytest`. Tests that compare against MATLAB references are skipped until `tests/ref/` is filled in
(see [legacy/README.md](legacy/README.md)).

## Glossary

- **EPI / EPTI**: EPI reads many k-space lines after a single excitation by alternating the readout gradient. EPTI zig-zags
  through ky over time, so every echo time gets its own image, and a model-based recon fills in the gaps.
- **PE (phase encoding)**: the ky direction. PE lines are the ky rows of k-space.
- **ESP (echo spacing)**: the time between EPI readout lines. GE forbids ESPs whose readout frequency hits mechanical
  resonances of the gradient coil.
- **SMS / MB (simultaneous multi-slice / multiband factor)**: one RF pulse excites MB slices at once. The recon separates
  them using coil sensitivities plus blipped-CAIPI kz blips, which shift one slice by FOV/2.
- **K**: the number of subspace basis vectors. The recon solves for K coefficient maps instead of every echo image.
  Each echo image is Φ(t)·coefficients, where Φ comes from an SVD of simulated T2 decay curves.
- **GRAPPA**: k-space parallel imaging. Missing ky lines are filled from neighbouring acquired lines across coils, using
  weights fit on the fully sampled center (ACS) lines.
- **Ghost (Nyquist) correction**: odd and even EPI lines are read with opposite gradient polarity. That gives them a phase
  mismatch, which shows up as an N/2 ghost. Navigator lines measure the mismatch so it can be removed.
- **Hybrid / x-space**: data after an inverse FFT along the readout direction only, i.e. (x, ky).
- **TRID**: a Pulseq label marking the start of a repeating segment. The GE interpreter requires it.
