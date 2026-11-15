"""Run configuration: one TOML file per run, ``configs/<run>.toml``, with sections

    [general]         which scripts to run: sequence = "<name>" (sequences/scripts/<name>.py), recon = "<name>"
    [sequence]        protocol (SI units) passed to the sequence script's ``build_sequence``, plus which slices to write
    [simulation]      phantom options for sim/run_sim.py (optional ``ncoils``: simulate fewer coils than [hardware];
                      optional ``threads``: cap KomaMRI's threads, default all)
    [reconstruction]  options for the recon script
    [hardware]        a ``mrtools.hardware.HardwareProfile``

File names follow the config's file name (the "run" name), so several configs can share the same scripts:

    sequences/<run>.seq          full sequence, for the scanner
    sequences/<run>_sim.seq      the slice subset in ``sequence.sim_slices``, for the simulator
    out/<seq stem>/raw.npz       raw data of a simulated (or scanned) .seq, plus truth.npz and recon outputs

Each path can be overridden: ``simulation.seq`` / ``simulation.raw`` and ``reconstruction.seq`` /
``reconstruction.raw`` (relative to the repo root). By default the recon reads what the simulator wrote.
"""

import tomllib
from dataclasses import dataclass
from pathlib import Path

from mrtools.hardware import HardwareProfile

REPO = Path(__file__).resolve().parents[1]
SECTIONS = ("general", "sequence", "simulation", "reconstruction", "hardware")


@dataclass
class Config:
    run: str                       # config file name without .toml
    general: dict
    sequence: dict
    simulation: dict
    reconstruction: dict
    hardware: HardwareProfile

    @property
    def seq_path(self):
        """The full sequence (for the scanner)."""
        return REPO / "sequences" / f"{self.run}.seq"

    @property
    def sim_seq_path(self):
        """The short sequence for the simulator; ``None`` if ``sequence.sim_slices`` is empty."""
        return REPO / "sequences" / f"{self.run}_sim.seq" if self.sequence.get("sim_slices") else None

    def _path(self, section, key, default):
        value = getattr(self, section).get(key)
        return REPO / value if value else default

    @property
    def sim_input(self):
        """.seq file the simulator runs: ``simulation.seq``, else the sim sequence, else the full one."""
        return self._path("simulation", "seq", self.sim_seq_path or self.seq_path)

    @property
    def sim_raw(self):
        return self._path("simulation", "raw", REPO / "out" / self.sim_input.stem / "raw.npz")

    @property
    def recon_seq(self):
        return self._path("reconstruction", "seq", self.sim_input)

    @property
    def recon_raw(self):
        return self._path("reconstruction", "raw", self.sim_raw)


def display_path(path):
    """``path`` relative to the repo root for messages, or unchanged if it lies outside the repo (e.g. a config's
    ``raw = "/tmp/..."`` override)."""
    path = Path(path)
    return path.relative_to(REPO) if path.is_relative_to(REPO) else path


def load(path):
    """Read ``configs/<run>.toml`` into a ``Config``. Missing sections are empty; unknown sections raise an error."""
    path = Path(path)
    with open(path, "rb") as f:
        d = tomllib.load(f)
    unknown = set(d) - set(SECTIONS)
    if unknown:
        raise ValueError(f"{path}: unknown sections {sorted(unknown)}; allowed: {list(SECTIONS)}")
    general = d.get("general", {})
    if "sequence" not in general:
        raise ValueError(f"{path}: [general] needs `sequence` (the script name in sequences/scripts/)")
    general.setdefault("recon", general["sequence"])
    return Config(run=path.stem, general=general, sequence=d.get("sequence", {}), simulation=d.get("simulation", {}),
                  reconstruction=d.get("reconstruction", {}), hardware=HardwareProfile.from_dict(d["hardware"]))
