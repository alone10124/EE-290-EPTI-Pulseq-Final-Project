"""Scan protocols (what to acquire) and YAML config loading.

The defaults are the MGH "P1" protocol (2D single-shot EPTI, gradient echo, 3 mm isotropic, 3T), copied from
``legacy/sequences/seq_EPTI_imaging_P1.m`` and ``legacy/sequences/seq_EPTI_calib_P1.m``. Line numbers below refer
to those files. Units are SI (m, s) unless the name says otherwise.
"""

from dataclasses import dataclass, fields
from pathlib import Path

import yaml

from epti.hardware import HardwareProfile, get_profile


@dataclass
class FatSat:
    ppm: float = -3.31            # L72
    flip_deg: float = 110.0       # L74
    duration: float = 5.12e-3     # L74
    spoil_area: float = -0.1 / 1e-4   # L77, 1/m ("spoil up to 0.1 mm")


@dataclass
class ImagingProtocol:
    """2D SMS single-shot EPTI imaging scan (seq_EPTI_imaging_P1.m)."""
    fov: tuple = (216e-3, 216e-3)     # L30 (x, y)
    nx: int = 72                      # L31
    ny: int = 72
    thickness: float = 3e-3           # L32
    mb: int = 2                       # L33, multiband factor (<= 4)
    n_slices: int = 42                # L34
    slice_gap: float = 0.0            # L35
    min_te: bool = True               # L37
    te1: float = 4e-3                 # L39, only used when min_te is False
    tr: float = 1.8                   # L40
    min_tr: bool = False              # L38
    ro_os: int = 2                    # L41, readout oversampling
    nav_etl: int = 3                  # L42, navigator lines for ghost correction
    etl: int = 125                    # L43, = 128 - nav_etl
    rseg: int = 72                    # L44
    rpe: int = 6                      # L45
    nseg: int = 1                     # L46
    pf_shift: int = 0                 # L47
    esp: float = 0.49e-3              # L48, echo spacing
    blip_dur: float = 0.12e-3         # L29
    rf_angle_deg: float = 60.0        # L49
    rf_dur: float = 4e-3              # L50
    rf_tbw: float = 4.0               # L51
    prephase_dur: float = 1.2e-3      # L102, duration of the slice rephaser (sets all prephasers)
    fatsat: FatSat = None
    name: str = "P1_EPTI_sms2.seq"    # L55
    # --- additions for simulation (not in the original; defaults reproduce it) ---
    slice_groups: tuple | None = None   # spatial slice-group indices to play (None = all); see imaging.py
    n_dynamics: int = 1               # volumes in the .seq (the scanner normally repeats the file instead)

    def __post_init__(self):
        if self.fatsat is None:
            self.fatsat = FatSat()
        elif isinstance(self.fatsat, dict):
            self.fatsat = FatSat(**self.fatsat)
        self.fov = tuple(self.fov)
        if self.slice_groups is not None:
            self.slice_groups = tuple(self.slice_groups)

    @property
    def n_groups(self) -> int:
        return self.n_slices // self.mb

    @property
    def slice_sep(self) -> float:
        """Center-to-center distance between simultaneously excited slices (L36)."""
        return (self.thickness + self.slice_gap) * (self.n_slices / self.mb)


@dataclass
class CalibProtocol:
    """Multi-echo GRE-EPI calibration scan for coil sensitivity, B0 and odd/even phase (seq_EPTI_calib_P1.m)."""
    fov: float = 216e-3               # L34
    nx: int = 72                      # L35
    ny: int = 72
    thickness: float = 3e-3           # L36
    slice_gap: float = 0.0            # L37
    pat_factor: int = 3               # L38, outer-k-space PE undersampling
    acs_npe: int = 18                 # L39, fully sampled center lines
    nseg: int = 49                    # L40, full PE matrix of the calibration scan
    n_slices: int = 42                # L41
    min_te: bool = True               # L42
    te: float = 10.1e-3               # L43
    min_tr: bool = False              # L44
    tr: float = 1.2                   # L45
    ro_os: int = 2                    # L46
    nav_etl: int = 3                  # L47
    etl: int = 14                     # L48
    esp: float = 0.49e-3              # L49
    blip_dur: float = 0.12e-3         # L33 (only used for the readout area correction)
    rf_dur: float = 2.05e-3           # L50
    rf_angle_deg: float = 50.0        # L51
    rf_tbw: float = 4.0               # L81
    rf_apodization: float = 0.3       # L81
    rf_spoil_inc_deg: float = 0.0     # L53
    fatsat: FatSat = None
    name: str = "P1_Calib.seq"        # L56
    # --- addition for simulation: only play these slice positions (m); None = all ---
    slice_positions: tuple | None = None

    def __post_init__(self):
        if self.fatsat is None:
            self.fatsat = FatSat()
        elif isinstance(self.fatsat, dict):
            self.fatsat = FatSat(**self.fatsat)
        if self.slice_positions is not None:
            self.slice_positions = tuple(self.slice_positions)


@dataclass
class Config:
    hardware: HardwareProfile
    imaging: ImagingProtocol
    calib: CalibProtocol


def _make(cls, values):
    known = {f.name: f for f in fields(cls)}
    unknown = set(values) - set(known)
    if unknown:
        raise ValueError(f"unknown {cls.__name__} options: {sorted(unknown)}")
    values = dict(values)
    for key, val in values.items():
        # PyYAML reads "1e-3" (no decimal point) as a string; coerce it for numeric fields.
        if isinstance(val, str) and known[key].type in (float, "float"):
            values[key] = float(val)
    return cls(**values)


def load_config(path) -> Config:
    """Load a YAML config::

        hardware: ge_mr750w            # or {profile: ge_mr750w, max_grad: 40}
        imaging: {esp: 0.6e-3, ...}    # overrides of ImagingProtocol defaults
        calib: {...}                   # overrides of CalibProtocol defaults
    """
    data = yaml.safe_load(Path(path).read_text()) or {}
    return Config(
        hardware=get_profile(data.get("hardware", "ge_mr750w")),
        imaging=_make(ImagingProtocol, data.get("imaging") or {}),
        calib=_make(CalibProtocol, data.get("calib") or {}),
    )
