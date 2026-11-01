"""Recon parameters read back from a ``.seq`` file's ``[DEFINITIONS]`` section.

Replaces ``legacy/libs/utils/EPTI_SMS_Ext_parameters_pulseq.m``. The MATLAB version loads the whole sequence with
``mr.Sequence.read``; only the definitions are needed, so this parses that section directly (fast, and no Pulseq
version warnings). It also reads the definitions our port adds (``AdcDwell``, ``SliceSep``, ``NumDynamics``,
``PlayedSliceGroups``), falling back to the original's conventions when they're missing (e.g. MATLAB-made files).
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


def read_definitions(path):
    """Return the ``[DEFINITIONS]`` section as a dict of floats, float arrays, or strings."""
    defs, inside = {}, False
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith("["):
            if inside:
                break
            inside = line == "[DEFINITIONS]"
            continue
        if not inside or not line or line.startswith("#"):
            continue
        key, *vals = line.split()
        try:
            nums = np.array([float(v) for v in vals])
            defs[key] = nums[0] if nums.size == 1 else nums
        except ValueError:
            defs[key] = " ".join(vals)
    return defs


@dataclass
class EptiParams:
    """Parameters shared by both scans; scan-specific ones are None for the other scan."""
    path: str
    fov: np.ndarray                 # (x, y, z) m
    resolution: np.ndarray          # (x, y, z) m
    nx: int                         # TargetGriddedSamples: readout points after regridding
    te: float                       # s; imaging: center of the echo train; calib: echo N/2
    rise: float                     # readout lobe timing (TrapezoidGriddingParameters), s
    flat: float
    fall: float
    adc_delay: float
    dwell: float                    # ADC dwell, s
    slice_positions: np.ndarray     # m, in acquisition order (SlicePositions)
    actual_slice_positions: np.ndarray
    thickness: float
    # imaging scan
    n_lobes: int | None = None      # GRE_Lines: navigators + echoes
    rseg: int | None = None
    nseg: int | None = None
    rpe: int | None = None
    pf_shift: int = 0
    mb: int = 1
    pattern_ge: np.ndarray | None = None   # 1-based ky line per lobe (as in the original)
    slice_sep: float | None = None
    n_dynamics: int = 1
    played_groups: np.ndarray | None = None
    # calibration scan
    etl: int | None = None          # ETL: navigators + echoes
    cnseg: int | None = None        # full PE matrix of the calibration scan
    pe_lines: int | None = None
    pat_factor: int | None = None
    extra: dict = field(default_factory=dict)

    @property
    def esp(self):
        return self.rise + self.flat + self.fall

    @property
    def is_imaging(self):
        return self.n_lobes is not None


def _int(d, key, default=None):
    return int(round(float(d[key]))) if key in d else default


def read_params(path) -> EptiParams:
    d = read_definitions(path)
    trap = np.atleast_1d(d["TrapezoidGriddingParameters"])
    p = EptiParams(
        path=str(path), fov=np.atleast_1d(d["FOV"]), resolution=np.atleast_1d(d["Resolution"]),
        nx=_int(d, "TargetGriddedSamples"), te=float(d["TE"]),
        rise=trap[0], flat=trap[1], fall=trap[2], adc_delay=trap[3],
        dwell=float(d.get("AdcDwell", d.get("AdcRasterTime"))),   # MATLAB files carry the dwell in AdcRasterTime
        slice_positions=np.atleast_1d(d["SlicePositions"]),
        actual_slice_positions=np.atleast_1d(d["actual_SlicePositions"]),
        thickness=float(d["SliceThickness"]), extra=d,
    )
    if "GRE_Lines" in d:
        p.n_lobes = _int(d, "GRE_Lines")
        p.rseg, p.nseg, p.rpe = _int(d, "Rseg"), _int(d, "Nseg"), _int(d, "Rpe")
        p.pf_shift, p.mb = _int(d, "PFshift", 0), _int(d, "SMS", 1)
        p.pattern_ge = np.atleast_1d(d["pattern_GE"])
        p.slice_sep = float(d["SliceSep"]) if "SliceSep" in d else None
        p.n_dynamics = _int(d, "NumDynamics", 1)
        p.played_groups = np.atleast_1d(d["PlayedSliceGroups"]).astype(int) if "PlayedSliceGroups" in d else None
    if "ETL" in d:
        p.etl, p.cnseg = _int(d, "ETL"), _int(d, "CNseg")
        p.pe_lines, p.pat_factor = _int(d, "pe_lines"), _int(d, "PAT_factor")
    return p
