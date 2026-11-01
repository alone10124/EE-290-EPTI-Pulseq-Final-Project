"""Build the calibration + imaging pair from a config, check them, and write the .seq files."""

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from epti.config import Config
from epti.seq.calib import build_calib
from epti.seq.checks import Report, check_sequence
from epti.seq.imaging import build_imaging, excited_slice_positions


def readout_consistency(img_info, cal_info) -> Report:
    """The ghost-correction coefficients come from the calibration navigators and are applied to the imaging data
    (legacy/recon/recon_EPTI_P1.m L81-83, L135), so both scans need the same readout lobe and ADC."""
    rep = Report("calibration vs imaging readout")
    pairs = {
        "readout amplitude": (img_info.gx.amplitude, cal_info.gx.amplitude),
        "readout rise time": (img_info.gx.rise_time, cal_info.gx.rise_time),
        "readout flat time": (img_info.gx.flat_time, cal_info.gx.flat_time),
        "ADC samples": (img_info.adc.num_samples, cal_info.adc.num_samples),
        "ADC dwell": (img_info.adc.dwell, cal_info.adc.dwell),
        "ADC delay": (img_info.adc.delay, cal_info.adc.delay),
    }
    for name, (a, b) in pairs.items():
        rep.add(name, np.isclose(a, b, rtol=1e-9, atol=1e-12), f"imaging {a:.6g}, calib {b:.6g}")
    return rep


def build_all(cfg: Config, check=True):
    """Build both sequences. If the imaging scan plays only some slice groups, the calibration scan plays the same
    slices (unless ``calib.slice_positions`` is set explicitly)."""
    calib_protocol = cfg.calib
    if cfg.imaging.slice_groups is not None and calib_protocol.slice_positions is None:
        calib_protocol = replace(calib_protocol,
                                 slice_positions=tuple(excited_slice_positions(cfg.imaging, cfg.imaging.slice_groups)))
    seq_cal, info_cal = build_calib(calib_protocol, cfg.hardware)
    seq_img, info_img = build_imaging(cfg.imaging, cfg.hardware)
    reports = []
    if check:
        reports = [
            check_sequence(seq_cal, cfg.hardware, info_cal.esp, f"calibration ({cfg.hardware.name})"),
            check_sequence(seq_img, cfg.hardware, info_img.esp, f"imaging ({cfg.hardware.name})"),
            readout_consistency(info_img, info_cal),
        ]
    return SimpleNamespace(calib=seq_cal, imaging=seq_img, calib_info=info_cal, imaging_info=info_img,
                           reports=reports)


def write_all(built, outdir, stem):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = SimpleNamespace(calib=outdir / f"{stem}_calib.seq", imaging=outdir / f"{stem}_imaging.seq")
    built.calib.write(str(paths.calib))
    built.imaging.write(str(paths.imaging))
    return paths
