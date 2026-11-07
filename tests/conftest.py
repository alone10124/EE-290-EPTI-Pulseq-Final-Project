"""Shared test settings."""

import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def load_script(relpath, module_name):
    """Import a script (e.g. sequences/scripts/simple_epi.py) as a module. Sequence and recon scripts share file names,
    so each one gets its own module name."""
    spec = importlib.util.spec_from_file_location(module_name, REPO / relpath)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
