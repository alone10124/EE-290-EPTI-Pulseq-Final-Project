"""Shared test fixtures.

Tests marked "Tier 2" compare against MATLAB outputs in tests/ref/, which are produced by running
legacy/make_reference_fixtures.m once (see legacy/README.md). They are skipped until those files exist.
"""

from pathlib import Path

import pytest
import scipy.io

REF = Path(__file__).parent / "ref"
REPO = Path(__file__).resolve().parents[1]



def _skip_missing(path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(REPO)} not generated yet: run legacy/make_reference_fixtures.m in MATLAB")


@pytest.fixture(scope="session")
def fixtures():
    path = REF / "fixtures.mat"
    _skip_missing(path)
    return scipy.io.loadmat(path, squeeze_me=True, struct_as_record=False)


@pytest.fixture(scope="session")
def ref_seq():
    """MATLAB-generated reference sequences, read with PyPulseq."""
    import pypulseq as pp

    def load(name):
        path = REF / name
        _skip_missing(path)
        seq = pp.Sequence()
        seq.read(str(path))
        return seq

    return load
