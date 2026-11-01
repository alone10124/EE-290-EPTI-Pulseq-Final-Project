"""k-space pattern helpers (epti.seq.patterns)."""

import numpy as np

from epti.seq.patterns import (cal_pe_index, calc_pattern_sms_blip, calc_pattern_wjg3, calc_pattern_wjg4,
                               interleaved_slice_order)


# ---------------------------------------------------------------- Tier 1: properties, no MATLAB needed
def test_epti_zigzag_properties():
    ky, ind = calc_pattern_wjg3(1, 72, 128, 6, 2, 0)
    assert ky.shape == (1, 128)
    assert np.all(ky[0, :3] == 36)                      # navigators at the k-space center
    assert 1 <= ky.min() and ky.max() <= 72             # stays inside the 72-line matrix
    steps = np.diff(ky[0])[2:]
    assert len(steps) == 125                            # one blip per imaging lobe
    assert set(np.abs(steps).astype(int)) <= {2, 3, 4, 6}
    assert len(ind) == 4 and all(3 <= i <= 128 for i in ind)
    assert np.array_equal(ky, calc_pattern_wjg4(1, 72, 128, 6, 2, 0)[0])   # same trajectory in WJG3 and WJG4


def test_pe_index_matches_known_values():
    expected = np.r_[np.arange(-24, -8, 3), np.arange(-8, 10), np.arange(12, 25, 3)]
    assert np.array_equal(cal_pe_index(18, 49, 3), expected)


def test_sms_blips():
    kz, dkz = calc_pattern_sms_blip(2, 6, 0.063)
    assert np.array_equal(kz, [1, -1, 1, -1, 1, -1]) and np.isclose(dkz, 1 / (2 * 0.063))
    kz4, dkz4 = calc_pattern_sms_blip(4, 6, 0.03)
    assert np.array_equal(kz4, [-1, -1, 2, -1, -1, 2]) and np.isclose(dkz4, 1 / (3 * 0.03))
    # net kz returns to a fixed set of phase states (CAIPI): cumulative sum stays bounded
    assert np.abs(np.cumsum(calc_pattern_sms_blip(2, 125, 0.063)[0])).max() == 1


def test_slice_order_is_a_permutation():
    for n in (1, 2, 7, 21, 42):
        assert sorted(interleaved_slice_order(n)) == list(range(n))


# ---------------------------------------------------------------- Tier 2: MATLAB fixtures
def test_zigzag_matches_matlab(fixtures):
    ky3, ind3 = calc_pattern_wjg3(1, 72, 128, 6, 2, 0)
    ky4, ind4 = calc_pattern_wjg4(1, 72, 128, 6, 2, 0)
    assert np.array_equal(ky3[0], np.ravel(fixtures["ky_WJG3"]))
    assert np.array_equal(ind3, np.ravel(fixtures["ind_WJG3"]))
    assert np.array_equal(ky4[0], np.ravel(fixtures["ky_WJG4"]))
    assert np.array_equal(ind4, np.ravel(fixtures["ind_WJG4"]))


def test_sms_blips_match_matlab(fixtures):
    for mb in (1, 2, 3, 4):
        kz, dkz = calc_pattern_sms_blip(mb, 125, 0.063)
        assert np.array_equal(kz, np.ravel(fixtures[f"kz_mb{mb}"]))
        assert np.isclose(dkz, fixtures[f"dkz_mb{mb}"])


def test_pe_index_matches_matlab(fixtures):
    assert np.array_equal(cal_pe_index(18, 49, 3), np.ravel(fixtures["PEindex_18_49_3"]))


def test_slice_order_matches_encrypted_tdr_sliceorder(fixtures):
    """interleaved_slice_order(n) must equal MATLAB's  [~, order] = sort(tdr_sliceorder(n, 2))  (0-based)."""
    for key in [k for k in fixtures if k.startswith("sliceorder_mode2_N")]:
        n = int(key.rsplit("N", 1)[1])
        tdr = np.atleast_1d(fixtures[key])
        assert np.array_equal(interleaved_slice_order(n), np.argsort(tdr, kind="stable")), key
