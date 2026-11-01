"""EPTI k-space sampling patterns, ported 1:1 from the MGH MATLAB helpers.

Ky positions keep the original 1-based PE-line numbering (e.g. the k-space center of a 72-line matrix is 36), because
these values are written into the ``.seq`` file as ``pattern_GE`` and read back by the recon.
"""

import numpy as np

from epti.util import mround


def _zigzag(nseg, rseg, etl, rpe, rpe2, pf_shift):
    """Shared core of calc_Pattern_WJG3/4: the EPTI zig-zag ky(t) trajectory and the center-crossing echo indices."""
    if nseg != 1:
        # The original compares a vector Start_point inside `&&` conditions, which errors in MATLAB for Nseg > 1.
        raise NotImplementedError("calc_Pattern_WJG3/4 only support nseg == 1 (as in the original)")
    start = (rseg + pf_shift) / 2
    ky = np.zeros((nseg, etl))
    center_index = []
    skip, skip2 = 1, 0
    pe = start
    sign = -1
    for i in range(1, etl + 1):
        ky[0, i - 1] = pe
        if i <= 2:
            continue
        if start - 3 <= pe <= start + 3 and sign == 1:
            pe += sign * rpe2
            center_index.append(i)
        elif start - 2 <= pe <= start + 4 and sign == -1:
            pe += sign * rpe2
            center_index.append(i)
        elif start - 9 <= pe <= start - 3 and sign == 1 and skip2 == 1:
            pe = start - rpe2 - 1
            skip2 = 0
        elif start + 4 <= pe <= start + 10 and sign == -1 and skip2 == 1:
            pe = start + 2 * rpe2
            skip2 = 0
        elif (pe >= rseg + pf_shift - rpe / 2 or pe <= rpe / 2 + pf_shift) and skip == 1:
            sign = -sign
            pe += sign * rpe / 2
            skip = 0
        elif (pe > rseg + pf_shift - rpe or pe <= rpe + pf_shift) and skip == 1:
            sign = -sign
            pe += sign * rpe / 2
            skip = 0
        else:
            pe += sign * rpe
            skip = 1
            skip2 = 1
    center_index = np.array(center_index, dtype=int)
    # MATLAB: ind = find(diff(diff(ceter_index)) < 0);  ceter_index(ind + 1)
    crossings = center_index[np.nonzero(np.diff(np.diff(center_index)) < 0)[0] + 1]
    return ky, list(crossings)


def calc_pattern_wjg3(nseg, rseg, etl, rpe, rpe2, pf_shift):
    """EPTI zig-zag ky schedule used by the sequence script.

    Ported from ``legacy/libs/libs_EPTI_Recon/library_EPTI_recon/calc_Pattern_WJG3.m`` (MGH, Wang/Dong/Wu).

    Returns
    -------
    ky_acq : (nseg, etl) array of 1-based PE line numbers, one per echo. The first entries (navigators) stay at center.
    ind : list of 1-based echo indices where the trajectory crosses the k-space center (padded to 4 entries).
    """
    ky, ind = _zigzag(nseg, rseg, etl, rpe, rpe2, pf_shift)
    ind = ind[:4]
    if len(ind) == 3:          # MATLAB: ind(4) = ind(3); ind(3) = ind(2)
        ind = [ind[0], ind[1], ind[1], ind[2]]
    elif len(ind) == 2:        # MATLAB: ind(4) = ind(2); ind(3) = ind(2)
        ind = [ind[0], ind[1], ind[1], ind[1]]
    return ky, ind


def calc_pattern_wjg4(nseg, rseg, etl, rpe, rpe2, pf_shift):
    """Same ky schedule as :func:`calc_pattern_wjg3`, with the recon's version of the center-crossing indices.

    Ported from ``legacy/libs/libs_EPTI_Recon/library_EPTI_recon/calc_Pattern_WJG4.m`` (used by recon_EPTI_P1.m L115).
    """
    ky, ind = _zigzag(nseg, rseg, etl, rpe, rpe2, pf_shift)
    # MATLAB: ind(4) = ind(3); ind(3) = ind(2)   (auto-extends when ind has 3 entries)
    if len(ind) < 3:
        raise ValueError("expected at least 3 center crossings")
    if len(ind) == 3:
        ind = ind + [ind[2]]
    else:
        ind[3] = ind[2]
    ind[2] = ind[1]
    return ky, ind


def calc_pattern_sms_blip(mb, etl, slice_sep):
    """Blipped-CAIPI kz steps (in units of ``deltakz``) for each echo, and ``deltakz`` in 1/m.

    Ported from ``legacy/libs/utils/calc_Pattern_SMS_blip.m``. MB=2 alternates +1/-1 (a FOV/2 shift between the two
    slices); MB=3 uses the same +-1 pattern with deltakz = 1/(2*slice_sep); MB=4 cycles (-1, -1, +2) / (3*slice_sep).
    """
    loop = np.arange(1, etl + 1)
    if mb == 1:
        return np.zeros(etl), 0.0
    if mb in (2, 3):
        steps = np.array([-1, 1])[loop % 2]
        return steps.astype(float), (1 / mb if mb == 2 else 1 / 2) / slice_sep
    if mb == 4:
        steps = np.array([2, -1, -1])[loop % 3]
        return steps.astype(float), (1 / 3) / slice_sep
    raise ValueError("calc_Pattern_SMS_blip only supports MB <= 4")


def cal_pe_index(acs_npe, nseg, pat_factor):
    """Phase-encoding line indices (centered on 0) for the calibration scan: fully sampled center, R-fold outside.

    Ported from ``legacy/libs/utils/cal_PEindex.m``. For (18, 49, 3) this gives 29 lines:
    ``[-24:3:-9, -8:9, 12:3:24]``.
    """
    n_outer = int(np.ceil((nseg - acs_npe) / pat_factor))
    n_pe = acs_npe + n_outer
    index = np.zeros(n_pe)
    for i in range(1, n_pe + 1):
        if i < n_outer / 2 + 1:
            index[i - 1] = -(nseg - 1) / 2 + pat_factor * (i - 1)
        elif i >= n_outer / 2 + acs_npe + 1:
            index[i - 1] = -(nseg - 1) / 2 + pat_factor * (i - acs_npe + mround(acs_npe / pat_factor) - 1)
        else:
            index[i - 1] = -acs_npe / 2 + (i - np.floor(n_outer / 2) - 1)
    return index


def interleaved_slice_order(n):
    """Acquisition order (0-based spatial indices, in time order) for interleaved 2D slice excitation.

    Replaces the encrypted ``tdr_sliceorder(n, 2)`` P-code
    (``legacy/libs/libs_EPTI_Recon/library_EPTI_rawprocessing/library_EPI/mrir_toolbox/tdr_sliceorder.p``).
    The scripts use it as ``SliceOrder = tdr_sliceorder(n, 2); [~, SliceOrder] = sort(SliceOrder)``, so this function
    corresponds to that *sorted* result. Assumed rule (Siemens-style interleaving): odd n -> 0, 2, 4, ..., 1, 3, ...;
    even n -> 1, 3, 5, ..., 0, 2, .... UNVERIFIED until tests/ref/fixtures.mat exists (see test_patterns.py).
    Either way, the recon reads the order actually used from the ``actual_SlicePositions`` definition.
    """
    idx = np.arange(n)
    if n % 2:
        return np.concatenate([idx[0::2], idx[1::2]])
    return np.concatenate([idx[1::2], idx[0::2]])
