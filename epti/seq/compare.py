"""Compare two sequences by what the scanner would play, not by their block lists.

Used for the MATLAB-parity tests. A block-by-block comparison isn't meaningful here, because the Python port attaches
counter labels to the next block instead of using zero-duration label blocks (see ``epti.seq.common.BlockWriter``).
"""

import numpy as np

from epti.seq.checks import gradient_waveforms


def compare_sequences(a, b, time_range=None):
    """Return a dict of differences between sequences ``a`` and ``b``.

    Keys: ``duration`` (s); ``n_adc_samples`` (pair); ``adc_time_maxdiff`` (s); ``grad_rel_maxdiff`` (max abs
    difference on a 1 us grid, relative to the largest gradient); ``rf_time_maxdiff`` (s) and ``rf_phase_maxdiff``
    (rad) for excitation pulse centers; ``definitions`` (keys of ``a`` whose values differ or are missing in ``b``);
    ``labels`` (labels whose per-ADC values differ).
    """
    out = {"duration": abs(a.duration()[0] - b.duration()[0])}

    wa = a.waveforms_and_times(time_range=time_range)
    wb = b.waveforms_and_times(time_range=time_range)
    adc_a, adc_b = np.asarray(wa[3]), np.asarray(wb[3])
    out["n_adc_samples"] = (adc_a.size, adc_b.size)
    out["adc_time_maxdiff"] = float(np.abs(adc_a - adc_b).max()) if adc_a.size == adc_b.size else np.inf

    t_end = time_range[1] if time_range else max(a.duration()[0], b.duration()[0])
    window = (time_range[0] if time_range else 0.0, t_end)
    _, ga = gradient_waveforms(a, 1e-6, window)
    _, gb = gradient_waveforms(b, 1e-6, window)
    n = min(ga.shape[1], gb.shape[1])
    scale = max(np.abs(ga).max(), np.abs(gb).max(), 1e-12)
    out["grad_rel_maxdiff"] = float(np.abs(ga[:, :n] - gb[:, :n]).max() / scale)

    ea, eb = np.asarray(wa[1]), np.asarray(wb[1])        # excitation: rows (time, frequency, phase)
    if ea.shape == eb.shape and ea.size:
        out["rf_time_maxdiff"] = float(np.abs(ea[0] - eb[0]).max())
        out["rf_phase_maxdiff"] = float(np.abs(np.angle(np.exp(1j * (ea[2] - eb[2])))).max())
    else:
        out["rf_time_maxdiff"] = out["rf_phase_maxdiff"] = np.inf

    diff_defs = []
    for key, va in a.definitions.items():
        vb = b.definitions.get(key)
        if vb is None:
            diff_defs.append(key)
        elif isinstance(va, str) or isinstance(vb, str):
            if str(va) != str(vb):
                diff_defs.append(key)
        elif np.shape(va) != np.shape(vb) or not np.allclose(va, vb, rtol=1e-6, atol=1e-9):
            diff_defs.append(key)
    out["definitions"] = diff_defs

    la, lb = a.evaluate_labels(evolution="adc"), b.evaluate_labels(evolution="adc")
    out["labels"] = [k for k in set(la) | set(lb)
                     if k not in la or k not in lb or not np.array_equal(la[k], lb[k])]
    return out
