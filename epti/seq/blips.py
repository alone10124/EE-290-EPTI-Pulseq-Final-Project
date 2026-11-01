"""EPTI phase-encoding (Gy) and blipped-CAIPI (Gz) blips that straddle consecutive readout lobes."""

import pypulseq as pp


def epti_blip(channel, start_step, end_step, deltak, blip_dur, gx, system, state):
    """Triangular blip(s) to play inside one readout lobe ``gx``.

    Each ky (or kz) step is a triangle of area ``step * deltak`` and duration ``blip_dur``, centered on the boundary
    between two readout lobes. It is split at its peak: the rising half sits at the *end* of one lobe and the falling
    half at the *start* of the next.

    state 1 -> rising half of ``start_step`` only (last navigator lobe)
    state 2 -> falling half of ``start_step`` + rising half of ``end_step`` (a middle lobe)
    state 3 -> falling half of ``start_step`` only (last lobe of the train)

    Returns ``None`` for a zero step (e.g. Gz when MB = 1), so the caller can leave it out of the block.

    Ported from ``legacy/libs/utils/calc_EPTIGyBlip_WJG.m`` / ``calc_EPTIGzBlip_WJG.m`` (MGH, Jian Wu); the two
    MATLAB functions differ only in the channel. Uses ``pp.split_gradient_at``, ``pp.align`` and ``pp.add_gradients``
    in place of ``mr.splitGradientAt``, ``mr.align`` and ``mr.addGradients``.
    """

    def halves(step):
        tri = pp.make_trapezoid(channel, system=system, amplitude=step * deltak / blip_dur * 2,
                                rise_time=blip_dur / 2, flat_time=0)
        up, down = pp.split_gradient_at(tri, blip_dur / 2, system)
        up, down, _ = pp.align(right=up, left=[down, gx])
        return up, down

    if state == 1:
        return halves(start_step)[0] if start_step else None
    if state == 3:
        return halves(start_step)[1] if start_step else None
    if state == 2:
        parts = []
        if start_step:
            parts.append(halves(start_step)[1])
        if end_step:
            parts.append(halves(end_step)[0])
        if not parts:
            return None
        return parts[0] if len(parts) == 1 else pp.add_gradients(parts, system=system)
    raise ValueError(f"unknown blip state {state}")
