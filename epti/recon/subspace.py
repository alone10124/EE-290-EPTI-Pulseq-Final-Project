"""Temporal subspace basis for EPTI: the echo images are modeled as ``K`` coefficient maps times basis curves.

Port of ``legacy/libs/libs_EPTI_Recon/library_EPTI_recon/EPTI_gen_Subspace_Basis_GE.m`` and ``gen_GE_basis_T2B0.m``
(Zijing Dong, MGH; adapted from Jon Tamir's T2-shuffling basis code). Simulates ``exp(-TE / T2)`` for 91 T2 values
from 10 to 300 ms, takes the SVD, and keeps the first ``K`` left singular vectors. recon_EPTI_P1.m L98-100 then scales
them by ``S^(1/3)``. Reference: Dong et al., MRM 2020;84:2442 (EPTI with subspace reconstruction).

Note: the original builds the dictionary on times *relative to the start of the echo train* (``gen_GE_basis_T2B0`` is
called with ``T0 = 0``, i.e. t = ESP, 2 ESP, ..., N ESP), not on the absolute TEs. The absolute TEs are used only
for the B0 phase. Absolute times would weight each T2 curve by ``exp(-t0/T2)`` and change the leading subspace; the
MATLAB fixture test caught this.
"""

import numpy as np

T2_VALUES = np.concatenate([np.arange(10, 51, 1), np.arange(52, 101, 2), np.arange(102, 201, 5),
                            np.arange(220, 301, 20)]) / 1000.0


def imaging_echo_times(p):
    """TE of every lobe of the imaging echo train: TE + (t - N/2) * ESP for t = 1..N (EPTI_gen_Subspace_Basis_GE.m)."""
    n = p.n_lobes
    return p.te + (np.arange(1, n + 1) - n / 2) * p.esp


def basis_times(p):
    """Times the basis is built on: ESP, 2 ESP, ..., N ESP (gen_GE_basis_T2B0 with T0 = 0)."""
    return p.esp * np.arange(1, p.n_lobes + 1)


def basis(tes, k=4, t2_values=T2_VALUES):
    """Return ``(Phi, S)``: ``Phi [n_te, k]`` (scaled by S^(1/3)) and all singular values ``S``.

    Pass relative times (``basis_times``) to match the original."""
    x = np.exp(-np.asarray(tes)[:, None] / t2_values[None, :])
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    return u[:, :k] * s[:k] ** (1 / 3), s
