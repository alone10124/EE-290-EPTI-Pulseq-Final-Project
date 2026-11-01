"""Locally low-rank (LLR) proximal operator for the subspace coefficients (SigPy has none).

Replaces the LLR regularizer in the EPTI BART fork (``bart picsS -R L:7:7:lambda -b 6``; block thresholding in
``legacy/bart-0.8.00_EPTI/src/lowrank/lrthresh.c``). For each block of ``block x block`` voxels in one slice, the
``(block^2) x K`` matrix of coefficients is replaced by its singular-value soft-thresholded version. Blocks never mix
slices: the EPTI fork hard-codes this too (lrthresh.c L516, "avoid LLR along SMS slices").

Difference from the original: BART was run with ``-N`` (fully overlapping blocks). Here the block grid is shifted by
a random offset at each call, as in BART's default and in Trzasko & Manduca (ISMRM 2011) / Zhang et al. (MRM 2015).
That gives the same shift-invariance on average at a fraction of the cost.
"""

import numpy as np
import sigpy as sp


class LLR(sp.prox.Prox):
    def __init__(self, shape, lamda, block=6, randshift=True, seed=0):
        k, ns, ny, nx = shape
        if ny % block or nx % block:
            raise ValueError(f"image size {ny}x{nx} must be a multiple of the LLR block size {block}")
        self.lamda, self.block, self.randshift = lamda, block, randshift
        self.rng = np.random.default_rng(seed)
        super().__init__(shape)

    def _prox(self, alpha, input):
        xp = sp.get_array_module(input)
        k, ns, ny, nx = self.shape
        b = self.block
        shift = tuple(int(v) for v in self.rng.integers(0, b, 2)) if self.randshift else (0, 0)
        x = xp.roll(input, shift, axis=(-2, -1))
        # [K, s, ny/b, b, nx/b, b] -> [s, ny/b, nx/b, b*b, K]
        blocks = x.reshape(k, ns, ny // b, b, nx // b, b).transpose(1, 2, 4, 3, 5, 0).reshape(-1, b * b, k)
        u, s, vh = xp.linalg.svd(blocks, full_matrices=False)
        s = xp.maximum(s - self.lamda * alpha, 0)
        blocks = (u * s[:, None, :]) @ vh
        x = blocks.reshape(ns, ny // b, nx // b, b, b, k).transpose(5, 0, 1, 3, 2, 4).reshape(k, ns, ny, nx)
        return xp.roll(x, (-shift[0], -shift[1]), axis=(-2, -1))
