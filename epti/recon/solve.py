"""Subspace reconstruction of one slice group: ADMM with locally-low-rank regularization (SigPy).

Replaces ``bart picsS -m -u 0.01 -i 70 -C 8 -R L:7:7:0.001 -b 6 -w <scaling> -S`` in
``legacy/libs/libs_EPTI_Recon/library_EPTI_recon/SubspaceRecon_Process_v5.m`` L68-69 with
``sigpy.app.LinearLeastSquares(solver="ADMM")`` and :class:`epti.recon.llr.LLR`.

As with BART's ``-w``, the data are divided by an intensity scale before solving, so ``lamda`` is relative to image
intensity, and the result is scaled back (``-S``). The scale follows ``calculate_imagescaling_EPTIv2.m``: the 90th
percentile of the coil-combined zero-filled image, or its maximum if the top of the histogram is unusually wide.
"""

from dataclasses import dataclass

import numpy as np
import sigpy as sp

from epti.recon.llr import LLR


@dataclass
class SolverOptions:
    lamda: float = 1e-3          # LLR weight (recon_EPTI_P1.m L52); 0 -> plain least squares (CG)
    rho: float = 1e-2            # ADMM penalty (L51)
    max_iter: int = 70           # ADMM iterations (SubspaceRecon_Process_v5.m L68)
    max_cg_iter: int = 8         # CG iterations per ADMM step
    block: int = 6               # LLR block size (L53)
    device: int = -1             # -1 = CPU; >= 0 = CUDA GPU via CuPy
    show_pbar: bool = False


def intensity_scale(y):
    """Scale of the zero-filled data ``y [t, coil, ky, kx]`` (calculate_imagescaling_EPTIv2.m)."""
    img = sp.ifft(y, axes=(-2, -1))
    sos = np.sqrt((np.abs(img) ** 2).sum(axis=1)).mean(axis=0).ravel()
    sos.sort()
    p100, p90, p50 = sos[-1], sos[int(round(0.9 * len(sos))) - 1], sos[int(round(0.5 * len(sos))) - 1]
    return p90 if (p100 - p90) < 2 * (p90 - p50) else p100


def solve(A, y, coeff_shape, opts: SolverOptions):
    """Return the coefficient maps ``[K, s, y, x]`` minimizing ``||A c - y||^2 + lamda * LLR(c)``."""
    device = sp.Device(opts.device)
    xp = device.xp
    scale = intensity_scale(y)
    with device:
        A_dev = A
        y_dev = sp.to_device(y / scale, device)
        if opts.lamda > 0:
            app = sp.app.LinearLeastSquares(A_dev, y_dev, proxg=LLR(coeff_shape, opts.lamda, opts.block),
                                            solver="ADMM", rho=opts.rho, max_iter=opts.max_iter,
                                            max_cg_iter=opts.max_cg_iter, show_pbar=opts.show_pbar)
        else:
            app = sp.app.LinearLeastSquares(A_dev, y_dev, solver="ConjugateGradient", max_iter=opts.max_iter,
                                            show_pbar=opts.show_pbar)
        coeff = app.run()
    return sp.to_device(coeff, sp.cpu_device) * scale
