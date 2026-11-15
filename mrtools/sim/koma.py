"""Run a .seq file through KomaMRI and return the raw signal as a NumPy array."""

import time

import numpy as np

import mrtools.sim  # noqa: F401  (sets the juliacall environment variables before komamripy loads)


def simulate(seq_path, koma_phantom, scanner, threads=None):
    """Simulate ``seq_path`` on the CPU and return ``(signal [n_adc_samples, n_coils], seconds)``.

    ``signal`` holds every ADC sample in acquisition order (complex). ``koma_phantom`` comes from ``phantom.to_koma``
    and ``scanner`` from ``HardwareProfile.make_scanner``.

    ``threads``: how many parts KomaMRI splits the phantom's spins into, one per thread (its ``Nthreads`` setting);
    ``None`` = all Julia threads (see ``mrtools.sim``). Each part keeps its own complex128 copy of the whole signal,
    so memory grows by ``n_adc_samples x n_coils x 16 bytes`` per thread (2.1M samples x 8 coils: 270 MB per thread).
    More parts than Julia threads costs memory without running faster.
    """
    import komamripy as km
    from juliacall import Main as jl

    seq = km.read_seq(str(seq_path))
    params = jl.seval('() -> Dict{String,Any}("return_type" => "mat")')()
    if threads is not None:
        jl.seval("(d, n) -> d[\"Nthreads\"] = n")(params, int(threads))
    t0 = time.time()
    signal = np.asarray(km.simulate(koma_phantom, seq, scanner, sim_params=params))
    return signal.reshape(signal.shape[0], -1), time.time() - t0


def julia_threads():
    """Number of Julia threads KomaMRI can use."""
    from juliacall import Main as jl

    return int(jl.seval("Threads.nthreads()"))
