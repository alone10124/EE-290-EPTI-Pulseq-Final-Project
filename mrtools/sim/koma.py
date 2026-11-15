"""Run a .seq file through KomaMRI and return the raw signal as a NumPy array."""

import time

import numpy as np

import mrtools.sim  # noqa: F401  (sets the juliacall environment variables before komamripy loads)


def simulate(seq_path, koma_phantom, scanner):
    """Simulate ``seq_path`` on the CPU and return ``(signal [n_adc_samples, n_coils], seconds)``.

    ``signal`` holds every ADC sample in acquisition order (complex). ``koma_phantom`` comes from ``phantom.to_koma``
    and ``scanner`` from ``HardwareProfile.make_scanner``. KomaMRI uses all Julia threads (see ``mrtools.sim``).
    """
    import komamripy as km
    from juliacall import Main as jl

    seq = km.read_seq(str(seq_path))
    params = jl.seval('() -> Dict{String,Any}("return_type" => "mat")')()
    t0 = time.time()
    signal = np.asarray(km.simulate(koma_phantom, seq, scanner, sim_params=params))
    return signal.reshape(signal.shape[0], -1), time.time() - t0


def julia_threads():
    """Number of Julia threads KomaMRI can use."""
    from juliacall import Main as jl

    return int(jl.seval("Threads.nthreads()"))
