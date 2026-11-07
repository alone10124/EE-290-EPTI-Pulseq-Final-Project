"""Run a .seq file through KomaMRI and return the raw signal as a NumPy array."""

import time

import numpy as np

from mrtools.sim import phantom  # noqa: F401  (sets the juliacall environment variables before komamripy loads)


def simulate(seq_path, koma_phantom, scanner, gpu=False):
    """Simulate ``seq_path`` and return ``(signal [n_adc_samples, n_coils], seconds)``.

    ``signal`` holds every ADC sample in acquisition order (complex). ``koma_phantom`` and ``scanner`` come from
    ``phantom.to_koma`` and ``phantom.make_scanner``. ``gpu=True`` runs on an NVIDIA GPU (CUDA); otherwise KomaMRI
    uses all CPU threads.
    """
    import komamripy as km
    from juliacall import Main as jl

    if gpu:
        km.load_cuda()
    seq = km.read_seq(str(seq_path))
    params = jl.seval('(gpu) -> Dict{String,Any}("return_type" => "mat", "gpu" => gpu)')(gpu)
    t0 = time.time()
    signal = np.asarray(km.simulate(koma_phantom, seq, scanner, sim_params=params))
    return signal.reshape(signal.shape[0], -1), time.time() - t0
