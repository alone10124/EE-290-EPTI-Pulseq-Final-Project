"""Numerical phantom and scanner for simulating multi-slice EPI scans with KomaMRI.

Builds on KomaMRI's ``brain_phantom2D`` (BrainWeb-based, 0.5 mm native grid). For multi-slice EPI the default phantom alone
falls short in three ways, which this module fixes:

* All of its spins sit in one plane at z = 0, but a multi-slice scan excites slices at other z positions. We place one copy
  at each slice position, every other one rotated by 90 degrees so leakage between slices is easy to see.
* A smooth synthetic B0 field map (Hz) per slice, so the simulation shows realistic EPI distortion.
* KomaMRI's simulator ignores the phantom's ``T2s`` field (KomaMRICore only uses T2). A gradient-echo signal decays
  with T2*, so we set ``T2 := T2*``, i.e. intravoxel dephasing is modeled as a mono-exponential decay.

Fat: the phantom's fat offset is hard-coded for 1.5 T. By default fat spins are treated as water (``fat="water"``),
which keeps validation simple; ``fat="3T"`` puts them at -3.4 ppm at the scanner's B0 so a fat-sat pulse has something
to suppress.

Sign convention: KomaMRI precesses an off-resonance of +df Hz as ``exp(-i 2 pi df t)``, so the phase of the simulated
signal *decreases* with TE where the phantom's B0 map is positive.

Returns plain NumPy arrays alongside the KomaMRI objects, so the ground truth can be compared with reconstructions.
(First written for the EPTI port; see the ``epti`` branch.)
"""

import os
from dataclasses import dataclass

import numpy as np

# KomaMRI/juliacall settings must be in place before the first `import komamripy` (see README).
os.environ.setdefault("PYTHON_JULIACALL_THREADS", "auto")
os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")

GAMMA_HZ_PER_T = 42.576e6


@dataclass
class SpinPhantom:
    """Spin-level description, in SI units, with Δw given in Hz (KomaMRI uses rad/s)."""
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    rho: np.ndarray
    t1: np.ndarray
    t2: np.ndarray
    t2s: np.ndarray
    df_hz: np.ndarray
    slice_id: np.ndarray        # index into slice_positions

    def __len__(self):
        return len(self.x)


def _field_map(x, y, slice_index):
    """Smooth synthetic B0 map (Hz): a frontal 'sinus' hot spot plus a weak linear gradient, different per slice."""
    if slice_index % 2 == 0:
        hot = 40.0 * np.exp(-((x - 0.0) ** 2 + (y - 0.07) ** 2) / (2 * 0.02 ** 2))
        lin = 60.0 * x
    else:
        hot = -30.0 * np.exp(-((x + 0.04) ** 2 + (y + 0.03) ** 2) / (2 * 0.025 ** 2))
        lin = 50.0 * y
    return hot + lin


def make_phantom(slice_positions, ss=4, fat="water", b0_scale=1.0, B0=3.0):
    """Build the multi-slab phantom.

    Parameters
    ----------
    slice_positions : z positions (m) of the slices to fill, e.g. ``(-0.0315, 0.0315)``.
    ss : KomaMRI subsampling of the 0.5 mm BrainWeb grid (4 -> 2 mm spin spacing, ~6.5k spins per slice).
    fat : "water" (fat spins get no chemical shift) or "3T" (-3.4 ppm at ``B0``).
    b0_scale : multiplies the synthetic field map (0 turns B0 inhomogeneity off).
    """
    import komamripy as km

    base = km.brain_phantom2D(ss=ss)
    x0, y0 = np.asarray(base.x, dtype=float), np.asarray(base.y, dtype=float)
    rho0, t10 = np.asarray(base.ρ, dtype=float), np.asarray(base.T1, dtype=float)
    t20, t2s0 = np.asarray(base.T2, dtype=float), np.asarray(base.T2s, dtype=float)
    is_fat = np.asarray(base.Δw, dtype=float) != 0
    fat_hz = -3.4e-6 * GAMMA_HZ_PER_T * B0 if fat == "3T" else 0.0

    parts = []
    for k, z in enumerate(slice_positions):
        # every other slice is rotated by 90 degrees: (x, y) -> (y, -x)
        x, y = (x0, y0) if k % 2 == 0 else (y0, -x0)
        df = b0_scale * _field_map(x, y, k) + np.where(is_fat, fat_hz, 0.0)
        parts.append(SpinPhantom(x=x, y=y, z=np.full_like(x, z), rho=rho0, t1=t10, t2=t20, t2s=t2s0,
                                 df_hz=df, slice_id=np.full(len(x), k)))
    return SpinPhantom(*(np.concatenate([getattr(p, f) for p in parts])
                         for f in ("x", "y", "z", "rho", "t1", "t2", "t2s", "df_hz", "slice_id")))


def to_koma(ph: SpinPhantom, name="multi-slab brain"):
    """KomaMRI ``Phantom``. T2 is set to T2* because KomaMRI's simulator does not use the T2s field (see module doc)."""
    import komamripy as km
    from juliacall import Main as jl

    vec = jl.seval("v -> collect(Float64, v)")
    return km.Phantom(name=name, x=vec(ph.x), y=vec(ph.y), z=vec(ph.z), ρ=vec(ph.rho), T1=vec(ph.t1),
                      T2=vec(ph.t2s), T2s=vec(ph.t2s), Δw=vec(2 * np.pi * ph.df_hz))


def make_scanner(hw, ncoils=8, coil_radius=0.15, coil_half_length=0.12):
    """KomaMRI ``Scanner`` with the hardware profile's limits and a multi-channel receive array.

    ``BirdcageCoilSens`` gives one receive channel per rung (so it's an array, despite the name). Several coils are
    needed for parallel imaging and for simultaneous multi-slice methods. Its sensitivity is symmetric in z, so two
    slices at +-z see the same coil maps.
    ``hw`` is a ``mrtools.hardware.HardwareProfile``.
    """
    from juliacall import Main as jl

    jl.seval("using KomaMRI")
    limits = jl.seval("(B0, Gmax, Smax, adc, dur, gr, rf, rd, rfd, adcd) -> HardwareLimits(B0=B0, Gmax=Gmax, Smax=Smax, "
                      "ADC_Δt=adc, DUR_Δt=dur, GR_Δt=gr, RF_Δt=rf, RF_ring_down_time=rd, RF_dead_time=rfd, "
                      "ADC_dead_time=adcd)")(
        hw.B0, hw.max_grad * 1e-3, hw.max_slew, hw.adc_raster_time, hw.block_duration_raster, hw.grad_raster_time,
        hw.rf_raster_time, hw.rf_ringdown_time, hw.rf_dead_time, hw.adc_dead_time)
    receiver = jl.seval("(n, r, L) -> BirdcageCoilSens(ncoils=n, radius=r, L=L)")(ncoils, coil_radius, coil_half_length)
    return jl.seval("(lim, rx) -> Scanner(limits=lim, receiver=rx)")(limits, receiver)


def ground_truth_maps(ph: SpinPhantom, fov, n, slice_positions):
    """Proton density, T2* (s) and B0 (Hz) on the reconstruction grid, one map per slice.

    Spins are binned into ``n x n`` voxels over ``fov`` (centered at 0). Proton density is the mean over the spins in
    each voxel; a plain sum would turn the spin-grid/voxel-grid mismatch into a checkerboard. T2* and B0 are
    rho-weighted voxel averages.
    Map axes are [slice, y, x] with x and y increasing along the array axes.
    """
    edges = np.linspace(-fov[0] / 2, fov[0] / 2, n[0] + 1), np.linspace(-fov[1] / 2, fov[1] / 2, n[1] + 1)
    out = {k: np.zeros((len(slice_positions), n[1], n[0])) for k in ("rho", "t2s", "b0_hz")}
    for s in range(len(slice_positions)):
        m = ph.slice_id == s
        bins = (edges[1], edges[0])
        count, _, _ = np.histogram2d(ph.y[m], ph.x[m], bins=bins)
        w, _, _ = np.histogram2d(ph.y[m], ph.x[m], bins=bins, weights=ph.rho[m])
        for key, val in (("t2s", ph.t2s[m]), ("b0_hz", ph.df_hz[m])):
            num, _, _ = np.histogram2d(ph.y[m], ph.x[m], bins=bins, weights=ph.rho[m] * val)
            out[key][s] = np.divide(num, w, out=np.zeros_like(num), where=w > 0)
        out["rho"][s] = np.divide(w, count, out=np.zeros_like(w), where=count > 0)
    out["x_edges"], out["y_edges"] = edges
    return out
