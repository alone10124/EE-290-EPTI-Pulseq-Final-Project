"""Numerical phantom for simulating multi-slice EPI scans with KomaMRI, and its ground truth.

Builds on KomaMRI's ``brain_phantom2D`` (BrainWeb-based, 0.5 mm native grid). For multi-slice EPI the default phantom alone
falls short in three ways, which this module fixes:

* All of its spins sit in one plane at z = 0, but a multi-slice scan excites slices at other z positions. We place one copy
  at each slice position (exactly at the slice center, so there are no through-slice effects). Slices with an odd
  global index are rotated by 90 degrees, so leakage between neighbouring slices is easy to see. A consequence: z
  gradient spoilers dephase nothing (every spin of a slice gets the same phase), so a fat-sat pulse doesn't suppress
  fat in the simulation (fat signal ratio with/without fat-sat: 1.00 here, 0.53 with spins spread over +-1 mm).
* A smooth synthetic B0 field map (Hz) per slice, so the simulation shows realistic EPI distortion.
* KomaMRI's simulator ignores the phantom's ``T2s`` field: KomaMRICore only uses T2 (checked in the KomaMRICore 0.14
  source: ``T2s`` appears only in its precompile example). A gradient-echo signal decays with T2*, so we set
  ``T2 := T2*``, i.e. intravoxel dephasing is modeled as a mono-exponential decay.

Fat: the phantom's fat offset is hard-coded for 1.5 T. By default fat spins are treated as water (``fat="water"``),
which keeps validation simple; ``fat="3T"`` puts them at -3.4 ppm at the scanner's B0 so a fat-sat pulse has something
to suppress.

Sign convention: KomaMRI precesses an off-resonance of +df Hz as ``exp(-i 2 pi df t)``, so the phase of the simulated
signal *decreases* with TE where the phantom's B0 map is positive.

Image-grid convention (shared with the recon): maps are [slice, y, x]; pixel ``i`` of an ``n``-pixel axis sits at
``(i - n // 2) * fov / n``, the pixel grid of a centered FFT (``sigpy.ifft``).

Ground truth on disk: sim/run_sim.py saves the spins themselves (``save_truth`` -> ``truth.npz``), with each spin's
receive-coil sensitivity, so the file doesn't depend on any sequence or recon grid. A recon script calls ``reference``
with its own FOV, matrix and per-spin signal weighting (its contrast: TE, TR, fat-sat, ...) to get the true maps and the
image an ideal Cartesian scan with that contrast would give.

(First written for the EPTI port; see the ``epti`` branch.)
"""

from dataclasses import dataclass

import numpy as np

import mrtools.sim  # noqa: F401  (sets the juliacall environment variables before komamripy loads)

GAMMA_HZ_PER_T = 42.576e6


SPIN_FIELDS = ("x", "y", "z", "rho", "t1", "t2", "t2s", "df_hz", "slice_id", "is_fat")


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
    slice_id: np.ndarray        # index into the slice_positions passed to make_phantom
    is_fat: np.ndarray = None   # bool per spin (None: no fat information)

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


def make_phantom(slice_positions, slice_indices=None, ss=4, fat="water", b0_scale=1.0, B0=3.0, relaxation=True):
    """Build the multi-slice phantom.

    Parameters
    ----------
    slice_positions : z positions (m) of the slices to fill, e.g. ``(-0.0315, 0.0315)``.
    slice_indices : the global index of each slice (the sequence's ``SliceIndices``). Odd indices get the rotated
        phantom and the second field map, so a slice looks the same whichever subset is simulated. Default: 0, 1, ...
    ss : KomaMRI subsampling of the 0.5 mm BrainWeb grid (4 -> 2 mm spin spacing, ~6.5k spins per slice).
    fat : "water" (fat spins get no chemical shift) or "3T" (-3.4 ppm at ``B0``).
    b0_scale : multiplies the synthetic field map (0 turns B0 inhomogeneity off).
    relaxation : False sets T1, T2 and T2* to 1000 s (no decay), for debugging.
    """
    import komamripy as km

    if slice_indices is None:
        slice_indices = range(len(slice_positions))
    base = km.brain_phantom2D(ss=ss)
    x0, y0 = np.asarray(base.x, dtype=float), np.asarray(base.y, dtype=float)
    rho0, t10 = np.asarray(base.ρ, dtype=float), np.asarray(base.T1, dtype=float)
    t20, t2s0 = np.asarray(base.T2, dtype=float), np.asarray(base.T2s, dtype=float)
    if not relaxation:
        t10, t20, t2s0 = (np.full_like(t10, 1e3) for _ in range(3))
    is_fat = np.asarray(base.Δw, dtype=float) != 0
    fat_hz = -3.4e-6 * GAMMA_HZ_PER_T * B0 if fat == "3T" else 0.0

    parts = []
    for k, (z, index) in enumerate(zip(slice_positions, slice_indices)):
        # odd slices are rotated by 90 degrees: (x, y) -> (y, -x)
        x, y = (x0, y0) if index % 2 == 0 else (y0, -x0)
        df = b0_scale * _field_map(x, y, index) + np.where(is_fat, fat_hz, 0.0)
        parts.append(SpinPhantom(x=x, y=y, z=np.full_like(x, z), rho=rho0, t1=t10, t2=t20, t2s=t2s0,
                                 df_hz=df, slice_id=np.full(len(x), k), is_fat=is_fat))
    return SpinPhantom(*(np.concatenate([getattr(p, f) for p in parts]) for f in SPIN_FIELDS))


def to_koma(ph: SpinPhantom, name="multi-slice brain"):
    """KomaMRI ``Phantom``. T2 is set to T2* because KomaMRI's simulator does not use the T2s field (see module doc)."""
    import komamripy as km
    from juliacall import Main as jl

    vec = jl.seval("v -> collect(Float64, v)")
    return km.Phantom(name=name, x=vec(ph.x), y=vec(ph.y), z=vec(ph.z), ρ=vec(ph.rho), T1=vec(ph.t1),
                      T2=vec(ph.t2s), T2s=vec(ph.t2s), Δw=vec(2 * np.pi * ph.df_hz))


def coil_rss(scanner, ph: SpinPhantom):
    """Root-sum-of-squares receive sensitivity ``sqrt(sum_c |S_c|^2)`` at each spin, scaled to a median of 1.

    Uses KomaMRI's own coil model (``get_sens`` of the scanner's receiver), i.e. exactly what the simulation used.
    Both RSS and ESPIRiT/SENSE-1 coil combination return the object multiplied by this (ESPIRiT's maps are normalized
    to unit length across coils), which is why images are brighter near the coils.
    """
    from juliacall import Main as jl

    sens = jl.seval("(rx, x, y, z) -> get_sens(rx, collect(Float64, x), collect(Float64, y), collect(Float64, z))")(
        scanner.receiver, ph.x, ph.y, ph.z)
    rss = np.sqrt((np.abs(np.asarray(sens)) ** 2).sum(axis=1))
    return rss / np.median(rss)


def pixel_edges(fov, n):
    """Bin edges whose centers are the centered-FFT pixel positions ``(i - n // 2) * fov / n``."""
    return (np.arange(n + 1) - n // 2 - 0.5) * fov / n


def ground_truth_maps(ph: SpinPhantom, fov, n, n_slices):
    """Proton density, T2* (s) and B0 (Hz) on the reconstruction grid, one map per slice (``slice_id``).

    Spins are binned into ``n x n`` voxels over ``fov`` on the FFT pixel grid (see module doc). Proton density is the
    mean over the spins in each voxel; a plain sum would turn the spin-grid/voxel-grid mismatch into a checkerboard.
    T2* and B0 are rho-weighted voxel averages. Map axes are [slice, y, x].
    """
    bins = pixel_edges(fov[1], n[1]), pixel_edges(fov[0], n[0])
    out = {k: np.zeros((n_slices, n[1], n[0])) for k in ("rho", "t2s", "b0_hz")}
    for s in range(n_slices):
        m = ph.slice_id == s
        count, _, _ = np.histogram2d(ph.y[m], ph.x[m], bins=bins)
        w, _, _ = np.histogram2d(ph.y[m], ph.x[m], bins=bins, weights=ph.rho[m])
        for key, val in (("t2s", ph.t2s[m]), ("b0_hz", ph.df_hz[m])):
            num, _, _ = np.histogram2d(ph.y[m], ph.x[m], bins=bins, weights=ph.rho[m] * val)
            out[key][s] = np.divide(num, w, out=np.zeros_like(num), where=w > 0)
        out["rho"][s] = np.divide(w, count, out=np.zeros_like(w), where=count > 0)
    return out


def bandlimited_image(x, y, weights, fov, n, chunk=20000):
    """Image of point spins as an ideal Cartesian scan would see it: the exact Fourier transform of the spins at the
    ``n[1] x n[0]`` k-space grid points ``k = (m - n // 2) / fov``, inverse-FFT'd onto the pixel grid.

    So it has the same resolution and Gibbs ringing as a fully sampled recon, unlike a box-average of the spins.
    Scaled so that a uniform object with spin weight ``w`` comes out as ``w`` (spins per voxel divided out).
    Returns a complex [y, x] image.
    """
    ky = (np.arange(n[1]) - n[1] // 2) / fov[1]
    kx = (np.arange(n[0]) - n[0] // 2) / fov[0]
    k = np.zeros((n[1], n[0]), complex)
    for i in range(0, len(x), chunk):
        s = slice(i, i + chunk)
        k += (np.exp(-2j * np.pi * np.outer(ky, y[s])) * weights[s]) @ np.exp(-2j * np.pi * np.outer(x[s], kx))
    img = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(k)))
    spacing = np.diff(np.unique(np.round(x, 7))).min()                  # spin grid spacing (m)
    return img / (fov[0] / n[0] * fov[1] / n[1] / spacing ** 2)


def save_truth(path, ph: SpinPhantom, slice_positions, coil_rss):
    """Write the ground truth of a simulation: every spin (``SpinPhantom`` fields), its RSS receive sensitivity
    ``coil_rss`` (see ``coil_rss``) and the z position of each ``slice_id``. Read it back with ``load_truth``."""
    is_fat = ph.is_fat if ph.is_fat is not None else np.zeros(len(ph), bool)
    fields = {f: getattr(ph, f) for f in SPIN_FIELDS if f != "is_fat"}
    np.savez(path, **fields, is_fat=is_fat, coil_rss=coil_rss, slice_positions=np.asarray(slice_positions))


def load_truth(path):
    """``(SpinPhantom, coil_rss, slice_positions)`` from a file written by ``save_truth``, or ``None`` if the file is
    missing or in the old format (maps on a grid instead of spins)."""
    try:
        d = np.load(path)
    except FileNotFoundError:
        return None
    if "x" not in d:
        return None
    return SpinPhantom(**{f: d[f] for f in SPIN_FIELDS}), d["coil_rss"], d["slice_positions"]


def reference(path, fov, matrix, weight):
    """Ground truth for a recon on a ``fov`` / ``matrix`` grid (``matrix`` = (nx, ny)), from ``truth.npz`` at ``path``.

    ``weight(spins, coil_rss)`` returns each spin's signal in an ideal scan with the recon's contrast, e.g.
    ``spins.rho * np.exp(-te / spins.t2s) * coil_rss`` for a gradient echo at ``te``. Returns ``None`` if there is no
    usable truth file, else a dict with ``slice_positions``, the voxel maps ``rho``, ``t2s`` (s), ``b0_hz``
    (``ground_truth_maps``) and ``expected``: |band-limited image| of the weights (``bandlimited_image``), all
    [slice, y, x] with slices in the order of ``slice_positions``.
    """
    loaded = load_truth(path)
    if loaded is None:
        return None
    spins, rss, positions = loaded
    out = ground_truth_maps(spins, fov, matrix, len(positions))
    w = weight(spins, rss)
    out["expected"] = np.stack([np.abs(bandlimited_image(spins.x[m], spins.y[m], w[m], fov, matrix))
                                for m in (spins.slice_id == k for k in range(len(positions)))])
    out["slice_positions"] = positions
    return out
