"""
make_cyclic_epi_ge.py -- cyclically shifted multi-shot EPI for a GE Discovery MR750w 3T
(Pulseq -> pulseg -> pge2), multi-slice, with recorded per-sample timing.

Protocol defaults
    FOV 23 cm in-plane, 64 x 64 (3.6 mm), 24 x 5 mm interleaved slices = 12 cm in z
    Short phase-encode blips; one long blip for the bottom->top wrap (non-uniform
    echo spacing).  The real time and k-space position of every ADC sample are saved,
    and recon_model_based.py uses them directly (no uniform-ESP assumption).

Shot s (0-based) starts at ky row s and wraps:  row(s, j) = (s + j) mod N.

GE notes
  * MR750w spec: 44 mT/m, 200 T/m/s.  Defaults are derated (40 mT/m, 150 T/m/s)
    for EPI/PNS headroom; for oblique slices divide by sqrt(3).  Run pge2.check()
    with your site's PNS model before scanning.
  * Every segment (one slice of one shot) begins with a TRID label.  Segment layout
    depends on where the wrap blip falls, so TRID = s + 1 for acquired shot s and
    TRID = N + 1 for the dummy (no-ADC) shots.  Slices reuse the same TRIDs.
  * Raster/dead-time settings follow the Pulseq-on-GE example (adcRasterTime 2 us,
    rfRingdownTime 60 us, adcDeadTime 40 us).  Check against your pge2 version.

Outputs
    <out>.seq               multi-slice sequence for the scanner (with TRID labels)
    <out>_koma1slice.seq    the same shots/TR for one slice at z = 0, no labels (KomaMRI)
    <out>_meta.npz          ordering, slice table, per-sample times and k-space positions
    <out>_slices.csv        slice positions for simulate_koma_brain3d.jl

Usage
    python make_cyclic_epi_ge.py
    python make_cyclic_epi_ge.py --N 32 --slices 4        # quick test
    python make_cyclic_epi_ge.py --grad-delay 2e-6         # measured gradient delay [s]
"""
import argparse
import copy

import numpy as np
import pypulseq as pp


def ge_system(max_grad=40, max_slew=150):
    return pp.Opts(max_grad=max_grad, grad_unit="mT/m", max_slew=max_slew, slew_unit="T/m/s",
                   rf_dead_time=100e-6, rf_ringdown_time=60e-6, adc_dead_time=40e-6,
                   adc_raster_time=2e-6, rf_raster_time=2e-6, grad_raster_time=10e-6,
                   block_duration_raster=10e-6, B0=3.0)


def scaled(grad, scale):
    g = copy.deepcopy(grad)
    g.amplitude *= scale
    g.area *= scale
    g.flat_area *= scale
    return g


def interleaved(n):
    return list(range(0, n, 2)) + list(range(1, n, 2))


def build(system, N, fov, slice_z, slice_thk, shots, n_dummy, TR, flip_deg, dwell, labels):
    """Return (seq, info).  slice_z: positions in acquisition order."""
    seq = pp.Sequence(system)
    n_sl = len(slice_z)

    rf, gz, _ = pp.make_sinc_pulse(flip_angle=np.deg2rad(flip_deg), system=system, duration=3e-3,
                                   slice_thickness=slice_thk, apodization=0.5,
                                   time_bw_product=4, delay=system.rf_dead_time, return_gz=True)
    rf_center = pp.calc_rf_center(rf)[0] + rf.delay

    dk = 1.0 / fov
    readout_time = N * dwell
    # flat time a multiple of 2 gradient rasters -> the centred ADC start lands on the 2 us ADC raster
    flat_time = np.ceil(readout_time / (2 * system.grad_raster_time)) * 2 * system.grad_raster_time
    gx = pp.make_trapezoid("x", system=system, amplitude=N * dk / readout_time, flat_time=flat_time)
    # samples centred on the flat top, at +/-(n + 1/2) dk -> odd/even echoes share kx positions
    adc_delay = gx.rise_time + flat_time / 2 - readout_time / 2
    adc_delay = np.round(adc_delay / system.adc_raster_time) * system.adc_raster_time
    adc = pp.make_adc(num_samples=N, dwell=dwell, delay=adc_delay, system=system)

    gy_step = pp.make_trapezoid("y", system=system, area=dk)
    gy_wrap = pp.make_trapezoid("y", system=system, area=-(N - 1) * dk)

    gx_pre0 = pp.make_trapezoid("x", system=system, area=-gx.area / 2)
    gy_pre0 = pp.make_trapezoid("y", system=system, area=(N / 2) * dk)
    gz_reph0 = pp.make_trapezoid("z", system=system, area=-gz.area / 2)
    pre_time = max(pp.calc_duration(g) for g in (gx_pre0, gy_pre0, gz_reph0))
    gx_pre = pp.make_trapezoid("x", system=system, area=-gx.area / 2, duration=pre_time)
    gz_reph = pp.make_trapezoid("z", system=system, area=-gz.area / 2, duration=pre_time)
    gy_pre_unit = pp.make_trapezoid("y", system=system, area=-(N / 2) * dk, duration=pre_time)

    gx_spoil = pp.make_trapezoid("x", system=system, area=2 * N * dk)
    gz_spoil = pp.make_trapezoid("z", system=system, area=4 / slice_thk)

    def segment_blocks(s, acquire):
        """List of block-argument tuples for shot s (one slice), excluding RF."""
        blocks = []
        scale = -(s - N / 2) / (N / 2)
        pre = [gx_pre, gz_reph] + ([scaled(gy_pre_unit, scale)] if scale != 0 else [])
        blocks.append(tuple(pre))
        for j in range(N):
            g = scaled(gx, (-1) ** j)                 # polarity depends on echo index only
            blocks.append((g, adc) if acquire else (g, pp.make_delay(pp.calc_duration(g, adc))))
            if j < N - 1:
                blocks.append((gy_wrap,) if (s + j) % N == N - 1 else (gy_step,))
        blocks.append((gx_spoil, gz_spoil))
        return blocks

    # longest segment sets the minimum per-slice TR
    seg_dur = {s: pp.calc_duration(rf, gz) + sum(pp.calc_duration(*b) for b in segment_blocks(s, True))
               for s in set(shots)}
    tr_slice_min = max(seg_dur.values())
    tr_slice = max(TR / n_sl, tr_slice_min)
    tr_slice = np.ceil(tr_slice / system.block_duration_raster) * system.block_duration_raster

    order = [(shots[0], False)] * n_dummy + [(s, True) for s in shots]
    for s, acq in order:
        trid = (s + 1) if acq else (N + 1)
        for z in slice_z:
            rf.freq_offset = gz.amplitude * z
            rf.phase_offset = -2 * np.pi * rf.freq_offset * pp.calc_rf_center(rf)[0]
            first = [rf, gz] + ([pp.make_label(label="TRID", type="SET", value=int(trid))] if labels else [])
            seq.add_block(*first)
            for b in segment_blocks(s, acq):
                seq.add_block(*b)
            pad = tr_slice - seg_dur[s]
            if pad > 0:
                seq.add_block(pp.make_delay(np.round(pad / system.block_duration_raster)
                                            * system.block_duration_raster))

    info = dict(tr_slice=tr_slice, TR=tr_slice * n_sl, d_step=pp.calc_duration(gy_step),
                d_wrap=pp.calc_duration(gy_wrap), rf_center=rf_center,
                gx_amp_mTm=gx.amplitude / system.gamma * 1e3)
    return seq, info


def main(a):
    system = ge_system(a.gmax, a.smax)
    N, fov = a.N, a.fov
    shots = list(range(N)) if a.shots is None else a.shots
    n_sl = a.slices
    # slice centres, spatial order, spanning a.coverage in z
    pitch = a.coverage / n_sl
    z_spatial = (np.arange(n_sl) - (n_sl - 1) / 2) * pitch
    acq_order = interleaved(n_sl)                       # spatial indices in acquisition order
    slice_z = z_spatial[acq_order]

    # 1) multi-slice scanner version
    seq, info = build(system, N, fov, slice_z, a.thk, shots, a.dummy, a.TR, a.flip, a.dwell, labels=True)
    ok, rep = seq.check_timing()
    if not ok:
        print("Timing check FAILED:", "".join(str(r) for r in rep[:10]))
    seq.set_definition("FOV", [fov, fov, a.coverage])
    seq.set_definition("Name", "cyclic_epi_ge")
    seq.set_definition("SliceThickness", a.thk)
    seq.write(a.out + ".seq")

    # 2) single slice at z = 0 with the SAME TR (same T1 history) for KomaMRI; also
    #    used to record per-sample timing, which is identical for every slice
    seq1, _ = build(system, N, fov, [0.0], a.thk, shots, a.dummy, info["TR"], a.flip, a.dwell, labels=False)
    seq1.set_definition("FOV", [fov, fov, a.thk])
    seq1.set_definition("Name", "cyclic_epi_koma")
    seq1.write(a.out + "_koma1slice.seq")

    # 3) record the real sampling times and k-space positions (one slice; identical per slice)
    k_adc, _, t_exc, _, t_adc = seq1.calculate_kspace(trajectory_delay=a.grad_delay)
    n_acq = len(shots)
    t_adc = np.asarray(t_adc).reshape(n_acq, N, N)
    t_exc = np.asarray(t_exc)
    t_rel = t_adc - t_exc[np.searchsorted(t_exc, t_adc, side="right") - 1]
    te = t_rel.mean(axis=2)
    rows = np.array([[(s + j) % N for j in range(N)] for s in shots])

    np.savez(a.out + "_meta.npz", N=N, fov=fov, shots=np.array(shots), rows=rows, te=te,
             t_rel=t_rel, k_adc=np.asarray(k_adc)[:, :n_acq * N * N], uniform_esp=False,
             n_slices=n_sl, slice_thk=a.thk, z_spatial=z_spatial, acq_order=np.array(acq_order),
             TR=info["TR"], grad_delay=a.grad_delay, d_step=info["d_step"], d_wrap=info["d_wrap"])

    # slice table for simulate_koma_brain3d.jl (spatial order)
    np.savetxt(a.out + "_slices.csv", np.c_[np.arange(n_sl), z_spatial], delimiter=",",
               header="spatial_index,z_m", comments="", fmt=["%d", "%.6f"])

    esp = np.diff(te, axis=1)
    dur = seq.duration()[0]
    print(f"GE MR750w cyclic EPI: N={N}, FOV={fov*100:.0f} cm ({fov/N*1e3:.2f} mm), "
          f"{n_sl} x {a.thk*1e3:.1f} mm slices over {a.coverage*100:.0f} cm, {len(shots)} shots")
    print(f"  readout {info['gx_amp_mTm']:.1f} mT/m, BW/px {1/(N*a.dwell):.0f} Hz, "
          f"normal blip {info['d_step']*1e3:.2f} ms, wrap blip {info['d_wrap']*1e3:.2f} ms")
    print(f"  ESP {esp.min()*1e3:.3f}-{esp.max()*1e3:.3f} ms, TE {te.min()*1e3:.2f}-{te.max()*1e3:.2f} ms, "
          f"TR {info['TR']:.3f} s, scan time {dur/60:.1f} min")
    print(f"  wrote {a.out}.seq, {a.out}_koma1slice.seq, {a.out}_meta.npz, {a.out}_slices.csv")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=64)
    p.add_argument("--fov", type=float, default=0.23)
    p.add_argument("--slices", type=int, default=24)
    p.add_argument("--coverage", type=float, default=0.12, help="z coverage [m]")
    p.add_argument("--thk", type=float, default=5e-3, help="slice thickness [m]")
    p.add_argument("--shots", type=int, nargs="*", default=None)
    p.add_argument("--TR", type=float, default=0.0, help="0 = shortest for the slice count")
    p.add_argument("--flip", type=float, default=90.0)
    p.add_argument("--dwell", type=float, default=4e-6)
    p.add_argument("--dummy", type=int, default=1)
    p.add_argument("--gmax", type=float, default=40.0, help="mT/m (MR750w max 44)")
    p.add_argument("--smax", type=float, default=150.0, help="T/m/s (MR750w max 200)")
    p.add_argument("--grad-delay", type=float, default=0.0, help="measured gradient delay [s]")
    p.add_argument("--out", default="cyclic_epi_ge")
    main(p.parse_args())
