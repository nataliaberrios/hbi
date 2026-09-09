#!/usr/bin/env python3
"""Shift june_clean.txt's clock to a new t = 0 and write a new injection file.

june_clean.txt is UNTOUCHED, always. md5 603285de6a5b1145d61bea05d97ebbf7.
This reads it and writes a new file; it never opens it for writing.

WHY A SECOND SHIFT. june_clean_from_d4300.txt starts at data-day 4.300, where
injection resumes in earnest, and therefore DISCARDS THE TRICKLE: the record
carries 1.6-2.7 L/s from 3.5575 d to 4.2999 d, 17.8 hours of it, before the
rate steps to 7.9 L/s.

That trickle is only 0.1417 ML, 0.42% of the 33.65 ML injected after 4.300, so
by volume it looks ignorable. It is not, because it arrives FIRST and has time
to spread: in 17.8 h at the near-well diffusivity kpmax/(eta phi beta) =
8.749 m^2/s, sqrt(4 D t) = 1498 m. The pressure field the main ramp starts from
is therefore already charged out to roughly the radius the front later grows
through, and the model as run starts from a uniform pfinit instead. That is a
candidate explanation for the model's early wellhead rise being too steep and
too early against the observed curved build-up.

CONSEQUENCE, and it cannot be avoided: including the trickle MOVES t = 0. There
is no way to prepend 17.8 hours of injection to a file whose clock starts at
4.300, because HBI's clock starts at zero. So a run using this file is on a
DIFFERENT CLOCK from 632950-632994, shifted by -0.7425 d, and everything read
against it -- observed front, wellhead, slip -- must be shifted by the new t0.
score_cycle2.py takes T0 as a module constant for that reason.

WHAT IS WRITTEN. Times in seconds, rates in m^2/s, in the format
m_diffusion.f90:58-66 reads: nwell / npoint / times / i j / rates. The
construction is the same as june_clean_from_d4300.txt's -- keep every point at
or after t0, shift by -t0, and prepend an INTERPOLATED point at t = 0 so the
file starts exactly at the new origin rather than at the first surviving
sample.

VERIFIED BEFORE WRITING: the volume from t0 to the end of the new file must
equal june_clean.txt's own volume over the same interval to better than 1e-6
relative. That catches a dropped point, a units slip, and an off-by-one in the
interpolated first point, which are the three ways this goes wrong.

Usage:
  python build_injection_shift.py --t0 3.5575            # dry run
  python build_injection_shift.py --t0 3.5575 --write
"""
import argparse
import hashlib
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
SRC = "june_clean.txt"
SRC_MD5 = "603285de6a5b1145d61bea05d97ebbf7"


def read_inj(name):
    """(nwell, i, j, t_days, q_raw) from an HBI injection file."""
    f = (IN / name).read_text().split("\n")
    nwell, npoint = int(f[0]), int(f[1])
    t = np.array(f[2].split(), float)
    ij = f[3].split()
    q = np.array(f[4].split(), float)
    assert nwell == 1, f"{name}: nwell = {nwell}, this tool assumes 1"
    assert len(t) == npoint == len(q), \
        f"{name}: npoint {npoint} vs {len(t)} times, {len(q)} rates"
    return nwell, int(ij[0]), int(ij[1]), t / 86400.0, q


def volume(t_days, q):
    """Trapezoidal integral in the file's own units x seconds."""
    return float(np.trapz(q, t_days * 86400.0))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--t0", type=float, required=True,
                    help="data-day to become t = 0")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args(argv)

    md5 = hashlib.md5((IN / SRC).read_bytes()).hexdigest()
    assert md5 == SRC_MD5, \
        f"{SRC} has md5 {md5}, expected {SRC_MD5} -- it must not be modified"
    nwell, iw, jw, t, q = read_inj(SRC)
    print(f"{SRC}: md5 verified, {len(t)} points, 0 to {t.max():.4f} d, "
          f"q up to {q.max():.6g} (file units)")

    keep = t >= a.t0
    q0 = float(np.interp(a.t0, t, q))
    tn = np.concatenate([[0.0], t[keep] - a.t0])
    qn = np.concatenate([[q0], q[keep]])
    # if a sample sits exactly on t0 the prepended point duplicates it
    if len(tn) > 1 and tn[1] < 1e-9:
        tn, qn = tn[1:], qn[1:]

    vsrc = volume(t[t >= a.t0], q[t >= a.t0])
    # the source's volume from t0 must include the partial interval from t0 to
    # the first surviving sample -- omitting it was a real error once, and it
    # produced a spurious 0.06% failure
    i0 = int(np.argmax(t >= a.t0))
    if i0 > 0:
        vsrc += 0.5 * (q0 + q[i0]) * (t[i0] - a.t0) * 86400.0
    vnew = volume(tn, qn)
    rel = abs(vnew - vsrc) / vsrc
    print(f"t0 = {a.t0} d -> {len(tn)} points, 0 to {tn.max():.4f} d")
    print(f"  interpolated first point q(t0) = {q0:.6e}")
    print(f"  volume  source {vsrc:.8e}   new {vnew:.8e}   rel {rel:.3e}")
    assert rel < 1e-6, f"volume not conserved: {100*rel:.6f}%"
    print(f"  volume conserved to {100*rel:.6f}%")

    # what the trickle adds, in the units a reader can check
    W = 6.0
    tri = (t >= a.t0) & (t < 4.300)
    if tri.any():
        vt = volume(np.concatenate([[a.t0], t[tri], [4.300]]),
                    np.concatenate([[q0], q[tri],
                                    [float(np.interp(4.300, t, q))]]))
        print(f"  the {(4.300-a.t0)*24:.1f} h before data-day 4.300 carries "
              f"{vt*W*1000/1e6:.4f} ML at {q[tri].min()*W*1000:.2f}-"
              f"{q[tri].max()*W*1000:.2f} L/s")

    name = f"june_clean_from_d{int(round(a.t0*1000)):04d}.txt"
    if not a.write:
        print(f"\ndry run -- would write {name}. re-run with --write")
        return
    assert not (IN / name).exists(), \
        f"{name} already exists; refusing to overwrite an input file"
    body = (f"{nwell}\n{len(tn)}\n"
            + " ".join(f"{x*86400.0:.6e}" for x in tn) + "\n"
            + f"{iw} {jw}\n"
            + " ".join(f"{x:.6e}" for x in qn) + "\n")
    (IN / name).write_text(body)
    print(f"\nwrote {IN / name}")
    # read it back through the same reader, so the file that HBI will parse is
    # the thing that gets checked, not the arrays in memory
    _, i2, j2, t2, q2 = read_inj(name)
    assert (i2, j2) == (iw, jw), "well index changed"
    assert abs(volume(t2, q2) - vnew) / vnew < 1e-9, "round-trip volume changed"
    print(f"  read back: {len(t2)} points, well {i2} {j2}, "
          f"0 to {t2.max():.4f} d, volume round-trips")


if __name__ == "__main__":
    main()
