"""Load HBI slip/pressure profiles for overlay on the GC analytical solution.

Import from the notebook:

    from sim_curves import load_slip, load_dp, sim_params, front_lambda

Everything is returned in SI (r in m, slip in m, dp in Pa) so it lines up with
the GC script's internals; multiply by 100 and divide by 1000 at plot time as
that script does.

Conventions worth knowing before comparing:

  * Profiles are AZIMUTHALLY AVERAGED by default. The GC solution is a function
    of r alone, so an azimuthal average is the right thing to compare, and the
    spread across angle doubles as an axisymmetry check (return_std=True).
  * The crack edge in the theory is where slip -> 0. Any finite threshold puts
    the measured edge INSIDE the true one and biases lambda low: for 632892 at
    17 d, lambda is 0.2954 at a 1e-3 m threshold, 0.2991 at 1e-4 and 0.3046 at
    1e-8. front_lambda() therefore takes the threshold explicitly rather than
    hiding one.
  * Requested times are matched to the nearest output snapshot, and the actual
    time is returned. Do not assume you got the time you asked for -- a run
    ending at 16.999999 d does not contain t = 17 exactly.
"""
import glob
import os

import numpy as np

SCRATCH = "/scratch/users/nberrios/3dhbi/output"
RUNS = "/scratch/users/nberrios/3dhbi/runs"
DECKS = "/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs"


def _deck(job):
    d = {}
    for line in open(f"{DECKS}/res{job}.in"):
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            d[w[0]] = " ".join(w[1:]).strip('"')
    return d


def _ff(x):
    return float(str(x).replace("d", "e").replace("D", "e"))


def _path(job, field):
    """Finished runs live under SCRATCH; a run in flight is still under RUNS."""
    for c in [f"{SCRATCH}/{job}/{field}{job}.dat"] + sorted(
            glob.glob(f"{RUNS}/*/output/{field}{job}.dat")):
        if os.path.exists(c) and os.path.getsize(c):
            return c
    raise FileNotFoundError(f"no {field}{job}.dat in {SCRATCH} or {RUNS}")


def sim_params(job):
    """The GC solution's physical inputs, as this deck actually sets them.

    alpha uses KPMAX when permev is T: with enhancement the slipping region
    saturates at kpmax, and that is the permeability the point-source solution
    should be evaluated at. Using kp instead gives a diffusivity that
    corresponds to nothing in the model.
    """
    d = _deck(job)
    eta, phi, beta = _ff(d["eta"]), _ff(d["phi"]), _ff(d["beta"])
    on = d.get("permev", "F").upper() == "T"
    k_eff = _ff(d["kpmax"]) if on else _ff(d["kp"])
    qfile = d.get("injection_file")
    q_over_w = float(open(f"{DECKS}/{qfile}").read().splitlines()[4].split()[0])
    return dict(
        f=_ff(d["f0"]),
        sigma_0=_ff(d["sigmainit"]) * 1e6,
        tau_0=_ff(d["muinit"]) * _ff(d["sigmainit"]) * 1e6,
        mu=_ff(d["rigid"]) * 1e9,
        alpha=k_eff / (eta * phi * beta),
        deltaP=q_over_w * eta / (4 * np.pi * k_eff),
        k_eff=k_eff, permev=d.get("permev", "F"),
        ds_m=_ff(d["ds"]) * 1000.0, imax=int(d["imax"]),
    )


def _load(job, field, t_days, azimuthal, return_std):
    d = _deck(job)
    IM, JM = int(d["imax"]), int(d["jmax"])
    ds_m = _ff(d["ds"]) * 1000.0
    NC = IM * JM
    p = _path(job, field)
    t = np.atleast_2d(np.loadtxt(p.replace(field, "time", 1)))[:, 1] / 86400.0
    nt = min(os.path.getsize(p) // (8 * NC), len(t))
    arr = np.memmap(p, np.float64, "r", shape=(nt, NC))
    k = int(np.argmin(np.abs(t[:nt] - t_days)))
    frame = np.asarray(arr[k]).reshape(IM, JM)
    c = (IM - 1) // 2

    if not azimuthal:
        prof = frame[c, c:]
        return np.arange(len(prof)) * ds_m, prof, float(t[k]), None

    ii, jj = np.mgrid[0:IM, 0:JM]
    idx = np.hypot(ii - c, jj - c).astype(int)
    nb = idx.max() + 1
    cnt = np.bincount(idx.ravel(), minlength=nb).astype(float)
    tot = np.bincount(idx.ravel(), weights=frame.ravel(), minlength=nb)
    mean = tot / np.maximum(cnt, 1)
    std = None
    if return_std:
        sq = np.bincount(idx.ravel(), weights=frame.ravel() ** 2, minlength=nb)
        std = np.sqrt(np.maximum(sq / np.maximum(cnt, 1) - mean ** 2, 0.0))
    return np.arange(nb) * ds_m, mean, float(t[k]), std


def load_slip(job, t_days, azimuthal=True, return_std=False):
    """(r_m, slip_m, t_actual_days[, std_m]) at the snapshot nearest t_days."""
    r, v, ta, sd = _load(job, "slip", t_days, azimuthal, return_std)
    return (r, v, ta, sd) if return_std else (r, v, ta)


def load_dp(job, t_days, azimuthal=True, return_std=False):
    """(r_m, dp_Pa, t_actual_days[, std_Pa]). pf is written in MPa; converted."""
    r, v, ta, sd = _load(job, "pf", t_days, azimuthal, return_std)
    return ((r, v * 1e6, ta, None if sd is None else sd * 1e6) if return_std
            else (r, v * 1e6, ta))


def front_lambda(job, t_days, alpha, thresh=1e-6):
    """(R_m, lambda, t_actual) with lambda = R/sqrt(4*alpha*t).

    thresh is explicit because it biases the answer: the theory's edge is where
    slip vanishes, so a larger threshold reports a smaller crack. See the module
    docstring for the measured sensitivity.
    """
    r, sl, ta = load_slip(job, t_days)
    b = np.where(sl < thresh)[0]
    R = r[b[0]] if len(b) else r[-1]
    return float(R), float(R / np.sqrt(4 * alpha * ta * 86400)), ta
