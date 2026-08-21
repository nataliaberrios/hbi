#!/usr/bin/env python3
"""Census of every run that has BOTH permeability enhancement and a nonuniform
initial permeability -- the combination Natalia's hypothesis rests on.

Selection is read from the decks, not from a hand-kept list:
    permev T
  AND parameterfromfile T with a parameter_file that is genuinely nonuniform
      (max/min > 1.001 -- NOT np.isclose, whose default atol=1e-8 dwarfs
       permeabilities of order 1e-15 and reported every map as uniform)

Two passes, because a full lambda fit reads every frame of a ~300 MB slip file and
most of these runs never slip at all:

  PASS 1 (cheap, final frame only) -- peak slip and the radius where slip crosses
    dc. A run whose peak slip never reaches dc has no front by definition, and
    that is the single most common outcome here. Reported as peak/dc so the margin
    is visible: 1.6 is a razor-thin front whose lambda is essentially noise, 175
    is a robust one.

  PASS 2 (full history) -- lambda, fit on a shared window, only for runs that
    cleared dc in pass 1.

Also flags decks that set permev T but leave kpmax and/or kpmin unset. HBI's
t_params does not initialise them and the defaults block does not either, so
permeability decays without a floor -- k reached 1.1e-16 in 632537. Those runs
are not evidence about enhancement and are listed separately.

Writes:
    figures/enhancement_census/census.csv        every run, both passes
    figures/enhancement_census/census.txt        readable table
    figures/enhancement_census/overview.png      the whole population, both targets

Usage:  python census_enhancement.py            # full census
        python census_enhancement.py --quick    # pass 1 only
"""
import argparse
import glob
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
OUT = H / "figures" / "enhancement_census"
IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
RHO, G, HW, DW, FD, P0, W = 1000.0, 9.81, 4077.0, 0.178, 0.015, 73.8, 6.0

INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
OBSC, OK, BAD = "#a8071a", "#009E73", "#0072BD"

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11.5, "axes.labelsize": 10,
    "axes.edgecolor": MUTED, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "legend.fontsize": 8.5,
})


def ffloat(x):
    try:
        return float(str(x).replace("d", "e").replace("D", "e"))
    except (TypeError, ValueError):
        return None


def deck(p):
    d = {}
    for ln in Path(p).read_text().splitlines():
        if ln.startswith("!"):
            continue
        w = ln.split()
        if len(w) >= 2:
            d[w[0]] = w[1].strip('"')
    return d


def outdir(n):
    for pat in (f"/scratch/users/nberrios/3dhbi/output/{n}",
                f"/oak/stanford/groups/edunham/nberrios/3doutput/{n}"):
        if os.path.exists(f"{pat}/time{n}.dat"):
            return pat
    g = glob.glob(f"/scratch/users/nberrios/3dhbi/runs/*/output/time{n}.dat")
    return os.path.dirname(g[0]) if g else None


def select():
    """Every deck with enhancement AND a nonuniform initial permeability map."""
    keep, invalid = [], []
    for p in sorted(IN.glob("res*.in")):
        st = p.stem[3:]
        if not st.isdigit():
            continue
        dk = deck(p)
        if not dk.get("permev", "F").upper().startswith("T"):
            continue
        if dk.get("parameterfromfile", "F").upper() not in ("T", "TRUE", ".TRUE."):
            continue
        pf = dk.get("parameter_file")
        if not pf or not (IN / pf).exists():
            continue
        k = np.loadtxt(IN / pf, skiprows=1)
        lo, hi = float(k.min()), float(k.max())
        if not (lo > 0 and hi / lo > 1.001):
            continue
        rec = dict(n=int(st), dk=dk, near=hi, far=lo, pf=pf)
        # kpmax/kpmin unset -> permeability decays with no floor; not evidence
        if dk.get("kpmax") and dk.get("kpmin"):
            keep.append(rec)
        else:
            invalid.append(rec)
    return keep, invalid


def pass1(rec):
    """Final frame only: peak slip, and the radius where slip crosses dc."""
    n, dk = rec["n"], rec["dk"]
    base = outdir(n)
    if base is None:
        return None
    sp = f"{base}/slip{n}.dat"
    tp = f"{base}/time{n}.dat"
    if not (os.path.exists(sp) and os.path.getsize(sp)):
        return None
    IM, JM = int(dk["imax"]), int(dk["jmax"])
    ds, NC = ffloat(dk["ds"]), IM * JM
    try:
        t = np.atleast_2d(np.loadtxt(tp))[:, 1] / 86400.0
    except Exception:
        return None
    nt = min(os.path.getsize(sp) // (8 * NC), len(t))
    if nt < 2:
        return None
    # read ONLY the last frame: offset straight to it rather than mapping the file
    with open(sp, "rb") as f:
        f.seek((nt - 1) * NC * 8)
        last = np.frombuffer(f.read(NC * 8), np.float64).reshape(IM, JM)
    dc = ffloat(dk["dc"])
    pk = float(last.max())
    prof = last[:, JM // 2]
    x = (np.arange(IM) - IM // 2) * ds
    over = prof > dc
    r_dc = float(np.abs(x[over]).max()) if over.any() else 0.0
    return dict(t_end=float(t[nt - 1]), nframes=int(nt), peak=pk, dc=dc,
                peak_over_dc=pk / dc, r_dc_km=r_dc, base=base)


def observed_lambda(tcut):
    from datetime import datetime as dt, timedelta as td
    from scipy.io import loadmat
    mat = loadmat(H / "Cooper_Basin_Catalog_HAB_4.mat", squeeze_me=True,
                  struct_as_record=False)
    cat = {e.field: e.val for e in mat["Catalog"]}
    tt = cat["Time"].astype(float)
    dts = np.array([dt.fromordinal(int(x)) + td(days=x % 1) - td(days=366) for x in tt])
    dur = (dts >= dt(2012, 11, 13)) & (dts <= dt(2012, 11, 30))
    la, lo = -27.8115, 140.7596
    d = np.sqrt(((cat["Lat"][dur] - la) * 111.0) ** 2
                + ((cat["Long"][dur] - lo) * 111.0 * np.cos(np.radians(la))) ** 2)
    t = tt[dur] - tt[dur][0]
    org = np.median(d[:10])
    o = np.argsort(t); st, sd = t[o], (d - org)[o]
    ft, fd = [], []
    for i in range(0, len(st), 100):
        bt, bd = st[i:i + 100], sd[i:i + 100]
        if not len(bt):
            continue
        l_, h_ = np.percentile(bd, 90), np.percentile(bd, 95)
        k = (bd >= l_) & (bd <= h_)
        ft.extend(bt[k]); fd.extend(bd[k])
    ft, fd = np.array(ft), np.array(fd)
    m = ft <= tcut
    b = np.sqrt(ft[m])
    return float(np.sum(b * fd[m]) / np.sum(b ** 2))


def pass2(rec, p1, tcut):
    """Full slip history -> lambda on the shared window."""
    n, dk = rec["n"], rec["dk"]
    IM, JM = int(dk["imax"]), int(dk["jmax"])
    ds, NC = ffloat(dk["ds"]), IM * JM
    sp = f"{p1['base']}/slip{n}.dat"
    t = np.atleast_2d(np.loadtxt(f"{p1['base']}/time{n}.dat"))[:, 1] / 86400.0
    nt = min(os.path.getsize(sp) // (8 * NC), len(t))
    mm = np.memmap(sp, np.float64, "r", shape=(nt, NC))
    cs = np.array([np.asarray(mm[k]).reshape(IM, JM)[:, JM // 2] for k in range(nt)]).T
    x = (np.arange(IM) - IM // 2) * ds
    dc = p1["dc"]
    R, T = [], []
    for i in range(IM):
        a = cs[i, :] > dc
        if a.any():
            R.append(abs(x[i])); T.append(t[np.argmax(a)])
    R, T = np.array(R), np.array(T)
    m = T <= tcut
    if m.sum() < 2:
        return np.nan
    b = np.sqrt(T[m])
    return float(np.sum(b * R[m]) / np.sum(b ** 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="pass 1 only")
    ap.add_argument("--tcut", type=float, default=8.0)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    keep, invalid = select()
    print(f"{len(keep)} valid decks (enhancement + nonuniform perm, bounds set)")
    print(f"{len(invalid)} excluded: kpmax/kpmin unset\n")

    rows = []
    for rec in keep:
        p1 = pass1(rec)
        if p1 is None:
            continue
        dk = rec["dk"]
        rows.append(dict(
            n=rec["n"], mu=ffloat(dk.get("muinit")), sig=ffloat(dk.get("sigmainit")),
            pb=(ffloat(dk.get("phi")) or 0) * (ffloat(dk.get("beta")) or 0),
            kpmax=dk.get("kpmax"), kpmin=dk.get("kpmin"), kL=dk.get("kL"),
            near=rec["near"], far=rec["far"], pf=rec["pf"], **p1))
        print(f"  {rec['n']}  peak {p1['peak']:.2e}  peak/dc {p1['peak_over_dc']:8.2f}  "
              f"r_dc {p1['r_dc_km']:.3f} km  {p1['t_end']:.2f} d", flush=True)

    slipped = [r for r in rows if r["peak_over_dc"] >= 1.0]
    print(f"\n{len(rows)} runs with output; {len(slipped)} ever reached dc")

    lam_obs = observed_lambda(a.tcut)
    if not a.quick:
        for r in slipped:
            if r["t_end"] < 0.9 * a.tcut:
                r["lam"] = np.nan
                continue
            r["lam"] = pass2({"n": r["n"], "dk": deck(IN / f"res{r['n']}.in")}, r, a.tcut)
            print(f"  {r['n']}  lambda {r['lam']:.4f}  = {r['lam']/lam_obs:.2f}x obs",
                  flush=True)

    hdr = ("run,muinit,sigmainit,tau0,phibeta,kpmax,kpmin,kL,k_near,k_far,"
           "t_end_days,nframes,peak_slip_m,dc,peak_over_dc,r_dc_km,lambda,lam_over_obs,"
           "perm_file")
    lines = [hdr]
    for r in sorted(rows, key=lambda z: -z["peak_over_dc"]):
        lam = r.get("lam", np.nan)
        lines.append(",".join(str(v) for v in [
            r["n"], r["mu"], r["sig"], round((r["mu"] or 0) * (r["sig"] or 0), 3),
            f"{r['pb']:.3e}", r["kpmax"], r["kpmin"], r["kL"],
            f"{r['near']:.3e}", f"{r['far']:.3e}", f"{r['t_end']:.3f}", r["nframes"],
            f"{r['peak']:.4e}", r["dc"], f"{r['peak_over_dc']:.3f}",
            f"{r['r_dc_km']:.4f}",
            "" if not np.isfinite(lam) else f"{lam:.4f}",
            "" if not np.isfinite(lam) else f"{lam/lam_obs:.3f}", r["pf"]]))
    (OUT / "census.csv").write_text("\n".join(lines) + "\n")

    L = [f"Runs with BOTH permeability enhancement and nonuniform initial perm",
         "=" * 104,
         f"observed lambda on 0-{a.tcut:.1f} d = {lam_obs:.4f}",
         f"{len(rows)} runs with output; {len(slipped)} ever reached dc "
         f"(peak slip >= dc, so they have a front at all)", "",
         f"{'run':>8s} {'mu':>5s} {'tau0':>6s} {'phibeta':>9s} {'kpmax':>8s} "
         f"{'kL':>6s} {'peak slip':>10s} {'peak/dc':>8s} {'lam/obs':>8s} {'days':>6s}",
         "-" * 104]
    for r in sorted(rows, key=lambda z: -z["peak_over_dc"]):
        lam = r.get("lam", np.nan)
        ls = f"{lam/lam_obs:.2f}" if np.isfinite(lam) else ("-" if r["peak_over_dc"] >= 1
                                                            else "no front")
        L.append(f"{r['n']:>8d} {r['mu']:>5.2f} {(r['mu'] or 0)*(r['sig'] or 0):>6.2f} "
                 f"{r['pb']:>9.1e} {str(r['kpmax']):>8s} {str(r['kL']):>6s} "
                 f"{r['peak']:>10.2e} {r['peak_over_dc']:>8.2f} {ls:>8s} "
                 f"{r['t_end']:>6.2f}")
    L += ["", "EXCLUDED -- permev T but kpmax and/or kpmin unset, so permeability "
          "decays with no floor (k reached 1.1e-16 in 632537):"]
    L += ["  " + ", ".join(str(r["n"]) for r in invalid)]
    (OUT / "census.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L[:8]))
    print(f"\nwrote {OUT}/census.csv and census.txt")

    # ---- overview: the whole population against both targets
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8), constrained_layout=True)
    ax = axes[0]
    for r in rows:
        c = OK if r["peak_over_dc"] >= 10 else (BAD if r["peak_over_dc"] >= 1 else MUTED)
        ax.scatter([r["pb"]], [r["peak_over_dc"]], s=44, color=c, alpha=0.85,
                   edgecolors=INK, linewidths=0.4)
    ax.axhline(1.0, color=OBSC, lw=1.6, ls="--")
    ax.annotate("peak slip = dc\nbelow this line there is no front at all",
                xy=(0.02, 1.0), xycoords=("axes fraction", "data"),
                xytext=(0, 6), textcoords="offset points", fontsize=9, color=OBSC)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set(xlabel=r"storage $\phi\beta$ (Pa$^{-1}$)",
           ylabel="peak slip / dc",
           title="Does it slip at all?")
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax = axes[1]
    got = [r for r in rows if np.isfinite(r.get("lam", np.nan))]
    if got:
        ax.axhspan(0.85, 1.15, color=OK, alpha=0.13, zorder=0)
        ax.axhline(1.0, color=OBSC, lw=1.6, ls="--")
        for r in got:
            ax.scatter([(r["mu"] or 0) * (r["sig"] or 0)], [r["lam"] / lam_obs],
                       s=52, color=OK if abs(r["lam"] / lam_obs - 1) <= 0.15 else BAD,
                       edgecolors=INK, linewidths=0.4, alpha=0.9)
            ax.annotate(str(r["n"]), xy=((r["mu"] or 0) * (r["sig"] or 0),
                                         r["lam"] / lam_obs),
                        xytext=(5, 3), textcoords="offset points", fontsize=7.5,
                        color=INK)
        ax.axvline(15.0, color=MUTED, lw=1.1, ls=":")
        ax.annotate("Taiyi $\\tau_0$=15.0 MPa\n(understressed is left of this)",
                    xy=(15.0, 0.02), xycoords=("data", "axes fraction"),
                    xytext=(-6, 0), textcoords="offset points", ha="right",
                    fontsize=9, color=MUTED)
        ax.set_yscale("log")
    else:
        ax.text(0.5, 0.5, "no run in this population produced a front",
                transform=ax.transAxes, ha="center", fontsize=13, color=OBSC)
        ax.set(xticks=[], yticks=[])
    ax.set(xlabel=r"initial shear stress $\tau_0$ (MPa)",
           ylabel=r"$\lambda\,/\,\lambda_{observed}$",
           title=f"Of those that slip, how good is the front?  "
                 f"(fits on 0–{a.tcut:.0f} d)")
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.suptitle("Every run with permeability enhancement AND nonuniform initial "
                 "permeability\n"
                 f"{len(rows)} runs with output; {len(slipped)} ever reach dc; "
                 f"{len(got)} have a fittable front", fontsize=12)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"overview.{e}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/overview.png")


if __name__ == "__main__":
    main()
