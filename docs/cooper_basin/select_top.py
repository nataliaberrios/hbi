#!/usr/bin/env python3
"""Pick the N best cycle-2 runs, DE-DUPLICATED, and write the figure set.

rank_grid.py ranks everything and is the record. This is the selector the
overnight loop uses to decide which runs get figures, and it differs from
rank_grid in two ways that matter:

DE-DUPLICATION. Taken straight, the global top 10 spends three slots on
633069/633070/633071 -- the Sw_fwid x10/x100/x1000 sweep -- which land at
lambda 1.03/1.03/1.01 and dp +5.5/+5.4/+3.7. Sw_fwid was measured to move dp by
2.5 points and lambda not at all, so those are three pictures of one thing while
genuinely different parameter combinations get pushed off the list. Two runs are
treated as the same match if they agree to TOL_LAM in lambda ratio AND TOL_DP in
dp percent; the better-scoring one survives. This is a rule about the TARGETS,
not about the parameters, so it does not need to know which knob was turned.

ARM 2 ONLY. eta = 1.27e-4 Pa s is this project's fluid; arm 1 is Wang & Dunham's
0.89e-3 and belongs only in runs reproducing them. Several arm-1 runs score well
enough to enter a top 10 (632900 at 0.182) and would be misleading there.

THE FIGURES ARE compare_pair.py's, not make_run_figures.py's. make_run_figures
plots ABSOLUTE wellhead pressure against "days since injection began" on the
unshifted observed clock, which is the pre-cycle-2 convention: for a run whose
t = 0 is data-day 4.300 it compares the wrong instants and ignores the Holl
datum. compare_pair gives the four panels asked for -- front vs time (RT), front
vs injected volume (RV), pressure change on the datum, and cumulative slip at
the injector -- with the datum and the front definition imported rather than
reimplemented.

SCORE is rank_grid's, unchanged: |lambda/lambda_obs - 1| + |dp|/100. Slip is
reported and not scored, because it sits at 0.35-0.65x in every run and would
rank on a constant.

Usage:
  python select_top.py                 # list the selection, no figures
  python select_top.py --figures       # and build them
  python select_top.py -n 6 --figures
"""
import argparse
import glob
import importlib.util as iu
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
H = "/home/users/nberrios/3dhbi/hbi_analysis"
_s = iu.spec_from_file_location("sf", H + "/make_sweep_figures.py")
sf = iu.module_from_spec(_s); _s.loader.exec_module(sf)
_a = iu.spec_from_file_location("ca", HERE / "compare_arms.py")
ca = iu.module_from_spec(_a); _a.loader.exec_module(ca)
_f = iu.spec_from_file_location("fp", HERE / "fault_pressure.py")
fp = iu.module_from_spec(_f); _f.loader.exec_module(fp)

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
OUTD = Path("/scratch/users/nberrios/3dhbi/output")
T0, T_EVAL = 4.300, 8.7
TOL_LAM, TOL_DP = 0.02, 2.0
ETA2 = 1.27e-4


def completed():
    """Every run with a deck, an output directory and the d4300 clock."""
    out = []
    for f in sorted(glob.glob(str(IN / "res632[89][0-9][0-9].in"))
                    + glob.glob(str(IN / "res633[0-1][0-9][0-9].in"))):
        n = int(re.search(r"res(\d+)", f).group(1))
        if not (OUTD / str(n)).is_dir():
            continue
        try:
            dk = sf.deck(n)
            if abs(sf.ffloat(dk["eta"]) - ETA2) > 1e-12:
                continue
            if "d4300" not in dk.get("injection_file", ""):
                continue
        except KeyError:
            continue
        out.append(n)
    return out


def score_one(n, obs, ref):
    dk = sf.deck(n)
    d = sf.run_data(n, dk)
    if d is None or d.get("tpw") is None or len(d["T"]) < 3:
        return None
    LOBS, to, dpo = ref
    o = np.argsort(d["T"])
    T, R = np.asarray(d["T"])[o], np.asarray(d["R"])[o] * 1000.0
    L = ca.lam(T, R)
    gr = np.linspace(0.05, min(d["tpw"][-1], to.max(), T_EVAL), 2000)
    ps = np.interp(gr, d["tpw"], d["ppw"]) - (sf.P0 - sf.RHO * sf.G * sf.HW / 1e6)
    ob = np.interp(gr, to, dpo)
    qg = np.interp(gr, obs["ti"] - T0, obs["q"])
    fl = (qg > 0.25 * np.nanmax(obs["q"])) \
        & (np.interp(gr, obs["tp"] - T0, obs["pm"]) > 5.0)
    if fl.sum() < 30:
        return None
    e = 100.0 * float(np.mean(ps[fl] - ob[fl])) / float(np.mean(ob[fl]))
    m = re.search(r"disc(\d+)", dk.get("parameter_file", ""))
    # PLAINNESS. When two runs land on the same match, the one to SHOW is the
    # canonical grid run, not a variant of it that happens to score 0.003
    # better. 633070 (Sw_fwid x100) beat 632961 on exactly that, even though
    # Sw_fwid was measured to move dp 2.5 points and lambda not at all -- so the
    # figure would have been captioned with a parameter doing no work. Counts
    # how far the deck sits from the arm-2 baseline; lower is plainer.
    exotic = 0
    exotic += abs(sf.ffloat(dk.get("Sw_fwid", "7.4e-9")) - 7.4e-9) > 1e-13
    exotic += abs(sf.ffloat(dk.get("phi", "0.01")) - 0.01) > 1e-12
    exotic += "skin" in dk
    exotic += abs(sf.ffloat(dk.get("rw", "8.9e-2")) - 8.9e-2) > 1e-12
    return dict(n=n, lam=L / LOBS, dp=e, t_end=d["t_end"],
                tau=sf.ffloat(dk["muinit"]) * sf.ffloat(dk["sigmainit"]),
                disc=int(m.group(1)) if m else None,
                kpmax=sf.ffloat(dk["kpmax"]),
                permev=dk.get("permev", "T"), exotic=int(exotic),
                score=abs(L / LOBS - 1.0) + abs(e) / 100.0)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--top", type=int, default=10)
    ap.add_argument("--figures", action="store_true")
    a = ap.parse_args(argv)

    obs = sf.observed()
    tc, rc, _ = ca.fb.catalogue()
    k = (tc - T0) > 0
    te, re_ = tc[k] - T0, rc[k]
    org = float(np.median(re_[np.argsort(te)][:10]))
    ft, fd = ca.seismicity_front(te, re_ - org)
    ref = (ca.lam(ft, fd), *fp.dp_observed(obs))
    print(f"observed lambda {ref[0]:.1f} m/sqrt(d); datum {fp.datum():.3f} MPa")

    runs = completed()
    rows = [r for r in (score_one(n, obs, ref) for n in runs) if r]
    rows.sort(key=lambda r: r["score"])
    print(f"{len(runs)} arm-2 runs on the d4300 clock, {len(rows)} scorable")

    sel, dropped = [], []
    for r in rows:
        dup = next((s for s in sel
                    if abs(s["lam"] - r["lam"]) < TOL_LAM
                    and abs(s["dp"] - r["dp"]) < TOL_DP), None)
        if dup:
            # same match -- keep the PLAINER deck even though it scored worse,
            # so the figure is captioned with parameters that did the work
            if r["exotic"] < dup["exotic"]:
                sel[sel.index(dup)] = r
                dropped.append((dup, r["n"]))
            else:
                dropped.append((r, dup["n"]))
            continue
        sel.append(r)
        if len(sel) == a.top:
            break
    sel.sort(key=lambda r: r["score"])

    print(f"\nTOP {len(sel)}, de-duplicated at "
          f"|dlambda| < {TOL_LAM} and |ddp| < {TOL_DP} points")
    print(f"{'#':>3} {'run':>7} {'tau_0':>6} {'disc':>5} {'kpmax':>9} "
          f"{'permev':>6} {'lam/obs':>8} {'dp %':>7} {'score':>6} {'t_end':>6}")
    for i, r in enumerate(sel, 1):
        print(f"{i:>3} {r['n']:>7} {r['tau']:>6.2f} {str(r['disc']):>5} "
              f"{r['kpmax']:>9.2e} {r['permev']:>6} {r['lam']:>8.2f} "
              f"{r['dp']:>+7.1f} {r['score']:>6.3f} {r['t_end']:>6.2f}")
    if dropped:
        print(f"\ncollapsed as duplicates of an already-selected match:")
        for r, of in dropped[:12]:
            print(f"    {r['n']}  lam {r['lam']:.2f}  dp {r['dp']:+.1f}  "
                  f"-> same match as {of}")

    if not a.figures:
        print("\nno figures written. re-run with --figures")
        return sel

    print()
    for i, r in enumerate(sel, 1):
        stem = f"top10/top{i:02d}_{r['n']}"
        (Path(H) / "figures" / "cycle2" / "top10").mkdir(parents=True,
                                                         exist_ok=True)
        cmd = [sys.executable, str(HERE / "compare_pair.py"), str(r["n"]),
               "--out", stem]
        p = subprocess.run(cmd, capture_output=True, text=True)
        tag = "ok " if p.returncode == 0 else "FAIL"
        print(f"  {tag} top{i:02d}_{r['n']}"
              + ("" if p.returncode == 0
                 else "  " + p.stderr.strip().splitlines()[-1][:90]))
    return sel


if __name__ == "__main__":
    main()
