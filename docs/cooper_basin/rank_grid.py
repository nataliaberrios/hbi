#!/usr/bin/env python3
"""Rank every completed cycle-2 run on ONE consistent pair of metrics.

WHY THIS EXISTS, and it is a correction. Two different front metrics were in
use at once and I quoted both as though they were the same number:

    lambda ratio      fitted to the PROJECT'S OWN percentile front (bins of
                      100 events, 90-95th percentile band, initial cloud radius
                      subtracted) over the whole curve, via
                      lambda = sum(sqrt(t) R)/sum(t). This is the definition in
                      cooper_basin_plots-30.cleaned.ipynb cells 49/51/55.
    front ratio       R_sim/R_obs at 8.7 d, with R_obs the RUNNING MAXIMUM.
                      This is what score_cycle2.py still computes.

They disagree, and not slightly: 632961 is lambda 1.03x but front ratio 0.78.
Worse, they disagree about the STRUCTURE of the grid. On the running-max ratio
the front is flat in disc radius, which led me to report that "the front
depends on tau_0 alone and the disc moves only dp" -- an orthogonal pair of
levers. On lambda that is false: at tau_0 11.53 lambda runs 229.7 to 188.0
(1.09x to 0.89x) across disc 150 to 450 m. The orthogonality was an artifact of
a metric that integrates out the disc's effect.

So everything here is on lambda, because that is the project's own definition
and the one every figure now uses. dp is on the Holl datum via
fault_pressure.dp_observed, which carries the corrected column length, the
measured flowing density and the segmented friction.

SCORE = |lambda/lambda_obs - 1| + |dp error| / 100. Equal weight on a
fractional front error and a fractional pressure error. Slip is reported but
NOT scored -- it sits at 0.35-0.65x throughout and no run has ever matched it,
so including it would rank on a constant.

Only arm-2 runs with the standard fluid (eta 1.27e-4, beta 2.25e-8, phi 0.01)
and no skin are compared, since those are the ones that differ in tau_0 and
disc alone. The phi and Sw_fwid sweeps are listed separately.

Usage:  python rank_grid.py
"""
import glob
import importlib.util as iu
import re
from pathlib import Path

import numpy as np

H = "/home/users/nberrios/3dhbi/hbi_analysis"
_s = iu.spec_from_file_location("sf", H + "/make_sweep_figures.py")
sf = iu.module_from_spec(_s); _s.loader.exec_module(sf)
_a = iu.spec_from_file_location(
    "ca", "/home/users/nberrios/3dhbi/hbi_git/docs/cooper_basin/compare_arms.py")
ca = iu.module_from_spec(_a); _a.loader.exec_module(ca)
_f = iu.spec_from_file_location(
    "fp", "/home/users/nberrios/3dhbi/hbi_git/docs/cooper_basin/fault_pressure.py")
fp = iu.module_from_spec(_f); _f.loader.exec_module(fp)

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
OUTD = Path("/scratch/users/nberrios/3dhbi/output")
MD = Path("/home/users/nberrios/3dhbi/hbi_git/docs/figs/cycle2/RANKING.md")
T0, TEVAL = 4.300, 8.7


def observed():
    obs = sf.observed()
    tc, rc, _ = ca.fb.catalogue()
    te = tc - T0
    k = te > 0
    te, re_ = te[k], rc[k]
    org = float(np.median(re_[np.argsort(te)][:10]))
    ft, fd = ca.seismicity_front(te, re_ - org)
    lo = ca.lam(ft, fd)
    to, dpo = fp.dp_observed(obs)
    return obs, lo, to, dpo, len(ft)


def score_one(n, obs, lam_obs, to, dpo):
    dk = sf.deck(n)
    d = sf.run_data(n, dk)
    if d is None or d.get("tpw") is None or len(d["T"]) < 3:
        return None
    o = np.argsort(d["T"])
    T, R = np.asarray(d["T"])[o], np.asarray(d["R"])[o] * 1000.0
    L = ca.lam(T, R)
    gr = np.linspace(0.05, min(d["tpw"][-1], to.max(), TEVAL), 2000)
    ps = np.interp(gr, d["tpw"], d["ppw"]) - (sf.P0 - sf.RHO*sf.G*sf.HW/1e6)
    ob = np.interp(gr, to, dpo)
    qg = np.interp(gr, obs["ti"] - T0, obs["q"])
    fl = (qg > 0.25*np.nanmax(obs["q"])) \
        & (np.interp(gr, obs["tp"] - T0, obs["pm"]) > 5.0)
    e = 100.0 * float(np.mean(ps[fl] - ob[fl])) / float(np.mean(ob[fl]))
    try:
        import sys
        sys.path.insert(0, H + "/notebooks")
        from sim_curves import load_slip
        _, sl, ta = load_slip(n, TEVAL, how="strike")
        slip = float(sl[0]) * 100.0 if abs(ta - TEVAL) < 0.3 else np.nan
    except Exception:
        slip = np.nan
    m = re.search(r"disc(\d+)", dk["parameter_file"])
    return dict(n=n, tau=sf.ffloat(dk["muinit"])*sf.ffloat(dk["sigmainit"]),
                disc=int(m.group(1)) if m else None, lam=L, ratio=L/lam_obs,
                dp=e, slip=slip,
                score=abs(L/lam_obs - 1) + abs(e)/100.0,
                phi=sf.ffloat(dk["phi"]), swf=sf.ffloat(dk["Sw_fwid"]),
                skin=dk.get("skin"), eta=sf.ffloat(dk["eta"]),
                beta=sf.ffloat(dk["beta"]),
                inj=dk["injection_file"].strip('"'))


def main(argv=None):
    obs, lam_obs, to, dpo, nft = observed()
    rows = []
    for f in sorted(glob.glob(str(IN / "res6329[0-9][0-9].in"))
                    + glob.glob(str(IN / "res6330[0-9][0-9].in"))):
        n = int(re.search(r"res(\d+)", f).group(1))
        if not (OUTD / str(n)).is_dir():
            continue
        r = score_one(n, obs, lam_obs, to, dpo)
        if r:
            rows.append(r)

    def std(r):
        return (abs(r["eta"] - 1.27e-4) < 1e-12 and abs(r["beta"] - 2.25e-8) < 1e-12
                and abs(r["phi"] - 0.01) < 1e-12 and r["skin"] is None
                and abs(r["swf"] - 7.4e-9) < 1e-13 and "d4300" in r["inj"])

    grid = sorted([r for r in rows if std(r)], key=lambda z: z["score"])
    other = [r for r in rows if not std(r)]

    L = []
    A = L.append
    A("# Cycle-2 ranking, on the project's own front definition\n")
    A("Generated by `docs/cooper_basin/rank_grid.py`. "
      f"{len(rows)} completed runs.\n")
    A(f"\nObserved **λ = {lam_obs:.1f} m/√d** from {nft} percentile-front "
      f"points. Δp on the Holl datum, **{fp.datum():.3f} MPa** at the wellhead.\n")
    A("\n`score = |λ/λ_obs − 1| + |Δp error|/100`. Slip is reported but not "
      "scored — it sits at 0.35–0.65× in every run, so scoring it would rank "
      "on a constant.\n")

    A("\n## The τ₀ × disc grid, ranked\n")
    A("| rank | run | τ₀ MPa | disc m | λ | λ/λ_obs | Δp % | slip | score |")
    A("|---|---|---|---|---|---|---|---|---|")
    for i, r in enumerate(grid, 1):
        A(f"| {i} | {r['n']} | {r['tau']:.2f} | {r['disc']} | {r['lam']:.1f} | "
          f"**{r['ratio']:.2f}** | **{r['dp']:+.1f}** | {r['slip']/2.598:.2f} | "
          f"{r['score']:.3f} |")

    A("\n## λ/λ_obs against (τ₀, disc) — the structure\n")
    taus = sorted({r["tau"] for r in grid})
    discs = sorted({r["disc"] for r in grid})
    A("| τ₀ | " + " | ".join(f"{d} m" for d in discs) + " |")
    A("|" + "---|" * (len(discs) + 1))
    for t in taus:
        cells = []
        for d in discs:
            m = [r for r in grid if r["tau"] == t and r["disc"] == d]
            cells.append(f"{m[0]['ratio']:.2f} / {m[0]['dp']:+.0f}%" if m else "—")
        A(f"| **{t:.2f}** | " + " | ".join(cells) + " |")
    A("\nEach cell is λ/λ_obs / Δp error. **λ varies with disc**, which the "
      "running-maximum front ratio hid: it is flat in disc, which is why an "
      "earlier note claimed the two levers were orthogonal. On λ they are not.\n")

    A("\n## The other sweeps, same metrics\n")
    A("| run | what | λ/λ_obs | Δp % | slip | score |")
    A("|---|---|---|---|---|---|")
    for r in sorted(other, key=lambda z: z["score"]):
        if r["skin"] is not None:
            w = f"skin {r['skin']}"
        elif abs(r["phi"] - 0.01) > 1e-12:
            w = f"phi {r['phi']:g} (φβ ×{r['phi']/0.01:.0f})"
        elif abs(r["swf"] - 7.4e-9) > 1e-13:
            w = f"Sw_fwid ×{r['swf']/7.4e-9:.0f}"
        elif "d3558" in r["inj"]:
            w = "trickle included"
        elif abs(r["eta"] - 0.89e-3) < 1e-9:
            w = "arm 1, eta 0.89e-3"
        elif r["beta"] > 1e-7:
            w = "arm 3, beta 1.58e-7"
        elif r["beta"] < 1e-8:
            w = "arm 4, beta 7.72e-9"
        else:
            w = f"disc {r['disc']} m"
        A(f"| {r['n']} | τ₀ {r['tau']:.2f}, {w} | {r['ratio']:.2f} | "
          f"{r['dp']:+.1f} | {r['slip']/2.598:.2f} | {r['score']:.3f} |")

    MD.parent.mkdir(parents=True, exist_ok=True)
    MD.write_text("\n".join(L) + "\n")
    print("\n".join(L[:8]))
    print(f"\ntop 8 of the grid:")
    for r in grid[:8]:
        print(f"  {r['n']}  tau_0 {r['tau']:5.2f}  disc {r['disc']:3d}  "
              f"lam {r['ratio']:.2f}x  dp {r['dp']:+6.1f}%  "
              f"slip {r['slip']/2.598:.2f}  score {r['score']:.3f}")
    print(f"\nwrote {MD}")


if __name__ == "__main__":
    main()
