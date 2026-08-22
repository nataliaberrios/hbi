#!/usr/bin/env python3
"""Regenerate the run-inventory table in the calibration README.

The table is built from the decks and the output, never typed by hand, so it
cannot drift from what was actually run. Rewrites everything between the
INVENTORY-START and INVENTORY-END markers in the README and leaves the prose
alone.

Columns, and why each is there rather than "read the deck":
  tau0        = muinit * sigmainit. The understressed/not judgement, against
                Wang & Dunham's 15.0 MPa.
  phi*beta    the storage product, which is the lever that moves the front.
  kpmax/kpmin only shown for permev T; meaningless when permeability is fixed.
  k near/far  read from the PERMEABILITY MAP FILE, not the deck. A deck reading
                "kp 4e-13" can be running a map with 1.1e-12 near the well.
  lam/obs     front coefficient over the observed one, both fit on the same
                window (lambda is window-dependent, so a single shared window is
                the only fair comparison).
  press %     wellhead error while FLOWING. Blank when the run is driven by an
                injection file other than june_clean.txt, since the measured
                record is the June 2012 stage and the comparison would be
                between two different stimulations.

Usage:  python make_readme_table.py
"""
import glob
import os
from pathlib import Path

import numpy as np
import importlib.util as iu

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
REPO = Path("/home/users/nberrios/3dhbi/hbi_git/docs/figs/cooper_grid")
TAU0_TAIYI = 15.0
TCUT = 8.0            # shared front-fit window, set by the 8 d runs

spec = iu.spec_from_file_location("sf", str(H / "make_sweep_figures.py"))
sf = iu.module_from_spec(spec)
spec.loader.exec_module(sf)

RHO, G, HW, DW, FD, P0, W = 1000.0, 9.81, 4077.0, 0.178, 0.015, 73.8, 6.0


def perm_map(dk):
    """The permeability actually in force, from the map file."""
    if dk.get("parameterfromfile", "F").upper() not in ("T", "TRUE", ".TRUE."):
        k = sf.ffloat(dk.get("kp"))
        return (k, k)
    pf = dk.get("parameter_file")
    if not pf or not (IN / pf).exists():
        return (None, None)
    k = np.loadtxt(IN / pf, skiprows=1)
    return (float(k.max()), float(k.min()))


def scores(n, dk, obs, lam_obs):
    d = sf.run_data(n, dk)
    if d is None:
        return None, None, None, None
    m = d["T"] <= TCUT
    lam = sf.fit(d["T"][m], d["R"][m]) if m.sum() >= 2 else np.nan
    pct = rms = None
    if d["tpw"] is not None and dk.get("injection_file") == "june_clean.txt":
        hi = min(d["tpw"][-1], obs["tp"].max(), TCUT)
        if hi > 0.2:
            gr = np.linspace(0.05, hi, 2000)
            ps = np.interp(gr, d["tpw"], d["ppw"])
            ob = np.interp(gr, obs["tp"], obs["pm"])
            qg = np.interp(gr, obs["ti"], obs["q"])
            fl = (qg > 0.25 * np.nanmax(obs["q"])) & (ob > 5.0)
            if fl.sum() >= 30:
                pct = 100.0 * (ps[fl].mean() / ob[fl].mean() - 1.0)
                rms = float(np.sqrt(np.mean((ps[fl] - ob[fl]) ** 2)))
    return lam, pct, rms, d["t_end"]


def _passes(rows, lam_obs):
    """Runs inside BOTH bands, ranked by combined distance from the target."""
    ok = []
    for r in rows:
        if r["lam"] is None or not np.isfinite(r["lam"]) or r["pct"] is None:
            continue
        if r["tend"] is not None and r["tend"] < 0.95 * TCUT:
            continue
        lr, pr = r["lam"] / lam_obs, r["pct"]
        if abs(lr - 1) <= 0.15 and abs(pr) <= 15:
            ok.append((abs(lr - 1) + abs(pr) / 100.0, r["n"], lr, pr, r))
    if not ok:
        return "**Passing both targets: none.**"
    ok.sort()
    L = [f"**Passing both targets: {len(ok)} run"
         + ("s" if len(ok) > 1 else "") + "**, best first:", ""]
    for _, n, lr, pr, r in ok:
        pev = "enhancement ON" if r["pev"] else "enhancement off"
        L.append(f"- **{n}** — front {lr:.2f}x, pressure {pr:+.1f}%, "
                 f"tau_0 {r['tau0']:.2f} MPa, phi*beta {r['pb']:.1e}, {pev}")
    return "\n".join(L)


def main():
    obs = sf.observed()
    mo = obs["ft"] <= TCUT
    lam_obs = sf.fit(obs["ft"][mo], obs["fd"][mo])

    jobs = sorted(int(p.name) for p in REPO.iterdir()
                  if p.is_dir() and p.name.isdigit())
    rows = []
    for n in jobs:
        if not (IN / f"res{n}.in").exists():
            continue
        dk = sf.deck(n)
        near, far = perm_map(dk)
        mu, sig = sf.ffloat(dk.get("muinit")), sf.ffloat(dk.get("sigmainit"))
        tau0 = (mu or 0) * (sig or 0)
        pb = (sf.ffloat(dk.get("phi")) or 0) * (sf.ffloat(dk.get("beta")) or 0)
        pev = dk.get("permev", "F").upper().startswith("T")
        lam, pct, rms, tend = scores(n, dk, obs, lam_obs)
        rows.append(dict(n=n, mu=mu, sig=sig, tau0=tau0, pb=pb, pev=pev,
                         kpmax=dk.get("kpmax"), kpmin=dk.get("kpmin"),
                         kL=dk.get("kL"), dc=sf.ffloat(dk.get("dc")),
                         near=near, far=far, lam=lam, pct=pct, rms=rms,
                         tend=tend, inj=dk.get("injection_file")))
        print(f"  {n}: tau0 {tau0:.2f}  lam {lam}  press {pct}")

    L = ["<!-- INVENTORY-START: generated by docs/cooper_basin/make_readme_table.py"
         " -- do not edit by hand -->",
         "",
         "## Run inventory",
         "",
         f"Every run with a folder here. Front fits use one shared window "
         f"**0–{TCUT:.0f} d** with the observed front refit on it "
         f"(λ_obs = {lam_obs:.4f}); λ is window-dependent, so a single shared "
         f"window is the only fair comparison. Pressure error is measured while "
         f"the well is **flowing**.",
         "",
         "τ₀ = μ₀ × σ̄₀. **Understressed** means τ₀ < 15.0 MPa, Wang & Dunham's "
         "value — the whole point of the exercise. `k near / far` is read from the "
         "permeability **map file**, not the deck, because a deck reading "
         "`kp 4e-13` can be running a map with 1.1e-12 near the well.",
         "",
         "| run | μ₀ | σ̄₀ MPa | τ₀ MPa | φβ Pa⁻¹ | permev | kpmax | kpmin | kL | dc | k near / far m² | λ/λ_obs | pressure | reached |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]

    for r in sorted(rows, key=lambda z: z["n"]):
        pev = "**T**" if r["pev"] else "F"
        kmx = r["kpmax"] if r["pev"] else "—"
        kmn = r["kpmin"] if r["pev"] else "—"
        kl = r["kL"] if r["pev"] else "—"
        # A run that has not reached the shared window is fit on less data than
        # the observed front it is compared against, so its ratio is not
        # comparable to the others. Flag it rather than letting it sit in the
        # table looking like a result.
        short = r["tend"] is not None and r["tend"] < 0.95 * TCUT
        if r["lam"] is None or not np.isfinite(r["lam"]):
            lam = "no slip"
        elif short:
            lam = f"{r['lam']/lam_obs:.2f} \u2020"
        elif abs(r["lam"] / lam_obs - 1) <= 0.15:
            lam = f"**{r['lam']/lam_obs:.2f}**"
        else:
            lam = f"{r['lam']/lam_obs:.2f}"
        if r["pct"] is None:
            pr = "n/a" if r["inj"] != "june_clean.txt" else "—"
        else:
            pr = (f"**{r['pct']:+.1f}%**" if abs(r["pct"]) <= 15
                  else f"{r['pct']:+.1f}%")
        km = (f"{r['near']:.1e} / {r['far']:.1e}"
              if r["near"] is not None else "?")
        tau = (f"**{r['tau0']:.2f}**" if r["tau0"] < TAU0_TAIYI
               else f"{r['tau0']:.2f}")
        L.append(f"| [{r['n']}]({r['n']}) | {r['mu']:.2f} | {r['sig']:.2f} | {tau} "
                 f"| {r['pb']:.1e} | {pev} | {kmx} | {kmn} | {kl} | {r['dc']:.1e} "
                 f"| {km} | {lam} | {pr} | {r['tend']:.2f} d |")

    L += ["",
          "**Bold** τ₀ = understressed. **Bold** λ/λ_obs or pressure = inside the "
          "±15% tolerance. `no slip` means peak slip never reached dc anywhere, so "
          "the run has no front to fit — a result, not missing data. `n/a` in the "
          "pressure column means the run is driven by an injection file other than "
          "`june_clean.txt`, so it cannot be compared against the June 2012 "
          "measurement.",
          "",
          "\u2020 fit on less than the shared window (the run has not reached "
          f"{TCUT:.0f} d yet), so this ratio is not comparable to the others.",
          "",
          _passes(rows, lam_obs),
          "",
          "Full decks are in "
          "`/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs/res<run>.in`, "
          "and each folder's `params.txt` has the deck verbatim plus the derived "
          "quantities (diffusivity, Δp needed to fail, tmax in days, how the run "
          "terminated).",
          "",
          "<!-- INVENTORY-END -->"]

    table = "\n".join(L)
    rd = REPO / "README.md"
    txt = rd.read_text()
    if "<!-- INVENTORY-START" in txt:
        pre = txt.split("<!-- INVENTORY-START")[0]
        post = txt.split("<!-- INVENTORY-END -->")[-1]
        txt = pre + table + post
    else:
        txt = txt.rstrip() + "\n\n" + table + "\n"
    rd.write_text(txt)
    print(f"\nwrote {rd} ({len(rows)} runs)")


if __name__ == "__main__":
    main()
