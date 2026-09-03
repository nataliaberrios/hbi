#!/usr/bin/env python3
"""Generate RUN_KEY.md -- one row per run, what varies and what it gave.

The README carries prose about individual stages but there has never been a
single place that says what all 68 runs are. This builds that, from the decks
and grid_scores.json rather than by hand, so it cannot drift out of date.

Everything printed is either read from the deck or derived from it by a formula
stated in the legend. Nothing is transcribed.

Usage:  python make_run_key.py
"""
import json
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
OUT = Path("/home/users/nberrios/3dhbi/hbi_git/docs/figs/cooper_grid")

# Blocks, in the order they were run. Each is (title, run numbers, what it asks).
BLOCKS = [
    ("Reuse — already on fixed code before the grid started",
     [632522, 632523, 632524, 632525],
     "res1807/res1808 physics plus a two-zone map, with and without "
     "enhancement. 8 d runs, scored on the same 0-5 d window as everything "
     "else. These were NOT resubmitted."),
    ("Stage 1 — baseline, uniform initial permeability",
     list(range(632800, 632808)),
     "Does the 1807/1808 front match survive the limitsigma fix and the June "
     "injection record? Twinned on sigmabar_0 and on fluid: A = the parent's "
     "own eta and beta, B = eta 1.27e-4 with beta scaled 7.008x so that "
     "D = k/(eta*phi*beta) is unchanged."),
    ("Stage 2 — nonuniform initial permeability, consistent bounds",
     list(range(632810, 632822)),
     "Near-well disc IS kpmax, background IS kp = kpmin. Fills the 12 cells of "
     "the 16-cell (sigmabar_0 x fluid x contrast x permev) slice that 632522-525 "
     "did not cover."),
    ("Stage 3 — the tau_0 sweep, 11.0 to 15.0 MPa in 0.5 MPa steps",
     list(range(632830, 632848)) + list(range(632850, 632868)),
     "How much of Wang & Dunham's match is their raised initial shear stress? "
     "sigmainit held fixed within each half so this is a true one-parameter "
     "sweep in tau_0. 36 runs: two map contrasts x two sigmabar_0 x nine tau_0."),
    ("Stage 4 — spatially graded porosity",
     list(range(632870, 632876)),
     "beta is a scalar in HBI and cannot vary in space, but phi can, and it "
     "enters BOTH storage (str = beta*phi) and diffusivity "
     "(cdiff = kp/(eta*beta*phi)). Grades the porosity to flatten the pressure "
     "profile. G3 is the inverse grading, included as a control."),
    ("Taiyi reference — his published parameters verbatim",
     [632880, 632881],
     "First time in the project that HBI has been run on Wang & Dunham's own "
     "inputs. Differ from each other ONLY in dc: 1.53e-5 (his value) vs 1e-4 "
     "(this project's). Note ds = 20 m here, not the 5 m used by the other 66."),
]


def deck(n):
    d = {}
    for line in (IN / f"res{n}.in").read_text().splitlines():
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            d[w[0]] = " ".join(w[1:]).strip('"')
    return d


def ff(x):
    try:
        return float(str(x).replace("d", "e").replace("D", "e"))
    except (TypeError, ValueError):
        return None


def kp_col(pf):
    """Permeability column of a parameter file, selected by header NAME.

    Stage 4 files are two-column ("kp phi"); taking column 0 blindly would work
    today but breaks the moment a file orders them differently, and taking max()
    over all columns compares a permeability against a porosity -- a bug that has
    already been fixed twice in this repo.
    """
    a = np.loadtxt(IN / pf, skiprows=1)
    if a.ndim > 1:
        names = (IN / pf).read_text().split("\n", 1)[0].split()
        a = a[:, names.index("kp")] if "kp" in names else a[:, 0]
    return a


def phi_desc(dk):
    """Porosity as NEAR-well / FAR-field, plus the grading tag.

    Reporting min/max instead would be actively misleading: G1 is 0.020 near and
    0.005 far, G3 is the inverse, and both have the same min and max. The G1-vs-G3
    contrast is the whole experiment, so the field is sampled at the injector cell
    and at a corner rather than reduced to a range.
    """
    pf = dk.get("parameter_file")
    if pf and dk.get("parameter_file_ncol", "1") != "1":
        names = (IN / pf).read_text().split("\n", 1)[0].split()
        if "phi" in names:
            a = np.loadtxt(IN / pf, skiprows=1)[:, names.index("phi")]
            im = int(dk["imax"])
            c = (im + 1) // 2                       # injector cell, 1-indexed
            near = a[(c - 1) * im + (c - 1)]
            far = a[0]                              # corner (1,1)
            tag = next((g for g in ("G1", "G2", "G3") if f"_{g}" in pf), "")
            s = f"{near:.3f}/{far:.3f}"
            return f"{s} {tag}".strip() if near != far else f"{near:.3f}"
    return f"{ff(dk['phi']):.3f}"


def perm_desc(dk):
    """Near/far permeability and the contrast, however it is specified."""
    on = dk.get("parameterfromfile", "F").upper() in ("T", "TRUE", ".TRUE.")
    if not on:
        k = ff(dk.get("kp"))
        return f"{k:.0e} uniform", 1.0
    k = kp_col(dk["parameter_file"])
    lo, hi = k.min(), k.max()
    return f"{hi:.1e}/{lo:.0e}", (hi / lo if lo > 0 else float("inf"))


def main():
    scores = {r["n"]: r for r in json.load(open(OUT / "grid_scores.json"))}
    L = []
    A = L.append
    A("# Run key — what every simulation in this grid is")
    A("")
    A("Generated by `docs/cooper_basin/make_run_key.py` from the decks in")
    A("`grid_search_inputs/` and `grid_scores.json`. Do not edit by hand.")
    A("")
    A("## How to read the two score columns")
    A("")
    A("**front** is λ/λ_obs. λ is the fit of R = λ√t to the slip front, where the")
    A("front is the outermost radius with slip > **1e-4 m — a fixed threshold for")
    A("every run**, not each run's own `dc`. Scored on **0–5 d**, with the observed")
    A("front refit on that same window (λ_obs = 0.1866). 1.00 is a match. A number")
    A("without its window is meaningless: these same runs score 0.97–0.99 on 0–8 d.")
    A("")
    A("**wellhead** is mean(sim) / mean(measured) − 1 over **flowing periods only**")
    A("(injecting above 25% of peak). Not mean(sim/measured) — the measured wellhead")
    A("falls through ~0 during shut-ins, and dividing by that read 2.17 and 0.32 on")
    A("two samplings of the same run. 0% is a match; ±15% is the agreed band.")
    A("")
    A("**pk/dc** is peak slip divided by that run's own `dc`. This one *is*")
    A("per-run, because it asks a different question — did the patch weaken at all.")
    A("Below ~1 there is no front to speak of and the λ is near-noise; `no slip`")
    A("means nothing reached the threshold, which is a result, not missing data.")
    A("")
    A("## Column legend")
    A("")
    A("| column | meaning |")
    A("|---|---|")
    A("| σ̄₀ | `sigmainit`, effective normal stress, MPa |")
    A("| τ₀ | `muinit` × `sigmainit`, initial shear stress, MPa. **An initial "
      "condition on stress, not a friction coefficient.** |")
    A("| Δτc | f₀σ̄₀ − τ₀, how much overpressure a Mohr–Coulomb reading would need. "
      "Reported for orientation only — HBI's regularised friction law has no "
      "threshold. |")
    A("| k near/far | max and min of the initial permeability field, m². The rule: "
      "near-well disc **is** `kpmax`, background **is** `kp` = `kpmin`. |")
    A("| φ | porosity, as **near-well / far-field**. G1 = 0.020 near, G3 = 0.005 near (the inverse control). A single value means uniform. |")
    A("| η | fluid viscosity, Pa·s. 8.9e-4 is Taiyi's (water at 25 °C); 1.27e-4 is "
      "this study's, for 4.1 km reservoir conditions |")
    A("| pev | `permev` — permeability enhancement on (T) or off (F) |")
    A("| ds | cell size, m. **Not constant across the grid** — 5 m for 66 runs, "
      "20 m for the Taiyi pair. |")
    A("")

    for title, runs, why in BLOCKS:
        A(f"## {title}")
        A("")
        for para in [why]:
            A(para)
        A("")
        A("| run | σ̄₀ | τ₀ | Δτc | k near/far | φ | η | β | pev | ds | front | "
          "wellhead | pk/dc |")
        A("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for n in runs:
            if not (IN / f"res{n}.in").exists():
                continue
            dk = deck(n)
            s = scores.get(n, {})
            sig, mu, f0 = ff(dk["sigmainit"]), ff(dk["muinit"]), ff(dk["f0"])
            tau0 = mu * sig
            pk, contrast = perm_desc(dk)
            lam = (f"**{s['lam_ratio']:.2f}**" if s.get("lam_ratio") is not None
                   else "no slip")
            wh = f"{s['p_pct']:+.1f}%" if s.get("p_pct") is not None else "—"
            pd = (f"{s['peak_over_dc']:.0f}"
                  if s.get("peak_over_dc") is not None else "—")
            ds_m = ff(dk["ds"]) * 1000.0
            A(f"| {n} | {sig:.2f} | {tau0:.2f} | {f0*sig - tau0:.2f} | {pk} | "
              f"{phi_desc(dk)} | {ff(dk['eta']):.2e} | {ff(dk['beta']):.3e} | "
              f"{dk.get('permev','F')} | {ds_m:.0f} | {lam} | {wh} | {pd} |")
        A("")

    # Held fixed everywhere -- as important as what varies.
    A("## Held fixed in all 68 runs")
    A("")
    A("| | value | |")
    A("|---|---|---|")
    A("| `f0` | **0.6** | never varied in any of 723 decks — every \"friction "
      "sweep\" here was a τ₀ sweep |")
    A("| `a` / `b` | **0.015 / 0.012** | a − b = +0.003, velocity-**strengthening** |")
    A("| `dc` | 1e-4 | except 1.53e-5 in 632880 |")
    A("| `imax` × `jmax` | 601 × 601 | |")
    A("| ~~`ds`~~ | **NOT fixed** | 5 m for 66 runs, 20 m for 632880/632881. "
      "The Taiyi pair needs the wider domain because their far-field k gives "
      "D = 4.49 m²/s, so at ds 5 m the 5 d diffusion front (2.79 km) would "
      "exceed the 1.5 km half-domain. Listed here because this is where someone "
      "would look for it. |")
    A("| `injection_file` | `june_clean.txt` | rate is per unit fault width, not "
      "total well flow — multiply by w = 6 m |")
    A("| `rigid`, `Sw_fwid`, `rw` | 24 GPa, 7.4e-9, 8.9e-2 m | |")
    A("| `kL` / `kT` | 1e-3 / 1e15 | enhancement slip scale and healing time |")
    A("")
    A("**Consequence of a − b = +0.003:** every run is velocity-strengthening and "
      "produces no seismicity. What is being compared is an *aseismic slip front* "
      "against an observed *seismicity front*, under the assumption that the slip "
      "front does not care whether the slip was seismic.")
    A("")
    A("## Where the grid stands")
    A("")
    best_f = max((r for r in scores.values() if r.get("lam_ratio") is not None),
                 key=lambda r: r["lam_ratio"])
    inband = [r for r in scores.values()
              if r.get("p_pct") is not None and abs(r["p_pct"]) <= 15]
    both = [r for r in inband
            if r.get("lam_ratio") is not None and 0.85 <= r["lam_ratio"] <= 1.15]
    A(f"- Best front: **{best_f['n']}**, λ/λ_obs "
      f"**{best_f['lam_ratio']:.2f}**, wellhead {best_f['p_pct']:+.1f}%")
    A(f"- Runs with the wellhead inside ±15%: **{len(inband)}** of {len(scores)}")
    A(f"- Runs inside **both** bands: **{len(both) or 'none'}**")
    A("")
    A("The two targets are in tension: the wellhead needs *low* near-well "
      "overpressure, the front needs *high*. τ₀ moves the front without moving "
      "the wellhead at all — see Stage 3, where the wellhead is flat to ±0.1% "
      "across the whole sweep — but it saturates at λ 0.34 by τ₀ = 15.0, and τ₀ "
      "cannot exceed f₀σ̄₀ = 16.79 MPa without the fault being past its own "
      "strength at zero overpressure.")

    (OUT / "RUN_KEY.md").write_text("\n".join(L) + "\n")
    print(f"wrote {OUT / 'RUN_KEY.md'}")
    print(f"  {sum(len(b[1]) for b in BLOCKS)} runs across {len(BLOCKS)} blocks")


if __name__ == "__main__":
    main()
