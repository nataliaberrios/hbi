#!/usr/bin/env python3
"""Label every pre-grid run folder with whether its permeability bounds are valid.

docs/figs/cooper_grid/ contains, alongside the 68 scored grid runs, ~90 folders
from earlier work (632501-632595). Some of those runs are misconfigured and some
are load-bearing, and nothing on disk said which was which. This labels them.

THE RULE BEING CHECKED. kpmax is the ceiling, kp = kpmin is the floor and the
initial value; when the initial permeability is nonuniform the near-well disc IS
kpmax and the background IS kpmin. A deck with kpmax above anything present in
its map is asking permeability to evolve toward a value the medium never
contains, and the enhancement term
    -vel/kL*(kp - kpmax)        main_LH.f90:2413
then drives kp at a rate set by a ceiling that has no physical basis in that run.
Across 112 archived runs, the 35 obeying this rule had 29 reach dc; the 77 that
did not had 2.

WHY NOT JUST DELETE THE BAD ONES. Because the split is not clean. Of the ~90:
permev F runs are not subject to the rule at all, and among them are 632503 and
632507 (the runs that validated HBI's pressure solver against Wang & Dunham's
code, bias -0.011 MPa and RMS 0.148 at ds 20 m) and 632520 (the deck res632880
and res632881 were built from, which are in turn the parents of 632884/632885).
Deleting by run number would take the provenance of current work with it.

Writes two things, so a label is visible whether you are reading the directory
or browsing one folder:
  BOUNDS_AUDIT.md          master table, all folders, sorted by verdict
  <run>/README.md          per-folder verdict, which GitHub renders inline

Usage:  python audit_bounds.py [--write]
"""
import argparse
import json
import re
import subprocess
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
REPO = Path("/home/users/nberrios/3dhbi/hbi_git")
FIGS = REPO / "docs" / "figs" / "cooper_grid"

OK, BAD, NA, NODECK = "VALID", "MISCONFIGURED", "NOT APPLICABLE", "NO DECK"


def deck(n):
    p = IN / f"res{n}.in"
    if not p.exists():
        return None, None
    d, hdr = {}, None
    for line in p.read_text().splitlines():
        if line.startswith("!"):
            if hdr is None and len(line) > 2:
                hdr = line.lstrip("! ").rstrip()
            continue
        if not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            d[w[0]] = " ".join(w[1:]).strip('"')
    return d, hdr


def ff(x):
    try:
        return float(str(x).replace("d", "e").replace("D", "e"))
    except (TypeError, ValueError):
        return None


def map_bounds(d):
    """(max, min) of the initial permeability field, however it is specified."""
    if d.get("parameterfromfile", "F").upper() in ("T", "TRUE", ".TRUE."):
        pf = d.get("parameter_file")
        if not pf or not (IN / pf).exists():
            return None, None
        a = np.loadtxt(IN / pf, skiprows=1)
        if a.ndim > 1:
            names = (IN / pf).read_text().split("\n", 1)[0].split()
            a = a[:, names.index("kp")] if "kp" in names else a[:, 0]
        return float(a.max()), float(a.min())
    k = ff(d.get("kp"))
    return k, k


def audit(n):
    d, hdr = deck(n)
    if d is None:
        return dict(n=n, verdict=NODECK, note="no deck in grid_search_inputs",
                    hdr=None)
    if d.get("permev", "F").upper() != "T":
        return dict(n=n, verdict=NA, hdr=hdr,
                    note="permev F — permeability is fixed, so kpmax/kpmin are "
                         "never used and the rule does not apply")
    mx, mn = map_bounds(d)
    kx, km = ff(d.get("kpmax")), ff(d.get("kpmin"))
    if mx is None:
        return dict(n=n, verdict=BAD, hdr=hdr,
                    note="permev T but the permeability field cannot be resolved")
    if kx is None or km is None:
        return dict(n=n, verdict=BAD, hdr=hdr,
                    note=f"permev T but the deck declares no "
                         f"{'kpmax' if kx is None else 'kpmin'}")
    okx, okm = abs(kx - mx) / mx < 1e-6, abs(km - mn) / mn < 1e-6
    if okx and okm:
        return dict(n=n, verdict=OK, hdr=hdr,
                    note=f"kpmax {kx:.3e} == map max, kpmin {km:.3e} == map min")
    parts = []
    if not okx:
        parts.append(f"**kpmax {kx:.2e} but the map never exceeds {mx:.2e}** "
                     f"({kx/mx:.0f}x above it)")
    if not okm:
        parts.append(f"kpmin {km:.2e} != map min {mn:.2e}")
    return dict(n=n, verdict=BAD, hdr=hdr, note="; ".join(parts))


def referenced(n):
    """Is this run cited in prose anywhere in the repo? Those are load-bearing."""
    try:
        out = subprocess.run(
            ["grep", "-rl", "--include=*.md", "--include=*.py", str(n), "docs"],
            cwd=REPO, capture_output=True, text=True, timeout=120).stdout
    except (subprocess.SubprocessError, OSError):
        return []
    # A run's own folder mentioning itself is not a citation, and neither is this
    # script's own docstring, which names the range "632501-632595" and therefore
    # matched the endpoints as if they were cited.
    return sorted({f for f in out.split()
                   if f and f"cooper_grid/{n}/" not in f
                   and not f.endswith(("audit_bounds.py", "BOUNDS_AUDIT.md"))})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    grid = {r["n"] for r in json.load(open(FIGS / "grid_scores.json"))}
    dirs = sorted(int(p.name) for p in FIGS.iterdir()
                  if p.is_dir() and re.fullmatch(r"\d{6}", p.name)
                  and int(p.name) not in grid)
    print(f"{len(dirs)} folders to label ({len(grid)} scored grid runs excluded)")

    rows = [audit(n) for n in dirs]
    for r in rows:
        r["refs"] = referenced(r["n"])
    by = {v: [r for r in rows if r["verdict"] == v] for v in (BAD, OK, NA, NODECK)}
    for v in (BAD, OK, NA, NODECK):
        print(f"  {v:<16s} {len(by[v]):>3d}")
    load = [r for r in rows if r["refs"]]
    print(f"  cited in repo prose (load-bearing): {[r['n'] for r in load]}")

    if not a.write:
        print("\ndry run — nothing written. re-run with --write")
        return

    L = ["# Bounds audit — the run folders that predate the grid", "",
         "Generated by `docs/cooper_basin/audit_bounds.py`. Do not edit by hand.",
         "",
         f"These **{len(dirs)}** folders are earlier work, kept for provenance. They",
         f"are **not** part of the {len(grid)} scored runs in "
         "[`RUN_KEY.md`](RUN_KEY.md) and none of them appear in `grid_scores.json`.",
         "",
         "## The rule",
         "",
         "`kpmax` is the ceiling, `kp` = `kpmin` is the floor and the initial value.",
         "When the initial permeability is nonuniform the near-well disc **is**",
         "`kpmax` and the background **is** `kpmin`. A deck with `kpmax` above",
         "anything in its map drives the enhancement term",
         "",
         "    tmp = -vel(i)/kL*(kp(i) - kpmax) - (kp(i) - kpmin)/kT    main_LH.f90:2413",
         "",
         "toward a value the medium never contains. This is not bookkeeping: across",
         "112 archived runs, the 35 obeying the rule had 29 reach `dc`; the 77 that",
         "did not had 2.",
         "",
         "## Verdicts", "",
         f"| verdict | count | meaning |",
         "|---|---|---|",
         f"| **{BAD}** | {len(by[BAD])} | `permev T` with bounds that do not match "
         "the map. **Do not read these as results.** |",
         f"| {OK} | {len(by[OK])} | `permev T`, `kpmax` == map max and `kpmin` == "
         "map min |",
         f"| {NA} | {len(by[NA])} | `permev F` — permeability fixed, so the rule "
         "does not apply |",
         f"| {NODECK} | {len(by[NODECK])} | no deck left in `grid_search_inputs/` |",
         "",
         "## Load-bearing folders", "",
         "Cited in prose elsewhere in this repo, so they are provenance for current",
         "work rather than dead ends:", ""]
    for r in load:
        L.append(f"- **{r['n']}** ({r['verdict']}) — cited in "
                 + ", ".join(f"`{f}`" for f in r["refs"]))
    L += ["",
          "Independently of citations, **632520** is the deck `res632880.in` and",
          "`res632881.in` were built from, and those are the parents of 632884 and",
          "632885. **632503** and **632507** are the runs that validated HBI's",
          "pressure solver against Wang & Dunham's code (bias −0.011 MPa, RMS 0.148",
          "at ds 20 m) — the basis for saying the two models' diffusion agrees.",
          "**632510** is the `limitsigma` cusp test.", ""]

    for v in (BAD, OK, NA, NODECK):
        if not by[v]:
            continue
        L += [f"## {v} — {len(by[v])} folders", "",
              "| run | detail | what it was |", "|---|---|---|"]
        for r in by[v]:
            h = (r["hdr"] or "—")
            h = h[:96] + "…" if len(h) > 96 else h
            L.append(f"| [{r['n']}]({r['n']}/) | {r['note']} | {h} |")
        L.append("")
    (FIGS / "BOUNDS_AUDIT.md").write_text("\n".join(L) + "\n")
    print(f"wrote {FIGS / 'BOUNDS_AUDIT.md'}")

    # Per-folder label, so the verdict is visible when browsing one folder.
    banner = {BAD: "> # ⚠ MISCONFIGURED — do not read as a result\n>\n"
                   "> `permev T` with permeability bounds that do not match this "
                   "run's own map.",
              OK: "> # Bounds valid\n>\n> `permev T`, `kpmax` == map max and "
                  "`kpmin` == map min.",
              NA: "> # Bounds rule not applicable\n>\n> `permev F` — permeability "
                  "is fixed, so `kpmax`/`kpmin` are never used.",
              NODECK: "> # No deck\n>\n> The input deck is no longer in "
                      "`grid_search_inputs/`, so this run cannot be audited or "
                      "reproduced."}
    for r in rows:
        t = [banner[r["verdict"]], "",
             f"# Run {r['n']}", "",
             f"**Verdict: {r['verdict']}** — {r['note']}", ""]
        if r["hdr"]:
            t += [f"What it was: {r['hdr']}", ""]
        t += [f"This folder predates the grid search and is **not** one of the "
              f"{len(grid)} scored runs. See [`../BOUNDS_AUDIT.md`]"
              f"(../BOUNDS_AUDIT.md) for the full audit and "
              f"[`../RUN_KEY.md`](../RUN_KEY.md) for the runs that are scored.", ""]
        if r["refs"]:
            t += ["This run **is** cited elsewhere in the repo, so it is "
                  "provenance for current work: "
                  + ", ".join(f"`{f}`" for f in r["refs"]), ""]
        (FIGS / str(r["n"]) / "README.md").write_text("\n".join(t))
    print(f"wrote {len(rows)} per-folder README.md labels")


if __name__ == "__main__":
    main()
