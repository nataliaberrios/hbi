#!/usr/bin/env python3
"""Every parameter of every cycle-2 run, in one table, read from the decks.

WHY THIS EXISTS. The arm-2 result is being shown to people, so the question
that follows is which of its numbers were measured, which are lab values, which
were swept, and which were simply chosen. That question is not answerable from
a figure, and answering it from memory is how a number's provenance quietly
becomes "we always used that".

TWO TABLES.

  1. THE VALUES. One row per run, one column per parameter, read out of
     res<n>.in. Columns where every run agrees are collapsed into a "common to
     all runs" list, so the table shows only what actually varies -- with 47
     keys per deck, a full matrix hides the five things that differ.

  2. THE PROVENANCE. One row per parameter: its value, where it came from, and
     whether it was swept. The three categories that matter are kept distinct:

       MEASURED   traceable to a field measurement or to the record
       LAB        a laboratory constraint on rock or fluid properties
       CHOSEN     set by us -- numerically, for resolution, or by assumption

     Anything whose justification is not documented in this repository is
     marked "not documented in-repo" rather than given a plausible-sounding
     source. There are several, and they are the honest answer to "where did
     this come from".

Usage:
  python params_table.py                 # print, and write PARAMS_cycle2.md
  python params_table.py --runs 632960 632963
"""
import argparse
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
OUT = Path("/home/users/nberrios/3dhbi/hbi_git/docs/figs/cycle2/PARAMS_cycle2.md")

RUNS = (list(range(632950, 632955)) + list(range(632960, 632965))
        + list(range(632970, 632975)) + list(range(632980, 632985))
        + [632990, 632991, 632992, 632993, 632994, 632995])

ARM = {**{n: "1" for n in range(632950, 632955)},
       **{n: "2" for n in range(632960, 632965)},
       **{n: "3" for n in range(632970, 632975)},
       **{n: "4" for n in range(632980, 632985)},
       **{n: "2" for n in (632990, 632991, 632992, 632993, 632994, 632995)}}

# One line per parameter: (value-as-set, category, source). Categories are
# MEASURED / LAB / CHOSEN / SWEPT. "not documented in-repo" is used literally
# where this repository does not record a justification -- those are the ones
# to be able to answer for.
PROV = [
    ("f0", "0.60", "LAB",
     "Lab friction for granite gouge is 0.69-0.74 (Zhang et al. 2022), so 0.60 "
     "already assumes some weakening. Held FIXED and not swept, at the user's "
     "instruction: well constrained by lab, uncomfortable to touch. "
     "CONCLUSIONS.md retracted an earlier f0 = 0.433 on these grounds."),
    ("sigmainit", "27.99 MPa", "MEASURED",
     "Effective normal stress on the fault from the Habanero stress "
     "measurements (Holl & Barton 2015). Cross-checked: the field-constrained "
     "Sv of ~95.3 MPa (Holl section 8.4) supports 27.99; a density-integrated "
     "Sv of 106.6 would not. Not swept -- tau_0 = muinit*sigmainit, so varying "
     "it cannot be done one key at a time."),
    ("muinit", "0.370 - 0.536", "SWEPT",
     "THE SWEEP VARIABLE. tau_0 = muinit*sigmainit spans 10.36 to 15.00 MPa, "
     "i.e. from the MEASURED shear stress to Wang & Dunham's assumed ratio. "
     "tau_0 is a result of this study, not an input to it."),
    ("eta", "1.27e-4 Pa s (arms 2-4)", "LAB",
     "Water viscosity at the reservoir's depth and temperature. Arm 1 uses "
     "0.89e-3, which is Wang & Dunham's value (setup_model.m:100; their Table 1 "
     "says 8e-4, a rounding) and is NOT appropriate for this reservoir. dp is "
     "linear in eta, so the 7.008x difference is the single largest effect in "
     "the whole study."),
    ("beta", "2.25e-8 1/Pa (arms 1-2)", "CHOSEN",
     "Fluid compressibility. Arm 3 uses 1.5768e-7 to hold diffusivity fixed "
     "against the eta correction -- that puts phi*beta 3.94x above the stated "
     "ceiling of 4e-10 and is why arm 3 is retired. Arm 4 uses 7.72e-9. "
     "NOT DOCUMENTED IN-REPO: where 2.25e-8 itself came from."),
    ("phi", "0.01", "CHOSEN",
     "Porosity. NOT DOCUMENTED IN-REPO. Enters only as the product phi*beta = "
     "2.25e-10 in the diffusivity and storage, so it is not separately "
     "identifiable from beta by any of these fits."),
    ("kpmax", "2.5e-13 m^2", "CHOSEN",
     "Near-well permeability, the value the whole disc starts at and that slip "
     "grows kp toward. CONCLUSIONS.md records that the front holds at "
     "0.98-1.03x across a 100x change in kpmax, so the front does not "
     "constrain it; the wellhead plateau does, falling from 11.90 MPa as kpmax "
     "rises. Every one of the previous 70 runs also sat at ~2.5e-13."),
    ("kpmin", "1e-15 m^2", "CHOSEN",
     "Background permeability, 250x below kpmax. NOT DOCUMENTED IN-REPO. It "
     "matters more than its obscurity suggests: the pressure drop through the "
     "kpmin background beyond the disc is what the disc-radius lever acts on."),
    ("kp", "1e-15 m^2", "CHOSEN",
     "Initial permeability away from the disc; set equal to kpmin by "
     "construction, so the map's minimum and this key agree."),
    ("kL, kT", "1e-3 m, 1e15 s", "CHOSEN",
     "Slip and time constants of the permeability evolution ODE "
     "(main_LH.f90:2426). kT = 1e15 s effectively disables healing, so kp only "
     "ever increases. NOT DOCUMENTED IN-REPO."),
    ("a, b", "0.015, 0.012", "CHOSEN",
     "Rate-and-state parameters, a-b = 0.003, velocity weakening. Stage 11 "
     "found a (not a-b) drives slip amplitude, but only ~2x and at the cost of "
     "the front. NOT DOCUMENTED IN-REPO as a lab or field value."),
    ("dc", "1e-4 m", "CHOSEN",
     "State evolution distance, and ALSO the slip threshold defining the model "
     "front. CONCLUSIONS.md flags this explicitly: 'the front is an arbitrary "
     "contour... defined as slip > dc = 1e-4 m'."),
    ("ds, imax, jmax", "10 m, 601, 601", "CHOSEN",
     "Cell size and grid. Set by the domain constraint, not by convergence: at "
     "ds 5 m the muinit sweep breaches the boundary by muinit 0.42. Half-domain "
     "3005 m; usable front ~2006 m after allowing for far-field diffusion."),
    ("pfinit", "0", "CHOSEN",
     "Initial pore-pressure perturbation. KNOWN TO BE WRONG and left "
     "uncorrected: the measured wellhead sits 1.14 MPa BELOW its pre-injection "
     "value at data-day 4.300, because the well was vented during the shut-in. "
     "A pfinit = -1.14 arm has been proposed and not run."),
    ("rw", "0.089 m", "MEASURED",
     "Wellbore radius, half the 0.178 m diameter used in the pipe-friction "
     "correction. Enters the Peaceman well-cell transmissivity."),
    ("skin", "0, or -0.5/-1.0/-2.0", "SWEPT",
     "Peaceman skin factor (m_diffusion.f90:669). Default 0 was never a "
     "measurement. Swept negative in 632990-632992 because Habanero 4 is a "
     "stimulated injector; the result was that it moves dp by only 9 points "
     "over the whole range and leaves the front untouched."),
    ("Sw_fwid", "7.4e-9", "CHOSEN",
     "Wellbore storage over fault width, in the gamma weighting at "
     "m_diffusion.f90:671. NOT DOCUMENTED IN-REPO."),
    ("injection_file", "june_clean_from_d4300.txt", "MEASURED",
     "The field rate record, shifted so t = 0 is data-day 4.300. Verified "
     "against the .mat record by compare_injection.py: residual rms 2.38 L/s, "
     "volume +0.32% at the scoring time. 632995 instead uses "
     "june_clean_from_d3558.txt, which INCLUDES the 17.8 h trickle the others "
     "discard."),
    ("W (analysis only)", "6.0 m", "CHOSEN",
     "Fault width used to turn the deck's m^2/s into L/s and the model's "
     "pressure into a wellhead value. NOT IN ANY DECK -- hard-coded at "
     "make_sweep_figures.py:55. A wrong W would rescale every rate and "
     "pressure comparison while leaving all the shapes correct."),
    ("P0, rho, H (analysis only)", "73.8 MPa, 1000 kg/m^3, 4077 m", "CHOSEN",
     "Used to convert the model's dp into an absolute wellhead pressure, "
     "P0 - rho g H = 33.805 MPa. The density assumption alone moves that by "
     "+-2 MPa against a mean observed dp of 7.9 MPa, which is why dp and not "
     "absolute pressure is scored."),
]


def read_deck(p):
    out = {}
    for line in Path(p).read_text().splitlines():
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            out[w[0]] = " ".join(w[1:])
    return out


def ff(x):
    return float(str(x).replace("d", "e").replace("D", "e"))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", type=int, default=RUNS)
    a = ap.parse_args(argv)

    decks, missing = {}, []
    for n in a.runs:
        p = IN / f"res{n}.in"
        if p.exists():
            decks[n] = read_deck(p)
        else:
            missing.append(n)
    keys = sorted({k for d in decks.values() for k in d})
    varying = [k for k in keys
               if len({d.get(k, "-") for d in decks.values()}) > 1]
    common = [k for k in keys if k not in varying and k != "filenumber"]

    L = []
    L.append("# Cycle-2 parameters, every run\n")
    L.append(f"Read from `res<n>.in` in `{IN}` by `params_table.py`. "
             f"{len(decks)} runs.\n")
    if missing:
        L.append(f"Decks not found: {missing}\n")

    L.append("\n## What varies between runs\n")
    hdr = ["run", "arm", "tau_0 MPa"] + [k for k in varying if k != "filenumber"]
    L.append("| " + " | ".join(hdr) + " |")
    L.append("|" + "---|" * len(hdr))
    for n, d in decks.items():
        tau = ff(d["muinit"]) * ff(d["sigmainit"])
        row = [str(n), ARM.get(n, "?"), f"{tau:.2f}"]
        row += [d.get(k, "—") for k in varying if k != "filenumber"]
        L.append("| " + " | ".join(row) + " |")

    L.append(f"\n## Common to all {len(decks)} runs\n")
    L.append("| key | value |")
    L.append("|---|---|")
    for k in common:
        L.append(f"| `{k}` | {next(iter(decks.values()))[k]} |")

    L.append("\n## Where each value comes from\n")
    L.append("`MEASURED` traceable to a field measurement or the record  ·  "
             "`LAB` a laboratory constraint  ·  `CHOSEN` set by us  ·  "
             "`SWEPT` varied as the study variable\n")
    L.append("| parameter | value | category | source |")
    L.append("|---|---|---|---|")
    for k, v, cat, src in PROV:
        L.append(f"| `{k}` | {v} | **{cat}** | {src} |")

    ndoc = sum(1 for _, _, _, s in PROV if "NOT DOCUMENTED IN-REPO" in s
               or "not documented in-repo" in s)
    L.append(f"\n### Summary of provenance\n")
    for cat in ("MEASURED", "LAB", "SWEPT", "CHOSEN"):
        got = [k for k, _, c, _ in PROV if c == cat]
        L.append(f"- **{cat}** ({len(got)}): {', '.join('`'+g+'`' for g in got)}")
    L.append(f"\n**{ndoc} parameters have no justification recorded anywhere in "
             f"this repository.** They are marked in the table above. That is "
             f"the honest answer to where they came from, and the list is the "
             f"thing to work through if the result is going to be defended: "
             f"`beta`, `phi`, `kpmin`, `kL`/`kT`, `a`/`b`, `Sw_fwid`.")

    txt = "\n".join(L) + "\n"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(txt)
    print(txt)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
