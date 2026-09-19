#!/usr/bin/env python3
"""Emit the paper's parameter table as LaTeX, READ FROM THE DECK.

A hand-typed parameter table drifts from the runs it claims to describe, and
that drift is invisible until a referee recomputes something and gets a
different answer. So this reads res<n>.in and writes the .tex, and every
derived quantity is computed here rather than copied.

LAYOUT is the conventional JGR form -- symbol, description, value, units,
grouped by physics, booktabs rules, no vertical lines. It is NOT a transcript
of Taiyi's table, which I do not have; if his grouping is preferred, the
GROUPS list below is the only thing that needs reordering.

TWO CHOICES WORTH KNOWING

1. A SOURCE COLUMN. Each row is marked measured / lab / chosen / swept. This is
   not standard, and it is the single most useful thing the table can do for
   this paper: six of the parameters have no justification recorded anywhere
   (beta, phi, kpmin, kL, kT, Sw_fwid), and a table that silently presents them
   beside the measured stress invites exactly the objection it should be
   pre-empting. Use --no-source to drop it.

2. DERIVED ROWS ARE INCLUDED AND FLAGGED. tau_0, dp_crit, the two
   diffusivities and the disc radius are what a reader actually reasons with,
   and none of them is a deck key: tau_0 = muinit*sigmainit,
   dp_crit = sigmainit(1 - muinit/f0), D = kp/(eta*phi*beta), and the disc
   radius is parsed out of the permeability map's filename. They sit in a
   separate block so nothing suggests they were set independently.

SWEPT RANGES. For the parameters the study varies, the reference value is given
with the swept range beneath it, since quoting only the preferred value hides
that the paper's argument is about the range.

Usage:
  python table_params_tex.py                       # 633001, the reference
  python table_params_tex.py --run 633001 --no-source
  python table_params_tex.py --out /path/table1.tex
"""
import argparse
import re
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
OUT_DEFAULT = Path("/home/users/nberrios/3dhbi/hbi_git/docs/table_params.tex")
REF = 633001

# The study variables, as (reference, range-as-typeset). Quoting only the
# preferred value would hide that the paper's claim is about the range.
SWEPT = {
    "muinit": r"0.370--0.536",
    "kpmax": r"$2.75\times10^{-15}$--$4.0\times10^{-13}$",
    "R_disc": r"150--450",
}


def read_deck(n):
    d = {}
    for line in (IN / f"res{n}.in").read_text().splitlines():
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            d[w[0]] = " ".join(w[1:])
    return d


def ff(x):
    return float(str(x).replace("d", "e").replace("D", "e"))


def sci(x, sig=2):
    """LaTeX math for a float: 2.5e-13 -> $2.5\\times10^{-13}$."""
    if x == 0:
        return "0"
    s = f"{x:.{sig}e}"
    m, e = s.split("e")
    m = m.rstrip("0").rstrip(".")
    e = int(e)
    if -2 <= e <= 3:
        v = x
        return f"${v:g}$"
    return rf"${m}\times10^{{{e}}}$"


def build(n, source=True):
    d = read_deck(n)
    ds_m = ff(d["ds"]) * 1000.0
    imax = int(d["imax"])
    eta, phi, beta = ff(d["eta"]), ff(d["phi"]), ff(d["beta"])
    kmin, kmax = ff(d["kpmin"]), ff(d["kpmax"])
    mu0, sig0, f0 = ff(d["muinit"]), ff(d["sigmainit"]), ff(d["f0"])
    den = eta * phi * beta
    m = re.search(r"disc(\d+)", d.get("parameter_file", ""))
    disc = int(m.group(1)) if m else None
    tau0 = mu0 * sig0
    dpc = sig0 * (1.0 - mu0 / f0)

    # (symbol, description, value, units, source)
    GROUPS = [
        ("Elasticity and geometry", [
            (r"$G$", "Shear modulus", f"${ff(d['rigid']):g}$", r"GPa", "chosen"),
            (r"$\Delta s$", "Cell size", f"${ds_m:g}$", r"m", "chosen"),
            (r"$N$", "Cells per side", rf"${imax}\times{imax}$", "--", "chosen"),
            (r"$L$", "Fault half-length", f"${imax*ds_m/2:.0f}$", r"m", "chosen"),
        ]),
        ("Initial stress", [
            (r"$\bar\sigma_0$", "Effective normal stress",
             f"${sig0:g}$", r"MPa", "measured"),
            (r"$\mu_0$", "Initial friction coefficient",
             f"${mu0:.3f}$", "--", "swept"),
        ]),
        ("Rate-and-state friction", [
            (r"$a$", "Direct-effect parameter", f"${ff(d['a']):g}$", "--", "chosen"),
            (r"$b$", "Evolution-effect parameter", f"${ff(d['b']):g}$", "--", "chosen"),
            (r"$D_c$", "State-evolution distance",
             sci(ff(d["dc"])), r"m", "chosen"),
            (r"$f_0$", "Reference friction coefficient",
             f"${f0:.2f}$", "--", "lab"),
            (r"$V_0$", "Initial slip rate",
             sci(ff(d["velinit"])), r"m\,s$^{-1}$", "chosen"),
        ]),
        ("Fluid and pore properties", [
            (r"$\eta$", "Fluid viscosity", sci(eta), r"Pa\,s", "measured"),
            (r"$\beta$", "Fluid compressibility",
             sci(beta), r"Pa$^{-1}$", "chosen"),
            (r"$\phi$", "Porosity", f"${phi:g}$", "--", "chosen"),
        ]),
        ("Fault-zone permeability", [
            (r"$k_{p,\min}$", "Background permeability",
             sci(kmin), r"m$^2$", "chosen"),
            (r"$k_{p,\max}$", "Enhanced permeability",
             sci(kmax), r"m$^2$", "swept"),
            (r"$R_{\rm disc}$", "Initial enhanced-zone radius",
             f"${disc}$" if disc else "--", r"m", "swept"),
            (r"$k_L$", "Permeability-enhancement slip scale",
             sci(ff(d["kL"])), r"m", "chosen"),
            (r"$k_T$", "Permeability-healing time scale",
             sci(ff(d["kT"])), r"s", "chosen"),
        ]),
        ("Wellbore", [
            (r"$r_w$", "Wellbore radius", f"${ff(d['rw']):.3f}$", r"m", "measured"),
            (r"$S_w$", "Wellbore storage per unit fault width",
             sci(ff(d["Sw_fwid"])), r"m$^2$\,Pa$^{-1}$", "chosen"),
            (r"$s$", "Skin factor", f"${ff(d['skin']):g}$" if "skin" in d
             else "$0$", "--", "chosen"),
        ]),
    ]
    DERIVED = [
        (r"$\tau_0$", r"Initial shear stress, $\mu_0\bar\sigma_0$",
         f"${tau0:.2f}$", r"MPa"),
        (r"$\Delta p_{\rm crit}$",
         r"Overpressure to reach failure, $\bar\sigma_0(1-\mu_0/f_0)$",
         f"${dpc:.2f}$", r"MPa"),
        (r"$D_{\rm disc}$", r"Diffusivity in the enhanced zone, "
         r"$k_{p,\max}/\eta\phi\beta$",
         f"${kmax/den:.2f}$", r"m$^2$\,s$^{-1}$"),
        (r"$D_{\rm far}$", r"Background diffusivity, $k_{p,\min}/\eta\phi\beta$",
         f"${kmin/den:.4f}$", r"m$^2$\,s$^{-1}$"),
        (r"$k_{p,\max}/k_{p,\min}$", "Permeability contrast",
         f"${kmax/kmin:.0f}$", "--"),
    ]

    ncol = 5 if source else 4
    spec = "llrl" + ("l" if source else "")
    head = [r"\begin{table}[t]", r"\centering",
            r"\caption{Model parameters. Values are those of the reference "
            r"simulation; ranges beneath the swept parameters give the span "
            r"explored. Derived quantities are not independent inputs.}",
            r"\label{tab:params}",
            rf"\begin{{tabular}}{{{spec}}}", r"\toprule"]
    cols = ["Symbol", "Description", "Value", "Units"]
    if source:
        cols.append("Source")
    head.append(" & ".join(cols) + r" \\")
    head.append(r"\midrule")

    body = []
    for title, rows in GROUPS:
        body.append(rf"\multicolumn{{{ncol}}}{{l}}{{\textit{{{title}}}}} \\")
        for sym, desc, val, un, src in rows:
            key = {"$\\mu_0$": "muinit", "$k_{p,\\max}$": "kpmax",
                   "$R_{\\rm disc}$": "R_disc"}.get(sym)
            cells = [sym, desc, val, un] + ([src] if source else [])
            body.append(" & ".join(cells) + r" \\")
            if key in SWEPT:
                rng = [" ", r"\quad swept range", SWEPT[key], un] \
                    + ([" "] if source else [])
                body.append(" & ".join(rng) + r" \\")
        body.append(r"\addlinespace")
    body.append(rf"\multicolumn{{{ncol}}}{{l}}{{\textit{{Derived}}}} \\")
    for sym, desc, val, un in DERIVED:
        cells = [sym, desc, val, un] + ([" "] if source else [])
        body.append(" & ".join(cells) + r" \\")

    tail = [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(head + body + tail) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=int, default=REF)
    ap.add_argument("--no-source", action="store_true")
    ap.add_argument("--out", default=str(OUT_DEFAULT))
    a = ap.parse_args(argv)
    tex = build(a.run, source=not a.no_source)
    Path(a.out).write_text(tex)
    # structural checks -- there is no pdflatex on the cluster
    assert tex.count(r"\begin{tabular}") == tex.count(r"\end{tabular}") == 1
    assert tex.count(r"\begin{table}") == tex.count(r"\end{table}") == 1
    assert tex.count("$") % 2 == 0, "unbalanced math delimiters"
    ncol = 4 if a.no_source else 5
    for ln in tex.splitlines():
        if ln.endswith(r"\\") and r"\multicolumn" not in ln:
            got = ln.count("&") + 1
            assert got == ncol, f"{got} cells, expected {ncol}: {ln}"
    print(tex)
    print(f"% wrote {a.out}  (run {a.run}, "
          f"{'no ' if a.no_source else ''}source column)")


if __name__ == "__main__":
    main()
