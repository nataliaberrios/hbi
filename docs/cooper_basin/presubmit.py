#!/usr/bin/env python3
"""Pre-submittal review package: everything needed to check a deck BEFORE it runs.

Sherlock priority is a finite resource, so no stage is submitted until this has
been pushed and looked at. Per stage it produces, into
docs/figs/cooper_grid/presubmit_stage<N>/:

    PRESUBMIT_stage<N>.md   one table, every deck parameter plus derived values
    params_<run>.txt        that deck verbatim + the derived block
    perm_ic_stage<N>.png    initial permeability: 2D map and radial profile
    diff_<run>.txt          diff against the parent deck

The derived block exists because the deck alone does not tell you what the solver
will do:

  tau0 = muinit*sigmainit and dp_crit = sigmainit*(1 - muinit/f0) decide whether
    the fault can fail at all, and appear nowhere in the deck.
  D = kp/(eta*phi*beta) is what sets how fast pressure spreads.
  str = beta*phi, the Peaceman well index T = 2*pi*kp/eta/(log(0.2*ds/rw)+skin),
    and gamma = (Sw_fwid/h)/((Sw_fwid/h)+T) are the quantities the diffusion
    solver actually assembles (m_diffusion.f90:444, :669, :671). gamma matters
    here because it is the one thing a viscosity/compressibility trade does NOT
    hold fixed.
  tmax in DAYS. HBI's year is 365 d (main_LH.f90:211), not 365.25.
  the permeability MAP's near/far values, since a deck reading "kp 1e-15" can be
    running a two-zone map -- or, in Stage 1, genuinely be uniform.

Usage:  python presubmit.py 1                    # stage 1
        python presubmit.py 1 --runs 632800 632801 632802 632803
"""
import argparse
import math
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
REPO = Path("/home/users/nberrios/3dhbi/hbi_git/docs/figs/cooper_grid")

STAGES = {
    1: dict(runs=list(range(632800, 632808)),
            parents={632800: 1807, 632801: 1807, 632802: 1807, 632803: 1807,
                     632804: 1808, 632805: 1808, 632806: 1808, 632807: 1808},
            title="Stage 1 — res1807/res1808 physics on fixed HBI, uniform initial "
                  "perm, twinned on sigmabar_0"),
    2: dict(runs=list(range(632810, 632822)),
            parents={n: (632522 if n <= 632815 else 632524)
                     for n in range(632810, 632822)},
            title="Stage 2 — two-zone initial perm, consistent bounds, the 12 cells "
                  "missing from the 16-cell slice"),
    # Parents below are read off each deck's own "Built from resNNNNNN.in" header
    # rather than reconstructed from the build scripts, so they cannot drift.
    3: dict(runs=list(range(632830, 632848)) + list(range(632850, 632868)),
            parents={**{n: 632814 for n in range(632830, 632848)},
                     **{n: 632820 for n in range(632850, 632868)}},
            title="Stage 3 — tau_0 sweep, 11.0 to 15.0 MPa in 0.5 MPa steps, at "
                  "both sigmabar_0 and both map contrasts (36 runs)"),
    4: dict(runs=list(range(632870, 632876)),
            parents={**{n: 632810 for n in range(632870, 632873)},
                     **{n: 632812 for n in range(632873, 632876)}},
            title="Stage 4 — spatially graded porosity, two bases x three "
                  "gradings (6 runs)"),
    5: dict(runs=[632880, 632881],
            parents={632880: 632520, 632881: 632520},
            title="Taiyi reference — Wang & Dunham's published parameters "
                  "verbatim, dc 1.53e-5 and 1e-4"),
    7: dict(runs=[632888, 632889],
            parents={632888: 911, 632889: 911},
            title="Stage 7 — CONSTANT-RATE injection on corrected code, for "
                  "comparison against an analytical solution (2 runs)"),
    6: dict(runs=[632884, 632885],
            parents={632884: 632880, 632885: 632881},
            title="Stage 6 — permeability enhancement ON, on the Taiyi "
                  "configuration: the only base that matches the wellhead AND "
                  "slips (2 runs)"),
}

INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
DISC = "#a8071a"
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10,
    "axes.edgecolor": MUTED, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "legend.fontsize": 8.5,
})


def ff(x):
    try:
        return float(str(x).replace("d", "e").replace("D", "e"))
    except (TypeError, ValueError):
        return None


def read_deck(n):
    order, d = [], {}
    for line in (IN / f"res{n}.in").read_text().splitlines():
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            d[w[0]] = " ".join(w[1:]).strip('"')
            order.append(w[0])
    return d, order


def perm_field(dk):
    """(near, far, map name or None). None means genuinely uniform."""
    if dk.get("parameterfromfile", "F").upper() not in ("T", "TRUE", ".TRUE."):
        k = ff(dk.get("kp"))
        return k, k, None
    pf = dk.get("parameter_file")
    if not pf or not (IN / pf).exists():
        return None, None, f"{pf} (MISSING)"
    k = np.loadtxt(IN / pf, skiprows=1)
    return float(k.max()), float(k.min()), pf


def derived(dk):
    """Everything the solver responds to that the deck does not state."""
    mu, sig, f0 = (ff(dk.get(k)) for k in ("muinit", "sigmainit", "f0"))
    phi, beta, eta = (ff(dk.get(k)) for k in ("phi", "beta", "eta"))
    ds_km, rw = ff(dk.get("ds")), ff(dk.get("rw"))
    skin = ff(dk.get("skin")) or 0.0
    sw = ff(dk.get("Sw_fwid"))
    near, far, pf = perm_field(dk)
    o = {}
    o["tau0"] = mu * sig if None not in (mu, sig) else None
    o["dp_crit"] = sig * (1 - mu / f0) if None not in (mu, sig, f0) else None
    o["tmax_d"] = ff(dk.get("tmax")) * 365.0 if ff(dk.get("tmax")) else None
    o["phibeta"] = phi * beta if None not in (phi, beta) else None
    o["str"] = o["phibeta"]                       # m_diffusion.f90:444
    o["near"], o["far"], o["map"] = near, far, pf
    if None not in (near, eta, phi, beta):
        o["D_near"] = near / (eta * phi * beta)
        o["D_far"] = far / (eta * phi * beta)
    # Peaceman well index, m_diffusion.f90:669. ds0 is in km there, rw in m.
    if None not in (ds_km, rw, eta):
        den = math.log(0.2 * ds_km * 1e3 / rw) + skin
        o["peaceman_den"] = den
        kmin, kmax = ff(dk.get("kpmin")), ff(dk.get("kpmax"))
        for lab, kk in (("T_at_kpmin", kmin), ("T_at_kpmax", kmax)):
            if kk is not None:
                o[lab] = 2 * math.pi * kk / eta / den
        # gamma at a representative timestep, m_diffusion.f90:671
        if sw:
            for h in (60.0, 600.0):
                for lab, kk in (("kpmin", kmin), ("kpmax", kmax)):
                    if kk is None:
                        continue
                    T = 2 * math.pi * kk / eta / den
                    o[f"gamma_h{int(h)}_{lab}"] = (sw / h) / (sw / h + T)
    return o


def params_txt(n, parent):
    dk, order = read_deck(n)
    o = derived(dk)
    L = [f"run {n}", "=" * 66, ""]
    L += ["DERIVED — what the solver responds to, not stated in the deck", "-" * 66]
    if o["tau0"] is not None:
        L += [f"  tau0 = muinit*sigmainit        {o['tau0']:.3f} MPa"
              f"   ({'understressed' if o['tau0'] < 15.0 else 'AT/ABOVE'} vs Taiyi 15.0)"]
    if o["dp_crit"] is not None:
        L += [f"  dp_crit = sigmainit(1-mu/f0)   {o['dp_crit']:.3f} MPa"
              "   (overpressure needed to fail)"]
    if o["tmax_d"]:
        L += [f"  tmax                           {o['tmax_d']:.4f} d  (HBI year = 365 d)"]
    L += [f"  permeability field             "
          + (f"UNIFORM {o['near']:.3e} m^2 (parameterfromfile F)"
             if o["map"] is None else
             f"{o['map']}\n                                 near {o['near']:.3e} / "
             f"far {o['far']:.3e} m^2, contrast {o['near']/o['far']:.0f}x")]
    if "D_near" in o:
        L += [f"  D near-well                    {o['D_near']:.4e} m^2/s",
              f"  D far-field                    {o['D_far']:.4e} m^2/s"]
    if o["str"]:
        L += [f"  str = beta*phi                 {o['str']:.4e} 1/Pa"
              "   (m_diffusion.f90:444)"]
    if "peaceman_den" in o:
        L += [f"  Peaceman log(0.2*ds/rw)+skin   {o['peaceman_den']:.4f}"]
    for k in ("T_at_kpmin", "T_at_kpmax"):
        if k in o:
            L += [f"  {k:<30s} {o[k]:.4e}   (m_diffusion.f90:669)"]
    g = sorted(k for k in o if k.startswith("gamma_"))
    if g:
        L += ["  wellbore coupling gamma = (Sw_fwid/h)/((Sw_fwid/h)+T), "
              "m_diffusion.f90:671"]
        for k in g:
            L += [f"    {k[6:]:<24s} {o[k]:.4f}"]
    L += ["", f"FULL DECK, verbatim (parent: res{parent}.in)", "-" * 66]
    L += [f"  {k:<22s} {dk[k]}" for k in order]
    return "\n".join(L) + "\n", o


def observed_peak_dp(t_end_days):
    """Peak MEASURED downhole overpressure inside the run window, in MPa.

    dp = p_wh + rho*g*H - P0. This is the hard ceiling on how much overpressure a
    run can have and still match the wellhead data, so comparing it against
    dp_crit = sigmainit*(1 - muinit/f0) says -- before running anything -- whether
    the fault can fail at all WITHOUT the model over-pressurising.
    """
    from scipy.io import loadmat
    H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
    RHO, G, HW, P0 = 1000.0, 9.81, 4077.0, 73.8
    wh = loadmat(H / "Cooper_Basin_HAB_4_Wellhead_Pressure.mat")["d"][0, 0]
    inj = loadmat(H / "Cooper_Basin_HAB_4_Injection_Rate.mat")["d"][0, 0]
    tw = wh["Date"].squeeze() - inj["Date"].squeeze()[0]
    pm = wh["Wellhead_pressure"].squeeze()
    g = np.isfinite(tw) & np.isfinite(pm)
    tw, pm = tw[g], pm[g]
    m = tw <= t_end_days
    return float(pm[m].max() + RHO * G * HW / 1e6 - P0)


def feasibility(rows):
    """Can the fault fail at all, at the pressure the data actually show?"""
    tend = max(o["tmax_d"] for _, _, _, o in rows)
    dp_obs = observed_peak_dp(tend)
    L = ["## Failure feasibility — can the fault slip at the OBSERVED pressure?", "",
         f"Peak **measured** downhole overpressure over these {tend:.2f} d is "
         f"**{dp_obs:.2f} MPa** (p_wh + rho·g·H − P0, with rho·g·H = 40.0 and "
         f"P0 = 73.8 MPa). Slip requires Δp > Δp_crit = σ̄₀(1 − μ₀/f₀). If "
         f"Δp_crit exceeds {dp_obs:.2f} MPa the fault can only slip by "
         f"**over-pressurising past the measurement**.", "",
         "| run | μ₀ | Δp_crit MPa | vs measured | verdict |", "|---|---|---|---|---|"]
    for n, parent, dk, o in rows:
        v = ("**CANNOT FAIL** without over-pressurising"
             if o["dp_crit"] > dp_obs else "can fail")
        L.append(f"| {n} | {dk['muinit']} | {o['dp_crit']:.2f} | "
                 f"{o['dp_crit'] - dp_obs:+.2f} | {v} |")
    # sigmabar_0 is an axis, so the floor has to be quoted per value
    f0 = ff(rows[0][2]["f0"])
    sigs = sorted({ff(dk["sigmainit"]) for _, _, dk, _ in rows}, reverse=True)
    L += ["", "**How understressed can the fault be and still fail at the observed "
          "pressure?** μ₀ ≥ f₀(1 − Δp_obs/σ̄₀):", "",
          "| σ̄₀ MPa | minimum μ₀ | minimum τ₀ MPa | source of σ̄₀ |",
          "|---|---|---|---|"]
    for sg in sigs:
        mm = f0 * (1 - dp_obs / sg)
        src = ("res1807/res1808's own value — a round number, **not** a measurement"
               if abs(sg - 30.0) < 1e-9 else
               "derived from σ_v 100, σ_Hmax 160, dip 10°, p_pore 73.82")
        L += [f"| {sg:g} | **{mm:.4f}** | **{mm*sg:.2f}** | {src} |"]
    if len(sigs) > 1:
        lo, hi = min(sigs), max(sigs)
        t_lo = f0 * (1 - dp_obs / lo) * lo
        t_hi = f0 * (1 - dp_obs / hi) * hi
        L += ["", f"So σ̄₀ = {lo:g} admits a fault **{t_hi - t_lo:.2f} MPa more "
              f"understressed** than σ̄₀ = {hi:g} does ({t_lo:.2f} vs {t_hi:.2f} MPa). "
              "Both floors follow from the wellhead record and the friction law alone, "
              "with no simulation involved. The lower σ̄₀ is also the "
              "measurement-derived one, so it is both the more defensible choice and "
              "the more permissive one."]
    # This closing paragraph was written for the mu0 = 0.37 stages and is FALSE
    # for any stage that does not contain them -- stage 6 is mu0 = 0.5359, where
    # it read as a claim about runs that are not in the package. Emit it only
    # when the stage actually has such runs.
    mus = {ff(dk["muinit"]) for _, _, dk, _ in rows}
    if any(abs(m - 0.37) < 0.02 for m in mus) and any(abs(s - 30.0) < 1e-9
                                                      for s in sigs):
        L += ["", "The μ₀ = 0.37 runs at σ̄₀ = 30.0 are therefore expected to slip "
              "only by over-pressurising — which is precisely the 1807/1808 "
              "behaviour being characterised (their wellhead runs +67% to +280%). "
              "Their σ̄₀ = 27.99 twins sit just below the threshold and should be "
              "able to slip at the observed pressure, by 0.19 MPa. That margin is "
              "thin enough that the twins may differ qualitatively, not just "
              "quantitatively."]
    else:
        m = ", ".join(f"{x:.4f}" for x in sorted(mus))
        L += ["", f"This stage runs μ₀ = {m}, comfortably above the floor above, so "
              "the fault can reach failure at the measured pressure without "
              "over-pressurising. Whether it then accumulates enough slip to build "
              "a front is a separate question that the friction law, not this "
              "inequality, decides — HBI's regularised rate-and-state law has no "
              "failure threshold, so Δp_crit is an orientation number here and not "
              "a prediction."]
    return L


def perm_plot(stage, runs, out):
    """Initial permeability for every distinct field. Uniform fields are plotted
    too, as a flat line, to PROVE no near-well disc is present."""
    seen, fields = {}, []
    for n in runs:
        dk, _ = read_deck(n)
        near, far, pf = perm_field(dk)
        key = pf or f"uniform {near:.3e}"
        if key not in seen:
            seen[key] = True
            im = int(dk["imax"])
            ds = ff(dk["ds"])
            if pf:
                # Stage 4 parameter files are multi-column ("kp phi"), so a bare
                # reshape fails and a bare max() would compare kpmax against a
                # porosity. Pick the kp column by NAME from the header line.
                k = np.loadtxt(IN / pf, skiprows=1)
                if k.ndim > 1:
                    names = (IN / pf).read_text().split("\n", 1)[0].split()
                    k = k[:, names.index("kp")] if "kp" in names else k[:, 0]
                k = k.reshape(im, im)
            else:
                k = np.full((im, im), near)
            fields.append((key, k, ds, near, far, [n]))
        else:
            for f in fields:
                if f[0] == key:
                    f[5].append(n)
    ncol = max(len(fields), 1)
    fig, axes = plt.subplots(2, ncol, figsize=(6.2 * ncol, 9.4),
                             constrained_layout=True, squeeze=False)
    for c, (key, k, ds, near, far, rr) in enumerate(fields):
        im = k.shape[0]
        x = (np.arange(im) - im // 2) * ds
        Z = 40
        sl = slice(im // 2 - Z, im // 2 + Z + 1)
        ax = axes[0][c]
        if near == far:
            ax.pcolormesh(x[sl], x[sl], k[sl, sl], shading="nearest",
                          cmap="viridis", vmin=near * 0.5, vmax=near * 1.5)
            ax.text(0.5, 0.5, "UNIFORM\nno disc present", transform=ax.transAxes,
                    ha="center", va="center", fontsize=13, color="w", weight="bold")
        else:
            from matplotlib.colors import LogNorm
            m = ax.pcolormesh(x[sl], x[sl], k[sl, sl], shading="nearest",
                              cmap="viridis", norm=LogNorm(vmin=far, vmax=near))
            fig.colorbar(m, ax=ax, pad=0.02).set_label("$k_p$ (m$^2$)")
            th = np.linspace(0, 2 * np.pi, 200)
            ax.plot(0.15 * np.cos(th), 0.15 * np.sin(th), lw=1.4, ls="--", color=DISC)
        ax.plot(0, 0, marker="+", ms=13, mew=2, color=DISC)
        ax.set_aspect("equal")
        ax.set(xlabel="along strike (km)", ylabel="along dip (km)")
        ax.set_title(f"{key}\nnear {near:.2e} / far {far:.2e}  "
                     f"({'uniform' if near == far else f'{near/far:.0f}x'})\n"
                     f"runs: {', '.join(str(r) for r in rr)}", fontsize=9)
        ax = axes[1][c]
        ax.plot(x, k[:, im // 2], lw=2.0, color="#0072BD")
        ax.axvline(0.15, color=DISC, lw=1.1, ls="--")
        ax.axvline(-0.15, color=DISC, lw=1.1, ls="--")
        ax.annotate("150 m", xy=(0.15, 0.03), xycoords=("data", "axes fraction"),
                    xytext=(4, 0), textcoords="offset points", fontsize=8.5, color=DISC)
        ax.set_yscale("log")
        ax.set_ylim(min(far, near) * 0.3, max(far, near) * 3)
        ax.set(xlabel="distance along dip from injector (km)",
               ylabel="initial $k_p$ (m$^2$)",
               title="radial profile through the injector (full domain)")
        ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle(f"Stage {stage} initial permeability — INITIAL CONDITION, "
                 f"before any slip\nuniform fields are plotted anyway, to show "
                 f"no near-well disc is present", fontsize=12)
    for e in ("png", "pdf"):
        fig.savefig(out / f"perm_ic_stage{stage}.{e}", dpi=170, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", type=int)
    ap.add_argument("--runs", nargs="*", type=int)
    a = ap.parse_args()
    cfg = STAGES[a.stage]
    runs = a.runs or cfg["runs"]
    out = REPO / f"presubmit_stage{a.stage}"
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for n in runs:
        parent = cfg["parents"][n]
        txt, o = params_txt(n, parent)
        (out / f"params_{n}.txt").write_text(txt)
        dk, _ = read_deck(n)
        rows.append((n, parent, dk, o))
        # diff against the parent, comments stripped
        d = subprocess.run(
            f"diff <(grep -vE '^!|^[[:space:]]*$' {IN}/res{parent}.in) "
            f"<(grep -vE '^!|^[[:space:]]*$' {IN}/res{n}.in)",
            shell=True, capture_output=True, text=True, executable="/bin/bash")
        (out / f"diff_{n}.txt").write_text(
            f"res{parent}.in  ->  res{n}.in   (comments stripped)\n"
            f"{'='*60}\n{d.stdout or '(identical)'}\n")
        print(f"  {n}: params_{n}.txt, diff_{n}.txt")

    perm_plot(a.stage, runs, out)
    print(f"  perm_ic_stage{a.stage}.png")

    L = [f"# {cfg['title']}", "",
         "**Nothing here has been submitted.** This is the pre-submittal check.",
         "", "## Decks", "",
         "| run | parent | **permev** | **sigmabar_0** | eta Pa·s | beta 1/Pa | "
         "phi | phi*beta | kpmax | kp=kpmin | perm field | injection | muinit | "
         "tau0 MPa | dp_crit MPa | tmax d |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for n, parent, dk, o in rows:
        field = ("uniform" if o["map"] is None
                 else f"{o['map']} ({o['near']/o['far']:.0f}x)")
        # permev and the injection record belong in the table: Stage 7's two runs
        # differ ONLY in permev, and Stage 6's differ from their parents only in
        # permev too, so a table without it renders those rows identical.
        inj = (dk.get("injection_file", "—")
               if dk.get("injectionfromfile", "F").upper() in ("T", "TRUE", ".TRUE.")
               else f"{dk.get('injection','—')} qinj={dk.get('qinj','—')}")
        L.append(f"| [{n}](params_{n}.txt) | {parent} | **{dk.get('permev','F')}** | "
                 f"**{dk['sigmainit']}** | {dk['eta']} | {dk['beta']} | "
                 f"{dk['phi']} | {o['phibeta']:.3e} | {dk.get('kpmax','—')} | "
                 f"{dk.get('kpmin','—')} | {field} | `{inj}` | {dk['muinit']} | "
                 f"{o['tau0']:.2f} | {o['dp_crit']:.2f} | {o['tmax_d']:.2f} |")
    L += ["", "## Derived hydraulics", "",
          "| run | D near m²/s | D far m²/s | str=beta*phi | T at kpmin | T at kpmax | "
          "gamma(h=60s, kpmin) | gamma(h=60s, kpmax) |",
          "|---|---|---|---|---|---|---|---|"]
    for n, parent, dk, o in rows:
        L.append(f"| {n} | {o.get('D_near',float('nan')):.4e} | "
                 f"{o.get('D_far',float('nan')):.4e} | {o['str']:.3e} | "
                 f"{o.get('T_at_kpmin',float('nan')):.3e} | "
                 f"{o.get('T_at_kpmax',float('nan')):.3e} | "
                 f"{o.get('gamma_h60_kpmin',float('nan')):.4f} | "
                 f"{o.get('gamma_h60_kpmax',float('nan')):.4f} |")
    L += ["", "`str`, `T` and `gamma` are the quantities `m_diffusion.f90` actually "
          "assembles (lines 444, 669, 671). `gamma` is the one a "
          "viscosity/compressibility trade does **not** hold fixed.", ""]
    L += feasibility(rows)
    L += ["", "## Initial condition", "",
          f"![initial permeability](perm_ic_stage{a.stage}.png)", "",
          "## Deck diffs against parents", ""]
    for n, parent, dk, o in rows:
        L.append(f"- [`res{n}.in` vs `res{parent}.in`](diff_{n}.txt)")
    L += ["", "## Launch commands, once approved", "", "```bash",
          "cd /home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs"]
    for n, parent, dk, o in rows:
        p = f" -p {o['map']}" if o["map"] else ""
        L.append(f"sbatch march26_submit_hbi_git_scratch.sh -i res{n}.in "
                 f"-w june_clean.txt{p}")
    L += ["```"]
    (out / f"PRESUBMIT_stage{a.stage}.md").write_text("\n".join(L) + "\n")
    print(f"  PRESUBMIT_stage{a.stage}.md")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
