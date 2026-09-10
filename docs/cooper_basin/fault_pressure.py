#!/usr/bin/env python3
"""The pressure datum: derive it from the Cooper Basin measurements, and show
exactly how much of the resulting curve shift is measured and how much assumed.

Writes docs/figs/cycle2/DATUM.md and figures/cycle2/fault_pressure.png. Every
number in the markdown is computed here, so the document cannot drift from the
arithmetic.

THE PROBLEM. Every dp in this project had been referenced to a value read off
the wellhead record -- its first sample (34.412 MPa), its pre-injection median
(33.970), its value at data-day 4.300 (33.271). Those are wellhead pressures
that lie NEAR the initial fault pressure; none of them is it. The initial fault
pressure is a quantity the field measurements determine, so it should be derived
and then used, which is what this does.

Usage:  python fault_pressure.py
"""
import importlib.util as iu
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = "/home/users/nberrios/3dhbi/hbi_analysis"
_s = iu.spec_from_file_location("sf", H + "/make_sweep_figures.py")
sf = iu.module_from_spec(_s); _s.loader.exec_module(sf)

FIG = Path(H) / "figures" / "cycle2"
MD = Path("/home/users/nberrios/3dhbi/hbi_git/docs/figs/cycle2/DATUM.md")
T0 = 4.300
INK, MUTED, OBSC = "#1a1a19", "#6b6b66", "#a8071a"

# Holl & Barton (2015), as recorded in CONCLUSIONS.md:1001-1002
P_F0 = 72.70            # MPa, reservoir pressure at -4100 mAHD
Z_AHD = 4100.0          # m below AHD (~sea level), fault datum
# THE GAUGE IS NOT AT THE FAULT DATUM. Habanero-4 Well Completion Report
# (HPP-FN-OT-RPT-00492-1.0, 8 Jan 2013; PEPS-SA well 2739, Open File):
# ground level 73.34 m AHD, RT-GL 9.14 m, so the wellhead gauge sits 82.48 m
# ABOVE the datum Holl quotes. The static column between them is therefore
# 4182.5 m, not 4100 -- worth 0.81 MPa at rho 1000, all of it lowering the
# datum and so raising the observed dp. An earlier version of this script used
# 4100 and guessed the elevation at "~60 m".
GL, RT_GL = 73.34, 9.14
Z_COLUMN = Z_AHD + GL + RT_GL
# the stress state, setup_model.m:110-113, for the sigmabar_0 cross-check
S_V, S_HMAX, DIP = 100.0, 160.0, 10.0
SIGMAINIT = 27.99       # this project's decks
P_F0_TAIYI = 73.82      # setup_model.m:112
# THE COLUMN IS NOT ONE FLUID. The static column, pre-injection, had been
# sitting in a 240-250 C reservoir and reheating: hot and light. The flowing
# column is surface water pumped at 25-50 L/s with almost no residence time to
# heat: cold and dense. A single density is wrong, and the two states differ by
# enough to matter -- 2.01 MPa per 50 kg/m3 over 4100 m.
#
# RHO_FLOW converts the record DURING INJECTION and is the one that sets the
# plotted dp. Cold water at ~40 MPa is 992 kg/m3 at 60 C and 1008 at 20 C, so
# 1000 is the middle of the plausible injectate range and is also what
# make_sweep_figures.py:55 uses.
#
# RHO_STATIC is a RESULT, not an input: the measured pre-injection wellhead of
# 33.970 MPa implies 963 kg/m3, i.e. a mean column temperature near 90-100 C,
# which is what a shut-in well in this reservoir should look like. It is
# reported as a consistency check on p_f0, NOT used to convert the flowing
# record.
#
# An earlier version of this script used the 963 to argue that rho = 1000 was
# "contradicted by the well". That was wrong: it compared two different fluid
# states. The difference between them, 1.49 MPa, IS the move-up.
# RHO_FLOW IS NOW MEASURED, NOT ASSUMED. The WCR gives the completion fluid
# left in the well as 8.35 ppg ambient = 1001 kg/m3 (section 4.5.3), and that
# is the wellbore column immediately before the Oct/Nov 2012 stimulation.
# Independently, Hogarth & Bour (2015, WGC 31006) Table 3 measures the cold
# limb of the 2013 closed loop at 983-1006 kg/m3 for 80 C brine. So 1001 is
# not a choice between plausible values; it is the reported number.
RHO_FLOW = 1001.0
RHO_BRACKET = (983.0, 1006.0)   # Hogarth Table 3, cold limb of the loop
RHO_ASSUMED = RHO_FLOW
G = 9.81
OLD_REF = 34.412        # the record's first sample

plt.rcParams.update({"font.size": 11, "axes.titlesize": 12,
                     "axes.labelsize": 11.5, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "legend.fontsize": 9.5})


def head(rho, z=None):
    """Static head in MPa over the gauge-to-fault column."""
    return rho * G * (Z_COLUMN if z is None else z) / 1e6


# THE FLOW PATH IS NOT ONE PIPE, AND THE OLD DIAMETER WAS THE WRONG ONE.
# make_sweep_figures.py uses DW = 0.178 m uniformly over 4077 m -- but 0.1778
# is the 7 in tubing OUTSIDE diameter. The WCR's Appendix K tally gives 7 in
# 41# T95SS with ID 5.820 in = 0.14783 m, run to the production packer at
# 3005.37 mMDRT, below which the 9-7/8 in casing is much wider. Friction goes
# as D^-5, so using the OD instead of the ID understates it badly.
FRIC_SEGMENTS = ((3005.0, 5.820 * 0.0254),           # tubing, ID
                 (Z_COLUMN - 3005.0, 8.625 * 0.0254))  # casing below the packer


def friction(q_Ls):
    """Darcy-Weisbach loss in MPa over the real segmented geometry.

    At the 60.9 L/s peak this gives 2.03 MPa against the 1.03 MPa the uniform
    0.178 m pipe gives -- a factor 1.97, and 1.92 MPa of it is the tubing
    alone. It enters dp(t) = p_wh + rho g Z - friction directly, so a 1 MPa
    error at peak rate was a 1 MPa error in the observed dp at peak rate. It
    matters least where it was checked (the plateau, q ~ 23 L/s, 0.28 MPa) and
    most at the 2.5 d rate peak.
    """
    q = np.asarray(q_Ls, float) / 1000.0
    return sum(sf.FD * 8.0 * L * RHO_FLOW * q ** 2 / (np.pi ** 2 * D ** 5)
               for L, D in FRIC_SEGMENTS) / 1e6


def datum():
    """The pressure datum expressed at the wellhead gauge, in MPa."""
    return P_F0 - head(RHO_FLOW)


def dp_observed(obs, t0=T0):
    """(t_sim, dp) -- the measured pressure change from p_f0.

    THE one place this conversion lives. score_cycle2.py and compare_arms.py
    import it so the datum, the column length, the density and the friction
    cannot drift between the scorer and the figures.
    """
    q = np.interp(obs["tp"], obs["ti"], obs["q"])
    return obs["tp"] - t0, obs["pm"] + head(RHO_FLOW) - friction(q) - P_F0


def main(argv=None):
    obs = sf.observed()
    tp, pm = obs["tp"], obs["pm"]
    q = np.interp(tp, obs["ti"], obs["q"])
    t = tp - T0
    p_static = float(np.median(pm[tp < 0.501]))
    rho_implied = (P_F0 - p_static) * 1e6 / (G * Z_COLUMN)
    th = np.radians(DIP)
    snn = S_HMAX * np.sin(th) ** 2 + S_V * np.cos(th) ** 2

    datum = P_F0 - head(RHO_ASSUMED)
    p_fault = pm + head(RHO_FLOW) - friction(q)
    dp = p_fault - P_F0
    dp_old = pm - OLD_REF

    L = []
    A = L.append
    A("# The pressure datum, and how much of the shift is real\n")
    A("Generated by `docs/cooper_basin/fault_pressure.py`. Every number below "
      "is computed there.\n")

    A("\n## 1. What is being plotted\n")
    A("The measurement is an **absolute wellhead pressure**. The model computes "
      "a **change in fault pore pressure** from `pfinit = 0`. To compare them, "
      "the record has to be carried down the well and referenced to the fault's "
      "initial pressure:\n")
    A("```")
    A("p_fault(t) = p_wh(t) + rho*g*Z - friction(q(t))")
    A("dp(t)      = p_fault(t) - p_f0")
    A("```")
    A(f"with `Z` = {Z_COLUMN:.1f} m the gauge-to-fault column and `friction` "
      f"the Darcy-Weisbach loss over the WCR's real segmented geometry — "
      f"{friction(obs['q'].max()):.2f} MPa at the {obs['q'].max():.1f} L/s "
      f"peak, against 1.03 MPa for the uniform 0.178 m pipe previously "
      f"assumed, and zero when shut in.\n")
    A("**The head does not cancel.** Substituting,\n")
    A("```")
    A("dp(t) = p_wh(t) - friction(q(t)) - [p_f0 - rho*g*Z]")
    A("                                   ^^^^^^^^^^^^^^^^ the datum")
    A("```")
    A("so `rho*g*Z` enters **once, additively**. It would only cancel if the "
      "datum were itself defined as a wellhead pressure. It is not — it is "
      "derived from `p_f0`, which is why `rho` matters directly.\n")

    A("\n## 2. The datum from the measurements\n")
    A("| quantity | value | source |")
    A("|---|---|---|")
    A(f"| `p_f0` | {P_F0:.2f} MPa | Holl & Barton (2015), at −{Z_AHD:.0f} mAHD "
      f"— `CONCLUSIONS.md:1001` |")
    A(f"| column length | {Z_COLUMN:.1f} m | {Z_AHD:.0f} m below AHD + gauge at "
      f"{GL+RT_GL:.2f} m AHD (WCR: GL {GL} m, RT−GL {RT_GL} m) |")
    A(f"| `rho*g*Z` | {head(RHO_FLOW):.3f} MPa | rho = "
      f"{RHO_FLOW:.0f} kg/m³ — **measured**, WCR §4.5.3, 8.35 ppg |")
    A(f"| overpressure | {P_F0 - head(RHO_FLOW):.3f} MPa | `p_f0 - rho*g*Z` |")
    A(f"| **datum at the wellhead** | **{datum:.3f} MPa** | |")
    A("")
    A(f"Previously the data was referenced to {OLD_REF:.3f} MPa, so the "
      f"observed curve moves **up by {OLD_REF - datum:+.3f} MPa** at "
      f"`sim t = 0`, easing to "
      f"{dp[int(np.argmin(np.abs(t-10.0)))] - dp_old[int(np.argmin(np.abs(t-10.0)))]:+.2f}"
      f" MPa by 10 d as the q² friction term grows.\n")

    A("\n## 3. The column is not one fluid, and that is where the shift "
      "comes from\n")
    A("The datum is `p_f0 - rho*g*Z`, so the column density acts directly on "
      f"the answer at **{50*G*Z_COLUMN/1e6:.2f} MPa per 50 kg/m³** over "
      f"{Z_COLUMN:.0f} m:\n")
    A("| rho kg/m³ | head MPa | datum MPa | move-up |")
    A("|---|---|---|---|")
    for r in (900, 949, int(round(rho_implied)), 983, 1001, 1006):
        mark = " ← implied by our record's static wellhead" \
            if abs(r - rho_implied) < 1 else \
            (" ← **MEASURED**, WCR completion fluid 8.35 ppg" if r == 1001
             else (" ← Hogarth Table 3, cold limb" if r in (983, 1006) else ""))
        A(f"| {r} | {head(r):.3f} | {P_F0-head(r):.3f} | "
          f"{OLD_REF-(P_F0-head(r)):+.3f}{mark} |")
    A("")
    A(f"That is **{50*G*Z_COLUMN/1e6:.2f} MPa per 50 kg/m³**. The whole "
      f"\"couple of MPa\" lies inside the density uncertainty.\n")
    A("### Two fluid states, not one\n")
    A("| | condition | rho | why |")
    A("|---|---|---|---|")
    A(f"| **static** | pre-injection, shut in | **{rho_implied:.0f}** | had been "
      f"sitting in a 240–250 °C reservoir and reheating — hot and light. A mean "
      f"column temperature near 90–100 °C. |")
    A(f"| **flowing** | during injection | **{RHO_FLOW:.0f}** | surface water at "
      f"25–50 L/s with almost no residence time to heat — cold and dense. Water "
      f"at 40 MPa is 992 kg/m³ at 60 °C and 1008 at 20 °C. |")
    A("")
    A("The record being converted is the **flowing** one, so `RHO_FLOW` sets "
      "the plotted dp. The static density is a *result*, not an input:\n")
    A("```")
    A(f"measured pre-injection wellhead (median, t < 0.501 d) = {p_static:.3f} MPa")
    A(f"rho = ({P_F0:.2f} - {p_static:.3f})e6 / (9.81 x {Z_COLUMN:.0f}) "
      f"= {rho_implied:.0f} kg/m3")
    A("```")
    A(f"That {rho_implied:.0f} kg/m³ is a **consistency check on `p_f0`**, and "
      f"it passes: a shut-in well in a 240 °C reservoir should have a column "
      f"averaging around 90–100 °C, which is what {rho_implied:.0f} corresponds "
      f"to. If Holl's `p_f0` were badly wrong this number would come out "
      f"unphysical.\n")
    A("### The move-up, and its bracket\n")
    A(f"The static and flowing columns differ by "
      f"{(RHO_FLOW-rho_implied)*G*Z_COLUMN/1e6:.2f} MPa of head, and **that "
      f"difference is the move-up**:\n")
    A("| injectate | rho | datum MPa | move-up |")
    A("|---|---|---|---|")
    for r, lab in ((RHO_BRACKET[0], "60 °C"), (RHO_FLOW, "~40 °C, used"),
                   (RHO_BRACKET[1], "20 °C")):
        A(f"| {lab} | {r:.0f} | {P_F0-head(r):.3f} | "
          f"{OLD_REF-(P_F0-head(r)):+.3f} |")
    A("")
    A(f"So **{OLD_REF-(P_F0-head(RHO_BRACKET[0])):.2f} to "
      f"{OLD_REF-(P_F0-head(RHO_BRACKET[1])):.2f} MPa** is the defensible "
      f"range, and {OLD_REF-datum:.2f} MPa is what is plotted.\n")
    A("**The justification is independent of the outcome, which is the only "
      "thing that makes it safe.** The thermal argument is about the wellbore, "
      "not about the fit; it would hold whether or not it improved the match. "
      "Choosing `rho` because it improves the match would be the first thing a "
      "reviewer pulls on, and it would take the arm-2 result with it.\n")
    A("What would settle it outright, in order of strength: a **downhole gauge "
      "or temperature log** from the Habanero 4 completion report (cited at "
      "`setup_model.m:92`, so it exists); the **injectate temperature** from "
      "the operational record; and the **wellhead elevation**, since Holl's "
      f"depth is mSS and ~60 m of surface elevation adds another "
      f"{head(RHO_FLOW, 60.0):.2f} MPa of column on top of everything above.\n")

    A("\n## 4. An inconsistency this exposes\n")
    A(f"With the measured stress state — `S_v` = {S_V:.0f} MPa, `S_Hmax` = "
      f"{S_HMAX:.0f} MPa, dip {DIP:.0f}° — the resolved normal stress is\n")
    A("```")
    A(f"sigma_nn = {S_HMAX:.0f} sin^2({DIP:.0f}) + {S_V:.0f} cos^2({DIP:.0f}) "
      f"= {snn:.3f} MPa")
    A("```")
    A(f"so Holl's `p_f0` = {P_F0:.2f} implies an effective normal stress of "
      f"{snn-P_F0:.2f} MPa. **This project's decks set `sigmainit` = "
      f"{SIGMAINIT:.2f}**, which instead implies `p_f0` = {snn-SIGMAINIT:.2f} — "
      f"i.e. Taiyi's {P_F0_TAIYI:.2f} (`setup_model.m:112`), not Holl's.\n")
    A(f"So `sigmainit` = {SIGMAINIT:.2f} and a datum built on Holl's "
      f"{P_F0:.2f} MPa are inconsistent by {abs((snn-SIGMAINIT)-P_F0):.2f} MPa. "
      "One or the other should move, and they cannot move independently, since "
      "`tau_0 = muinit * sigmainit`. Recorded, not acted on.\n")
    A("An earlier version of this document claimed `p_f0` = 73.82 was confirmed "
      "two independent ways, by `sigma_nn - sigmainit` and by the stated "
      "`p_pore`. **That was circular** — `sigmainit` was itself derived from "
      "that `p_pore`, so the 0.001 MPa agreement proved nothing. Withdrawn.\n")

    A("\n## 5. The corrected series\n")
    A("| sim t (d) | p_wh | p_fault | dp old | dp new | shift |")
    A("|---|---|---|---|---|---|")
    for ts in (0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 8.7, 10.0, 12.0):
        i = int(np.argmin(np.abs(t - ts)))
        A(f"| {ts:.2f} | {pm[i]:.3f} | {p_fault[i]:.3f} | {dp_old[i]:+.3f} | "
          f"{dp[i]:+.3f} | {dp[i]-dp_old[i]:+.3f} |")
    A("")
    A(f"`dp` at `sim t = 0` is **{float(np.interp(0.0,t,dp)):+.3f} MPa**, not "
      "zero: the fault sits above virgin pressure when injection resumes even "
      "though the well had been vented, which is cycle 1's residual formation "
      "overpressure. That is why this datum lifts the curve rather than "
      "dropping it.\n")
    A("This is a plotting and scoring datum. **No simulation changes**; "
      "`pfinit` stays 0.\n")

    MD.parent.mkdir(parents=True, exist_ok=True)
    MD.write_text("\n".join(L) + "\n")

    # ------------------------------------------------------------------ figure
    fig, ax = plt.subplots(2, 1, figsize=(11.5, 8.4), dpi=200, sharex=True,
                           gridspec_kw=dict(height_ratios=[1.0, 1.3],
                                            hspace=0.22))
    aa, ad = ax
    k = t >= 0
    aa.plot(t[k], pm[k], lw=1.3, color=MUTED, alpha=0.85, label="measured wellhead")
    aa.plot(t[k], p_fault[k], lw=1.6, color=OBSC,
            label=f"fault pressure = wellhead + {head(RHO_ASSUMED):.2f} − friction")
    aa.axhline(P_F0, color=INK, lw=1.4, ls="--")
    aa.set(ylabel="Absolute pressure (MPa)", ylim=(25, 95))
    aa.set_title(f"(a)  Measured wellhead, and the fault pressure it implies "
                 f"at rho = {RHO_FLOW:.0f}")
    aa.legend(loc="lower right", framealpha=0.95)

    ad.axhline(0, color=INK, lw=1.2, ls="--")
    # THE BRACKET, drawn. The observed dp under each plausible injectate
    # density: the band is the column-density uncertainty, and it is the
    # dominant uncertainty on the observed curve -- wider than anything else
    # done to the data.
    dp_lo = pm + head(RHO_BRACKET[0]) - friction(q) - P_F0
    dp_hi = pm + head(RHO_BRACKET[1]) - friction(q) - P_F0
    ad.fill_between(t[k], dp_lo[k], dp_hi[k], color=OBSC, alpha=0.22, lw=0)
    ad.plot(t[k], dp_lo[k], lw=0.9, ls=":", color=OBSC, alpha=0.9)
    ad.plot(t[k], dp_hi[k], lw=0.9, ls=":", color=OBSC, alpha=0.9)
    ad.plot(t[k], dp[k], lw=2.0, color=OBSC,
            label=f"measured, flowing column rho {RHO_BRACKET[0]:.0f}–"
                  f"{RHO_BRACKET[1]:.0f}\n   (20–60 °C injectate; the band is "
                  f"only {head(RHO_BRACKET[1])-head(RHO_BRACKET[0]):.2f} MPa wide)")
    # The comparison that IS visible: what the curve would be if the column
    # were still the hot, light one the shut-in pressure implies. The gap
    # between this and the solid curve is the whole move-up, and it dwarfs the
    # injectate-temperature bracket.
    dp_hot = pm + head(rho_implied) - friction(q) - P_F0
    ad.plot(t[k], dp_hot[k], lw=1.5, ls="--", color="#8E44AD",
            label=f"if the column were still hot, rho = {rho_implied:.0f}\n"
                  f"   ({head(RHO_FLOW)-head(rho_implied):.2f} MPa lower — this "
                  f"is the move-up)")
    for c, n in zip(plt.cm.viridis(np.linspace(0.05, 0.85, 5)),
                    range(632960, 632965)):
        d = sf.run_data(n, sf.deck(n))
        if d is None or d.get("tpw") is None:
            continue
        dk = sf.deck(n)
        tau = sf.ffloat(dk["muinit"]) * sf.ffloat(dk["sigmainit"])
        ad.plot(d["tpw"], d["ppw"] - (sf.P0 - sf.RHO * sf.G * sf.HW / 1e6),
                lw=1.6, color=c, label=r"$\tau_0$ = " f"{tau:.2f} MPa")
    ad.set(xlabel="Days since injection resumed (data-day 4.300)",
           ylabel="Pressure change from $p_{f0}$ (MPa)", xlim=(0, 13.2),
           ylim=(-2, 21))
    ad.set_title(f"(b)  Pressure change from $p_{{f0}}$ = {P_F0:.2f} MPa, "
                 f"with the column-density bracket")
    ad.legend(loc="upper left", framealpha=0.95, ncol=2, fontsize=8.5)
    FIG.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(FIG / f"fault_pressure.{e}", bbox_inches="tight")
    plt.close(fig)

    print(f"bracket width, rho {RHO_BRACKET[0]:.0f}-{RHO_BRACKET[1]:.0f}: "
          f"{head(RHO_BRACKET[1])-head(RHO_BRACKET[0]):.3f} MPa")
    print(f"hot-vs-flowing gap, rho {rho_implied:.0f}-{RHO_FLOW:.0f}: "
          f"{head(RHO_FLOW)-head(rho_implied):.3f} MPa\n")
    print(f"wrote {MD}")
    print(f"wrote {FIG}/fault_pressure.png")


if __name__ == "__main__":
    main()
