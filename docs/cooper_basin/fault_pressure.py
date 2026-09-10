#!/usr/bin/env python3
"""Solve for the initial fault pressure from the Cooper Basin measurements, then
plot the measured pressure change relative to it.

THE POINT. Every dp in this project has been referenced to some value picked
off the wellhead record -- its first sample, its pre-injection median, its value
at data-day 4.300. None of those is the initial fault pressure; they are
wellhead pressures that happen to be near it. The fault pressure is a quantity
the field measurements DETERMINE, and once solved for it fixes the reference
without any choice being made.

STEP 1, SOLVE FOR p_f0. Two independent routes, and they must agree or the
stress state and the pore pressure in this project are inconsistent:

  A. from the measured stress state and this project's own sigmainit.
     Holl & Barton (2015): s_v = 100 MPa, s_Hmax = 160 MPa, fault dip 10 deg
     (encoded at setup_model.m:110-113). Resolving onto the fault,

         sigma_nn = s_Hmax sin^2(theta) + s_v cos^2(theta) = 101.809 MPa

     and our decks set sigmainit = 27.99 MPa, which is the EFFECTIVE normal
     stress. So

         p_f0 = sigma_nn - sigmainit = 101.809 - 27.99 = 73.819 MPa

  B. stated directly. setup_model.m:112, p_pore = 73.82 MPa, commented
     "hydrostatic pressure + over pressure".

They agree to 0.001 MPa. That is not a coincidence and it is worth stating
plainly: sigmainit = 27.99 IS the measured stress state minus the measured pore
pressure. It was never an independent parameter, which also means p_f0 cannot
be adjusted without moving sigmainit and therefore tau_0.

At the fault's 4100 m median depth (setup_model.m:108), hydrostatic is
40.221 MPa, so the reservoir is OVERPRESSURED BY 33.599 MPa. That overpressure
is why a shut-in well reads ~34 MPa at surface, and it is what Holl reports.

STEP 2, CONVERT THE WELLHEAD RECORD TO FAULT PRESSURE.

    p_fault(t) = p_wh(t) + rho g z_fault - friction(q(t))

The static head is the dominant term, 40.221 MPa. The friction loss is the
Darcy-Weisbach pipe term already used in make_sweep_figures.run_data, worth up
to 0.64 MPa at the peak 48 L/s and zero when shut in -- so it matters for the
shape at high rate, not for the reference.

STEP 3, PLOT dp = p_fault(t) - p_f0. This is directly what HBI computes: its pf
starts at pfinit = 0 and sigmainit already contains p_f0, so the model's dp and
this dp are the same quantity measured from the same zero. No baseline mismatch
term survives.

WHAT COMES OUT. At sim t = 0 the measured fault pressure is 0.345 MPa BELOW
p_f0 -- the well was vented during the shut-in and had not fully recovered. That
0.345 MPa is the physically meaningful residual depletion, and it is much
smaller than the 1.141 MPa I had been applying by referencing to the record's
first sample instead.

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

OUT = Path(H) / "figures" / "cycle2"
T0 = 4.300
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
OBSC = "#a8071a"

# Cooper Basin measurements, Holl & Barton (2015) as encoded in
# taiyi-wang-seis3D/source_code/setup_model.m:108-113
S_V, S_HMAX, DIP_DEG = 100.0, 160.0, 10.0
Z_FAULT, Z_WELL = 4100.0, 4077.0        # fault median depth; well depth
SIGMAINIT = 27.99                       # our decks' effective normal stress
P_PORE_STATED = 73.82                   # setup_model.m:112

plt.rcParams.update({"font.size": 11, "axes.titlesize": 12,
                     "axes.labelsize": 11.5, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "legend.fontsize": 9.5})


def solve_p_f0():
    """(p_f0, sigma_nn) in MPa, with the two routes cross-checked."""
    th = np.radians(DIP_DEG)
    snn = S_HMAX * np.sin(th) ** 2 + S_V * np.cos(th) ** 2
    p_a = snn - SIGMAINIT
    assert abs(p_a - P_PORE_STATED) < 0.01, (
        f"the two routes to p_f0 disagree: sigma_nn - sigmainit = {p_a:.3f} "
        f"but setup_model.m states {P_PORE_STATED}. Either the stress state or "
        f"sigmainit is inconsistent and the reference cannot be trusted.")
    return p_a, snn


def friction(q_Ls):
    """Darcy-Weisbach pipe loss in MPa, same form as run_data uses."""
    q = np.asarray(q_Ls, float) / 1000.0        # m^3/s
    return (sf.FD * 8.0 * sf.HW * sf.RHO * q ** 2
            / (np.pi ** 2 * sf.DW ** 5) / 1e6)


def fault_pressure(obs, p_f0):
    """(t_sim, p_fault MPa, dp MPa) from the wellhead record."""
    tp, pm = obs["tp"], obs["pm"]
    q = np.interp(tp, obs["ti"], obs["q"])
    head = sf.RHO * sf.G * Z_FAULT / 1e6
    p_fault = pm + head - friction(q)
    return tp - T0, p_fault, p_fault - p_f0


def main(argv=None):
    p_f0, snn = solve_p_f0()
    obs = sf.observed()
    t, p_fault, dp = fault_pressure(obs, p_f0)
    head = sf.RHO * sf.G * Z_FAULT / 1e6

    print(f"sigma_nn = {S_HMAX} sin^2({DIP_DEG}) + {S_V} cos^2({DIP_DEG}) "
          f"= {snn:.3f} MPa")
    print(f"p_f0 = sigma_nn - sigmainit = {snn:.3f} - {SIGMAINIT} = {p_f0:.3f} MPa")
    print(f"     cross-check, setup_model.m:112 states {P_PORE_STATED} MPa\n")
    print(f"hydrostatic at {Z_FAULT:.0f} m = {head:.3f} MPa   ->  "
          f"OVERPRESSURE {p_f0 - head:.3f} MPa")
    print(f"wellhead-equivalent reference, p_f0 - rho g z = "
          f"{p_f0 - head:.3f} MPa\n")
    print(f"peak pipe friction at {obs['q'].max():.1f} L/s = "
          f"{friction(obs['q'].max()):.3f} MPa\n")

    k = t >= 0
    print(f"{'sim t':>7} {'p_wh':>8} {'p_fault':>9} {'dp':>8}")
    for ts in (0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 8.7, 10.0, 12.0):
        i = int(np.argmin(np.abs(t - ts)))
        print(f"{ts:>7.2f} {obs['pm'][i]:>8.3f} {p_fault[i]:>9.3f} {dp[i]:>+8.3f}")
    d0 = float(np.interp(0.0, t[k], dp[k]))
    print(f"\nat sim t = 0 the fault sits {d0:+.3f} MPa relative to p_f0 -- the")
    print(f"residual depletion from the vented shut-in. Referencing to the "
          f"record's\nfirst sample instead implied {-1.141:+.3f} MPa, i.e. "
          f"{abs(-1.141-d0):.3f} MPa too much.")

    fig, ax = plt.subplots(2, 1, figsize=(11.5, 8.4), dpi=200, sharex=True,
                           gridspec_kw=dict(height_ratios=[1.0, 1.3],
                                            hspace=0.22))
    aa, ad = ax

    aa.plot(t[k], obs["pm"][k], lw=1.3, color=MUTED, alpha=0.85,
            label="measured wellhead")
    aa.plot(t[k], p_fault[k], lw=1.6, color=OBSC,
            label=f"fault pressure = wellhead + {head:.2f} MPa − friction")
    aa.axhline(p_f0, color=INK, lw=1.4, ls="--")
    aa.set(ylabel="Absolute pressure (MPa)", ylim=(25, 95))
    aa.set_title("(a)  Measured wellhead, and the fault pressure it implies")
    aa.legend(loc="lower right", framealpha=0.95)

    ad.axhline(0, color=INK, lw=1.2, ls="--")
    ad.plot(t[k], dp[k], lw=1.6, color=OBSC, label="measured")
    cols = plt.cm.viridis(np.linspace(0.05, 0.85, 5))
    for c, n in zip(cols, range(632960, 632965)):
        d = sf.run_data(n, sf.deck(n))
        if d is None or d.get("tpw") is None:
            continue
        dk = sf.deck(n)
        tau = sf.ffloat(dk["muinit"]) * sf.ffloat(dk["sigmainit"])
        # the model's wellhead already carries P0 - rho g H_well; its dp is
        # what remains once that is removed
        dpm = d["ppw"] - (sf.P0 - sf.RHO * sf.G * sf.HW / 1e6)
        ad.plot(d["tpw"], dpm, lw=1.6, color=c,
                label=r"$\tau_0$ = " f"{tau:.2f} MPa")
    ad.set(xlabel="Days since injection resumed (data-day 4.300)",
           ylabel="Pressure change from $p_{f0}$ (MPa)", xlim=(0, 13.2),
           ylim=(-3, 25))
    ad.set_title(f"(b)  Pressure change relative to the initial fault "
                 f"pressure, $p_{{f0}}$ = {p_f0:.2f} MPa")
    ad.legend(loc="upper left", framealpha=0.95, ncol=2)

    OUT.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"fault_pressure.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}/fault_pressure.png")
    return p_f0


if __name__ == "__main__":
    main()
