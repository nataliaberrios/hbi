# Cooper Basin joint calibration — figures

Goal: match the observed slip front **and** the wellhead pressure history with an
**understressed** fault (τ₀ < 15.0 MPa, Wang & Dunham's value), using a nonuniform
initial permeability and permeability enhancement.

Every figure is one parameter sweep: one parameter varies, everything else is held
fixed and stated on the figure. Three panels each — wellhead pressure, front vs
time (R–T), front vs injected volume (R–V). The `.txt` beside each figure has the
numbers so the plot can be checked.

## Read these in order

| figure | what it shows |
|---|---|
| `sweeps/sweep_muinit.png` | **Pressure cannot constrain friction.** All four μ₀ give an identical pressure curve (bias −0.3 MPa, RMS 2.0 MPa, indistinguishable) while λ/λ_obs runs 0 → 0.46. With fixed permeability the pressure problem is decoupled from the mechanics, so only the front discriminates. |
| `sweeps/sweep_storage.png` | **Storage is the lever that works.** At understressed μ₀=0.46, φβ 3e-11 → 5e-12 moves λ/λ_obs 0.18 → **1.15**. 632721 reaches the front target with τ₀ = 12.9 MPa. |
| `sweeps/sweep_storage_mu050.png` | Same sweep at μ₀=0.50 (τ₀ 14.0 MPa) — λ/λ_obs 0.74 → 1.37, so the target is crossed at a *larger* φβ. Strength and storage trade off against each other. |
| `sweeps/sweep_enhancement.png` | **Enhancement is a negative feedback, not a positive one.** In three matched pairs, permev T *lowers* the wellhead pressure (bias +3.1 → −2.7 MPa) and the fault then never slips at all, where the identical permev F run does. |
| `sweeps/sweep_kpmax.png` | The same effect as a dose-response: kpmax 1.1e-12 → 1e-10 drives the pressure bias −0.27 → −0.94 MPa. No slip at any ceiling. |

## Per-run folders

`<jobnumber>/` holds that run's own three figures plus `params.txt` — the full deck
verbatim, plus what the deck does not state: the permeability map's near/far values
and disc radius, hydraulic diffusivity, τ₀, the overpressure needed to fail
(σ̄₀(1 − μ₀/f₀)), tmax in days, and how the run actually terminated.

Input decks for every run are in
`/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs/res<jobnumber>.in`.

## Two conventions that matter

**λ is window-dependent.** The same run gives 0.1847 over 30.7 d and 0.1730 over
16.4 d. Every sweep fits all its runs on one shared window and refits the observed
front on that same window; the window is printed on the figure. Runs too short to
reach it are named as excluded rather than silently dropped.

**Pressure is scored while flowing, in MPa.** HBI has no wellbore bleed-off, so
during a shut-in the simulated wellhead stays near its flowing value while the
measurement falls ~20 MPa (sim 34.6 vs measured 16.1). That gap is identical for
every parameter set and so cannot discriminate between them. `mean(sim/measured)`
over the full span is not reported: the measured wellhead passes through ~0 in
shut-ins, and dividing by it made the same run read 2.17 on one sampling and 0.32
on another.

## Caveat that no parameter choice fixes

These runs have a > b everywhere, so they are velocity-strengthening and produce no
seismicity. The comparison is an *aseismic slip front* against an observed
*seismicity front*.
