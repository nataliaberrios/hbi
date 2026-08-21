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
| `sweeps/sweep_enhancement.png` | **SUPERSEDED — these runs are misconfigured.** See "The consistency requirement" below. permev T lowers the wellhead pressure and the fault never slips, but every permev T run here has `kpmax` 2e-11 against a map whose near-well value is 1.1e-12, so slip drives permeability 18× above anything the initial condition contains. |
| `sweeps/sweep_kpmax.png` | **SUPERSEDED, same reason.** kpmax swept 1.1e-12 → 1e-10 against a map topping out at 1.1e-12. |

## The consistency requirement

With a nonuniform initial permeability, the map's **near-well value must equal
`kpmax`** and its **background must equal `kpmin`**. Then the near-well disc starts
already at the enhanced ceiling and slip propagates that conduit outward from the
floor — the initial condition and the evolution law describe the same rock. Setting
`kpmax` above anything present in the map instead lets slip manufacture permeability
out of nowhere, which drains pressure away from the well and shuts the fault down.

That single property sorts the entire archive of 112 runs that have both enhancement
and a nonuniform map:

| | reached dc | median λ/λ_obs | within ±15% of the front |
|---|---|---|---|
| `kpmax` = k_near **and** `kpmin` = k_far — 35 runs | **29 / 35 (83%)** | **0.89** | **12** |
| bounds do not match the map — 77 runs | 2 / 77 (3%) | 0.42 | 0 |

**The earlier conclusion that "enhancement suppresses slip" was drawn from the
misconfigured 77 and does not hold.** Correctly configured, enhancement plus a
nonuniform map reaches λ/λ_obs 0.97–0.99 on a fault at τ₀ = 11.10 MPa — deeply
understressed. What it has not yet done is reach that front *and* the pressure: the
runs that match the front all have a far-field `kpmin` of 1e-15, some 400× below the
4e-13 the paper calibrates, and that is what lets wellhead pressure run to +280%.

`enhancement_census/` has the full census, the trade-off figure, and the per-run
numbers. Runs 632770–632777 are the untested combination: consistent bounds, a
**physical** 4e-13 far field, contrasts of 50× and 250× matching the front-matching
family, understressed at τ₀ 12.88 MPa, with storage swept to its physical ceiling.

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

<!-- INVENTORY-START: generated by docs/cooper_basin/make_readme_table.py -- do not edit by hand -->

## Run inventory

Every run with a folder here. Front fits use one shared window **0–8 d** with the observed front refit on it (λ_obs = 0.1412); λ is window-dependent, so a single shared window is the only fair comparison. Pressure error is measured while the well is **flowing**.

τ₀ = μ₀ × σ̄₀. **Understressed** means τ₀ < 15.0 MPa, Wang & Dunham's value — the whole point of the exercise. `k near / far` is read from the permeability **map file**, not the deck, because a deck reading `kp 4e-13` can be running a map with 1.1e-12 near the well.

| run | μ₀ | σ̄₀ MPa | τ₀ MPa | φβ Pa⁻¹ | permev | kpmax | kpmin | kL | dc | k near / far m² | λ/λ_obs | pressure | reached |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [632507](632507) | 0.50 | 30.00 | 15.00 | 1.0e-10 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | 0.13 | **+5.3%** | 17.20 d |
| [632532](632532) | 0.37 | 30.00 | **11.10** | 1.0e-10 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **+5.3%** | 17.20 d |
| [632534](632534) | 0.42 | 30.00 | **12.60** | 1.0e-10 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **+5.3%** | 17.20 d |
| [632536](632536) | 0.46 | 30.00 | **13.80** | 1.0e-10 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **+5.3%** | 17.20 d |
| [632548](632548) | 0.37 | 27.99 | **10.36** | 1.0e-10 | **T** | 5e-12 | 4e-13 | 1e-3 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **+5.3%** | 17.20 d |
| [632549](632549) | 0.37 | 27.99 | **10.36** | 1.0e-10 | **T** | 2e-11 | 4e-13 | 1e-3 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **+5.3%** | 17.20 d |
| [632550](632550) | 0.37 | 27.99 | **10.36** | 1.0e-10 | **T** | 1e-10 | 4e-13 | 1e-3 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **+5.2%** | 17.20 d |
| [632551](632551) | 0.37 | 27.99 | **10.36** | 1.0e-10 | **T** | 1.1e-12 | 4e-13 | 1e-3 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **+5.3%** | 17.20 d |
| [632568](632568) | 0.37 | 27.99 | **10.36** | 1.0e-10 | **T** | 2e-11 | 4e-13 | 1e-5 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **+3.9%** | 8.00 d |
| [632700](632700) | 0.46 | 27.99 | **12.88** | 3.0e-11 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | 0.18 | **+7.9%** | 8.00 d |
| [632701](632701) | 0.46 | 27.99 | **12.88** | 3.0e-11 | **T** | 2e-11 | 4e-13 | 1e-5 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **-6.8%** | 8.00 d |
| [632702](632702) | 0.46 | 27.99 | **12.88** | 2.0e-11 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | 0.35 | **+8.7%** | 8.00 d |
| [632703](632703) | 0.46 | 27.99 | **12.88** | 2.0e-11 | **T** | 2e-11 | 4e-13 | 1e-5 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **-6.9%** | 8.00 d |
| [632704](632704) | 0.46 | 27.99 | **12.88** | 1.5e-11 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | 0.48 | **+9.3%** | 8.00 d |
| [632705](632705) | 0.46 | 27.99 | **12.88** | 1.5e-11 | **T** | 2e-11 | 4e-13 | 1e-5 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **-6.9%** | 8.00 d |
| [632706](632706) | 0.50 | 27.99 | **13.99** | 3.0e-11 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | 0.74 | **+7.9%** | 8.00 d |
| [632707](632707) | 0.50 | 27.99 | **13.99** | 3.0e-11 | **T** | 2e-11 | 4e-13 | 1e-5 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **-9.5%** | 8.00 d |
| [632708](632708) | 0.50 | 27.99 | **13.99** | 2.0e-11 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | **0.93** | **+8.7%** | 8.00 d |
| [632709](632709) | 0.50 | 27.99 | **13.99** | 2.0e-11 | **T** | 2e-11 | 4e-13 | 1e-5 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **-9.5%** | 8.00 d |
| [632710](632710) | 0.50 | 27.99 | **13.99** | 1.5e-11 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | **1.11** | **+9.3%** | 8.00 d |
| [632711](632711) | 0.50 | 27.99 | **13.99** | 1.5e-11 | **T** | 2e-11 | 4e-13 | 1e-5 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **-9.6%** | 8.00 d |
| [632720](632720) | 0.46 | 27.99 | **12.88** | 1.0e-11 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | 0.68 | **+10.0%** | 8.00 d |
| [632721](632721) | 0.46 | 27.99 | **12.88** | 5.0e-12 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | **1.15** | **+11.2%** | 8.00 d |
| [632722](632722) | 0.50 | 27.99 | **13.99** | 1.0e-11 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | 1.37 | **+10.0%** | 8.00 d |
| [632723](632723) | 0.50 | 27.99 | **13.99** | 5.0e-12 | F | — | — | — | 1.0e-04 | 1.1e-12 / 4.0e-13 | 2.34 † | +16.2% | 5.95 d |
| [632750](632750) | 0.46 | 27.99 | **12.88** | 5.0e-12 | **T** | 2e-11 | 4e-13 | 1e-5 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **-7.2%** | 8.00 d |
| [632751](632751) | 0.46 | 27.99 | **12.88** | 5.0e-12 | **T** | 2e-11 | 4e-13 | 1e-4 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **-3.9%** | 8.00 d |
| [632752](632752) | 0.46 | 27.99 | **12.88** | 5.0e-12 | **T** | 5e-12 | 4e-13 | 1e-4 | 1.0e-04 | 1.1e-12 / 4.0e-13 | 0.42 | **-0.8%** | 8.00 d |
| [632753](632753) | 0.46 | 27.99 | **12.88** | 5.0e-12 | **T** | 2e-11 | 4e-13 | 1e-5 | 1.0e-04 | 1.1e-12 / 4.0e-13 | no slip | **-7.2%** | 30.00 d |

**Bold** τ₀ = understressed. **Bold** λ/λ_obs or pressure = inside the ±15% tolerance. `no slip` means peak slip never reached dc anywhere, so the run has no front to fit — a result, not missing data. `n/a` in the pressure column means the run is driven by an injection file other than `june_clean.txt`, so it cannot be compared against the June 2012 measurement.

† fit on less than the shared window (the run has not reached 8 d yet), so this ratio is not comparable to the others.

**Passing both targets: 3 runs**, best first:

- **632708** — front 0.93x, pressure +8.7%, tau_0 13.99 MPa, phi*beta 2.0e-11, enhancement off
- **632710** — front 1.11x, pressure +9.3%, tau_0 13.99 MPa, phi*beta 1.5e-11, enhancement off
- **632721** — front 1.15x, pressure +11.2%, tau_0 12.88 MPa, phi*beta 5.0e-12, enhancement off

Full decks are in `/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs/res<run>.in`, and each folder's `params.txt` has the deck verbatim plus the derived quantities (diffusivity, Δp needed to fail, tmax in days, how the run terminated).

<!-- INVENTORY-END -->
