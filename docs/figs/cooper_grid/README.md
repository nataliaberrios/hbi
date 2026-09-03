# Cooper Basin grid search — clean branch

Branch `cooper-grid-search`, off `limitsigma-bug-notes` (which is where the
`limitsigma` fix lives — it is **not** on master).

Everything in *this* directory has been checked. The older
`docs/figs/cooper_basin_calibration/` is left in place but is **not** to be
trusted: 77 of its 112 enhancement runs had `kpmax` set above anything present in
their permeability map, and four of its sweep figures are marked SUPERSEDED for
that reason.

## Start here

**[`BOUNDS_AUDIT.md`](BOUNDS_AUDIT.md)** — the 90 folders here that predate
the grid, each labelled VALID / MISCONFIGURED / NOT APPLICABLE. 48 are
misconfigured and must not be read as results; 26 are `permev F` so the rule
does not apply; several are provenance for current work. Every one of those
folders also carries its own `README.md` with its verdict.

**[`RUN_KEY.md`](RUN_KEY.md)** — one row per run for all 68 simulations:
what each block was testing, what varies, and the two scores. Generated
from the decks by `docs/cooper_basin/make_run_key.py`, so it cannot drift.

## Goal

Match the Cooper Basin wellhead pressure **and** the slip front with an
understressed fault, without discarding the stress measurements. Hypothesis:
permeability enhancement (which Wang & Dunham did not have) plus a nonuniform
initial permeability gets there at lower τ₀ than their 15.0 MPa.

## The rule that governs every deck

`kpmax` is the ceiling, `kp` = `kpmin` is the floor and the initial value. When
the initial permeability is nonuniform, the near-well disc **is** `kpmax` and the
background **is** `kpmin`. No other permeability value may appear anywhere. This
is not bookkeeping — across 112 archived runs, the 35 obeying it had 29 reach dc;
the 77 that did not had 2.

## Verified against Taiyi's source

`taiyi-wang-seis3D/source_code/setup_model.m`. Every value checks out:

| | source | |
|---|---|---|
| σ_v, σ_Hmax, p_pore, dip | 100e6 (:111), 160e6 (:113), 73.82e6 (:112), 10/180·π (:110) | ✓ |
| **σ̄₀ = 27.99** | `rotate_stress(-160e6,0,-100e6,10°)` → −101.8092, −73.82 → **27.9892** | ✓ |
| τ₀ = 15.0 | `taux_as_0 = 15e6` (:118) | ✓ |
| η = 0.89e-3 | `eta_v = 8.90e-4` — *"assume water, at 25 degree celsius"* | ✓ |
| `rw` = 0.089 | `R_w = 0.1778/2` | ✓ |
| `Sw_fwid` = 7.4e-9 | V_w·β_w/w = **7.4401e-9** | ✓ derivation |

Three corrections that follow:

1. **Taiyi's implied μ₀ is 0.5359**, not 0.50 (15.0 / 27.9892). "0.5 is
   Taiyi-equivalent" only holds if σ̄₀ is also moved to 30.
2. `setup_model.m` ships `k = 7e-13` **uniform** (:95, "from trial and error").
   The 1.1e-12 / 4e-13 / 150 m values are from **Table 1 of the paper**, cited at
   `injection2.m:70-71`. The code default and the paper differ.
3. His fluid viscosity is water at 25 °C for a 4.1 km geothermal well, which is
   why this study uses **η = 1.27e-4** — a 7.008× change in every diffusivity.

## A bound that needs no simulation

Slip needs Δp > Δp_crit = σ̄₀(1 − μ₀/f₀). Peak **measured** downhole overpressure
in the first 5 d is **10.92 MPa** (p_wh + ρgH − P0). So:

| σ̄₀ | minimum μ₀ | minimum τ₀ | |
|---|---|---|---|
| 30.0 | 0.3816 | **11.45 MPa** | 1807/1808's round number, not a measurement |
| 27.99 | 0.3659 | **10.24 MPa** | measurement-derived |

The measurement-derived stress admits a fault **1.21 MPa more understressed**. It
is both the more defensible choice and the more permissive one, so σ̄₀ is treated
as an axis and every run is twinned on it.

## What is here

**`presubmit_stage1/`** — 8 decks built and checked before submission:
res1807/res1808 physics on fixed code, uniform initial perm, twinned on σ̄₀ and
on fluid combination (A = parent's η and β; B = η 1.27e-4 with β × 7.007874 so D
is unchanged). Submitted as 632800–632807.

Caveat on B: holding D fixed does **not** reproduce A. The Peaceman well index
`T = 2πk/η` (`m_diffusion.f90:669`) depends on η alone, so the wellbore coupling
`gamma` is not preserved — 0.977 → 0.858 at the initial `kp`, and 0.145 → 0.024
once enhancement reaches `kpmax`. Also φβ = 1.577e-9 for B is ~3.9× above what
φ ≤ 0.02, β ≤ 2e-8 allow; deliberate, for diffusivity consistency.

**`632522/ 632523/ 632524/ 632525/`** — already run on fixed code, 8 d, June
injection, consistent bounds. These are res1807/res1808 plus a two-zone map, with
and without enhancement. Scored on **0–5 d** (λ_obs = 0.1866):

| run | parent | permev | λ/λ_obs | wellhead | peak slip / dc |
|---|---|---|---|---|---|
| 632522 | 1807 | **T** | 0.66 | **+72%** | 401 |
| 632523 | 1807 | F | 0.57 | +75% | 1065 |
| 632524 | 1808 | **T** | 0.73 | **+310%** | 1105 |
| 632525 | 1808 | F | 0.62 | +312% | 1092 |

Two things to read off this:

- **Enhancement helps the front by ~17%** and leaves the pressure essentially
  untouched (632522 vs 632523, 632524 vs 632525). That is a real effect in the
  right direction, and small.
- **Window dependence is severe.** These same runs score λ/λ_obs 0.97–0.99 on a
  0–8 d window, because λ_obs is 0.1866 over 5 d and 0.1412 over 8 d — the
  observed front rises steeply early. Any front number quoted without its window
  is meaningless.

**`perm_maps/`** — initial permeability for the map families, 2D and radial.

**`GRID_coverage.txt`** — generated by `docs/cooper_basin/grid.py`: the axes, and
which cells have a finished run on fixed code that reached the window.

## Tooling

`docs/cooper_basin/` — `build_stage1.py` (decks), `presubmit.py` (pre-submittal
package), `grid.py` (axes + coverage), plus the figure scripts.

## Caveat no parameter choice fixes

a > b everywhere, so these runs are velocity-strengthening and produce no
seismicity. The comparison is an *aseismic slip front* against an observed
*seismicity front*.
