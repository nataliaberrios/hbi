# Is there literature justification for varying `f0`?

Short answer: **yes, and the strongest argument comes from Wang & Dunham's own
paper.** They invoke phyllosilicates in the Habanero fault zone to justify
velocity-strengthening behaviour, then keep the bare-granite friction value.
Those two choices are inconsistent, and resolving the inconsistency is what
licenses a sweep.

## What `f0` is, and what it is not

In HBI these are different objects and must not be conflated:

- **`muinit`** is an initial condition on stress. `main_LH.f90:861`,
  `tau = sigma*muinit`, executed once. Taiyi's move — raising initial shear
  stress from the resolved 10.3 MPa to a calibrated 15 MPa — lives here.
- **`f0`** is the nominal friction coefficient in the rate-and-state law, passed
  per-cell into `deriv()` at every timestep (`main_LH.f90:1908`). It is a
  material property of the fault rock.

Varying `f0` is therefore **not** an alternative way of doing what Taiyi did. It
is a claim about what the fault gouge is made of.

## The case that `f0 = 0.6` is the wrong material value here

**Wang & Dunham already argue the fault zone is phyllosilicate-bearing.** From
their discussion, justifying a velocity-strengthening main fault:

> "in hydrothermal environments, the preferential flow of hydrothermal fluid
> through the fault zone can form phyllosilicates, a mineral group known for
> promoting stable fault sliding. Hydrothermal fluid can form phyllosilicates
> through either weathering of the granitic host rock or precipitation."

They use phyllosilicates to set the *velocity dependence* (a > b) and bare
granite to set the *friction level* (f0 = 0.6). The same mineral group controls
both. If phyllosilicates are present in sufficient quantity to make the fault
velocity-strengthening, they are present in sufficient quantity to lower its
frictional strength.

**Their own stated reason for rejecting the measured stress is a statement about
f0.** They compute the resolved shear stress as 10.3 MPa on a 101.8 MPa normal
stress, note the ratio is "much less than the typical granite friction
coefficient of 0.6", conclude the fault is too far from failure to produce the
observed seismicity, and calibrate the initial shear stress upward instead. The
incompatibility they identify is between the resolved stress and an *assumed*
f0. It is not a stress measurement problem.

## Published values

| material | friction | source |
|---|---|---|
| most rocks except clays, high normal stress | **0.6 ± 0.05** | Aharonov & Scholz-type compilation via Barbot (2024) |
| grain-boundary reference friction at 1 µm/s | **0.60** | Chen & Spiers, via Chen et al. (2017) |
| K-feldspar / albite gouge, 120–400 °C | **0.71** (0.70–0.85 with pore pressure) | Hu et al. (2024) |
| granite/chlorite mixed gouge | decreases **monotonically** with chlorite content | Zhang et al. (2022) |
| pure chlorite gouge | **0.37** | Zhang et al. (2022) |
| illite | **0.28** | via Zhang et al. (2022) |
| talc | **0.20** | via Zhang et al. (2022) |
| montmorillonite | **0.13** | via Zhang et al. (2022) |

Two things this table settles:

1. **0.6 is well grounded for unaltered granite** — if anything conservative,
   since the primary minerals (feldspar 0.71) are stronger. So the sweep must be
   justified by alteration, not by claiming 0.6 is loose.
2. **The route below 0.6 runs through phyllosilicates specifically.** Chlorite
   caps the plausible weakening at 0.37; illite and clays go lower but are less
   likely in a 250 °C granite.

## Relevance of the conditions

Habanero is granite at ~250 °C and ~4.1 km, repeatedly injected.

- **Zhang et al. (2022)** ran granite/chlorite gouges explicitly for
  injection-induced seismicity in EGS, and report measured chlorite contents in
  real reservoirs: **Pohang < 20 wt.%, Gonghe up to 35 wt.%**. So intermediate
  chlorite fractions — and hence intermediate friction — are what geothermal
  granite reservoirs actually contain.
- **Jeppson et al. (2023)** slid Westerly granite at up to 250 °C and found
  time-dependent **weakening** for holds beyond 14 h, requiring "a second,
  strongly negative, state variable," which they attribute to "significant
  hydrothermal alteration ... consistent with microstructural observations of
  dissolution and secondary mineral precipitation." Habanero sits at exactly the
  temperature where they see this.
- **Zhang et al. (2022)** also report that chloritization *increases* (a − b),
  i.e. pushes gouge toward velocity-strengthening. That is the same direction as
  the a > b assumption already in every HBI deck here, so lowering `f0` on
  chloritization grounds is *consistent* with the existing friction setup rather
  than in tension with it.

## Proposed sweep

Hold every measured quantity fixed: `sigmainit 27.99`, `muinit 0.3666`
(τ₀ = 10.26 MPa, the resolved shear stress), and vary only `f0`.

| `f0` | interpretation | Δp_crit MPa | margin vs measured 10.92 |
|---|---|---|---|
| 0.60 | unaltered granite — Taiyi's assumption | 10.89 | +0.03 |
| 0.55 | lightly chloritized | 9.33 | +1.59 |
| 0.50 | ~20–35 wt.% chlorite, i.e. Pohang/Gonghe range | 7.47 | +3.45 |
| 0.45 | strongly chloritized | 5.19 | +5.73 |
| 0.40 | approaching pure chlorite (0.37 floor) | 2.34 | +8.58 |

**Hard floor: `f0` must exceed `muinit` = 0.3666.** Below that the fault is
already at failure with zero overpressure and the initial condition is
meaningless. That is why the sweep stops at 0.40 and not at chlorite's 0.37.

## What this does and does not establish

**Does:** it makes the model internally consistent. If the fault zone is
phyllosilicate-bearing enough to slide stably, its friction is below bare-granite
values, and the stress paradox that forced Taiyi to calibrate τ₀ upward may not
exist.

**Does not:** prove the Habanero fault is chloritized. **I could find no
published mineralogy for the Habanero fault zone** — no core, cuttings or
XRD study. Searches for Habanero/Cooper Basin/Innamincka combined with
chlorite/illite/clay/alteration returned nothing on the granite fault zone; the
one Cooper Basin hit concerns sedimentary diagenesis. Wang & Dunham's own
phyllosilicate argument is likewise generic — "in hydrothermal environments" —
not site-specific. So the strongest honest framing is:

> The phyllosilicate assumption is already load-bearing in the published model.
> This sweep tests what follows if it is applied consistently to friction as well
> as to velocity dependence.

**Competing explanation that must be acknowledged.** The paper states that
"estimates of horizontal stress are contingent on assumptions of rock strength
and therefore have significant uncertainties." The failure margin moves 1 MPa per
MPa of normal stress, so a 4 MPa error in σ_n does the same work as f0 = 0.487.
That is a cheaper hypothesis than alteration and cannot be excluded.

**Provenance note.** The paper reports S_Hmax = 150 MPa and S_hmin = 120 MPa;
`setup_model.m` uses `s_Hmax = 160e6`. Both give σ_n ≈ 101.8 MPa, so conclusions
are unaffected, but code and paper disagree.

## References

Chen, J., Niemeijer, A. R., & Spiers, C. J. (2017). Microphysically derived
expressions for rate-and-state friction parameters, a, b, and Dc. *Journal of
Geophysical Research: Solid Earth*. https://doi.org/10.1002/2017jb014226

Hu, Z., Zhang, C., & Zhang, L. (2024). Frictional properties of feldspar-chlorite
gouges and implications for fault reactivation in hydrothermal systems. *Earth
and Space Science*, 11(7). https://doi.org/10.1029/2023ea003492

Jeppson, T. N., Lockner, D. A., & Beeler, N. M. (2023). Time-dependent weakening
of granite at hydrothermal conditions. *Geophysical Research Letters*.
https://doi.org/10.1029/2023gl105517

Wang, T., & Dunham, E. M. (2022). Hindcasting injection-induced aseismic slip and
microseismicity at the Cooper Basin Enhanced Geothermal Systems Project.
*Scientific Reports*, 12(1). https://doi.org/10.1038/s41598-022-23812-7

Zhang, F., Huang, R., & An, M. (2022). Competing controls of effective stress
variation and chloritization on friction and stability of faults in granite:
Implications for seismicity triggered by fluid injection. *Journal of Geophysical
Research: Solid Earth*, 127(8). https://doi.org/10.1029/2022jb024310

Additional context consulted: Transient and steady-state friction in non-isobaric
conditions, https://doi.org/10.1029/2023gc011279 (the 0.6 ± 0.05 compilation);
Frictional stability of laumontite under hydrothermal conditions,
https://doi.org/10.1029/2023gl108103 (Gonghe geothermal reservoir).
