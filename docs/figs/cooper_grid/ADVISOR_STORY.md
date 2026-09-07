# Cooper Basin parameter search — the ten figures, and what they say

Written for a 20-minute walkthrough. Every number here is reproducible from
`grid_scores.json` and the scripts in `docs/cooper_basin/`.

**Two scoring windows appear below and they are not interchangeable.** Scores on
**0–5 d** are used when all 82 runs must be comparable, because most runs stop
there. Scores on **0–18 d** are used where longer runs exist, and that is the
honest window — §3 is entirely about why. Never quote one without saying which.

---

## The question

The Cooper Basin stress measurements (Holl & Barton 2015) put the resolved shear
stress on the Habanero Fault at τ₀ ≈ 10.4 MPa. Wang & Dunham (2022) reproduce
the data with τ₀ = 15.0 MPa — a critically stressed fault. **Can an
understressed fault, at the measured stress, reproduce the same observations if
it has permeability enhancement, which their model does not?**

Three observables: the **slip front** radius R(t), the **wellhead pressure**
history, and the **cumulative slip** magnitude.

---

## 1. The search was exhaustive — and the target corner is empty

**`story/story_search_summary.png`**

82 runs; 55 have all three observables measurable. Swept: τ₀ (11.0–15.0 MPa,
36 runs), `kpmax` (4 decades), permeability structure (uniform vs two-zone,
50× and 250×), enhancement on/off, viscosity and storage (two fluids,
diffusivity-matched), porosity grading, and friction `a`/`b`.

Panel B plots every run in (front, wellhead) with slip as colour. **Nothing
lands in the target corner.** Pairwise, at ±15%:

| pair | runs matching both |
|---|---|
| front + wellhead | 7 |
| front + slip | **0** |
| wellhead + slip | **0** |

All three at once: **0 runs, at ±15%, ±25% and ±50%.**

## 2. The invariant — the single most robust result

**Same figure, panel C.**

Six runs reproduce the observed slip to within 3%. **All six sit at +70 to +80%
wellhead** — spanning τ₀ = 10.4–11.1 MPa, uniform and two-zone permeability
maps, and every one with enhancement on:

| run | τ₀ | front | wellhead | slip/obs |
|---|---|---|---|---|
| 632874 | 10.36 | 0.68 | +70.2% | 0.97× |
| 632812 | 10.36 | 0.70 | +70.7% | 1.03× |
| 632913 | 10.36 | 0.70 | +70.7% | 1.03× |
| 632522 | 11.10 | 0.66 | +71.8% | 1.01× |
| 632801 | 10.36 | 0.79 | +77.8% | 1.04× |
| 632800 | 11.10 | 0.76 | +79.6% | 1.02× |

**The model cannot produce the observed slip on less than ~1.75× the measured
wellhead pressure.** That is a property of the model, not of a tuning choice.

## 3. The front "matches" were a 5-day artifact

**`stage11/stage11_front_window.png`**

Every earlier claim of a front match used a 0–5 d window. Refit on matched
windows with the observed front refit too:

| run | front 0–5 d | front 0–18 d |
|---|---|---|
| 632915 | 1.04 | **1.66** |
| 632917 | 1.03 | 1.49 |
| 632916 | 0.97 | 1.45 |
| 632913 | 0.70 | **1.08** |

λ_obs is 0.1866 over 0–5 d and 0.1709 over 0–18 d, so this is not a convention
artifact — the simulated front keeps growing while the observed one does not.
**The run that looked worst at 5 d is the only one in band at 18 d.** On the
honest window, zero runs match even the original two targets.

The left panel carries a methodological point worth making: **neither the
simulated nor the observed front is really √t.** Both step, with plateaus
through the shut-ins. λ from R = λ√t is this project's headline metric and it is
fitting a shape neither curve has.

## 4–6. The comparison set, on the honest window

Five runs spanning the distinct behaviours, all 18 d, all fits 0–18 d.

**4. `compare_story5_18d/pressure_compare_story5_18d.png`** — wellhead.
Note the measured trace crashing to zero at 1.58 d: the well is bled off, HBI
has no wellbore bleed-off, so shut-ins cannot be reproduced and the metric masks
them. On 0–18 d the flowing mask is 73% of the record; on 0–5 d it is 14%.

**5. `compare_story5_18d/RT_compare_story5_18d.png`** — front vs time.
**6. `compare_story5_18d/RV_compare_story5_18d.png`** — front vs injected volume.

| run | `kpmax` | front 0–18 d | slip 17 d | wellhead 0–18 d |
|---|---|---|---|---|
| **632913** | 2.5e-13 | **1.08** | **0.95×** | +58.0% |
| 632911 | 2.5e-12 | 1.57 | 0.55× | +9.9% |
| 632916 | 2.5e-11 | 1.45 | 0.17× | +8.8% |
| 632917 | 2.5e-11 | 1.49 | 0.09× | +8.3% |
| 632910 | 2.5e-11 | 1.70 | 0.16× | −1.0% |

Monotonic in `kpmax` on every column. One run gets the front and the slip; every
run that fixes the wellhead breaks both.

## 7. Slip distribution against the data

**`compare_story5_18d/slip_compare_story5_18d.png`** and
**`slip_poster_style/slip_poster_story5_18d.png`** (eight times, shared axis).

632913 against the observation at all eight times, cm:

| t | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 |
|---|---|---|---|---|---|---|---|---|
| observed | 2.81 | 2.81 | 3.02 | 4.37 | 4.99 | 5.40 | 7.56 | 9.17 |
| 632913 | 2.91 | 2.91 | 4.05 | 4.11 | 4.60 | 5.03 | 8.72 | 8.73 |

Within 10–35% everywhere, with the right growth rate. **Nothing else is within a
factor of two at 17 d** — and 632913 is the *original* configuration, not a
tuned one. Ten stages of tuning moved away from the answer, and the 5-day window
is why.

## 8. Why `kpmax` cannot deliver both — the mechanism

**`stage9/stage9_mechanism.png`**

The enhanced zone behaves as a sealed disc: a pressure plateau with a
logarithmic peak at the well and a cliff at the disc edge, coinciding with the
slip front to within 15–20 m. The peak follows Δp = qη/4π`kpmax` over three
decades. But the **plateau falls with `kpmax` too**, and slip tracks
plateau − Δp_crit — so `kpmax` buys the wellhead by draining the pressure that
drives the slip. That is the trade-off in §2, mechanistically.

## 9. It was never the stress state

**`story/story_activation_pressure.png`**, panel C, and
**`strength_margin/strength_margin.png`**

Holl & Barton give stress *ratios* (SHmax/Shmin/Sv ≈ 1.35–1.45/1.10–1.25/1.0,
Pp = 72.7 MPa), not magnitudes, and they are modelled from wellbore breakouts
with an assumed rock strength of 130–150 MPa; Shmin is bounded only by an
inequality. With Sv ≈ 95.3 MPa (from Holl §8.4: first seismicity 2.6 MPa above
Pp and ~20 MPa below the overburden), resolving onto an 18–20° thrust gives
τ = 9.8–12.6 MPa and σ̄ₙ = 25.8–26.7 MPa.

**τ₀ = 10.36 and σ̄₀ = 27.99 sit inside that band — and so does 15.0.** The
stress data cannot distinguish the two models. "They ignored the stress
measurements" is not a defensible claim, and neither is a search that only
varies τ₀.

The `strength_margin` figure shows why: their own front criterion
(`cmp_seis_extent.m`) pins the **strength margin** Δτc = f₀σ̄₀ − τ₀ ≈ 1.86 MPa,
not τ₀. Reproduced here from their Table 1, τ₀ = 15.0 at f₀ = 0.60 gives
196/340/438 m against an observed 187/323/417 — and τ₀ = 10.26 at f₀ = 0.433
gives 185/321/414. Two points on one locus; the data cannot separate them.

## 10. What the search actually diagnosed

**`story/story_activation_pressure.png`**, panels A and B.

The Habanero Field Development Plan (§4, Fig 4-3) and Holl (2015, §8.4) both
report the overpressure at which seismicity *began*: **Δp = 2.6 MPa** for the
H04 local stimulation and **Δp ≈ 0.4 MPa** for the extended November 2012
stimulation — the record this project models. That is Δp_crit, measured: no
stress tensor, no assumed friction, no fault dip.

| | Δp_crit |
|---|---|
| **measured, H04 2012** | **0.4–2.6 MPa** |
| Wang & Dunham (τ₀ 15.0, f₀ 0.60) | 3.0 |
| this model (τ₀ 10.36, f₀ 0.60) | **10.72** |
| 632913's measured plateau | 11.9 |

**The model needs 4–27× the overpressure at which the fault actually failed.**
That is the +58% wellhead restated as a measurement, and it explains all 82 runs
in one line: anything that reproduced the observed slip had to over-pressurise,
because the threshold it had to overcome was set an order of magnitude too high.

Requiring the measured threshold at the *measured* τ₀ and σ̄₀ gives

**f₀ = 0.375 – 0.408**

This is not the f₀ = 0.433 argument retracted earlier in `CONCLUSIONS.md`. That
retraction was correct — 0.433 would need ~82 wt.% chlorite against a documented
maximum of 35 — but it answered "what gouge composition gives f₀ = 0.433?" The
field data asks "what f₀ is consistent with a fault that slipped at 0.4 MPa?",
and that question does not go through friction databases at all.

It constrains Wang & Dunham too: their Δp_crit = 3.0 MPa is closer but still
above the 0.4 MPa of the very stimulation they model, and they reach it by
raising τ₀ rather than lowering f₀. The activation pressure constrains the
combination, so their solution is strained by it as well.

---

## The one-paragraph version

Eighty-two simulations spanning permeability magnitude and structure,
permeability enhancement, fluid properties, porosity grading, initial shear
stress and rate-and-state friction fail to reproduce the Cooper Basin slip
front, wellhead pressure and slip magnitude simultaneously — at any tolerance.
The failure is not a tuning gap but a single invariant: every configuration that
reproduces the observed slip requires ~1.75× the measured wellhead pressure. The
Cooper Basin stress measurements turn out to admit both the understressed state
used here and the critically stressed state of Wang & Dunham (2022), so stress
is not what separates the models. What does is the frictional strength: the
field record shows the Habanero Fault slipping at 0.4–2.6 MPa of overpressure,
while a fault at the measured stress with laboratory granite friction
(f₀ = 0.60) requires 10.7 MPa. **The parameter search's value is negative but
specific — it excludes the whole permeability/stress space and localises the
discrepancy in f₀, for which the field activation pressure gives an independent
value of 0.375–0.408.**

---

## What to expect to be asked, and the honest answers

**"Is f₀ ≈ 0.4 physically justifiable for granite?"** Not from laboratory
granite gouge, which is 0.69–0.74 (Zhang et al. 2022), and f₀ = 0.60 is already
below that. The argument is observational, not mineralogical: the fault
demonstrably failed at ~1 MPa, and f₀ = 0.60 at the measured stresses cannot
produce that. Wang & Dunham already invoke phyllosilicates in this fault zone to
justify velocity-strengthening behaviour.

**"Have you tested f₀ ≈ 0.4 in HBI?"** In progress — 632924–632927, f₀ = 0.40
against f₀ = 0.60 controls at two mesh resolutions. Not yet reported. Prediction
on record: the pressure falls while the slip survives.

**"Could the discrepancy be in the comparison rather than the model?"** Partly,
and it is bounded. HBI has no wellbore bleed-off so shut-ins are excluded from
the pressure metric (73% of the 0–18 d record survives). But +58% is far too
large to attribute to that, and the activation pressure is a formation
measurement that bypasses the wellhead comparison entirely.

**"What physics is missing?"** Dilatancy (a dead stub in the code), fault-zone
viscous flow (implemented, never used), and anisotropic enhancement — Holl §8.4
reports the enhancement *is* anisotropic, oriented **perpendicular** to slip,
with Llanos et al. (2015) using 2:1. Of these, only anisotropy has observational
support here; dilatancy has the wrong sign, since it makes slip *harder* and the
model already needs too much pressure.

**"What was never varied?"** `dc`, initial velocity, rigidity, well storage,
skin, and the three physics extensions above. `dc` is the most exposed: the front
metric is a slip contour, and Wang & Dunham use 1.53e-5 against this project's
1e-4.
