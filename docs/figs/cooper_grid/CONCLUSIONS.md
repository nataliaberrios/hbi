# What was tried after the code fix, and why both targets were not matched

All runs below used the fixed `limitsigma` code, the November 2012 injection
record, and a single scoring convention: front and pressure both on **0–5 d**,
with the observed front refit on that window (λ_obs = 0.1866), and pressure
scored while the well is flowing. 60 runs, all verified on fixed code with
`kpmax` = map maximum and `kpmin` = map minimum.

## What was tried

**Stage 1 — 632800–632807.** res1807/res1808 physics, uniform initial
permeability (kp = kpmin = 1e-15, grown toward kpmax by slip). Twinned on
σ̄₀ (30.0 / 27.99 MPa) and on fluid: **A** = η 8.9e-4 with the parent's β,
**B** = η 1.27e-4 (correct reservoir viscosity) with β × 7.007874 so diffusivity
is unchanged.

| | front λ/λ_obs | wellhead |
|---|---|---|
| fluid A | 0.76 – 0.81 | +78 to +317% |
| fluid B | 0.24 – 0.28 | +24 to +57% |

**Stage 2 — 632810–632821, plus 632522–632525.** Two-zone initial permeability
with the near-well 150 m disc set to `kpmax` and the background to `kpmin`, at
both contrasts (250× and 50×), both σ̄₀, both fluids, and enhancement on and off.

| | front λ/λ_obs | wellhead |
|---|---|---|
| fluid A | 0.57 – 0.76 | +71 to +312% |
| fluid B, 250× map | **no slip** | **−0.5%** (RMS 1.4 MPa) |
| fluid B, 50× map | 0.02 | +34% |

Enhancement on versus off is worth about **+17% on the front** and leaves the
pressure untouched.

**Stage 3 — 632830–632867.** τ₀ swept 11.0 → 15.0 MPa in 0.5 MPa steps (μ₀
derived as τ₀/σ̄₀), at both σ̄₀ and both map contrasts, fluid B, enhancement on.

- 250× map: **no slip at any τ₀**, including Taiyi's own 15.0 MPa.
- 50× map: λ/λ_obs climbs monotonically **0.02 → 0.34**.
- Pressure is identical to three figures across the whole sweep — with
  negligible slip the permeability never evolves and the hydraulics decouple from
  friction entirely.

**Also tested earlier on fixed code** (632548–632552): the physically calibrated
far field (4e-13, so pressure reaches ~2.8 km) with enhancement at μ₀ 0.37 and
`kpmax` from 1.1e-12 to 1e-10 — no slip in any case.

**Result: none of the 60 runs lands inside both ±15% bands.** Best front with an
acceptable wellhead is λ/λ_obs = 0.34, at τ₀ = 15.0 MPa — i.e. already at Wang &
Dunham's stress state, and still three times short. Best wellhead is −0.5%, with
no slip at all.

## Why

The wellhead pressure and the slip front are both controlled by the same
quantity — the near-well pore pressure — and they demand opposite values of it.

The clearest measurement is a matched pair: 632810 and 632812 share the same
permeability map, the same μ₀ = 0.37, the same f₀ = 0.6. Only the fluid differs.

| | wellhead | max pf | pf at 150 m | slip front at 5 d | max slip |
|---|---|---|---|---|---|
| **632812** fluid A | 82.5 MPa | 35.96 MPa | 11.34 MPa | 255 m | 2.9e-2 m |
| **632810** fluid B | 42.4 MPa | 5.99 MPa | 3.06 MPa | — | 6.3e-9 m |
| *measured / observed* | **44.73 MPa** | | | **~417 m** | |

632810 reproduces the measured wellhead almost exactly and generates only 6 MPa
of pore pressure — slip stays seven orders of magnitude below dc. 632812 slips,
but only by pushing the wellhead to 82.5 MPa, nearly double the measurement.

And the slip front **trails** the pressure front: at 5 d, slip > dc reaches 255 m
while pressure above 0.1 MPa reaches 385 m, with 8–10 MPa of pore pressure
sitting at the slip front throughout. So the front is driven by pressure arriving
and lowering effective normal stress, and getting slip out to 417 m requires
pressure comfortably beyond that — which requires more injection pressure than
the wellhead record permits.

**Every parameter tried moves the two targets together rather than
independently:**

| knob | effect on front | effect on wellhead |
|---|---|---|
| viscosity + storage (fluid A → B) | 0.76 → 0.24, or to no slip | +317% → −0.5% |
| permeability contrast (250× vs 50×) | trades one against the other | trades the other way |
| permeability enhancement on/off | +17% | none |
| τ₀ (11 → 15 MPa) | 0.02 → 0.34 | none |
| σ̄₀ (30.0 vs 27.99) | small | small |
| initial perm structure (uniform vs two-zone) | small | small |

The two knobs that move the front appreciably — the fluid properties and the
permeability structure — are exactly the two that set the wellhead pressure. The
knobs that leave the pressure alone (τ₀, enhancement) are too weak to close a
factor of three in the front.

## Caveats that could change this conclusion

These are limits of the comparison, not of the parameter search, and any of them
could matter more than everything above.

1. **The model produces no seismicity.** Every deck has a > b, so the fault is
   velocity-strengthening throughout. The comparison is an *aseismic slip front*
   against an observed *seismicity front*. They need not coincide.
2. **Wang & Dunham did not match a slip front.** They calibrated τ₀ against
   *normalized cumulative seismic moment*, and modelled the seismicity with
   separate off-fault spring-sliders. So "matching the front" may not be the
   target their parameters were ever tuned to hit, and λ may not be the right
   metric for comparison with their result.
3. **The front is an arbitrary contour.** It is defined as slip > dc = 1e-4 m.
   An earlier test showed the threshold alone is worth about 1.6× in λ, and
   Wang & Dunham use dc = 1.53e-5.
4. **λ is strongly window-dependent.** λ_obs is 0.1866 over 0–5 d and 0.1412
   over 0–8 d. The same runs score 0.57–0.76 on the 5 d window and 0.97–0.99 on
   the 8 d window. Every number here is 5 d.
5. **Shut-ins cannot be matched at all.** HBI has no wellbore bleed-off, so the
   simulated wellhead stays near its flowing value while the measurement drops
   ~20 MPa. Pressure is therefore scored only on flowing periods, which are 14%
   of the 5 d window.
6. **No Habanero fault-zone mineralogy exists** that I could find, so the
   frictional parameters (`f0 = 0.6`, a, b, dc) rest on generic granite values.

## Honest statement of what this establishes

The claim supported by the data is **empirical, not mechanistic**: across 60 runs
spanning fluid viscosity, storage, permeability magnitude and structure,
permeability enhancement, initial shear stress and initial effective normal
stress, the slip front and the wellhead pressure move together, and no
combination places both inside ±15%.

It does **not** establish that a joint match is impossible. Untested directions
that could plausibly change the answer:

- a and b, dc, and the state evolution law — never varied; every deck uses
  a 0.015 / b 0.012 / dc 1e-4.
- `Sw_fwid`, `skin`, `pwinit` — never varied.
- the near-well disc radius (fixed at 150 m) and the fault zone width.
- a velocity-weakening patch, which would let the model produce seismicity and
  make the comparison like-for-like.
- an intermediate β at η = 1.27e-4, between the uncompensated 2.25e-8 and the
  D-matched 1.577e-7. This is the single gap most directly on the line between
  632810 (pressure right, no slip) and 632812 (slip, pressure wrong), and it was
  designed but never run.

---

# ADDENDUM: the actual reason, from Wang & Dunham's own front formula

Their repository contains a closed-form prediction for the seismicity front
(`source_code/cmp_seis_extent.m`, used at `projects/mk_plots.m:596-645`). It needs
no simulation, so it can be evaluated directly:

    Dp(r,t) = (Q0*eta/(4*pi*k*w)) * E1(phi*eta*beta*r^2/(4*k*t)) * exp(-dz^2/(2*std^2))
    front radius r solves  Dp(r,t) = Dtauc,   Dtauc = f0*sigmabar_0 - tau_0

Evaluated with their Table 1 values (Q0 20e-3, eta 8.9e-4, k 4e-13, w 6 m,
phi 0.01, beta 1e-8, sigmabar_0 28.0 MPa, dz 0.5 m, std 2 m):

| t (d) | observed front | their tau_0 = 15.0 | resolved tau_0 = 10.26 |
|---|---|---|---|
| 1 | 187 m | 196 m | 3 m |
| 3 | 323 m | 340 m | 5 m |
| 5 | **417 m** | **438 m** | **7 m** |

Their calibrated stress reproduces the observed front to about 5%. The resolved
stress predicts 7 m instead of 417 m — sixty times too small.

## Why the front is such a sharp constraint

Dp falls off **logarithmically** with radius (E1(x) ~ -ln x for small x), so the
radius at which Dp crosses a fixed threshold depends **exponentially** on that
threshold. The prefactor is A = Q0*eta/(4*pi*k*w) = 0.572 MPa, so:

| tau_0 MPa | Dtauc MPa | front at 5 d | vs observed |
|---|---|---|---|
| 10.26 (resolved) | 6.54 | 7 m | 0.02x |
| 12.00 | 4.80 | 31 m | 0.08x |
| 13.00 | 3.80 | 75 m | 0.18x |
| 14.00 | 2.80 | 181 m | 0.43x |
| **14.94** | **1.86** | **417 m** | **1.00x** |
| 15.50 | 1.30 | 691 m | 1.66x |

A 0.2 MPa change in tau_0 moves the front by 20%. This is not a soft target that
can be traded against pressure.

## This is why the parameter search could not succeed

The two observations constrain two different things, and between them leave no
freedom:

- **the wellhead record fixes the pressure amplitude A** — it is what "matching
  the pressure" means;
- **the front fixes the ratio Dtauc/A**, and therefore, with A pinned, fixes
  **Dtauc itself** to 1.86 MPa.

At f0 = 0.6 and sigmabar_0 = 28.0 MPa that forces tau_0 = 14.94 MPa, essentially
exactly the 15.00 MPa Wang & Dunham calibrated. So their value is not merely a
fit to seismic moment; the front requires it independently.

It also explains every HBI result above quantitatively. Our understressed runs
sat at tau_0 = 10.4-13.0 MPa with f0 = 0.6, i.e. Dtauc = 3.8-6.5 MPa, for which
the formula predicts fronts of 7-75 m. The only runs that reached 255 m did so by
raising A — over-pressurising to 82.5 MPa wellhead against a measured 44.73.
That is the trade-off, and it is not something a permeability, storage or
enhancement choice can escape, because those choices act on A, which the wellhead
already pins.

## What this leaves for the understressed hypothesis

The front constrains only the **combination** Dtauc = f0*sigmabar_0 - tau_0.
There are two ways to reach the required 1.86 MPa:

| route | tau_0 | f0 | keeps the stress measurement? |
|---|---|---|---|
| Wang & Dunham | **14.94** (calibrated) | 0.60 (lab granite) | no |
| alternative | 10.26 (resolved) | **0.433** | **yes** |

f0 = 0.433 sits inside the range published for chlorite-bearing granitic gouge
(pure chlorite 0.37; unaltered granite/feldspar 0.60-0.71; mixed gouge decreasing
monotonically with chlorite content). And Wang & Dunham already argue for
phyllosilicates in this fault zone, in order to justify velocity-strengthening
behaviour.

So the honest conclusion is not that an understressed fault is impossible. It is:

> **The observed front requires a strength margin of 1.86 MPa. An understressed
> fault can supply that only if f0 is about 0.43 rather than 0.60. Permeability,
> enhancement, storage and viscosity cannot supply it, because they act on the
> pressure amplitude, which the wellhead record already fixes.**

That is a single-parameter, falsifiable claim, and it is testable with a handful
of HBI runs at f0 near 0.43 rather than another grid.

## Limits of this addendum

- The formula is Wang & Dunham's **Coulomb triggering criterion for the secondary
  faults** (their spring-sliders), not HBI's rate-and-state main fault. f0 = 0.433
  is therefore the value required *in their framework*. Whether HBI's aseismic
  slip front responds identically has to be checked by running it.
- It assumes constant-rate injection at Q0 = 20e-3 m^3/s, whereas the real record
  is strongly variable, and uses their far-field k = 4e-13 with their viscosity
  8.9e-4 rather than the corrected 1.27e-4.
- dz = 0.5 m (the closest secondary-fault offset they plot) was used; larger dz
  reduces the predicted radius.

---

# CORRECTIONS to the addendum above, and the final state

## 1. f0 = 0.433 is NOT supported by the literature — retract that suggestion

The addendum proposed f0 = 0.433 as the way to keep the resolved shear stress.
Checking it against measured values kills it:

| | friction |
|---|---|
| 100% granite gouge (Zhang et al. 2022, their own measurement) | **0.69 – 0.74** |
| 100% chlorite gouge | 0.37 |
| decrease with chlorite content | monotonic |

Chlorite contents actually measured in granite EGS reservoirs: Pohang cores
9 wt.%, cuttings max 19 wt.% ("generally below 20"); Gonghe max ~35 wt.%. Their
own remark: "all much lower than the 50 wt.% content."

Interpolating between the measured endpoints, **f0 = 0.433 requires ~82 wt.%
chlorite.** The most ever documented is 35 wt.%, which gives f0 ~ 0.59 — worth
Dtauc = 6.3 MPa and a front of about 10 m.

Two consequences:

- The defensible range for f0 is roughly **0.59 – 0.74**, not the 0.40 – 0.60
  written in `F0_LITERATURE_REVIEW.md`. That document's proposed sweep is too low
  at its bottom end and should be read with this correction.
- **f0 = 0.6 is already below the unaltered-granite value of 0.69–0.74.** Every
  deck in this project, and Wang & Dunham's, already assumes some weakening.

So the f0 route cannot rescue the resolved stress. The internal-consistency
argument in the review — that invoking phyllosilicates for a > b while keeping
bare-granite friction is inconsistent — still stands, but the size of the
available effect is far too small to matter here.

## 2. No physical permeability and compressibility can do it, and this is exact

Not a search result — a 2x2 solve. Impose both requirements:

    (a) Dp(r_w)   = 10.93 MPa    reproduce the measured wellhead overpressure
    (b) Dp(417 m) = Dtauc        put the front at the observed radius

and solve for the k and beta that satisfy them:

| scenario | Dtauc | k needed | beta needed | physical? |
|---|---|---|---|---|
| resolved tau_0, f0 0.6 | 6.54 MPa | 1.3e-13 | **6.6e-18** | no |
| f0 0.59 (35 wt.% chlorite, the maximum measured) | 6.30 MPa | 1.2e-13 | **5.5e-17** | no |
| Wang & Dunham tau_0 15.0 | 1.80 MPa | 6.2e-14 | 1.0e-8 | **yes** |

beta = 6.6e-18 Pa^-1 is **eight orders of magnitude below the compressibility of
water** (4.4e-10). The required value does not exist in any rock or fluid.

Reason: matching the wellhead fixes Dp at the well to ~11 MPa. With the resolved
stress the front then needs 6.54 MPa at 417 m, i.e. the pressure may fall only
40% over 4700 well radii — a nearly flat profile. Dp falls logarithmically and no
physical beta flattens a logarithm that much. With Wang & Dunham's Dtauc = 1.80
MPa the profile may fall 84%, which a logarithm does naturally, and the solution
lands at k = 6.2e-14, beta = 1.0e-8 — essentially their published values. Their
parameters are close to the only ones that work.

Correction to the earlier wording: I wrote that permeability, storage and
viscosity "act on A, which the wellhead already pins." That was too glib — k
appears in the prefactor AND in the E1 argument, beta in the argument. The
correct statement is about profile SHAPE: for a uniform medium the whole profile
is one logarithm, so Dp(r_w)/Dp(417 m) is fixed by geometry regardless of k and
beta.

## 3. Injecting harder is not available, because the rate is also measured

HBI is rate-controlled: the measured injection rate is the input and the wellhead
pressure is an output. So "inject at higher pressure" means raising the rate,
which contradicts the record.

To get Dp = 6.54 MPa at 417 m with their published k, beta and phi you would need
**Q0 = 68 L/s, 3.4x the measured 20 L/s**, and the wellhead would then read
**74.4 MPa against a measured 44.73 — +66%.**

This is what run 632812 did in effect: it reached 255 m by running the wellhead to
82.5 MPa. The front is producible; it just costs the pressure match.

Both the injection RATE and the wellhead PRESSURE are measured. Two constraints on
the hydraulics leaves the fault's strength as the only freedom, which is why every
route ends up at Dtauc = f0*sigmabar_0 - tau_0.

## 4. Leakoff and hydraulic fracturing, assessed

**Leakoff** — worth doing for correctness, since Wang & Dunham explicitly neglect
it ("Leak-off outside the fault zone is neglected"), but the sign is predictable
and unhelpful: fluid leaving the fault zone means less Dp at large r for a fixed
measured wellhead, so the profile gets STEEPER. Wrong direction for the front.

**Hydraulic fracturing** — the right kind of mechanism, since an opening fracture
has high transmissivity and naturally flattens the profile. But the pressures do
not reach it:

| | | short by |
|---|---|---|
| peak absolute downhole pressure | 84.73 MPa | |
| sigma_n on the 10 deg fault | 101.8 | **17.1 MPa** |
| sigma_v (minimum principal, overthrust regime) | 100.0 | 15.3 MPa |
| sigma_hmin (paper) | 120.0 | 35.3 MPa |

Opening the fault plane needs a wellhead of 61.8 MPa against a measured 44.73
(+38%). Invoking it requires revising the stress state by ~17 MPa, which is more
than the ~6 MPa that would fix the front on its own — so it is a redundant
hypothesis rather than an independent one. Note HBI already has an `opening` flag
(main_LH.f90:2767), so no collaborator code is needed to test the sensitivity if
a stress state is ever adopted where jacking is reachable.

## 5. Stage 4, running: spatially graded porosity

The one remaining direction that requires no measurement to be overridden.

`beta` is a scalar in `t_params` and cannot vary in space, but `phi` can:
`phiG(:)` is a field, `case('phi')` in the parameter-file reader assigns it, and
the solver uses it in both `str = beta*phiG` (storage) and
`cdiff = kpG/(eta*beta*phiG)` (diffusivity). Lowering phi far from the well
therefore lowers storage AND raises diffusivity there — both carry pressure
further with less drop. They cannot be separated; phi enters both.

Runs 632870-632875, two bases x three gradings, phi held inside 0.005-0.02:

| base | what it is | gradings |
|---|---|---|
| 632810 | best PRESSURE match (wellhead -0.5%, no slip) | G1 0.020/0.005, G2 0.010/0.005, G3 inverse |
| 632812 | best FRONT (255 m, wellhead 82.5 MPa) | same three |

G3 inverts the grading and is a **control**: it should steepen the profile and
worsen the front. If it does not, the mechanism is not what is assumed here.

The kp column in every Stage 4 parameter file is copied byte-identically from the
base run's map, verified, so porosity is the only difference.

**Why this is not ruled out by section 2:** that solve assumed a uniform medium,
where the profile is a single logarithm. With phi(r) graded the ratio
Dp(r_w)/Dp(417 m) becomes a functional of the grading rather than fixed geometry,
so it is a genuinely different problem. Whether a factor of 4 in phi — the whole
physical range — buys enough flattening cannot be settled analytically, because
the graded problem has no closed form. That is what these runs test.

## Remaining directions, ranked

1. **Graded porosity** — running. No measurement overridden, code path exists.
2. **a, b, dc and the state evolution law** — never varied; every deck uses
   a 0.015 / b 0.012 / dc 1e-4. dc alone is worth ~1.6x in lambda.
3. **A velocity-weakening patch** — would let the model produce seismicity and
   make the comparison like-for-like instead of aseismic-front-vs-seismicity-front.
4. **Leakoff** — for correctness, expecting a worse front.
5. **Hydraulic fracturing** — only meaningful alongside a large stress revision.

---

# STAGE 4 RESULT: the control disproved the mechanism, and the search closes

All six graded-porosity runs completed, 5.00 d, fixed code, bounds verified.

| run | grading (phi near / far) | front lam/lam_obs | wellhead |
|---|---|---|---|
| **632875** | **G3 INVERSE 0.005 / 0.020** | **1.03** | +80.0% |
| 632874 | G2 0.010 / 0.005 | 0.68 | +70.2% |
| 632873 | G1 0.020 / 0.005 | 0.29 | +60.5% |
| *632812 base* | *uniform 0.01* | *0.70* | *+70.7%* |
| 632870-872 | all three, on the 632810 base | no slip | -2.6 to +3.9% |

## The control failed, and that is the useful part

G3 was included as a control that should make the front WORSE. It produced
**lam/lam_obs = 1.03, the best front in the entire 66-run project.** G1, the
grading argued to help, made it worse (0.29 against a base of 0.70).

The mechanism reasoning behind Stage 4 was therefore wrong. The claim was that
the front is limited by pressure REACHING far out, so lowering far-field porosity
would raise diffusivity there and extend it. In fact the front is limited by
pressure AMPLITUDE near the well: G3 lowers porosity near the well, which lowers
near-well storage and raises the near-well pressure, and that drives the slip.

This was already visible in data collected earlier and not connected: the slip
front TRAILS the pressure front (255 m against 385 m at 5 d). Pressure reach was
never the binding constraint.

Recorded plainly because it matters for how much weight the rest of this document
should carry: the mechanistic reasoning in this project failed three separate
tests (a Coulomb threshold that does not exist in rate-and-state, a
diffusion-length argument, and an elastic-amplification story), and each time it
was a measurement rather than an argument that caught it. The control is the only
reason the Stage 4 error was found.

## Why the search now closes, quantitatively

Graded porosity is a real lever -- it moved the front from 0.70 to 1.03 at fixed
permeability, fixed stress and fixed friction, and it is the FIRST parameter to
reach the front target at the measured stress state. But it buys the front the
same way everything else does, by over-pressurising.

Measured exchange rates, front gained per point of wellhead error:

| lever | lam per % |
|---|---|
| porosity grading | 0.0361 |
| map contrast (50x -> 250x) | 0.0182 |
| fluid / storage | 0.0098 |

632875 sits at lam 1.03, +80.0%. To enter both bands it must shed 65 points of
pressure while losing no more than 0.18 in lam -- a **required rate of 0.0028**.

The cheapest lever available is **0.0098, three and a half times too expensive.**
All three trade in the same direction (more pressure gives more front) and none is
flat enough to cut pressure without handing the front back. So no combination of
them reaches the target; they span a line, not the plane.

That is a stronger closing statement than "60-plus runs failed to find a match":
the levers are characterised, their exchange rates are measured, and the gap to
the target is a factor of 3.5 in a quantity that can be stated.

## What would still be worth doing, and why it is not more of the same

Nothing in the remaining list is another point in this family. Each changes the
structure of the problem rather than moving along the trade-off:

1. **a, b, dc and the state evolution law** -- never varied in any of the 66 runs
   (all use a 0.015 / b 0.012 / dc 1e-4). dc alone is worth ~1.6x in lambda and,
   unlike every lever above, it changes the front WITHOUT touching the pressure,
   so its exchange rate is not on the line at all. This is the most promising
   remaining direction and it is cheap.
2. **A velocity-weakening patch** -- would let the model generate seismicity and
   make the comparison like-for-like, instead of an aseismic slip front against a
   seismicity front.
3. **Leakoff** -- for correctness, expecting a worse front.
4. **Hydraulic fracturing** -- only meaningful alongside a stress revision of
   ~17 MPa, larger than the ~6 MPa that would fix the front unaided.

Point 1 deserves emphasis: dc sets the front CONTOUR as well as the friction
length scale, so it moves lambda without moving the wellhead at all. Every lever
measured above has a positive exchange rate; dc's is effectively infinite. It was
on the "not yet swept" list from the very beginning and never got swept.

---

# SESSION 2026-09-03: the dc confound, the Taiyi reference, and what his model actually is

Five things changed today. The first invalidates the paragraph immediately above.

## 1. CORRECTION: the dc sweep, as argued above, was not interpretable

The paragraph ending the Stage 4 section says dc "sets the front CONTOUR as well
as the friction length scale, so it moves lambda without moving the wellhead at
all... dc's exchange rate is effectively infinite." That treats the double role
as a *feature*. It is a **confound**, and the sweep as designed could not have
been read.

The front was measured as the contour where slip exceeds **that run's own dc**
(`make_sweep_figures.py:254`, `thr = ffloat(dk["dc"])`, and a duplicate in
`make_run_figures.py:157`). So dc entered the metric twice, **with the same
sign**:

- contour level — lower dc means more cells counted as slipped, bigger front
- nucleation size — lower dc means slip propagates further, bigger front

Nothing can separate those by inspection. A dc sweep scored that way would have
produced a clean monotonic trend that was part measurement artifact and part
physics in unknown proportion.

This was live, not hypothetical: 632880 uses Taiyi's dc = 1.53e-5 while the other
64 grid decks use 1e-4, so it would have been scored against a threshold **6.5x
lower** than everything it was being compared to.

**Fixed.** `FRONT_THR = 1e-4 m`, fixed for every run. 64 of 66 decks already had
dc = 1e-4, so every existing lambda is numerically unchanged — verified
bit-identical for all 41 runs that have one. `peak_slip/dc` still uses each run's
own dc, because that asks a different question (did the patch weaken at all).
Runs whose dc differs also record `lam_at_own_dc`, so the artifact is measured.

Threshold sensitivity, measured over a 100x range rather than assumed:

| run | peak/dc | lambda/lambda_obs, 1e-5 -> 1e-3 m | spread |
|---|---|---|---|
| 632875 | 605 | 1.04 -> 0.99 | 5% |
| 632874 | 271 | 0.70 -> 0.63 | 10% |
| 632800 | 47 | 0.76 -> 0.67 | 12% |
| 632873 | 116 | 0.35 -> 0.24 | 31% |

So for a well-developed front the metric is nearly threshold-free — the slip
profile is steep there. The threshold does real work only on marginal runs whose
front is a thin ring. A dc sweep is still worth running, but now with a known
noise floor: anything under ~12% is not physics.

## 2. Friction has never been swept, in any run

| parameter | across all 723 decks |
|---|---|
| `f0` | **0.6 in every single one** |
| `a` / `b` | 0.015 / 0.012 in 721 of 723 (the 2 others are 5-YEAR test decks) |
| `dc` | 1e-4 in 64 of the grid; 1.53e-5 only in the Taiyi pair |
| `muinit` | 19 distinct values |

Every "friction sweep" in this project was a sweep of initial shear stress.
`muinit` is an initial condition applied once (`main_LH.f90:861`); `f0` is a
material parameter passed into `deriv()` every timestep (`main_LH.f90:1908`).

An attempt to choose f0 values analytically **failed**, and the failure is worth
recording. Under Mohr-Coulomb, f0_req = (tau_0 + dp(R_obs))/sigmabar_0. But
**dp(417 m) is below 0.5 MPa in all 66 runs**; the maximum anywhere is 0.340 MPa
(632875). So f0_req collapses onto `muinit`, exceeding it by at most +0.012
across the entire grid, and the screen degenerates to "the fault must already be
at failure with zero overpressure." That is a statement about the screen, not
about f0 — HBI's regularised law has no threshold, and importing one has already
produced a wrong conclusion in this project once. See `f0_required.py`.

What the screen *did* establish: pressure falls by a factor of **56** between the
well and 417 m in the best run. So in 632875, which fits the front, the front
sits beyond both the Delta-tau_c contour (365 m) and any meaningful pressure. It
is carried by **elastic stress transfer from the slipping patch**, not by
pressure at the front.

## 3. The Taiyi reference runs: HBI reproduces his wellhead, not his front

First time in the project that HBI has been run on Wang & Dunham's own inputs.

| run | dc | reached | wellhead | slip front | Delta-p contour |
|---|---|---|---|---|---|
| 632880 | 1.53e-5 | **0.88 d** | +11.1% | 0.14 | 0.75x @ 0.88 d |
| 632881 | 1e-4 | 5.00 d | **+8.2%** | **0.22** | **0.68x** |

632880 stopped on HBI's own termination — `Slip rate below vmin at time step 348`
— not a crash or a step limit. The fault arrested at 0.88 d, just before the
first wellhead peak at 0.95 d.

632881's wellhead is **inside the +/-15% band**. Its slip front is 80 m against
an observed 417 m. Measured as a Delta-p = Delta-tau_c contour instead, the same
run reaches **380 m, or 91% of the observed radius**.

Both numbers are on a common window with lambda_obs refit on it; comparing a
0.88 d fit against a 0-5 d lambda_obs inflates 632880's contour score from 0.75
to 1.27, which is a window artifact and not a result.

## 4. What Wang & Dunham's model actually is, and why the comparison is delicate

His code is **two layers**, which is not documented anywhere else in this repo:

| | a - b | role |
|---|---|---|
| main fault, 215x215 over +/-10 km (**93.5 m cells**) | 0.015 - 0.012 = **+0.003** | velocity-strengthening, slips aseismically, cannot nucleate |
| **1000 off-fault spring sliders**, 20x20 m | 0.015 - 0.018 = **-0.003** | velocity-**weakening**; these are the earthquakes |

`seismicity.m` gives each slider its own rate-and-state solve driven by pore
pressure interpolated from the main fault **plus** elastic stress transfer from
main-fault slip (`dtaux = M_as_ss{1} * Dx_as`), with radiation damping.
`find_quakes.m` calls an event a slip-velocity peak of prominence >= 0.1 m/s. He
then classifies each event by which mechanism dominated — pressure weakening
(`wk_trigger`) versus aseismic loading (`as_trigger`) — and that classification
is the substance of his seismicity figure. `cmp_seis_extent.m` is an analytical
envelope plotted alongside, not the model output.

**HBI's fault is his main fault** — a - b = +0.003, identical, velocity-
strengthening. HBI has no spring-slider layer.

This project's working assumption is that the slip front does not care whether
the slip was seismic, so HBI's total slip front is the right thing to compare
against the catalogue. That assumption is what makes the comparison legitimate;
it should be stated explicitly wherever the comparison is made, because in his
model the seismicity front and the main-fault aseismic front are different
objects joined by a layer HBI does not have.

## 5. Resolution is NOT constant across the grid, and it is not a small difference

| | cell size | 417 m front = |
|---|---|---|
| 62 grid runs | 5 m | 83 cells |
| Taiyi pair (632880/632881) | **20 m** | 21 cells |
| Wang & Dunham's own main fault | **93.5 m** | 4.5 cells |

The Taiyi pair needs the wider box: their far-field k gives D = 4.49 m^2/s, so
the 5 d diffusion length is 2.79 km, against a half-domain of only 1.5 km at
ds 5 m on a 601 grid.

**Consequence: 632875 (lambda 1.03, ds 5 m) and 632881 (lambda 0.22, ds 20 m) are
not resolution-comparable at face value.** Part of that gap may be a 4x cell-size
difference rather than parameters. This is untested and is the cheapest
outstanding check: rerun 632875's configuration at ds 10 and 20 m — coarser means
a *bigger* domain for that config, since its kpmin is 1e-15.

## Where the grid stands, as a census rather than a claim

| | count of 68 |
|---|---|
| wellhead inside +/-15% | **27** — but **25 of them produce no slip at all** (pk/dc exactly 0) |
| wellhead matched AND slipping | **2** — 632880, 632881, the Taiyi pair, both `permev F` |
| slip front inside 0.85-1.15 | **1** — 632875, wellhead +80% |
| **both** | **none**; closest is 632867 at 4.4x the band width |

The muinit sweep is exhausted. tau_0 moves the front **without moving the
wellhead at all** (flat to +/-0.1% across all 36 Stage-3 runs), so it is a
genuinely orthogonal lever — but it saturates at lambda 0.34 by tau_0 = 15.0, and
tau_0 cannot exceed **f0*sigmabar_0 = 16.79 MPa** without the fault being past its
own strength at zero overpressure. At that hard ceiling the front reaches ~0.48.

## Stage 6, submitted today

632884 <- 632880 and 632885 <- 632881, each differing from its parent in exactly
four keys: `filenumber`, `permev` F -> T, `kpmax` 1.1e-12, `kpmin` 4e-13.

The Taiyi pair is the **only** base that matches the wellhead and also slips, and
both have enhancement **off**. So the hypothesis this project rests on — that
permeability enhancement plus a nonuniform initial permeability reaches the match
where Wang & Dunham could not — has never been tested on the one configuration
that matches the wellhead. 632881 needs its front to grow 4.5x; the alternative
base, 632875, needs its overshoot cut 5.3x against a direct measurement.

Bounds were read from the map, not typed: `perm_taiyi_601_ds20.txt` verified as
177 cells at 1.1e-12 (equivalent radius 150 m) and 361024 at 4e-13 — the paper's
Table 1 values — with the parents already carrying `kp 4e-13`. eta stays at
Taiyi's 0.89e-3 deliberately, so exactly one thing changes and the result is
attributable; an eta = 1.27e-4 variant is the follow-up, not a substitute.

Falsification conditions stated in advance: if the front stays near 0.22,
enhancement is not the missing ingredient at a wellhead-matching pressure. If the
front grows but the wellhead leaves the band, enhancement is on the same
trade-off line as porosity grading, map contrast and fluid storage, and does not
span the plane either.

## Two unflagged bad figures found while auditing

`BOUNDS_AUDIT.md` now labels the 90 pre-grid folders in this directory: 48
MISCONFIGURED, 16 VALID, 26 NOT APPLICABLE (`permev F`). Cross-referencing that
against `make_sweep_figures.py` turned up three sweep figures built on
misconfigured runs, of which only one was flagged:

| figure | bad runs | was it flagged |
|---|---|---|
| `sweep_kpmax.png` | 3 of 4 | SUPERSEDED |
| `sweep_muinit_permevT.png` | **3 of 3**, kpmax 18x the map max | **no — and not even listed in its README table** |
| `sweep_muinit_permevT_kmax1e-10.png` | **3 of 3**, kpmax 91x | **no** |

Both unflagged figures show mu_0 with enhancement ON producing almost no slip,
which is the misconfiguration rather than a property of enhancement — with
consistent bounds the same lineage reaches lambda/lambda_obs 0.97-0.99. They are
exactly the figures that would support the retracted claim that enhancement kills
slip. Now marked SUPERSEDED with their ratios and an explicit instruction not to
cite them as evidence about enhancement.

## SESSION 2026-09-05 — Stage 9: the first double match, and the constraint it exposed

**A run is inside both bands.** 632897 (`kpmax` 2.5e-11) scores front 0.98x and
wellhead +12.2% at **tau_0 = 10.36 MPa** — 31% below Wang & Dunham's 15.0. The
two targets are not intrinsically coupled, and matching them does not require
abandoning the stress measurements. That was the question this whole grid was
built to answer.

**What separates them is `kpmax`, and the reason is worth stating precisely.**
The enhanced near-well zone behaves as a sealed disc: overpressure forms a
plateau across it with a logarithmic peak at the well and a cliff at the disc
edge. `deltaP = q*eta/(4*pi*kpmax)` predicts the peak should scale as 1/`kpmax`,
and it does, over three decades and to within the linewidth (see
`figures/stage9/stage9_mechanism.png`, middle panel):

     run     kpmax   plateau   peak above   dp(r_w)   front   wellhead   peak slip
  632875   2.5e-13     11.90        7.05     18.96    1.03      +80.0%     6.03 cm
  632896   2.5e-12     10.66        0.70     11.36    0.99      +19.1%     1.27
  632897   2.5e-11      9.78        0.07      9.85    0.98      +12.2%     0.45
  632898   2.5e-10      8.63        0.01      8.64    0.59      +10.7%     0.03

The front holds at 0.98-1.03 across a 100x change in `kpmax`, breaking only at
2.5e-10 where slip nearly dies — so the front is set by volume balance into the
disc (phi*beta and dp_crit), not by `kpmax`. Every one of the previous 70 runs sat
at `kpmax` ~2.5e-13, which is why the two targets had always appeared locked.

**Two corrections to what this file said before.**

1. The plateau is NOT pinned at dp_crit. It falls with `kpmax` too, 11.90 to
   8.63 MPa. The "pinned" reading came from measuring over a single decade.
2. The wellhead bias has a FLOOR near +10.7%, reached by 632898 even at
   dp(r_w) = 8.64 MPa against a measured peak of 10.92, so pushing `kpmax`
   higher cannot close it and costs slip.

   **CORRECTION, added after inspecting the curve rather than the score.** An
   earlier version of this section attributed that floor to the shut-in periods,
   which HBI cannot follow without wellbore bleed-off. That is wrong: the
   pressure metric MASKS shut-ins out (`q > 25%` of peak and `p_obs > 5 MPa`,
   `score_grid.py:188`). The residual is a ramp-SHAPE mismatch inside the
   flowing window. Decomposed for 632901:

         window        sim    measured     bias
       0.50-0.80 d   40.99      35.94    +14.1%
       0.80-1.10 d   44.64      41.50     +7.6%
       1.10-1.40 d   42.93      37.08    +15.8%
       1.40-1.65 d   44.55      43.09     +3.4%
       PEAK          45.42      44.68     +1.6%

   The simulation pressurises too EARLY on the ramp and does not dip enough
   during the mid-stage rate reduction near 1.1-1.4 d. The headline +9.7% is
   those two errors partly cancelling, not a uniform offset. The defensible
   claim is that the PEAK wellhead pressure matches to 1.6% while the history
   matches only to about +/-16%, and closing that is a near-well
   storage/permeability question (`Sw_fwid`, or the disc permeability being too
   low early so pressure builds faster than it spreads) — NOT a tau_0 or
   `kpmax` question.

3. SCOPE OF THE PRESSURE SCORE, which should be quoted whenever the score is.
   The flowing mask covers t = 0.62 to 1.58 d, **14% of the 0–5 d window** —
   about 23 hours. Outside it the measured surface gauge reads ~0 because the
   well is bled off, and no formation model reproduces that. "Matched the
   wellhead" means matched over that 23-hour flowing window; saying it without
   the qualifier overstates the result.

**THE BINDING CONSTRAINT IS NOW SLIP AMPLITUDE, which was not one of the two
original targets.** Observed cumulative slip at the injector at 5 d is 2.81 cm.
632897 gives 0.45 cm, 6x low; 632896 gives 1.27 cm, 2.2x low; 632875 gives
6.03 cm, 2.2x high. The front is matching by RADIUS while the slip producing it
is too small — a thin ring. Slip tracks `plateau - dp_crit`, the excess
overpressure over the failure threshold, and that excess goes negative exactly
where the wellhead comes right.

**Stage 10 follows from this and is the natural next test.** At `kpmax` 2.5e-11
the plateau is fixed at 9.78 MPa, so raise the excess by lowering
dp_crit = sigma_0 - tau_0/f instead — a small `muinit` increase. Stage 3 measured
the wellhead flat to +/-0.1% across a 36-run tau_0 sweep, so this should buy slip
without spending the wellhead match:

    muinit   tau_0    dp_crit   excess   expected slip
     0.370   10.36      10.73    -0.95    0.45 cm (632897, measured)
     0.385   10.78      10.03    -0.25    ~1 cm      632900
     0.397   11.11       9.47    +0.31    ~2-3 cm    632901
     0.410   11.48       8.86    +0.92    ~5 cm      632902

The right panel of the Stage 9 figure crosses the observed 2.81 cm at an excess
of about +0.35 MPa, which is what 632901 targets. tau_0 = 11.11 MPa is still 26%
below 15.0, so this stays an understressed fault.

Falsifiable both ways: if the wellhead moves with `muinit` here, the Stage-3
independence does not survive at high `kpmax` and the targets are coupled after
all. If slip rises but the front leaves its band, the front was not
volume-controlled and Stage 9's agreement was luck.

## SESSION 2026-09-05, later — Stage 10: the tau_0 route to slip amplitude FAILS

**The hypothesis is falsified.** Raising `muinit` at fixed `kpmax` 2.5e-11 was
meant to lower dp_crit and so raise the excess `plateau - dp_crit` that slip
amplitude tracks. Predicted ~2-3 cm for 632901. Measured:

     run   muinit   tau_0   dp_crit   plateau   excess     slip   front   wellhead
  632897    0.370   10.36     10.73      9.78    -0.95   0.450cm   0.98     +12.2%
  632900    0.385   10.78     10.03      8.92    -1.11   0.472     1.02     +10.8%
  632901    0.397   11.11      9.47      8.24    -1.23   0.488     1.06      +9.7%

Slip moved 8%, not 6x. **The plateau tracked dp_crit downward at slope 1.22**, so
the excess got MORE negative rather than turning positive. In the tau_0
direction the plateau really is pinned to the failure threshold -- Stage 9's fall
was `kpmax` draining the disc, a separate mechanism. Both are true, and together
they mean **the excess is nearly invariant at fixed `kpmax`**. tau_0 cannot buy
amplitude, and no pressure route can: positive excess needs low `kpmax`, which
is exactly what breaks the wellhead.

**632902 (muinit 0.410) DID NOT COMPLETE** -- "Maximum iteration" at step 6102,
t = 4.51 d against a 5.00 d tmax, so it is below the 4.75 d usability gate and
is not scored. The RK solver stiffens as tau_0 rises, which is itself a limit on
this direction.

**What Stage 10 did establish, and it is worth keeping.** The front stayed in
band across the sweep (0.98 -> 1.02 -> 1.06) and the wellhead IMPROVED
(+12.2 -> +10.8 -> +9.7%), its best value in the project. So Stage 3's
tau_0/wellhead independence survives at high `kpmax`, and the double match is
robust across tau_0 = 10.4-11.1 MPa rather than being one lucky point.

**Where the project actually stands** (632901, and quoting the scope):

  * slip front: MATCHED. lambda 0.1971 vs lambda_obs 0.1866, 1.06x. Peak slip is
    49x dc, so it is a real front, not the razor-thin ring the plan warns about.
  * wellhead PEAK: MATCHED, +1.6% (45.42 vs 44.68 MPa).
  * wellhead HISTORY: close, not matched. Right magnitude, wrong ramp shape,
    +/-16% within the flowing window. See the correction above.
  * slip AMOUNT: 5.7x low, 0.49 cm against 2.81 observed at 5 d.

all at tau_0 = 11.11 MPa, 26% below Wang & Dunham's 15.0.

**Two independent things remain.** (a) The ramp shape, which is a near-well
storage/permeability question. (b) The slip amount, where the only remaining
knob is the friction: `a` = 0.015 and `b` = 0.012 are identical in all 75 scored
runs, giving a-b = +0.003. Since the fault never reaches failure and creeps the
whole time at v = vref*exp[(tau/sigmabar - f0)/(a-b)], a-b sets how far below f0
the fault keeps creeping before it stalls -- measured at 9.3 e-folds, i.e.
d(tau/sigmabar) = 0.0278, for 632901. That window is proportional to a-b, so
slip should be roughly LINEAR in a-b, needing a-b ~ 0.017 for 2.81 cm.

Checked, not assumed: raising the amplitude 5.7x moves the front only 5%
(380 -> 400 m), because the profile has a cliff at the disc edge -- slip falls
from 10% to 0.1% of peak between 370 and 420 m. And the exterior stays dead
(tau/sigmabar - f0 = -0.203 there, so exp(-11.9) even at a-b = 0.017). So the
amplitude and the front are separable. Untested: a and b have never been varied.

## SESSION 2026-09-06 — Stage 11: the 5 d front match does NOT survive to 18 d

**THE HEADLINE, AND IT REVERSES THE PREVIOUS TWO SECTIONS.** Every claim that a
run "matches the front" was made on a 0-5 d window. Refitting on matched
windows, the tuned runs overrun badly:

     run    front 0-5 d   front 0-18 d
  632915        1.04           1.66
  632917        1.03           1.49
  632916        0.97           1.45
  632913        0.70           1.08

lambda_obs is 0.1866 over 0-5 d and 0.1709 over 0-18 d, so this is not a
window-convention artifact -- the SIMULATED front keeps growing while the
observed one does not. The 5 d agreement was a crossing, not a match.

**THE RANKING INVERTS, and 632913 is now the interesting run.** It is 632812
(kpmax 2.5e-13, tau_0 10.36, a 0.015 / b 0.012 -- the ORIGINAL configuration,
not a tuned one) extended to 18 d:

                        5 d      17 d    growth
     observed slip     2.81cm   9.17cm     3.27x
     632913            2.90     8.72       3.01x     <- amplitude AND evolution
     632916            0.81     1.53       1.88x
     632917            0.60     0.83       1.38x
     632915            0.27     0.59       2.24x

632913 matches the observed slip magnitude at BOTH times and very nearly its
growth rate, and its 0-18 d front is 1.08. It fails on one thing only: the
wellhead, at +70.7%. The runs that match the wellhead (+7.4 to +8.8%) are 6-15x
low on slip at 17 d and overrun the front by 45-66%.

**THE a-b HYPOTHESIS IS FALSIFIED, and the control run is what proved it.**
632915 and 632917 share a-b = 0.010 with different `a`, and give 0.27 vs
0.60 cm -- a factor of 2.3. So a-b is NOT the controlling parameter:

     run       a       b     a-b   slip 5 d
  632915   0.015   0.005   0.010    0.27cm
  632901   0.015   0.012   0.003    0.49
  632917   0.022   0.012   0.010    0.60
  632916   0.022   0.005   0.017    0.81

Amplitude tracks `a` (0.015 -> 0.022 roughly doubles it) and b's effect even
changes sign with `a`. The prediction of linear-in-(a-b) growth came from
inverting the STEADY-STATE relation tau/sigmabar = f0 + (a-b)*ln(v/vref), but
this fault is in a transient, where `a` enters the regularised law directly
through 2*vref*exp(-psi/a)*sinh(tau/sigmabar/a). Wrong frame. Also worth
recording: larger `a` is LESS stiff, not more -- 632916 took 4118 steps against
632915's 10261, the opposite of the risk flagged before launch.

**THE TRADE-OFF, NOW QUANTIFIED AT THE INJECTOR.** At 5 d, 632913 carries
16.80 MPa of formation overpressure at the injector cell to produce the observed
slip, while the measured wellhead peak allows about 10.9 MPa. That is the whole
disagreement: a factor of ~1.5 in near-well overpressure. Everything else --
kpmax, tau_0, a, b -- trades along it.

`skin` was checked and DOES NOT fix it. The Peaceman well-to-formation drop in
632913 is only 2.94 MPa of the 19.74 MPa well overpressure, so even skin = -2.2
(T x11) leaves pw near 17 MPa against a measured ~10.9 peak. The excess is
FORMATION pressure, not well coupling.

**WHERE THE FACTOR OF ~2 MIGHT COME FROM: the asymmetry.** The observed lobe is
one-sided -- it spans -400 to +80 m and peaks at -150 m -- while every
simulation is symmetric about the injector. The same pressurised volume
concentrated on one side produces roughly twice the local slip, which would
close much of the amplitude gap at unchanged pressure. This is implementable
without a code change: main_LH.f90:447-451 reads per-cell `tau` and `sigma` from
the parameter file (and `a`, `b`, `dc`, `f0` besides). What is NOT available is a
spatially varying `kpmax` -- the parameter-file dispatch at :453 handles `kp` and
`phi` only, so kpmax is a scalar and a three-zone kpmax would need a code change.

**RUN STATUS.** 632910/632911/632912 still running at 10 h. 632914 stopped early
at 6.9 d -- "Slip rate below vmin", not tmax -- so Wang & Dunham's config gets
only the 3 and 5 d observed times, fewer than the 14 d domain limit allowed.
632918, the permev F control, is queued.

### Correction to the asymmetry claim, from the 17 d data

The section above says the observed lobe is one-sided and treats that as a
property of the record. It is not -- it is an EARLY-TIME feature, and the claim
was read off the 5 d column:

      t   west edge   east edge   peak at    W/E   peak
     3 d     -481m         96m     -136m    5.01   2.81cm
     5 d     -481m         96m     -136m    5.01   2.81
     7 d     -513m        144m      -96m    3.56   3.02
     9 d     -617m        240m      -96m    2.57   4.37
    11 d     -665m        441m      -96m    1.51   4.99
    13 d     -729m        441m      -96m    1.66   5.40
    15 d     -978m        609m      -96m    1.61   7.56
    17 d     -978m       1090m     -289m    0.90   9.17

W/E falls from 5.0 at 5 d to 0.90 by 17 d. So the observed slip patch becomes
essentially symmetric, and the symmetric model is the right shape for the full
record even though it is the wrong shape for the first week. That weakens the
"concentrate the slip on one side to gain a factor of 2" idea as a route to the
amplitude gap -- it would help at 5 d and hurt at 17 d.

It also removes the main objection to 632913, whose slip_poster_632913.png shows
a symmetric cone tracking the observed amplitude at all eight times.

### Figures

  stage11/stage11_front_window.png  left: front radius vs sqrt(t) to 18 d with
    the observed cloud, so the divergence is visible rather than inferred.
    Neither the simulated nor the observed front is really sqrt(t) -- both step
    with the injection cycles, which is worth knowing given that lambda from
    R = lambda*sqrt(t) is this project's headline metric. Right:
    lambda_sim/lambda_obs against fitting window, both refit per window. 632913
    sits inside +/-15% from 6 to 18 d; the others leave the band immediately
    after 5 d and settle at 1.4-1.7.
  stage11/stage11_ab.png            the a-b falsification. Left: slip at 5 d vs
    a-b, with 632915 and 632917 at the same a-b and 2.3x apart. Right: the same
    four runs on an (a, b) grid, where the pattern actually lives.
  slip_poster_style/slip_poster_632913.png   the first eight-curve version of
    the poster figure. This is the one to look at.
  per-run suites for 632913-632917, all linear.

### The top-5 set, rebuilt on 18 d runs (compare_story5_18d/, slip_poster_story5_18d)

Six of the eight Stage 11 runs have finished, so the comparison set is now on
the honest window. Everything below is a 0-18 d fit with the observed front
refit on the same window.

     run     kpmax   front 0-5d   front 0-18d   slip 5d   slip 17d   /obs17   wellhead
  632913   2.5e-13        0.70          1.08    2.90cm     8.72cm     0.95     +70.7%
  632911   2.5e-12        0.99          1.57    1.27       5.09       0.55     +19.1%
  632916   2.5e-11        0.97          1.45    0.81       1.53       0.17      +8.8%
  632917   2.5e-11        1.03          1.49    0.60       0.83       0.09      +8.3%
  632910   2.5e-11        1.06          1.70    0.49       1.49       0.16      +9.7%

Monotonic in kpmax on every column, and every run that looked matched on 0-5 d
overruns by 45-70% on 0-18 d. The one that looked WORST at 5 d (632913, 0.70) is
the only one inside the band at 18 d.

632913's slip against the observation, all eight times:

      t        3d     5d     7d     9d    11d    13d    15d    17d
   observed  2.81   2.81   3.02   4.37   4.99   5.40   7.56   9.17
   632913    2.91   2.91   4.05   4.11   4.60   5.03   8.72   8.73

Within 10-35% at every time, with the same growth shape. No other run is within
a factor of 2 at 17 d. And 632913 is the ORIGINAL configuration -- 632812's
physics at kpmax 2.5e-13, tau_0 10.36 MPa, a 0.015 / b 0.012 -- not a tuned one.
Ten stages of tuning moved away from the answer, and the 5 d window is why.

METHOD NOTE, second occurrence of the same trap. Building this set with 632914
in it collapsed the common fitting window to 6.898 d, because that is where
632914 stopped ("Slip rate below vmin"). Every lambda was then refit on a third
of the intended span. This is exactly the 632880 incident recorded earlier, and
the script does print the window -- it has to be read. 632914 is excluded from
the comparison; it also has no 17 d slip, so it does not belong in a 17 d
figure regardless.

STILL OUT: 632912 (= 632875 at 18 d) and 632918 (the permev F control).
