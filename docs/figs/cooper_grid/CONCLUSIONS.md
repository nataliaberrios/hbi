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
