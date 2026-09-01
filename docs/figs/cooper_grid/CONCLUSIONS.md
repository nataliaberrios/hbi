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
