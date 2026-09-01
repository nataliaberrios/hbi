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
