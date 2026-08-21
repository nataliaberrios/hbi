# Goal: jointly match the Cooper Basin slip front AND wellhead pressure history

**Hypothesis (Natalia's):** this is achievable with an *understressed* fault, using
permeability enhancement and a nonuniform initial permeability — neither of which
Wang & Dunham (2022) used.

**The goal may be marked done only if one of these is established:**

1. **Parameter combinations are found that match both targets** while staying inside
   the physically realistic ranges below, **or**
2. **Taiyi's parameters are shown to be the only ones that work** — i.e. the
   feasible region collapses onto his stress state.

A null result is only admissible as (2) if the sweep actually covered the space. A
failure to find a match in a corner I never explored is *my* incompleteness, not
evidence for (2). The "not yet swept" list below exists to keep that honest.

## Acceptance criteria, fixed in advance

Both must hold simultaneously, on a **common time window** shared by every run
being compared, with the observed front **refit on that same window** (λ is
window-dependent: the same run gives 0.1847 over 30.7 d and 0.1730 over 16.4 d).

Tolerance is **±15%**, agreed with Natalia: Taiyi did not match the wellhead
pressure exactly either, so demanding better of this model would be holding it to
a standard the reference does not meet. ±10% is recorded separately as "strong".

| target | metric | PASS (±15%) | STRONG (±10%) |
|---|---|---|---|
| wellhead pressure | mean(sim/measured) over the run's span | **0.85 – 1.15** | 0.90 – 1.10 |
| slip front | λ / λ_observed, both refit on the common window | **0.85 – 1.15** | 0.90 – 1.10 |
| wellhead pressure | RMS against measured, as a cross-check | < 1.6 MPa | < 1.1 MPa |

The RMS row is a guard, not a gate: a run can sit at ratio 1.0 while oscillating
wildly about the data, and the mean ratio alone would not catch it.

"Understressed" means **τ₀ < 15.0 MPa**, Taiyi's value. τ₀ = muinit × sigmainit.

Peak-versus-peak is **not** an acceptable pressure metric. It flattered an early
result badly — 632522 reads 1.9× on peaks and 4.8× on the mean.

## Physically realistic ranges being swept

| parameter | range | justification |
|---|---|---|
| `sigmainit` | 26 – 30 MPa | Taiyi derives 27.99 from s_v 100, s_Hmax 160, 10° dip, p_pore 73.82 |
| `muinit` | 0.30 – 0.50 | 0.37 is the measurement-respecting value; 0.5 is Taiyi-equivalent at σ̄₀=30 |
| `f0` | 0.6 | fixed, standard for granite |
| `phi` | 0.005 – 0.02 | fractured granite; Taiyi uses 0.01 |
| `beta` | 1e-9 – 2e-8 Pa⁻¹ | water alone is 4.4e-10; bulk with pore compressibility reaches ~1e-8. Below 1e-9 is unphysically stiff |
| `kp` far field | 1e-14 – 4e-13 m² | Taiyi 4e-13 |
| `kp` near well (initial) | 4e-13 – 5e-12 m² | Taiyi 1.1e-12 over a 150 m disc |
| `kpmax` (enhanced) | 1e-12 – 1e-10 m² | fault-zone conduit; 1e-10 ≈ 100 D is the upper edge of defensible |
| `kL` | 1e-5 – 1e-3 m | slip scale for enhancement; compare dc = 1e-4 |
| `kT` | 1e12 – 1e15 s | healing time; 1e15 is effectively no healing |
| `eta` | 0.89e-3 Pa·s | fixed, water at reservoir temperature |
| `dc` | 1e-5 – 1e-4 m | Taiyi uses 1.53e-5, HBI decks 1e-4 |

## The constraint that makes this non-trivial

The measured wellhead pressure is not just a target, it **bounds the strength
margin**. Measured wellhead ≈ 45 MPa, hydrostatic ρgH = 40.0 MPa, p₀ = 73.8, so
the downhole overpressure is **≈ 10.8 MPa** — and 632507's peak pf is 10.84 MPa,
which is why its pressure matches.

Slip needs Δp > Δp_crit = σ̄₀(1 − μ/f₀), and Δp is largest at the well. So:

- **Fixed permeability:** Δp decays logarithmically, so Δp at the front is well
  below 10.8, forcing Δp_crit ≈ 4–6 MPa, i.e. μ ≈ 0.48 — essentially Taiyi's
  stress state. This is the wall.
- **With enhancement:** a slip-formed high-k channel has little pressure drop
  along it, so tip pressure stays near the wellbore value and Δp_crit can approach
  10.8. **This is the mechanism the hypothesis rests on, and it is why Taiyi's
  fixed-k model cannot access this region.**

## Not yet swept — (2) cannot be claimed until these are covered

- near-well initial permeability value, and the disc radius (fixed at 1.1e-12 / 150 m)
- far-field permeability / `kpmin` (fixed at 4e-13)
- `phi` (only `beta` has been varied; φβ is the product that matters, but they
  enter the pressure transient differently)
- `dc` (fixed 1e-4; it sets both the friction length scale and the front definition)
- `kT` healing (fixed 1e15 = none)
- `a`, `b` (fixed 0.015 / 0.012, velocity-strengthening throughout)
- well model `rw`, `Sw_fwid`, skin

## Known caveat that no parameter choice fixes

`a` > `b` everywhere means these runs are velocity-strengthening and generate **no
seismicity at all**, while the observed front is a *seismicity* front. The
comparison is an aseismic slip front against a seismicity front. Worth stating
plainly rather than letting it be discovered.
