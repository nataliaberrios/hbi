# CYCLE 2 as a new project, Round 1 — t = 0 is data-day 4.300, cycle 1 DISCARDED. muinit swept 0.370 -> 0.536 (tau_0 10.36 -> 15.00 MPa) with f0 = 0.60 FIXED at the lab value and sigmainit fixed. tau_0 is a result, not an assumption (5 runs)

**Nothing here has been submitted.** This is the pre-submittal check.

## Decks

| run | parent | **permev** | **sigmabar_0** | eta Pa·s | beta 1/Pa | phi | phi*beta | kpmax | kp=kpmin | perm field | injection | muinit | tau0 MPa | dp_crit MPa | tmax d |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [632950](params_632950.txt) | 632925 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | perm_2zone_601_ds10_disc300_kmax2.5e-13.txt (250x) | `june_clean_from_d4300.txt` | 0.3700 | 10.36 | 10.73 | 13.10 |
| [632951](params_632951.txt) | 632950 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | perm_2zone_601_ds10_disc300_kmax2.5e-13.txt (250x) | `june_clean_from_d4300.txt` | 0.4120 | 11.53 | 8.77 | 13.10 |
| [632952](params_632952.txt) | 632950 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | perm_2zone_601_ds10_disc300_kmax2.5e-13.txt (250x) | `june_clean_from_d4300.txt` | 0.4540 | 12.71 | 6.81 | 13.10 |
| [632953](params_632953.txt) | 632950 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | perm_2zone_601_ds10_disc300_kmax2.5e-13.txt (250x) | `june_clean_from_d4300.txt` | 0.4950 | 13.86 | 4.90 | 13.10 |
| [632954](params_632954.txt) | 632950 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | perm_2zone_601_ds10_disc300_kmax2.5e-13.txt (250x) | `june_clean_from_d4300.txt` | 0.5360 | 15.00 | 2.99 | 13.10 |

## Derived hydraulics

| run | D near m²/s | D far m²/s | str=beta*phi | T at kpmin | T at kpmax | gamma(h=60s, kpmin) | gamma(h=60s, kpmax) |
|---|---|---|---|---|---|---|---|
| 632950 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.268e-12 | 5.671e-10 | 0.9819 | 0.1786 |
| 632951 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.268e-12 | 5.671e-10 | 0.9819 | 0.1786 |
| 632952 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.268e-12 | 5.671e-10 | 0.9819 | 0.1786 |
| 632953 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.268e-12 | 5.671e-10 | 0.9819 | 0.1786 |
| 632954 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.268e-12 | 5.671e-10 | 0.9819 | 0.1786 |

`str`, `T` and `gamma` are the quantities `m_diffusion.f90` actually assembles (lines 444, 669, 671). `gamma` is the one a viscosity/compressibility trade does **not** hold fixed.

## Failure feasibility — can the fault slip at the OBSERVED pressure?

Peak **measured** downhole overpressure over these 13.10 d is **11.84 MPa** (p_wh + rho·g·H − P0, with rho·g·H = 40.0 and P0 = 73.8 MPa). Slip requires Δp > Δp_crit = σ̄₀(1 − μ₀/f₀). If Δp_crit exceeds 11.84 MPa the fault can only slip by **over-pressurising past the measurement**.

| run | μ₀ | Δp_crit MPa | vs measured | verdict |
|---|---|---|---|---|
| 632950 | 0.3700 | 10.73 | -1.11 | can fail |
| 632951 | 0.4120 | 8.77 | -3.07 | can fail |
| 632952 | 0.4540 | 6.81 | -5.03 | can fail |
| 632953 | 0.4950 | 4.90 | -6.94 | can fail |
| 632954 | 0.5360 | 2.99 | -8.85 | can fail |

**How understressed can the fault be and still fail at the observed pressure?** μ₀ ≥ f₀(1 − Δp_obs/σ̄₀):

| σ̄₀ MPa | minimum μ₀ | minimum τ₀ MPa | source of σ̄₀ |
|---|---|---|---|
| 27.99 | **0.3462** | **9.69** | derived from σ_v 100, σ_Hmax 160, dip 10°, p_pore 73.82 |

This stage runs μ₀ = 0.3700, 0.4120, 0.4540, 0.4950, 0.5360, comfortably above the floor above, so the fault can reach failure at the measured pressure without over-pressurising. Whether it then accumulates enough slip to build a front is a separate question that the friction law, not this inequality, decides — HBI's regularised rate-and-state law has no failure threshold, so Δp_crit is an orientation number here and not a prediction.

## Initial condition

![initial permeability](perm_ic_stage20.png)

## Deck diffs against parents

- [`res632950.in` vs `res632925.in`](diff_632950.txt)
- [`res632951.in` vs `res632950.in`](diff_632951.txt)
- [`res632952.in` vs `res632950.in`](diff_632952.txt)
- [`res632953.in` vs `res632950.in`](diff_632953.txt)
- [`res632954.in` vs `res632950.in`](diff_632954.txt)

## Launch commands, once approved

```bash
cd /home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs
sbatch march26_submit_hbi_git_scratch.sh -i res632950.in -w june_clean.txt -p perm_2zone_601_ds10_disc300_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632951.in -w june_clean.txt -p perm_2zone_601_ds10_disc300_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632952.in -w june_clean.txt -p perm_2zone_601_ds10_disc300_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632953.in -w june_clean.txt -p perm_2zone_601_ds10_disc300_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632954.in -w june_clean.txt -p perm_2zone_601_ds10_disc300_kmax2.5e-13.txt
```
