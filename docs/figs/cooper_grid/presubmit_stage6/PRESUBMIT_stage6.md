# Stage 6 — permeability enhancement ON, on the Taiyi configuration: the only base that matches the wellhead AND slips (2 runs)

**Nothing here has been submitted.** This is the pre-submittal check.

## Decks

| run | parent | **permev** | **sigmabar_0** | eta Pa·s | beta 1/Pa | phi | phi*beta | kpmax | kp=kpmin | perm field | injection | muinit | tau0 MPa | dp_crit MPa | tmax d |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [632884](params_632884.txt) | 632880 | **T** | **27.99** | 0.89e-3 | 1e-8 | 0.01 | 1.000e-10 | 1.1000e-12 | 4.0000e-13 | perm_taiyi_601_ds20.txt (3x) | `june_clean.txt` | 0.5359 | 15.00 | 2.99 | 5.00 |
| [632885](params_632885.txt) | 632881 | **T** | **27.99** | 0.89e-3 | 1e-8 | 0.01 | 1.000e-10 | 1.1000e-12 | 4.0000e-13 | perm_taiyi_601_ds20.txt (3x) | `june_clean.txt` | 0.5359 | 15.00 | 2.99 | 5.00 |

## Derived hydraulics

| run | D near m²/s | D far m²/s | str=beta*phi | T at kpmin | T at kpmax | gamma(h=60s, kpmin) | gamma(h=60s, kpmax) |
|---|---|---|---|---|---|---|---|
| 632884 | 1.2360e+01 | 4.4944e+00 | 1.000e-10 | 7.421e-10 | 2.041e-09 | 0.1425 | 0.0570 |
| 632885 | 1.2360e+01 | 4.4944e+00 | 1.000e-10 | 7.421e-10 | 2.041e-09 | 0.1425 | 0.0570 |

`str`, `T` and `gamma` are the quantities `m_diffusion.f90` actually assembles (lines 444, 669, 671). `gamma` is the one a viscosity/compressibility trade does **not** hold fixed.

## Failure feasibility — can the fault slip at the OBSERVED pressure?

Peak **measured** downhole overpressure over these 5.00 d is **10.92 MPa** (p_wh + rho·g·H − P0, with rho·g·H = 40.0 and P0 = 73.8 MPa). Slip requires Δp > Δp_crit = σ̄₀(1 − μ₀/f₀). If Δp_crit exceeds 10.92 MPa the fault can only slip by **over-pressurising past the measurement**.

| run | μ₀ | Δp_crit MPa | vs measured | verdict |
|---|---|---|---|---|
| 632884 | 0.5359 | 2.99 | -7.93 | can fail |
| 632885 | 0.5359 | 2.99 | -7.93 | can fail |

**How understressed can the fault be and still fail at the observed pressure?** μ₀ ≥ f₀(1 − Δp_obs/σ̄₀):

| σ̄₀ MPa | minimum μ₀ | minimum τ₀ MPa | source of σ̄₀ |
|---|---|---|---|
| 27.99 | **0.3659** | **10.24** | derived from σ_v 100, σ_Hmax 160, dip 10°, p_pore 73.82 |

This stage runs μ₀ = 0.5359, comfortably above the floor above, so the fault can reach failure at the measured pressure without over-pressurising. Whether it then accumulates enough slip to build a front is a separate question that the friction law, not this inequality, decides — HBI's regularised rate-and-state law has no failure threshold, so Δp_crit is an orientation number here and not a prediction.

## Initial condition

![initial permeability](perm_ic_stage6.png)

## Deck diffs against parents

- [`res632884.in` vs `res632880.in`](diff_632884.txt)
- [`res632885.in` vs `res632881.in`](diff_632885.txt)

## Launch commands, once approved

```bash
cd /home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs
sbatch march26_submit_hbi_git_scratch.sh -i res632884.in -w june_clean.txt -p perm_taiyi_601_ds20.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632885.in -w june_clean.txt -p perm_taiyi_601_ds20.txt
```
