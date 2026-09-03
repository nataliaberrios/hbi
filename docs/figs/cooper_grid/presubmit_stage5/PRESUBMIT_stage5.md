# Taiyi reference — Wang & Dunham's published parameters verbatim, dc 1.53e-5 and 1e-4

**Nothing here has been submitted.** This is the pre-submittal check.

## Decks

| run | parent | **permev** | **sigmabar_0** | eta Pa·s | beta 1/Pa | phi | phi*beta | kpmax | kp=kpmin | perm field | injection | muinit | tau0 MPa | dp_crit MPa | tmax d |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [632880](params_632880.txt) | 632520 | **F** | **27.99** | 0.89e-3 | 1e-8 | 0.01 | 1.000e-10 | — | — | perm_taiyi_601_ds20.txt (3x) | `june_clean.txt` | 0.5359 | 15.00 | 2.99 | 5.00 |
| [632881](params_632881.txt) | 632520 | **F** | **27.99** | 0.89e-3 | 1e-8 | 0.01 | 1.000e-10 | — | — | perm_taiyi_601_ds20.txt (3x) | `june_clean.txt` | 0.5359 | 15.00 | 2.99 | 5.00 |

## Derived hydraulics

| run | D near m²/s | D far m²/s | str=beta*phi | T at kpmin | T at kpmax | gamma(h=60s, kpmin) | gamma(h=60s, kpmax) |
|---|---|---|---|---|---|---|---|
| 632880 | 1.2360e+01 | 4.4944e+00 | 1.000e-10 | nan | nan | nan | nan |
| 632881 | 1.2360e+01 | 4.4944e+00 | 1.000e-10 | nan | nan | nan | nan |

`str`, `T` and `gamma` are the quantities `m_diffusion.f90` actually assembles (lines 444, 669, 671). `gamma` is the one a viscosity/compressibility trade does **not** hold fixed.

## Failure feasibility — can the fault slip at the OBSERVED pressure?

Peak **measured** downhole overpressure over these 5.00 d is **10.92 MPa** (p_wh + rho·g·H − P0, with rho·g·H = 40.0 and P0 = 73.8 MPa). Slip requires Δp > Δp_crit = σ̄₀(1 − μ₀/f₀). If Δp_crit exceeds 10.92 MPa the fault can only slip by **over-pressurising past the measurement**.

| run | μ₀ | Δp_crit MPa | vs measured | verdict |
|---|---|---|---|---|
| 632880 | 0.5359 | 2.99 | -7.93 | can fail |
| 632881 | 0.5359 | 2.99 | -7.93 | can fail |

**How understressed can the fault be and still fail at the observed pressure?** μ₀ ≥ f₀(1 − Δp_obs/σ̄₀):

| σ̄₀ MPa | minimum μ₀ | minimum τ₀ MPa | source of σ̄₀ |
|---|---|---|---|
| 27.99 | **0.3659** | **10.24** | derived from σ_v 100, σ_Hmax 160, dip 10°, p_pore 73.82 |

This stage runs μ₀ = 0.5359, comfortably above the floor above, so the fault can reach failure at the measured pressure without over-pressurising. Whether it then accumulates enough slip to build a front is a separate question that the friction law, not this inequality, decides — HBI's regularised rate-and-state law has no failure threshold, so Δp_crit is an orientation number here and not a prediction.

## Initial condition

![initial permeability](perm_ic_stage5.png)

## Deck diffs against parents

- [`res632880.in` vs `res632520.in`](diff_632880.txt)
- [`res632881.in` vs `res632520.in`](diff_632881.txt)

## Launch commands, once approved

```bash
cd /home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs
sbatch march26_submit_hbi_git_scratch.sh -i res632880.in -w june_clean.txt -p perm_taiyi_601_ds20.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632881.in -w june_clean.txt -p perm_taiyi_601_ds20.txt
```
