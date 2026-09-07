# Stage 13 — f0 sweep on 632913 at the MEASURED tau_0, each run in a domain sized for its front. 632922 (f0 0.4366, the value the front criterion requires) is on HOLD: imax 1601 is past the partition limit (3 decks, 2 to launch)

**Nothing here has been submitted.** This is the pre-submittal check.

## Decks

| run | parent | **permev** | **sigmabar_0** | eta Pa·s | beta 1/Pa | phi | phi*beta | kpmax | kp=kpmin | perm field | injection | muinit | tau0 MPa | dp_crit MPa | tmax d |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [632920](params_632920.txt) | 632913 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | perm_2zone_701_ds5_kmax2.5e-13.txt (250x) | `june_clean.txt` | 0.37 | 10.36 | 9.16 | 18.00 |
| [632921](params_632921.txt) | 632913 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | perm_2zone_901_ds5_kmax2.5e-13.txt (250x) | `june_clean.txt` | 0.37 | 10.36 | 7.28 | 18.00 |
| [632922](params_632922.txt) | 632913 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | perm_2zone_1601_ds5_kmax2.5e-13.txt (250x) | `june_clean.txt` | 0.37 | 10.36 | 4.27 | 18.00 |

## Derived hydraulics

| run | D near m²/s | D far m²/s | str=beta*phi | T at kpmin | T at kpmax | gamma(h=60s, kpmin) | gamma(h=60s, kpmax) |
|---|---|---|---|---|---|---|---|
| 632920 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-10 | 0.9769 | 0.1446 |
| 632921 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-10 | 0.9769 | 0.1446 |
| 632922 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-10 | 0.9769 | 0.1446 |

`str`, `T` and `gamma` are the quantities `m_diffusion.f90` actually assembles (lines 444, 669, 671). `gamma` is the one a viscosity/compressibility trade does **not** hold fixed.

## Failure feasibility — can the fault slip at the OBSERVED pressure?

Peak **measured** downhole overpressure over these 18.00 d is **53.95 MPa** (p_wh + rho·g·H − P0, with rho·g·H = 40.0 and P0 = 73.8 MPa). Slip requires Δp > Δp_crit = σ̄₀(1 − μ₀/f₀). If Δp_crit exceeds 53.95 MPa the fault can only slip by **over-pressurising past the measurement**.

| run | μ₀ | Δp_crit MPa | vs measured | verdict |
|---|---|---|---|---|
| 632920 | 0.37 | 9.16 | -44.79 | can fail |
| 632921 | 0.37 | 7.28 | -46.67 | can fail |
| 632922 | 0.37 | 4.27 | -49.68 | can fail |

**How understressed can the fault be and still fail at the observed pressure?** μ₀ ≥ f₀(1 − Δp_obs/σ̄₀):

| σ̄₀ MPa | minimum μ₀ | minimum τ₀ MPa | source of σ̄₀ |
|---|---|---|---|
| 27.99 | **-0.5101** | **-14.28** | derived from σ_v 100, σ_Hmax 160, dip 10°, p_pore 73.82 |

This stage runs μ₀ = 0.3700, comfortably above the floor above, so the fault can reach failure at the measured pressure without over-pressurising. Whether it then accumulates enough slip to build a front is a separate question that the friction law, not this inequality, decides — HBI's regularised rate-and-state law has no failure threshold, so Δp_crit is an orientation number here and not a prediction.

## Initial condition

![initial permeability](perm_ic_stage13.png)

## Deck diffs against parents

- [`res632920.in` vs `res632913.in`](diff_632920.txt)
- [`res632921.in` vs `res632913.in`](diff_632921.txt)
- [`res632922.in` vs `res632913.in`](diff_632922.txt)

## Launch commands, once approved

```bash
cd /home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs
sbatch march26_submit_hbi_git_scratch.sh -i res632920.in -w june_clean.txt -p perm_2zone_701_ds5_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632921.in -w june_clean.txt -p perm_2zone_901_ds5_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632922.in -w june_clean.txt -p perm_2zone_1601_ds5_kmax2.5e-13.txt
```
