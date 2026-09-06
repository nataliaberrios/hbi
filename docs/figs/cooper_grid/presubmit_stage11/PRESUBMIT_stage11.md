# Stage 11 — 18 d extensions of the five plotted runs (632914 is 14 d, domain-limited) plus the a-b sweep, the only untouched knob left for slip amplitude (8 runs)

**Nothing here has been submitted.** This is the pre-submittal check.

## Decks

| run | parent | **permev** | **sigmabar_0** | eta Pa·s | beta 1/Pa | phi | phi*beta | kpmax | kp=kpmin | perm field | injection | muinit | tau0 MPa | dp_crit MPa | tmax d |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [632910](params_632910.txt) | 632901 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-11 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-11_G3.txt (25000x) | `june_clean.txt` | 0.3970 | 11.11 | 9.47 | 18.00 |
| [632911](params_632911.txt) | 632896 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-12 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-12_G3.txt (2500x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 18.00 |
| [632912](params_632912.txt) | 632875 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-13_G3.txt (250x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 18.00 |
| [632913](params_632913.txt) | 632812 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | perm_2zone_601_ds5_kmax2.5e-13.txt (250x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 18.00 |
| [632914](params_632914.txt) | 632881 | **F** | **27.99** | 0.89e-3 | 1e-8 | 0.01 | 1.000e-10 | — | — | perm_taiyi_601_ds20.txt (3x) | `june_clean.txt` | 0.5359 | 15.00 | 2.99 | 14.00 |
| [632915](params_632915.txt) | 632901 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-11 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-11_G3.txt (25000x) | `june_clean.txt` | 0.3970 | 11.11 | 9.47 | 18.00 |
| [632916](params_632916.txt) | 632901 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-11 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-11_G3.txt (25000x) | `june_clean.txt` | 0.3970 | 11.11 | 9.47 | 18.00 |
| [632917](params_632917.txt) | 632901 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-11 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-11_G3.txt (25000x) | `june_clean.txt` | 0.3970 | 11.11 | 9.47 | 18.00 |

## Derived hydraulics

| run | D near m²/s | D far m²/s | str=beta*phi | T at kpmin | T at kpmax | gamma(h=60s, kpmin) | gamma(h=60s, kpmax) |
|---|---|---|---|---|---|---|---|
| 632910 | 1.2484e+02 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-08 | 0.9769 | 0.0017 |
| 632911 | 1.2484e+01 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-09 | 0.9769 | 0.0166 |
| 632912 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-10 | 0.9769 | 0.1446 |
| 632913 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-10 | 0.9769 | 0.1446 |
| 632914 | 1.2360e+01 | 4.4944e+00 | 1.000e-10 | nan | nan | nan | nan |
| 632915 | 1.2484e+02 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-08 | 0.9769 | 0.0017 |
| 632916 | 1.2484e+02 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-08 | 0.9769 | 0.0017 |
| 632917 | 1.2484e+02 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-08 | 0.9769 | 0.0017 |

`str`, `T` and `gamma` are the quantities `m_diffusion.f90` actually assembles (lines 444, 669, 671). `gamma` is the one a viscosity/compressibility trade does **not** hold fixed.

## Failure feasibility — can the fault slip at the OBSERVED pressure?

Peak **measured** downhole overpressure over these 18.00 d is **53.95 MPa** (p_wh + rho·g·H − P0, with rho·g·H = 40.0 and P0 = 73.8 MPa). Slip requires Δp > Δp_crit = σ̄₀(1 − μ₀/f₀). If Δp_crit exceeds 53.95 MPa the fault can only slip by **over-pressurising past the measurement**.

| run | μ₀ | Δp_crit MPa | vs measured | verdict |
|---|---|---|---|---|
| 632910 | 0.3970 | 9.47 | -44.48 | can fail |
| 632911 | 0.37 | 10.73 | -43.22 | can fail |
| 632912 | 0.37 | 10.73 | -43.22 | can fail |
| 632913 | 0.37 | 10.73 | -43.22 | can fail |
| 632914 | 0.5359 | 2.99 | -50.96 | can fail |
| 632915 | 0.3970 | 9.47 | -44.48 | can fail |
| 632916 | 0.3970 | 9.47 | -44.48 | can fail |
| 632917 | 0.3970 | 9.47 | -44.48 | can fail |

**How understressed can the fault be and still fail at the observed pressure?** μ₀ ≥ f₀(1 − Δp_obs/σ̄₀):

| σ̄₀ MPa | minimum μ₀ | minimum τ₀ MPa | source of σ̄₀ |
|---|---|---|---|
| 27.99 | **-0.5565** | **-15.58** | derived from σ_v 100, σ_Hmax 160, dip 10°, p_pore 73.82 |

This stage runs μ₀ = 0.3700, 0.3970, 0.5359, comfortably above the floor above, so the fault can reach failure at the measured pressure without over-pressurising. Whether it then accumulates enough slip to build a front is a separate question that the friction law, not this inequality, decides — HBI's regularised rate-and-state law has no failure threshold, so Δp_crit is an orientation number here and not a prediction.

## Initial condition

![initial permeability](perm_ic_stage11.png)

## Deck diffs against parents

- [`res632910.in` vs `res632901.in`](diff_632910.txt)
- [`res632911.in` vs `res632896.in`](diff_632911.txt)
- [`res632912.in` vs `res632875.in`](diff_632912.txt)
- [`res632913.in` vs `res632812.in`](diff_632913.txt)
- [`res632914.in` vs `res632881.in`](diff_632914.txt)
- [`res632915.in` vs `res632901.in`](diff_632915.txt)
- [`res632916.in` vs `res632901.in`](diff_632916.txt)
- [`res632917.in` vs `res632901.in`](diff_632917.txt)

## Launch commands, once approved

```bash
cd /home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs
sbatch march26_submit_hbi_git_scratch.sh -i res632910.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-11_G3.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632911.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-12_G3.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632912.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-13_G3.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632913.in -w june_clean.txt -p perm_2zone_601_ds5_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632914.in -w june_clean.txt -p perm_taiyi_601_ds20.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632915.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-11_G3.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632916.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-11_G3.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632917.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-11_G3.txt
```
