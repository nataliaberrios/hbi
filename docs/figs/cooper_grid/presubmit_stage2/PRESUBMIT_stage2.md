# Stage 2 — two-zone initial perm, consistent bounds, the 12 cells missing from the 16-cell slice

**Nothing here has been submitted.** This is the pre-submittal check.

## Decks

| run | parent | **permev** | **sigmabar_0** | eta Pa·s | beta 1/Pa | phi | phi*beta | kpmax | kp=kpmin | perm field | injection | muinit | tau0 MPa | dp_crit MPa | tmax d |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [632810](params_632810.txt) | 632522 | **T** | **30.0** | 1.27e-4 | 1.5768e-07 | 0.01 | 1.577e-09 | 2.5e-13 | 1e-15 | perm_2zone_601_ds5_kmax2.5e-13.txt (250x) | `june_clean.txt` | 0.37 | 11.10 | 11.50 | 5.00 |
| [632811](params_632811.txt) | 632522 | **F** | **30.0** | 1.27e-4 | 1.5768e-07 | 0.01 | 1.577e-09 | 2.5e-13 | 1e-15 | perm_2zone_601_ds5_kmax2.5e-13.txt (250x) | `june_clean.txt` | 0.37 | 11.10 | 11.50 | 5.00 |
| [632812](params_632812.txt) | 632522 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | perm_2zone_601_ds5_kmax2.5e-13.txt (250x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 5.00 |
| [632813](params_632813.txt) | 632522 | **F** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | perm_2zone_601_ds5_kmax2.5e-13.txt (250x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 5.00 |
| [632814](params_632814.txt) | 632522 | **T** | **27.99** | 1.27e-4 | 1.5768e-07 | 0.01 | 1.577e-09 | 2.5e-13 | 1e-15 | perm_2zone_601_ds5_kmax2.5e-13.txt (250x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 5.00 |
| [632815](params_632815.txt) | 632522 | **F** | **27.99** | 1.27e-4 | 1.5768e-07 | 0.01 | 1.577e-09 | 2.5e-13 | 1e-15 | perm_2zone_601_ds5_kmax2.5e-13.txt (250x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 5.00 |
| [632816](params_632816.txt) | 632524 | **T** | **30.0** | 1.27e-4 | 1.4051e-07 | 0.01 | 1.405e-09 | 5e-14 | 1e-15 | perm_2zone_601_ds5_kmax5e-14.txt (50x) | `june_clean.txt` | 0.37 | 11.10 | 11.50 | 5.00 |
| [632817](params_632817.txt) | 632524 | **F** | **30.0** | 1.27e-4 | 1.4051e-07 | 0.01 | 1.405e-09 | 5e-14 | 1e-15 | perm_2zone_601_ds5_kmax5e-14.txt (50x) | `june_clean.txt` | 0.37 | 11.10 | 11.50 | 5.00 |
| [632818](params_632818.txt) | 632524 | **T** | **27.99** | 0.89e-3 | 2.005e-8 | 0.01 | 2.005e-10 | 5e-14 | 1e-15 | perm_2zone_601_ds5_kmax5e-14.txt (50x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 5.00 |
| [632819](params_632819.txt) | 632524 | **F** | **27.99** | 0.89e-3 | 2.005e-8 | 0.01 | 2.005e-10 | 5e-14 | 1e-15 | perm_2zone_601_ds5_kmax5e-14.txt (50x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 5.00 |
| [632820](params_632820.txt) | 632524 | **T** | **27.99** | 1.27e-4 | 1.4051e-07 | 0.01 | 1.405e-09 | 5e-14 | 1e-15 | perm_2zone_601_ds5_kmax5e-14.txt (50x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 5.00 |
| [632821](params_632821.txt) | 632524 | **F** | **27.99** | 1.27e-4 | 1.4051e-07 | 0.01 | 1.405e-09 | 5e-14 | 1e-15 | perm_2zone_601_ds5_kmax5e-14.txt (50x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 5.00 |

## Derived hydraulics

| run | D near m²/s | D far m²/s | str=beta*phi | T at kpmin | T at kpmax | gamma(h=60s, kpmin) | gamma(h=60s, kpmax) |
|---|---|---|---|---|---|---|---|
| 632810 | 1.2484e+00 | 4.9937e-03 | 1.577e-09 | 2.045e-11 | 5.113e-09 | 0.8578 | 0.0236 |
| 632811 | 1.2484e+00 | 4.9937e-03 | 1.577e-09 | 2.045e-11 | 5.113e-09 | 0.8578 | 0.0236 |
| 632812 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-10 | 0.9769 | 0.1446 |
| 632813 | 1.2484e+00 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-10 | 0.9769 | 0.1446 |
| 632814 | 1.2484e+00 | 4.9937e-03 | 1.577e-09 | 2.045e-11 | 5.113e-09 | 0.8578 | 0.0236 |
| 632815 | 1.2484e+00 | 4.9937e-03 | 1.577e-09 | 2.045e-11 | 5.113e-09 | 0.8578 | 0.0236 |
| 632816 | 2.8019e-01 | 5.6039e-03 | 1.405e-09 | 2.045e-11 | 1.023e-09 | 0.8578 | 0.1076 |
| 632817 | 2.8019e-01 | 5.6039e-03 | 1.405e-09 | 2.045e-11 | 1.023e-09 | 0.8578 | 0.1076 |
| 632818 | 2.8020e-01 | 5.6040e-03 | 2.005e-10 | 2.918e-12 | 1.459e-10 | 0.9769 | 0.4581 |
| 632819 | 2.8020e-01 | 5.6040e-03 | 2.005e-10 | 2.918e-12 | 1.459e-10 | 0.9769 | 0.4581 |
| 632820 | 2.8019e-01 | 5.6039e-03 | 1.405e-09 | 2.045e-11 | 1.023e-09 | 0.8578 | 0.1076 |
| 632821 | 2.8019e-01 | 5.6039e-03 | 1.405e-09 | 2.045e-11 | 1.023e-09 | 0.8578 | 0.1076 |

`str`, `T` and `gamma` are the quantities `m_diffusion.f90` actually assembles (lines 444, 669, 671). `gamma` is the one a viscosity/compressibility trade does **not** hold fixed.

## Failure feasibility — can the fault slip at the OBSERVED pressure?

Peak **measured** downhole overpressure over these 5.00 d is **10.92 MPa** (p_wh + rho·g·H − P0, with rho·g·H = 40.0 and P0 = 73.8 MPa). Slip requires Δp > Δp_crit = σ̄₀(1 − μ₀/f₀). If Δp_crit exceeds 10.92 MPa the fault can only slip by **over-pressurising past the measurement**.

| run | μ₀ | Δp_crit MPa | vs measured | verdict |
|---|---|---|---|---|
| 632810 | 0.37 | 11.50 | +0.58 | **CANNOT FAIL** without over-pressurising |
| 632811 | 0.37 | 11.50 | +0.58 | **CANNOT FAIL** without over-pressurising |
| 632812 | 0.37 | 10.73 | -0.19 | can fail |
| 632813 | 0.37 | 10.73 | -0.19 | can fail |
| 632814 | 0.37 | 10.73 | -0.19 | can fail |
| 632815 | 0.37 | 10.73 | -0.19 | can fail |
| 632816 | 0.37 | 11.50 | +0.58 | **CANNOT FAIL** without over-pressurising |
| 632817 | 0.37 | 11.50 | +0.58 | **CANNOT FAIL** without over-pressurising |
| 632818 | 0.37 | 10.73 | -0.19 | can fail |
| 632819 | 0.37 | 10.73 | -0.19 | can fail |
| 632820 | 0.37 | 10.73 | -0.19 | can fail |
| 632821 | 0.37 | 10.73 | -0.19 | can fail |

**How understressed can the fault be and still fail at the observed pressure?** μ₀ ≥ f₀(1 − Δp_obs/σ̄₀):

| σ̄₀ MPa | minimum μ₀ | minimum τ₀ MPa | source of σ̄₀ |
|---|---|---|---|
| 30 | **0.3816** | **11.45** | res1807/res1808's own value — a round number, **not** a measurement |
| 27.99 | **0.3659** | **10.24** | derived from σ_v 100, σ_Hmax 160, dip 10°, p_pore 73.82 |

So σ̄₀ = 27.99 admits a fault **1.21 MPa more understressed** than σ̄₀ = 30 does (10.24 vs 11.45 MPa). Both floors follow from the wellhead record and the friction law alone, with no simulation involved. The lower σ̄₀ is also the measurement-derived one, so it is both the more defensible choice and the more permissive one.

The μ₀ = 0.37 runs at σ̄₀ = 30.0 are therefore expected to slip only by over-pressurising — which is precisely the 1807/1808 behaviour being characterised (their wellhead runs +67% to +280%). Their σ̄₀ = 27.99 twins sit just below the threshold and should be able to slip at the observed pressure, by 0.19 MPa. That margin is thin enough that the twins may differ qualitatively, not just quantitatively.

## Initial condition

![initial permeability](perm_ic_stage2.png)

## Deck diffs against parents

- [`res632810.in` vs `res632522.in`](diff_632810.txt)
- [`res632811.in` vs `res632522.in`](diff_632811.txt)
- [`res632812.in` vs `res632522.in`](diff_632812.txt)
- [`res632813.in` vs `res632522.in`](diff_632813.txt)
- [`res632814.in` vs `res632522.in`](diff_632814.txt)
- [`res632815.in` vs `res632522.in`](diff_632815.txt)
- [`res632816.in` vs `res632524.in`](diff_632816.txt)
- [`res632817.in` vs `res632524.in`](diff_632817.txt)
- [`res632818.in` vs `res632524.in`](diff_632818.txt)
- [`res632819.in` vs `res632524.in`](diff_632819.txt)
- [`res632820.in` vs `res632524.in`](diff_632820.txt)
- [`res632821.in` vs `res632524.in`](diff_632821.txt)

## Launch commands, once approved

```bash
cd /home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs
sbatch march26_submit_hbi_git_scratch.sh -i res632810.in -w june_clean.txt -p perm_2zone_601_ds5_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632811.in -w june_clean.txt -p perm_2zone_601_ds5_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632812.in -w june_clean.txt -p perm_2zone_601_ds5_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632813.in -w june_clean.txt -p perm_2zone_601_ds5_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632814.in -w june_clean.txt -p perm_2zone_601_ds5_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632815.in -w june_clean.txt -p perm_2zone_601_ds5_kmax2.5e-13.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632816.in -w june_clean.txt -p perm_2zone_601_ds5_kmax5e-14.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632817.in -w june_clean.txt -p perm_2zone_601_ds5_kmax5e-14.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632818.in -w june_clean.txt -p perm_2zone_601_ds5_kmax5e-14.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632819.in -w june_clean.txt -p perm_2zone_601_ds5_kmax5e-14.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632820.in -w june_clean.txt -p perm_2zone_601_ds5_kmax5e-14.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632821.in -w june_clean.txt -p perm_2zone_601_ds5_kmax5e-14.txt
```
