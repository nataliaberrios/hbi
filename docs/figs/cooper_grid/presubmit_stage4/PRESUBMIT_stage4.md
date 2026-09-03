# Stage 4 — spatially graded porosity, two bases x three gradings (6 runs)

**Nothing here has been submitted.** This is the pre-submittal check.

## Decks

| run | parent | **permev** | **sigmabar_0** | eta Pa·s | beta 1/Pa | phi | phi*beta | kpmax | kp=kpmin | perm field | injection | muinit | tau0 MPa | dp_crit MPa | tmax d |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [632870](params_632870.txt) | 632810 | **T** | **30.0** | 1.27e-4 | 1.5768e-07 | 0.01 | 1.577e-09 | 2.5e-13 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-13_G1.txt (20000000000000x) | `june_clean.txt` | 0.37 | 11.10 | 11.50 | 5.00 |
| [632871](params_632871.txt) | 632810 | **T** | **30.0** | 1.27e-4 | 1.5768e-07 | 0.01 | 1.577e-09 | 2.5e-13 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-13_G2.txt (10000000000000x) | `june_clean.txt` | 0.37 | 11.10 | 11.50 | 5.00 |
| [632872](params_632872.txt) | 632810 | **T** | **30.0** | 1.27e-4 | 1.5768e-07 | 0.01 | 1.577e-09 | 2.5e-13 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-13_G3.txt (20000000000000x) | `june_clean.txt` | 0.37 | 11.10 | 11.50 | 5.00 |
| [632873](params_632873.txt) | 632812 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-13_G1.txt (20000000000000x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 5.00 |
| [632874](params_632874.txt) | 632812 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-13_G2.txt (10000000000000x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 5.00 |
| [632875](params_632875.txt) | 632812 | **T** | **27.99** | 0.89e-3 | 2.25e-8 | 0.01 | 2.250e-10 | 2.5e-13 | 1e-15 | permphi_2zone_601_ds5_kmax2.5e-13_G3.txt (20000000000000x) | `june_clean.txt` | 0.37 | 10.36 | 10.73 | 5.00 |

## Derived hydraulics

| run | D near m²/s | D far m²/s | str=beta*phi | T at kpmin | T at kpmax | gamma(h=60s, kpmin) | gamma(h=60s, kpmax) |
|---|---|---|---|---|---|---|---|
| 632870 | 9.9873e+10 | 4.9937e-03 | 1.577e-09 | 2.045e-11 | 5.113e-09 | 0.8578 | 0.0236 |
| 632871 | 4.9937e+10 | 4.9937e-03 | 1.577e-09 | 2.045e-11 | 5.113e-09 | 0.8578 | 0.0236 |
| 632872 | 9.9873e+10 | 4.9937e-03 | 1.577e-09 | 2.045e-11 | 5.113e-09 | 0.8578 | 0.0236 |
| 632873 | 9.9875e+10 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-10 | 0.9769 | 0.1446 |
| 632874 | 4.9938e+10 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-10 | 0.9769 | 0.1446 |
| 632875 | 9.9875e+10 | 4.9938e-03 | 2.250e-10 | 2.918e-12 | 7.296e-10 | 0.9769 | 0.1446 |

`str`, `T` and `gamma` are the quantities `m_diffusion.f90` actually assembles (lines 444, 669, 671). `gamma` is the one a viscosity/compressibility trade does **not** hold fixed.

## Failure feasibility — can the fault slip at the OBSERVED pressure?

Peak **measured** downhole overpressure over these 5.00 d is **10.92 MPa** (p_wh + rho·g·H − P0, with rho·g·H = 40.0 and P0 = 73.8 MPa). Slip requires Δp > Δp_crit = σ̄₀(1 − μ₀/f₀). If Δp_crit exceeds 10.92 MPa the fault can only slip by **over-pressurising past the measurement**.

| run | μ₀ | Δp_crit MPa | vs measured | verdict |
|---|---|---|---|---|
| 632870 | 0.37 | 11.50 | +0.58 | **CANNOT FAIL** without over-pressurising |
| 632871 | 0.37 | 11.50 | +0.58 | **CANNOT FAIL** without over-pressurising |
| 632872 | 0.37 | 11.50 | +0.58 | **CANNOT FAIL** without over-pressurising |
| 632873 | 0.37 | 10.73 | -0.19 | can fail |
| 632874 | 0.37 | 10.73 | -0.19 | can fail |
| 632875 | 0.37 | 10.73 | -0.19 | can fail |

**How understressed can the fault be and still fail at the observed pressure?** μ₀ ≥ f₀(1 − Δp_obs/σ̄₀):

| σ̄₀ MPa | minimum μ₀ | minimum τ₀ MPa | source of σ̄₀ |
|---|---|---|---|
| 30 | **0.3816** | **11.45** | res1807/res1808's own value — a round number, **not** a measurement |
| 27.99 | **0.3659** | **10.24** | derived from σ_v 100, σ_Hmax 160, dip 10°, p_pore 73.82 |

So σ̄₀ = 27.99 admits a fault **1.21 MPa more understressed** than σ̄₀ = 30 does (10.24 vs 11.45 MPa). Both floors follow from the wellhead record and the friction law alone, with no simulation involved. The lower σ̄₀ is also the measurement-derived one, so it is both the more defensible choice and the more permissive one.

The μ₀ = 0.37 runs at σ̄₀ = 30.0 are therefore expected to slip only by over-pressurising — which is precisely the 1807/1808 behaviour being characterised (their wellhead runs +67% to +280%). Their σ̄₀ = 27.99 twins sit just below the threshold and should be able to slip at the observed pressure, by 0.19 MPa. That margin is thin enough that the twins may differ qualitatively, not just quantitatively.

## Initial condition

![initial permeability](perm_ic_stage4.png)

## Deck diffs against parents

- [`res632870.in` vs `res632810.in`](diff_632870.txt)
- [`res632871.in` vs `res632810.in`](diff_632871.txt)
- [`res632872.in` vs `res632810.in`](diff_632872.txt)
- [`res632873.in` vs `res632812.in`](diff_632873.txt)
- [`res632874.in` vs `res632812.in`](diff_632874.txt)
- [`res632875.in` vs `res632812.in`](diff_632875.txt)

## Launch commands, once approved

```bash
cd /home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs
sbatch march26_submit_hbi_git_scratch.sh -i res632870.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-13_G1.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632871.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-13_G2.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632872.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-13_G3.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632873.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-13_G1.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632874.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-13_G2.txt
sbatch march26_submit_hbi_git_scratch.sh -i res632875.in -w june_clean.txt -p permphi_2zone_601_ds5_kmax2.5e-13_G3.txt
```
