# Which run validates the Gauss-Chebyshev slip solution

**Use 632892. It is the only run here that the solution applies to.**

The GC solution's two fitted parameters are both evaluated at the **enhanced**
permeability `kpmax`, not at the background `kp`:

    alpha  = kpmax/(eta*phi*beta) = 2.5e-13/(0.89e-3*0.01*2.25e-8) = 1.24844
    deltaP = qinj*eta/(4*pi*kpmax) = 4.2e-3*0.89e-3/(4*pi*2.5e-13) = 1.18984e6

against the script's hardcoded `alpha = 1.25` and `deltaP = 1.2e6`. Inverting
`deltaP` for the rate returns 4.2359e-3 against Job 911's `qinj 4.2e-3`. So the
solution is parameterised on 911 at `kpmax` throughout, and a run without
permeability enhancement has no `kpmax` for the medium to reach.

| file | run | permev | kp -> kpmax | ds | use for this? |
|---|---|---|---|---|---|
| `radial_profiles_632892.csv` | 632892 | **T** | 1e-15 -> 2.5e-13, **250x** | 5 m | **YES** |
| `job911style_632892.png` | 632892 | **T** | 250x | 5 m | **YES** |
| `job911style_632893.png` | 632893 | T | 250x | 5 m | no -- see below |
| `radial_profiles_632888.csv` | 632888 | **F** | none, kp 1e-13 fixed | 20 m | **NO** |
| `job911style_632888.png` | 632888 | **F** | none | 20 m | **NO** |
| `job911style_632889.png` | 632889 | T | 1e-13 -> 2.5e-13, only **2.5x** | 20 m | **NO** |

The 632888/632889 files are kept only as the record of a wrong turn: they were
built on the premise that a fixed-alpha solution needs constant permeability,
which is backwards. `radial_profiles_632888.csv` was at one point recommended
for the overlay. It should not be used -- its `kp` is 1e-13 with no enhancement,
so its diffusivity of 0.4994 m^2/s corresponds to nothing in the solution.

632893 is 632892 with `Sw_fwid` 7.4e-5 instead of 7.4e-7. It is the
too-compliant arm of a deliberate 100x bracket: the wellbore buffers the fluid
instead of delivering it, so the run produces no slip before ~4 d and its front
ratio DRIFTS from 0.106 at 5 d to 0.236 at 17 d rather than holding constant.
That drift is the diagnostic -- a self-similar crack has a constant ratio.

## Why the wellbore was touched at all

Job 911 ran on the pre-2026-03-20 code, where constant-rate injection for a 3D
fault was a direct source term into the injector cell (`m_diffusion.f90`, the
branch commented out by `24cf728`/`d00d1db` and marked "No longer correct!").
The current code routes injection through a Peaceman well model instead, and at
`kp = 1e-15` that model is stiff: `gamma = (Sw_fwid/h)/((Sw_fwid/h)+T) = 0.977`,
so the formation barely damps the well and the pressure ramps at essentially the
full `q/Sw_fwid`. A first attempt at 911's own `Sw_fwid` of 7.4e-9 drove the
timestep to 0.01 s and would have needed 252 million steps to finish.

Raising `Sw_fwid` makes the wellbore compliant and the ramp is `q/Sw_fwid`, so
the fix is linear and needs **no code change**. 632892 uses 7.4e-7 and completed
the full 17 d in 4.7 h.

## The comparison, at t = 17 d, alpha = 1.24844

    analytical  lambda_from_T(9.583)                        0.3170

    measured    along-strike line, threshold 1e-4 m         0.2991
                along-strike line, threshold 1e-8 m         0.3046
                azimuthal mean,    threshold 1e-4 m         0.3194
                azimuthal mean,    threshold 1e-8 m         0.3249

The prediction sits inside the measured range. Two conventions move it at the
few-percent level and neither is wrong: the threshold is one-sided, because the
theory's crack edge is where slip vanishes and any finite cutoff reports a
smaller crack; and an azimuthal mean is the right comparison in principle, since
the solution is a function of r alone, while a single line is what the reference
figure plots. Quoting one number with one percentage disagreement overstates the
precision.

What is solid is that lambda is constant in time: along-strike it holds
0.2948-0.2991 across 3-17 d, a 1.5% drift over a 5.7x span in time. The
simulated crack is genuinely self-similar in `sqrt(4*alpha*t)`, which is the
form the solution assumes, so lambda is a well-defined quantity rather than an
artefact of when it was sampled.

## The assumption underneath all of this, checked

The solution treats the medium as uniformly at `kpmax`. HBI starts everything at
`kpmin` and raises k only where slip occurs, so the two agree only if the
enhanced zone covers the crack. Measured from `kp632892.dat`, the radius where k
reaches 50% of `kpmax` tracks the slip front to within 5-10 m -- one to two
cells -- at every output time, and k falls to near `kpmin` right at the tip. The
medium really is uniform `kpmax` inside a growing disc.

`radial_profiles_632892.csv` carries `k_over_kpmax` alongside slip and dp at each
time so this is checkable directly rather than on trust.

## Columns

`radial_profiles_632892.csv`, azimuthally averaged, r <= 1500 m, SI units:

    r_m
    slip_m_t{3,5,7,9,11,13,15,17}d          m
    dp_Pa_t{...}d                           Pa
    k_over_kpmax_t{...}d                    dimensionless, 1.0 = fully enhanced
