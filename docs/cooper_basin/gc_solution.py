"""Gauss-Chebyshev point-source slip solution, as a function.

Same algebra as the notebook's `case == "point"` branch: same g(), same
T(lambda) inversion by quadrature, same n-point Gauss-Chebyshev grid, same
origin condition q0 = -2, same dim_factor = R_t*f*deltaP/mu.

This exists so a figure cell does not have to inherit r_full / slip_full / lam /
R_t from whichever cell ran last. The notebook's GC cell has
`case = "constant"` set, which leaves those names holding the Eshelby
constant-stress-drop solution at R_t = 50 m and lam = None -- and a figure cell
that reads them plots that, silently, against a simulation whose crack is 810 m.
Calling gc_point() instead makes the figure depend only on its arguments.

    from gc_solution import gc_point
    S = gc_point(t_days=17)
    S["r_full"], S["slip_full"], S["lam"], S["R_t"], S["T_final"]
"""
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.special import ellipe as E, ellipk as K, exp1

# 632892's deck, and the notebook's own values for the two that are not deck
# entries. alpha and deltaP are evaluated at KPMAX, the enhanced permeability:
#   alpha  = kpmax/(eta*phi*beta) = 2.5e-13/(0.89e-3*0.01*2.25e-8) = 1.2484
#   deltaP = qinj*eta/(4*pi*kpmax) = 4.2e-3*0.89e-3/(4*pi*2.5e-13) = 1.1898e6
DEFAULTS = dict(f=0.6, sigma_0=30e6, tau_0=11.1e6, deltaP=1.2e6,
                alpha=1.25, mu=24e9)


def g(xi, lam, epsilon=0.0):
    return (exp1(lam ** 2 * xi ** 2) - exp1(lam ** 2)
            + lam ** -2 * np.exp(-lam ** 2) - epsilon)


def compute_T(lam, epsilon=0.0):
    return quad(lambda xi: g(xi, lam, epsilon) * xi / np.sqrt(1 - xi ** 2),
                0, 1)[0]


_CACHE = {}


def lambda_from_T(T_pull, npts=600):
    """Invert T(lambda). Cached, since the quadrature sweep is the slow part."""
    if npts not in _CACHE:
        lams = np.logspace(-2, 2, npts)
        Ts = np.array([compute_T(l) for l in lams])
        o = np.argsort(Ts)
        _CACHE[npts] = interp1d(np.log(Ts[o]), np.log(lams[o]),
                                bounds_error=True)
    return float(np.exp(_CACHE[npts](np.log(T_pull))))


def gc_point(t_days, n=200, **over):
    """Point-source GC solution at t_days. Keyword overrides: f, sigma_0,
    tau_0, deltaP, alpha, mu."""
    P = dict(DEFAULTS, **over)
    f, s0, t0 = P["f"], P["sigma_0"], P["tau_0"]
    dP, alpha, mu = P["deltaP"], P["alpha"], P["mu"]
    t_sec = t_days * 24 * 3600

    a_idx, b_idx = np.arange(1, n), np.arange(1, n + 1)
    rt = np.cos(np.pi * a_idx / n)
    st = np.cos(np.pi * (b_idx - 0.5) / n)
    rbar, sbar = 0.5 * (rt + 1.0), 0.5 * (st + 1.0)

    ros = (rt[:, None] + 1.0) / (st[None, :] + 1.0)
    # form m = k^2 directly; the sqrt-then-square round trip can exceed 1, where
    # K and E return nan. Clamp just below 1, where K diverges only
    # logarithmically so K(1-eps) ~ 19.9 is the right magnitude.
    m = np.minimum(4.0 * ros / (1.0 + ros) ** 2, np.nextafter(1.0, 0.0))
    A = (-1.0 / n) * (K(m) / (st[None, :] + rt[:, None] + 2.0)
                      + E(m) / (st[None, :] - rt[:, None]))

    T_final = (f * s0 - t0) / (f * dP)
    lam = lambda_from_T(T_final)
    gamma_i = g(rbar, lam, 0.0) - T_final
    R_t = lam * np.sqrt(4.0 * alpha * t_sec)
    dim_factor = R_t * f * dP / mu

    orow = np.zeros((1, n))
    orow[0, -1] = 1.0 / np.sqrt(1.0 - st[-1] ** 2)   # origin is stilde = -1
    A_full = np.vstack([A, orow])
    rhs = np.concatenate([gamma_i, [-2.0]])
    assert np.isfinite(A_full).all(), "non-finite kernel"
    assert np.isfinite(rhs).all(), "non-finite rhs"
    phi = np.linalg.solve(A_full, rhs)
    assert np.isfinite(phi).all(), "solve returned non-finite phi"

    delta_tilde = -(np.pi / n) * np.concatenate([[0.0], np.cumsum(phi)[:-1]])
    r_phys, slip = sbar * R_t, dim_factor * delta_tilde
    return dict(
        r_full=np.concatenate([-r_phys, r_phys[::-1]]),
        slip_full=np.concatenate([slip, slip[::-1]]),
        r=r_phys, slip=slip, phi=phi, delta_tilde=delta_tilde,
        lam=lam, R_t=R_t, T_final=T_final, n=n, **P)
