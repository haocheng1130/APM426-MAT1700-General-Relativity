"""Independent reference integration of Kerr geodesics.

Hamilton's equations for H = ½ g^{μν} p_μ p_ν are integrated in the affine parameter λ, with g^{-1} and ∂H/∂x^μ
obtained symbolically from the Boyer–Lindquist metric. The separated functions R(r) and Θ(θ) of the report are not
used, so agreement with the Mino-time integrations checks both the numerics and those formulas.
Units G = c = M = 1.
"""
from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


@lru_cache(maxsize=1)
def _functions():
    r, th, a = sp.symbols("r theta a", real=True)
    pt, pr, pth, pph = sp.symbols("p_t p_r p_theta p_phi", real=True)
    S, D, s2 = r**2 + a**2 * sp.cos(th)**2, r**2 + a**2 - 2 * r, sp.sin(th)**2
    g = sp.zeros(4)
    g[0, 0] = -(1 - 2 * r / S)
    g[0, 3] = g[3, 0] = -2 * a * r * s2 / S
    g[1, 1] = S / D
    g[2, 2] = S
    g[3, 3] = (r**2 + a**2 + 2 * a**2 * r * s2 / S) * s2
    gi = sp.simplify(g.inv())
    p = sp.Matrix([pt, pr, pth, pph])
    H = (p.T * gi * p)[0] / 2
    args = (r, th, pt, pr, pth, pph, a)
    return (sp.lambdify(args, list(gi * p), "numpy"),
            sp.lambdify(args, sp.diff(H, r), "numpy"),
            sp.lambdify(args, sp.diff(H, th), "numpy"))


def positions(a, E, L, y0, lam):
    """Cartesian positions (x, y, z) and coordinate time t at the affine parameters lam.

    y0 = [r, dr/dσ, θ, dθ/dσ] is the Mino-time initial state used by the simulator and the figure script;
    it corresponds to p_r = (dr/dσ)/Δ and p_θ = dθ/dσ.
    """
    xdot, dHdr, dHdth = _functions()
    r0, ur0, th0, uth0 = y0[:4]

    def rhs(_, y):
        _, r, th, _, pr, pth = y
        v = xdot(r, th, -E, pr, pth, L, a)
        return [v[0], v[1], v[2], v[3], -dHdr(r, th, -E, pr, pth, L, a), -dHdth(r, th, -E, pr, pth, L, a)]

    start = [0.0, r0, th0, 0.0, ur0 / (r0**2 + a**2 - 2 * r0), uth0]
    sol = solve_ivp(rhs, (0, lam[-1]), start, method="DOP853", rtol=1e-12, atol=1e-12, t_eval=lam)
    t, r, th, ph = sol.y[:4]
    rho = np.sqrt(r**2 + a**2)
    xyz = np.stack([rho * np.sin(th) * np.cos(ph), rho * np.sin(th) * np.sin(ph), r * np.cos(th)], axis=1)
    return xyz, t


def relative_difference(xyz, reference):
    return np.linalg.norm(xyz - reference, axis=1) / np.linalg.norm(reference, axis=1)
