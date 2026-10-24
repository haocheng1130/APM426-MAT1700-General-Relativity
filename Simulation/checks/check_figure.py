"""Check the orbit of the report's Section 4 figure against an independent integration of Hamilton's equations.

The report quotes an agreement of 3e-9 in relative position for this orbit.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent), str(HERE)]
import hamilton  # noqa: E402
import kerr_figure as kf  # noqa: E402

a, rp, ra, inc, tmax = 0.9, 6.0, 12.0, 45.0, 2000.0  # the figure's parameters (kerr_figure.py defaults)

c = kf.solve_bound(a, rp, ra, np.cos(np.radians(inc)))
E, L = c[1], c[2]
y0 = [rp, 0.0, np.pi / 2, np.sqrt(kf.Theta(c, np.pi / 2)), 0.0, 0.0, 0.0]
reach_tmax = lambda s, y, c: y[5] - tmax  # noqa: E731
reach_tmax.terminal = True
sol = solve_ivp(kf.rhs, (0, 1e4), y0, method="DOP853", rtol=1e-11, atol=1e-12, args=(c,), events=reach_tmax,
                dense_output=True)
r, _, th, _, ph, t, lam = sol.sol(np.linspace(0, sol.t[-1], 200))
rho = np.sqrt(r**2 + a**2)
xyz = np.stack([rho * np.sin(th) * np.cos(ph), rho * np.sin(th) * np.sin(ph), r * np.cos(th)], axis=1)

ref_xyz, ref_t = hamilton.positions(a, E, L, y0, lam)
dx = hamilton.relative_difference(xyz, ref_xyz).max()
dt = (np.abs(t - ref_t) / np.maximum(1, ref_t)).max()
print(f"E = {E:.6f}, L = {L:.6f}, Q = {c[3]:.6f}")
print(f"figure orbit vs Hamilton's equations: max relative difference {dx:.1e} in position, {dt:.1e} in t")
if dx > 1e-7:
    sys.exit("FAILED: the figure orbit disagrees with Hamilton's equations")
print("OK")
