"""Figure for Section 4 of the APM426 report: a bound geodesic around a Kerr black hole.

Integrates, in Mino time sigma (d lambda = Sigma d sigma),
    d^2 r / d sigma^2 = R'(r) / 2,   d^2 theta / d sigma^2 = Theta'(theta) / 2,
together with d phi / d sigma, d t / d sigma and d lambda / d sigma from Theorem 3.4.
These are the same equations as in the interactive explorer (kerr_geodesics.html).
Units G = c = M = 1.

Usage: python kerr_figure.py OUTPUT.pdf [--a 0.9 --rp 6 --ra 12 --inc 45 --tmax 2000]
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


def R(c, r):
    a, E, L, Q, k = c
    D = r * r + a * a - 2 * r
    P = E * (r * r + a * a) - a * L
    return D * (k * r * r - (L - a * E) ** 2 - Q) + P * P


def dR(c, r):
    a, E, L, Q, k = c
    D = r * r + a * a - 2 * r
    P = E * (r * r + a * a) - a * L
    return (2 * r - 2) * (k * r * r - (L - a * E) ** 2 - Q) + 2 * k * r * D + 4 * E * r * P


def Theta(c, th):
    a, E, L, Q, k = c
    s, co = np.sin(th), np.cos(th)
    return Q + co * co * ((k + E * E) * a * a - L * L / (s * s))


def dTheta(c, th):
    a, E, L, Q, k = c
    s, co = np.sin(th), np.cos(th)
    return -2 * s * co * (k + E * E) * a * a + 2 * L * L * co / s ** 3


def rhs(sigma, y, c):
    a, E, L, Q, k = c
    r, ur, th, uth = y[:4]
    s2 = np.sin(th) ** 2
    D = r * r + a * a - 2 * r
    P = E * (r * r + a * a) - a * L
    return [
        ur,
        0.5 * dR(c, r),
        uth,
        0.5 * dTheta(c, th),
        -(a * E - L / s2) + a * P / D,
        -a * (a * E * s2 - L) + (r * r + a * a) * P / D,
        r * r + a * a * np.cos(th) ** 2,
    ]


def solve_bound(a, rp, ra, x):
    """E, L, Q of the bound orbit with turning points rp < ra and x = cos(iota) = L / sqrt(L^2 + Q).

    Unknowns E and C = L^2 + Q; the solution is followed from a wide orbit down to the target.
    """
    e = (ra - rp) / (ra + rp)
    p = 2 * ra * rp / (ra + rp)

    def consts(E, C, k=-1):
        return (a, E, x * np.sqrt(C), (1 - x * x) * C, k)

    p0 = max(p, 60.0)
    guess = [np.sqrt(((p0 - 2) ** 2 - 4 * e * e) / (p0 * (p0 - 3 - e * e))), p0 * p0 / (p0 - 3 - e * e)]
    for pp in p0 * (p / p0) ** np.linspace(0, 1, 41)[1:]:
        r1, r2 = pp / (1 + e), pp / (1 - e)
        guess = fsolve(lambda v: [R(consts(*v), r1) / r1 ** 4, R(consts(*v), r2) / r2 ** 4], guess, xtol=1e-12)
    return consts(*guess)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output")
    ap.add_argument("--a", type=float, default=0.9)
    ap.add_argument("--rp", type=float, default=6.0)
    ap.add_argument("--ra", type=float, default=12.0)
    ap.add_argument("--inc", type=float, default=45.0, help="inclination in degrees")
    ap.add_argument("--tmax", type=float, default=2000.0, help="duration in coordinate time t/M")
    args = ap.parse_args()

    a = args.a
    c = solve_bound(a, args.rp, args.ra, np.cos(np.radians(args.inc)))
    E, L, Q = c[1:4]
    y0 = [args.rp, 0.0, np.pi / 2, np.sqrt(Theta(c, np.pi / 2)), 0.0, 0.0, 0.0]
    reach_tmax = lambda s, y, c: y[5] - args.tmax
    reach_tmax.terminal = True
    sol = solve_ivp(rhs, (0, 1e4), y0, method="DOP853", rtol=1e-11, atol=1e-12,
                    args=(c,), events=reach_tmax, dense_output=True)
    sig = np.linspace(0, sol.t[-1], 40000)
    r, ur, th, uth, ph, t, lam = sol.sol(sig)
    drift = max(np.max(np.abs(ur ** 2 - R(c, r)) / (r * r + a * a) ** 2),
                np.max(np.abs(uth ** 2 - Theta(c, th)) / (1 + Q + L * L + a * a)))
    rh = 1 + np.sqrt(1 - a * a)
    print(f"a={a} rp={args.rp} ra={args.ra} inc={args.inc} deg: E={E:.6f} L={L:.6f} Q={Q:.6f}")
    print(f"r in [{r.min():.4f}, {r.max():.4f}], theta in [{np.degrees(th.min()):.2f}, {np.degrees(th.max()):.2f}] deg, "
          f"t_end={t[-1]:.1f}, proper time={lam[-1]:.1f}, azimuthal turns={ph[-1] / (2 * np.pi):.2f}, "
          f"sigma_end={sig[-1]:.3f}, drift={drift:.1e}")

    rho = np.sqrt(r * r + a * a)
    X, Y, Z = rho * np.sin(th) * np.cos(ph), rho * np.sin(th) * np.sin(ph), r * np.cos(th)

    plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 10})
    orbit, hole, ergo = "#1f4e8c", "black", "#2a8c82"
    fig = plt.figure(figsize=(6.4, 6.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.55], hspace=0.38, wspace=0.28)
    ax1, ax2, ax3 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
    lim = 1.1 * max(np.abs(X).max(), np.abs(Y).max(), np.abs(Z).max())
    u = np.linspace(0, 2 * np.pi, 400)
    th_s = np.linspace(0, np.pi, 200)
    r_ergo = 1 + np.sqrt(1 - a * a * np.cos(th_s) ** 2)

    ax1.plot(X, Y, color=orbit, lw=0.45)
    ax1.fill(np.sqrt(rh ** 2 + a * a) * np.cos(u), np.sqrt(rh ** 2 + a * a) * np.sin(u), color=hole)
    ax1.plot(np.sqrt(4 + a * a) * np.cos(u), np.sqrt(4 + a * a) * np.sin(u), color=ergo, lw=0.8, ls="--")
    ax1.set(xlabel=r"$x/M$", ylabel=r"$y/M$", title=r"(a) equatorial plane, $x$--$y$")

    ax2.plot(X, Z, color=orbit, lw=0.45)
    ax2.fill(np.sqrt(rh ** 2 + a * a) * np.cos(u), rh * np.sin(u), color=hole)
    xe = np.sqrt(r_ergo ** 2 + a * a) * np.sin(th_s)
    ax2.plot(np.concatenate([xe, -xe[::-1]]), np.concatenate([r_ergo * np.cos(th_s), r_ergo[::-1] * np.cos(th_s[::-1])]),
             color=ergo, lw=0.8, ls="--")
    ax2.set(xlabel=r"$x/M$", ylabel=r"$z/M$", title=r"(b) meridional plane, $x$--$z$")
    for ax in (ax1, ax2):
        ax.plot(X[0], [Y, Z][ax is ax2][0], "o", color=orbit, ms=3.5)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")

    rr = np.linspace(rh, 1.25 * args.ra, 600)
    Rr = R(c, rr)
    top = 1.15 * Rr.max()
    ax3.fill_between(rr, -top, top, where=Rr >= 0, color="#f2c46a", alpha=0.35, lw=0, label=r"$R(r) \geq 0$")
    ax3.axhline(0, color="0.4", lw=0.6)
    ax3.plot(rr, Rr, color=orbit, lw=1.2)
    ax3.axvline(rh, color="0.3", lw=0.8, ls=":")
    ax3.text(rh + 0.15, 0.8 * top, r"$r_+$")
    for rt, name in ((args.rp, r"$r_{\mathrm p}$"), (args.ra, r"$r_{\mathrm a}$")):
        ax3.plot(rt, 0, "o", color=orbit, ms=3.5)
        ax3.annotate(name, (rt, 0), xytext=(0, 7), textcoords="offset points", ha="center")
    ax3.set(xlim=(rh, rr[-1]), ylim=(-top, top), xlabel=r"$r/M$", ylabel=r"$R(r)$",
            title=r"(c) radial function $R(r)$ with turning points $R(r_{\mathrm p}) = R(r_{\mathrm a}) = 0$")
    fig.savefig(args.output, bbox_inches="tight")


if __name__ == "__main__":
    main()
