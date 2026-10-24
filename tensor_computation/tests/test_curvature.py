"""Curvature of metrics whose answers are known in closed form."""
import sympy as sp

from grtensor import EXAMPLES, Metric, from_line_element
from grtensor.metric import simplify_component

t, r, th, ph, M, a, Q, k = sp.symbols("t r theta phi M a Q k", positive=True)


def example(name):
    ex = EXAMPLES[name]
    return from_line_element(ex["line_element"], ex["coords"], ex["functions"], ex["definitions"])


def is_zero(expr):
    return sp.simplify(expr) == 0


def test_sphere_curvature_is_positive():
    # regression: the old notebook contracted the wrong Riemann index and got R = -2/R_0^2
    m = example("2-sphere of radius R_0")
    R0 = sp.Symbol("R_{0}", positive=True)
    assert is_zero(m.ricci_scalar - 2 / R0**2)
    assert all(is_zero(e) for e in m.ricci - m.g / R0**2)


def test_scalar_uses_off_diagonal_inverse_metric():
    # regression: the old notebook summed only g^{ii} R_{ii}; this is the unit sphere with v = phi + theta
    v = sp.Symbol("v", positive=True)
    g = sp.Matrix([[1 + sp.sin(th)**2, -sp.sin(th)**2], [-sp.sin(th)**2, sp.sin(th)**2]])
    assert is_zero(Metric([th, v], g).ricci_scalar - 2)


def test_christoffel_uses_exact_rationals():
    # regression: the old notebook used the float 1/2 and printed -0.5*sin(2*theta)
    m = example("2-sphere of radius R_0")
    gamma = m.christoffel[0][1][1]
    assert not gamma.atoms(sp.Float)
    assert is_zero(gamma + sp.sin(th) * sp.cos(th))


def test_schwarzschild():
    m = example("Schwarzschild")
    nonzero = [(i, j, l) for i in range(4) for j in range(4) for l in range(j, 4) if m.christoffel[i][j][l] != 0]
    assert len(nonzero) == 9
    assert is_zero(m.christoffel[1][0][0] - M * (r - 2 * M) / r**3)
    assert m.ricci == sp.zeros(4)
    assert m.ricci_scalar == 0
    assert is_zero(m.kretschmann - 48 * M**2 / r**6)


def test_static_spherical_einstein_tensor():
    # Schutz, A First Course in General Relativity, eqs. (10.14)-(10.15)
    m = example("Static spherically symmetric (Schutz §10.1)")
    Phi, Lam = sp.Function("Phi")(r), sp.Function("Lambda")(r)
    G_tt = sp.exp(2 * Phi) / r**2 * sp.diff(r * (1 - sp.exp(-2 * Lam)), r)
    G_rr = -sp.exp(2 * Lam) / r**2 * (1 - sp.exp(-2 * Lam)) + 2 / r * sp.diff(Phi, r)
    assert is_zero(m.einstein[0, 0] - G_tt)
    assert is_zero(m.einstein[1, 1] - G_rr)


def test_reissner_nordstrom():
    m = example("Reissner–Nordström")
    f = 1 - 2 * M / r + Q**2 / r**2
    assert m.ricci_scalar == 0
    assert is_zero(m.einstein[0, 0] - Q**2 * f / r**4)


def test_flrw_ricci_scalar():
    m = example("FLRW (curvature k)")
    A = sp.Function("a")(t)
    expected = 6 * (A * sp.diff(A, t, 2) + sp.diff(A, t)**2 + k) / A**2
    assert is_zero(m.ricci_scalar - expected)


def test_kerr_is_ricci_flat_with_known_kretschmann():
    m = example("Kerr (Boyer–Lindquist, as in the report)")
    assert m.ricci == sp.zeros(4)
    S = r**2 + a**2 * sp.cos(th)**2
    c = sp.cos(th)
    K = 48 * M**2 * (r**6 - 15 * a**2 * r**4 * c**2 + 15 * a**4 * r**2 * c**4 - a**6 * c**6) / S**6
    assert simplify_component(m.kretschmann - K) == 0
