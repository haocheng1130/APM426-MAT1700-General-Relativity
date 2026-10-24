"""Curvature of a metric given in coordinates.

Conventions (as in Schutz, Carroll and MTW):
    Γ^a_{bc}   = ½ g^{ad} (∂_b g_{dc} + ∂_c g_{db} − ∂_d g_{bc})
    R^a_{bcd}  = ∂_c Γ^a_{db} − ∂_d Γ^a_{cb} + Γ^a_{ce} Γ^e_{db} − Γ^a_{de} Γ^e_{cb}
    R_{bd}     = R^a_{bad}
    R          = g^{bd} R_{bd}
    G_{ab}     = R_{ab} − ½ g_{ab} R
    K          = R_{abcd} R^{abcd}
"""
import random
from functools import cached_property

import sympy as sp
from sympy.core.function import AppliedUndef


def _numerically_zero(expr, trials=3):
    """True if expr vanishes at a few random points.

    Used only to recognise zero components that symbolic simplification misses.
    Expressions containing undefined functions such as Phi(r) are never declared zero this way.
    """
    if expr.has(AppliedUndef) or expr.has(sp.Derivative):
        return False
    symbols = sorted(expr.free_symbols, key=str)
    rng = random.Random(426)
    for _ in range(trials):
        point = {s: sp.Rational(rng.randint(310, 1730), 1000) for s in symbols}
        try:
            value = complex(expr.evalf(40, subs=point))
        except (TypeError, ValueError, ZeroDivisionError):
            return False
        scale = max(1.0, max((abs(complex(t.evalf(20, subs=point))) for t in sp.Add.make_args(expr)), default=1.0))
        if not abs(value) < 1e-25 * scale:
            return False
    return True


_RECIPROCAL = {sp.tan: lambda x: sp.sin(x) / sp.cos(x), sp.cot: lambda x: sp.cos(x) / sp.sin(x),
               sp.sec: lambda x: 1 / sp.cos(x), sp.csc: lambda x: 1 / sp.sin(x),
               sp.tanh: lambda x: sp.sinh(x) / sp.cosh(x), sp.coth: lambda x: sp.cosh(x) / sp.sinh(x),
               sp.sech: lambda x: 1 / sp.cosh(x), sp.csch: lambda x: 1 / sp.sinh(x)}


def _lower_degree(poly_expr, var, square):
    """Replace var**k by var**(k % 2) * square**(k // 2), using var**2 = square."""
    terms = sp.Poly(poly_expr, var).as_dict(native=False)
    return sp.expand(sum((coeff * var**(k % 2) * square**(k // 2) for (k,), coeff in terms.items()), sp.Integer(0)))


def _trig_normal_form(expr, keep_cos):
    """Rational normal form modulo sin² + cos² = 1 and cosh² − sinh² = 1.

    With keep_cos, sin and sinh appear at most linearly; otherwise cos and cosh do. The result is zero
    exactly when expr is zero, because a polynomial A + sin·B in cos vanishes only if A = B = 0.
    """
    reps, rules = {}, []
    for f_s, f_c in ((sp.sin, sp.cos), (sp.sinh, sp.cosh)):
        for arg in {f.args[0] for f in expr.atoms(f_s, f_c)}:
            s, c = sp.Dummy("s"), sp.Dummy("c")
            reps[f_s(arg)], reps[f_c(arg)] = s, c
            if f_s is sp.sin:
                rules.append((s, 1 - c**2) if keep_cos else (c, 1 - s**2))
            else:
                rules.append((s, c**2 - 1) if keep_cos else (c, 1 + s**2))
    e = sp.cancel(expr.xreplace(reps))
    if rules and e != 0:
        num, den = sp.fraction(e)
        for var, square in rules:
            num, den = _lower_degree(num, var, square), _lower_degree(den, var, square)
        e = sp.cancel(num / den)
    return e.xreplace({v: k for k, v in reps.items()})


def simplify_component(expr):
    """Simplify one tensor component to the factored normal form in which sin appears at most linearly.

    The normal form decides exactly whether the component vanishes. Every component is kept in this one form
    so that common factors such as Σ = r² + a² cos²θ look the same everywhere; mixing forms makes later
    cancellations many times slower. sympy's general simplify is not used: on Kerr it can take minutes for a
    single component. See display_form for the shorter form used when printing.
    """
    expr = sp.sympify(expr)
    if expr == 0:
        return sp.Integer(0)
    for f, rewrite in _RECIPROCAL.items():
        if expr.has(f):
            expr = expr.replace(f, rewrite)
    normal = _trig_normal_form(expr, keep_cos=True)
    if normal == 0 or _numerically_zero(normal):
        return sp.Integer(0)
    return sp.factor(normal)


def display_form(expr):
    """A form for printing: numerator and denominator each take the shorter of their cos- and sin-forms.

    For Kerr this gives -2Mar sin²θ / (r² + a² cos²θ) rather than 2Mar (cosθ - 1)(cosθ + 1) / (r² + a² cos²θ).
    """
    expr = sp.sympify(expr)
    if not expr.has(sp.sin, sp.cos, sp.sinh, sp.cosh):
        return expr

    def shorter(part):
        if not part.has(sp.sin, sp.cos, sp.sinh, sp.cosh):
            return part
        return min((sp.factor(part), sp.factor(_trig_normal_form(part, keep_cos=False))), key=sp.count_ops)

    num, den = sp.fraction(sp.factor(expr))
    return shorter(num) / shorter(den)


class Metric:
    """A metric g_{ab} in coordinates x^a, with lazily computed curvature tensors."""

    def __init__(self, coords, g, simplify=simplify_component):
        self.coords = list(coords)
        self.g = sp.Matrix(g)
        self.dim = len(self.coords)
        self.simplify = simplify
        if self.g.shape != (self.dim, self.dim):
            raise ValueError(f"the metric is {self.g.shape[0]}x{self.g.shape[1]} but there are {self.dim} coordinates")
        if self.g != self.g.T:
            raise ValueError("the metric must be symmetric")
        if self.simplify(self.g.det()) == 0:
            raise ValueError("the metric is degenerate (det g = 0)")

    @cached_property
    def ginv(self):
        n, g = self.dim, self.g
        if g.is_diagonal():
            return sp.diag(*[self.simplify(1 / g[i, i]) for i in range(n)])
        inv = g.inv(method="LU")
        return sp.Matrix(n, n, lambda i, j: self.simplify(inv[i, j]))

    @cached_property
    def christoffel(self):
        """Γ[a][b][c] = Γ^a_{bc}."""
        n, x, g, gi = self.dim, self.coords, self.g, self.ginv
        dg = [[[sp.diff(g[a, b], x[c]) for c in range(n)] for b in range(n)] for a in range(n)]
        G = [[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)]
        for a in range(n):
            for b in range(n):
                for c in range(b, n):
                    s = sum((gi[a, d] * (dg[d][c][b] + dg[d][b][c] - dg[b][c][d]) for d in range(n) if gi[a, d] != 0),
                            sp.Integer(0))
                    G[a][b][c] = G[a][c][b] = self.simplify(s / 2)
        return G

    @cached_property
    def riemann(self):
        """R[a][b][c][d] = R^a_{bcd}."""
        n, x, G = self.dim, self.coords, self.christoffel
        R = [[[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(c + 1, n):
                        s = sp.diff(G[a][d][b], x[c]) - sp.diff(G[a][c][b], x[d])
                        s += sum((G[a][c][e] * G[e][d][b] - G[a][d][e] * G[e][c][b] for e in range(n)), sp.Integer(0))
                        R[a][b][c][d] = self.simplify(s)
                        R[a][b][d][c] = -R[a][b][c][d]
        return R

    @cached_property
    def riemann_lower(self):
        """R[a][b][c][d] = R_{abcd} = g_{ae} R^e_{bcd}."""
        n, g, R = self.dim, self.g, self.riemann
        out = [[[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(c + 1, n):
                        s = sum((g[a, e] * R[e][b][c][d] for e in range(n) if g[a, e] != 0), sp.Integer(0))
                        out[a][b][c][d] = self.simplify(s)
                        out[a][b][d][c] = -out[a][b][c][d]
        return out

    @cached_property
    def ricci(self):
        """Ricci tensor R_{bd} = R^a_{bad}."""
        n, R = self.dim, self.riemann
        Ric = sp.zeros(n, n)
        for b in range(n):
            for d in range(b, n):
                Ric[b, d] = Ric[d, b] = self.simplify(sum((R[a][b][a][d] for a in range(n)), sp.Integer(0)))
        return Ric

    @cached_property
    def ricci_scalar(self):
        n, gi, Ric = self.dim, self.ginv, self.ricci
        return self.simplify(sum((gi[a, b] * Ric[a, b] for a in range(n) for b in range(n) if gi[a, b] != 0), sp.Integer(0)))

    @cached_property
    def einstein(self):
        """Einstein tensor G_{ab} = R_{ab} − ½ g_{ab} R."""
        n = self.dim
        return sp.Matrix(n, n, lambda a, b: self.simplify(self.ricci[a, b] - self.g[a, b] * self.ricci_scalar / 2))

    @cached_property
    def kretschmann(self):
        """K = R_{abcd} R^{abcd} = R^{ab}_{cd} R^{cd}_{ab}.

        Each R^{ab}_{cd} is simplified before the products are summed; simplifying only the final sum is about
        a hundred times slower for Kerr.
        """
        n, gi, R = self.dim, self.ginv, self.riemann
        up = {}
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(c + 1, n):
                        v = self.simplify(sum((gi[b, f] * R[a][f][c][d] for f in range(n) if gi[b, f] != 0), sp.Integer(0)))
                        up[a, b, c, d], up[a, b, d, c] = v, -v
        return self.simplify(sum((v * up[c, d, a, b] for (a, b, c, d), v in up.items() if v != 0), sp.Integer(0)))

    @cached_property
    def cyclic_coordinates(self):
        """Coordinates the metric does not depend on; each gives a conserved momentum p_k = g_{ka} dx^a/dλ."""
        return [k for k, xk in enumerate(self.coords) if all(sp.diff(gij, xk) == 0 for gij in self.g)]
