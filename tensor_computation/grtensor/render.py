"""LaTeX for a Metric and its curvature, listing only nonzero independent components."""
import sympy as sp

from .metric import display_form


def _idx(m, *indices):
    return " ".join(sp.latex(m.coords[i]) for i in indices)


def _eq(lhs, rhs, pretty=True):
    if pretty:
        rhs = rhs.applyfunc(display_form) if isinstance(rhs, sp.MatrixBase) else display_form(rhs)
    return f"{lhs} = {sp.latex(rhs)}"


def _signed_terms(terms):
    """LaTeX terms of Σ c_k f_k in the given order, each coefficient written before its factor f_k."""
    out = []
    for c, factor in terms:
        c = display_form(c)
        if c == 0:
            continue
        negative = c.could_extract_minus_sign()
        c = -c if negative else c
        coeff = "" if c == 1 else (rf"\left({sp.latex(c)}\right)" if c.is_Add else sp.latex(c)) + r"\,"
        out.append(("- " if negative else ("+ " if out else "")) + coeff + factor)
    return out or ["0"]


def _layout(lhs, terms, tail=""):
    """lhs = terms on one line, or one term per line when that would be too wide to read."""
    if sum(len(t) for t in terms) < 250:
        return f"{lhs} = " + " ".join(terms) + tail
    return r"\begin{aligned} " + lhs + " = {} & " + r" \\ & ".join(terms) + tail + r" \end{aligned}"


def _quadratic_terms(m, symbols):
    """Terms of g_{ab} v^a v^b with v^a written as the given LaTeX symbols."""
    return _signed_terms((m.g[i, j] * (1 if i == j else 2),
                          f"{symbols[i]}^{{2}}" if i == j else rf"{symbols[i]}\,{symbols[j]}")
                         for i in range(m.dim) for j in range(i, m.dim))


def _dots(m):
    return [rf"\dot{{{sp.latex(x)}}}" for x in m.coords]


def line_element(m):
    return _layout(r"\mathrm{d}s^2", _quadratic_terms(m, [rf"\mathrm{{d}}{sp.latex(x)}" for x in m.coords]))


def _matrix_or_components(m, M, upper):
    """A matrix when it fits on screen, otherwise its nonzero independent components one per line."""
    M = M.applyfunc(display_form)
    full = _eq(r"g^{\mu\nu}" if upper else r"g_{\mu\nu}", M, pretty=False)
    if len(full) < 500:
        return [full]
    return [_eq(("g^" if upper else "g_") + f"{{{_idx(m, a, b)}}}", M[a, b], pretty=False)
            for a in range(m.dim) for b in range(a, m.dim) if M[a, b] != 0]


def metric(m):
    return _matrix_or_components(m, m.g, upper=False)


def inverse_metric(m):
    return _matrix_or_components(m, m.ginv, upper=True)


def christoffel(m):
    G = m.christoffel
    return [_eq(rf"\Gamma^{{{_idx(m, a)}}}_{{{_idx(m, b, c)}}}", G[a][b][c])
            for a in range(m.dim) for b in range(m.dim) for c in range(b, m.dim) if G[a][b][c] != 0]


def riemann(m, lower=False):
    """R^a_{bcd} with c < d, or the independent R_{abcd} with a < b, c < d and (a, b) ≤ (c, d)."""
    n = m.dim
    if not lower:
        R = m.riemann
        return [_eq(rf"R^{{{_idx(m, a)}}}{{}}_{{{_idx(m, b, c, d)}}}", R[a][b][c][d])
                for a in range(n) for b in range(n) for c in range(n) for d in range(c + 1, n) if R[a][b][c][d] != 0]
    R = m.riemann_lower
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    return [_eq(rf"R_{{{_idx(m, a, b, c, d)}}}", R[a][b][c][d])
            for i, (a, b) in enumerate(pairs) for (c, d) in pairs[i:] if R[a][b][c][d] != 0]


def symmetric_tensor(m, M, symbol):
    return [_eq(rf"{symbol}_{{{_idx(m, a, b)}}}", M[a, b]) for a in range(m.dim) for b in range(a, m.dim) if M[a, b] != 0]


def ricci(m):
    return symmetric_tensor(m, m.ricci, "R")


def einstein(m):
    return symmetric_tensor(m, m.einstein, "G")


def ricci_scalar(m):
    return [_eq("R", m.ricci_scalar)] if m.ricci_scalar != 0 else []


def kretschmann(m):
    return [_eq(r"R_{abcd}R^{abcd}", m.kretschmann)] if m.kretschmann != 0 else []


def geodesics(m):
    """Geodesic equations d²x^a/dλ² = −Γ^a_{bc} (dx^b/dλ)(dx^c/dλ)."""
    n, G, dot = m.dim, m.christoffel, _dots(m)
    return [_layout(rf"\ddot{{{sp.latex(m.coords[a])}}}",
                    _signed_terms((-(1 if b == c else 2) * G[a][b][c],
                                   f"{dot[b]}^{{2}}" if b == c else rf"{dot[b]}\,{dot[c]}")
                                  for b in range(n) for c in range(b, n)))
            for a in range(n)]


def conserved(m):
    """Momenta p_k = g_{ka} dx^a/dλ of cyclic coordinates, plus the normalization κ = g_{ab} ẋ^a ẋ^b."""
    dot = _dots(m)
    out = [_layout(rf"p_{{{_idx(m, k)}}}", _signed_terms((m.g[k, a], dot[a]) for a in range(m.dim)), r" = \text{const}")
           for k in m.cyclic_coordinates]
    out.append(_layout(r"\kappa", _quadratic_terms(m, dot)))
    return out


def document(m, sections):
    """A standalone LaTeX document; sections is a list of (title, lines)."""
    body = []
    for title, lines in sections:
        body.append(rf"\section*{{{title}}}")
        body.append(r"\begin{align*}" + "\n" + " \\\\\n".join(lines or ["0"]) + "\n" + r"\end{align*}")
    return "\n".join([r"\documentclass{article}", r"\usepackage{amsmath}", r"\usepackage[margin=1in]{geometry}",
                      r"\allowdisplaybreaks", r"\begin{document}", rf"\[{line_element(m)}\]", *body, r"\end{document}"])
