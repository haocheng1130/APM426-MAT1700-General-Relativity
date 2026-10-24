"""Turn a line element ds² written in LaTeX (or SymPy syntax) into a metric matrix g_{ab}.

The line element is a quadratic form in the coordinate differentials, ds² = g_{ab} dx^a dx^b, so
g_{ab} is read off from the coefficients of dx^a dx^b (cross terms are split equally between g_{ab} and g_{ba}).
"""
import re

import sympy as sp
from sympy.parsing.latex import parse_latex
from sympy.parsing.sympy_parser import convert_xor, implicit_multiplication, parse_expr, standard_transformations

from .metric import simplify_component

_FUNCS = r"arcsin|arccos|arctan|sinh|cosh|tanh|coth|sech|csch|sin|cos|tan|cot|sec|csc|exp|ln|log"
_FUNC_RE = re.compile(r"\\(?:" + _FUNCS + r")(?![A-Za-z])")
_TOKEN = r"\\[A-Za-z]+(?:_\{[^{}]*\})?|[A-Za-z](?:_\{[^{}]*\})?"
_SYMPY_TRANSFORMS = standard_transformations + (implicit_multiplication, convert_xor)


class ParseError(ValueError):
    """The input could not be turned into a metric; the message says what to fix."""


# ---------- LaTeX preprocessing ----------

def _clean(s):
    s = re.sub(r"\\[,;:!> ]|\\q?quad(?![A-Za-z])|~", " ", s)
    s = re.sub(r"\\(?:mathrm|text|operatorname|mathit)\{d\}|\{\\rm\s*d\}|\\rm\s*d", "d", s)
    s = re.sub(r"\\dd\{([^{}]*)\}", r"d \1", s)
    s = re.sub(r"\\(?:displaystyle|textstyle|bigl|bigr|Bigl|Bigr|big|Big)", "", s)
    s = re.sub(r"_([A-Za-z0-9])", r"_{\1}", s)
    return s


def _group_end(s, k):
    """Index just after the TeX group starting at s[k]: a {...} block, a \\command, or one character."""
    if k >= len(s):
        return k
    if s[k] == "{":
        depth = 0
        for i in range(k, len(s)):
            depth += {"{": 1, "}": -1}.get(s[i], 0)
            if depth == 0:
                return i + 1
        raise ParseError("unbalanced braces { }")
    m = re.match(r"\\[A-Za-z]+", s[k:])
    return k + (m.end() if m else 1)


def _skip_spaces(s, k):
    while k < len(s) and s[k].isspace():
        k += 1
    return k


def _bracket_function_args(s):
    r"""Write every argument of \sin, \cos, ... as \left(...\right).

    Without this, the LaTeX parser reads \sin^2\theta\, d\phi^2 as sin(θ dφ²)².
    """
    out, i = [], 0
    arg_re = re.compile(r"\s*(\d*\.?\d*)\s*(" + _TOKEN + ")")
    while True:
        m = _FUNC_RE.search(s, i)
        if not m:
            out.append(s[i:])
            return "".join(out)
        out.append(s[i:m.end()])
        j = m.end()
        k = _skip_spaces(s, j)
        if k < len(s) and s[k] == "^":
            end = _group_end(s, _skip_spaces(s, k + 1))
            out.append(s[j:end])
            j = end
        k = _skip_spaces(s, j)
        if s.startswith(r"\left(", k) or s.startswith("(", k):
            i = j
        elif k < len(s) and s[k] == "{":
            end = _group_end(s, k)
            out.append(r"\left(" + _bracket_function_args(s[k + 1:end - 1]) + r"\right)")
            i = end
        elif (a := arg_re.match(s, j)) and a.group(2):
            out.append(r"\left(" + a.group(1) + a.group(2) + r"\right)")
            i = a.end()
        else:
            i = j


def _strip_function_args(s, fn_tokens):
    r"""Drop explicit arguments of declared functions, e.g. e^{2\Phi(r)} -> e^{2\Phi}; they are restored later."""
    for tok in fn_tokens:
        s = re.sub(r"(?<![A-Za-z\\])(" + re.escape(tok) + r")(?![A-Za-z])(\^\{[^{}]*\}|\^[0-9A-Za-z])?\s*(?:\\left)?\([^()]*(?:\\right)?\)",
                   r"\1\2", s)
    return s


def _latex(s, fn_tokens=()):
    s = _strip_function_args(s, fn_tokens)
    s = _bracket_function_args(s)
    try:
        return parse_latex(s, backend="antlr")
    except Exception as exc:  # the parser raises several exception types
        raise ParseError(f"could not read the LaTeX: {str(exc).splitlines()[0]}") from exc


def _latex_tokens(text):
    return re.findall(_TOKEN, _clean(text))


def _latex_symbol_name(token):
    expr = parse_latex(token, backend="antlr")
    if not isinstance(expr, sp.Symbol):
        raise ParseError(f"{token!r} is not a single symbol")
    return expr.name


# ---------- main entry point ----------

def parse_line_element(text, coords, functions="", definitions="", syntax="latex", positive=True):
    """Read a line element and return (coordinate symbols, metric matrix g_{ab}).

    text        the line element, e.g. r"ds^2 = -\\left(1-\\frac{2M}{r}\\right)dt^2 + ..."
    coords      the coordinates in order, e.g. r"t, r, \\theta, \\phi" (SymPy syntax: "t, r, theta, phi")
    functions   unknown functions of the coordinates, e.g. r"\\Phi(r), \\Lambda(r)"
    definitions abbreviations used in text, one per line, e.g. r"\\Sigma = r^2 + a^2\\cos^2\\theta"
    syntax      "latex" or "sympy"
    positive    treat every symbol as positive (helps simplification of square roots); otherwise real
    """
    if syntax not in ("latex", "sympy"):
        raise ParseError('syntax must be "latex" or "sympy"')
    if "=" in text and re.search(r"d\s*s", text.split("=")[0]):
        text = text.split("=", 1)[1]
    if syntax == "latex":
        coord_tokens = _latex_tokens(coords)
        coord_names = [_latex_symbol_name(t) for t in coord_tokens]
        fn_specs = [(tok, _latex_symbol_name(tok), _latex_tokens(args))
                    for tok, args in re.findall(r"(" + _TOKEN + r")\s*\(([^()]*)\)", _clean(functions))]
        fn_specs = [(tok, name, [_latex_symbol_name(a) for a in args]) for tok, name, args in fn_specs]
    else:
        coord_tokens = coord_names = re.findall(r"[A-Za-z_][A-Za-z0-9_]*", coords)
        fn_specs = [(name, name, re.findall(r"[A-Za-z_][A-Za-z0-9_]*", args))
                    for name, args in re.findall(r"([A-Za-z_][A-Za-z0-9_]*)\s*\(([^()]*)\)", functions)]
    n = len(coord_names)
    if n < 1:
        raise ParseError("enter the coordinates, for example: t, r, \\theta, \\phi")
    if len(set(coord_names)) != n:
        raise ParseError("each coordinate must appear once")
    for _, name, args in fn_specs:
        if unknown := [a for a in args if a not in coord_names]:
            raise ParseError(f"function {name} depends on {', '.join(unknown)}, which is not a coordinate")

    xi = [sp.Symbol(f"__xi{i}") for i in range(n)]
    if syntax == "latex":
        Q = _parse_latex_quadratic(text, coord_tokens, [tok for tok, _, _ in fn_specs])
        defs = _latex_definitions(definitions, [tok for tok, _, _ in fn_specs])
    else:
        Q = _parse_sympy_quadratic(text, coord_names)
        defs = _sympy_definitions(definitions)

    # abbreviations, then unknown functions, then assumptions on every symbol
    for _ in range(10):
        new = Q.xreplace(defs)
        if new == Q:
            break
        Q = new
    fn_map = {sp.Symbol(name): sp.Function(name)(*[sp.Symbol(a) for a in args]) for _, name, args in fn_specs}
    Q = Q.xreplace(fn_map)
    if "e" not in coord_names and "e" not in {name for _, name, _ in fn_specs}:
        Q = Q.xreplace({sp.Symbol("e"): sp.E})
    assume = {"positive": True} if positive else {"real": True}
    symbol_map = {s: sp.Symbol(s.name, **assume) for s in Q.free_symbols if s not in xi}
    symbol_map.update({sp.Symbol(c): sp.Symbol(c, **assume) for c in coord_names})
    Q = Q.xreplace(symbol_map)
    x = [sp.Symbol(c, **assume) for c in coord_names]
    return x, _quadratic_form_to_matrix(Q, xi, coord_tokens)


def _parse_latex_quadratic(text, coord_tokens, fn_tokens):
    s = _clean(text)
    if r"\Xi" in s:
        raise ParseError(r"\Xi is reserved; please rename that symbol")
    omega = re.compile(r"(?<![A-Za-z\\])d\s*\\Omega(?:_\{2\})?\s*\^\s*\{?\s*2\s*\}?")
    if omega.search(s):
        if len(coord_tokens) < 2:
            raise ParseError(r"d\Omega^2 needs two angular coordinates (the last two coordinates)")
        th, ph = coord_tokens[-2], coord_tokens[-1]
        s = omega.sub(lambda _: rf"\left(d{th}^{{2}}+\sin^{{2}}{th}\,d{ph}^{{2}}\right)", s)
        s = _clean(s)
    for i, tok in enumerate(coord_tokens):
        pattern = r"(?<![A-Za-z\\])d\s*" + re.escape(tok) + (r"(?![A-Za-z])" if tok[-1].isalpha() else "")
        s = re.sub(pattern, lambda _: rf" \Xi_{{{i}}} ", s)
    expr = _latex(s, fn_tokens)
    return expr.xreplace({sp.Symbol(f"Xi_{{{i}}}"): sp.Symbol(f"__xi{i}") for i in range(len(coord_tokens))})


def _parse_sympy_quadratic(text, coord_names):
    s = text
    if len(coord_names) >= 2:
        th, ph = coord_names[-2], coord_names[-1]
        s = re.sub(r"\bdOmega\s*(?:\*\*|\^)\s*2\b", f"(d{th}**2 + sin({th})**2*d{ph}**2)", s)
    for i, name in enumerate(coord_names):
        s = re.sub(rf"\bd{re.escape(name)}\b", f"__xi{i}", s)
    try:
        return parse_expr(s, local_dict={"E": sp.E}, transformations=_SYMPY_TRANSFORMS)
    except Exception as exc:
        raise ParseError(f"could not read the expression: {exc}") from exc


def _definition_lines(definitions):
    for line in re.split(r"[\n;]", definitions):
        if line.strip():
            if "=" not in line:
                raise ParseError(f"a definition needs '=': {line.strip()}")
            yield [part.strip() for part in line.split("=", 1)]


def _latex_definitions(definitions, fn_tokens):
    out = {}
    for lhs, rhs in _definition_lines(definitions):
        out[sp.Symbol(_latex_symbol_name(_clean(lhs)))] = _latex(_clean(rhs), fn_tokens)
    return out


def _sympy_definitions(definitions):
    out = {}
    for lhs, rhs in _definition_lines(definitions):
        out[sp.Symbol(lhs)] = parse_expr(rhs, local_dict={"E": sp.E}, transformations=_SYMPY_TRANSFORMS)
    return out


def _quadratic_form_to_matrix(Q, xi, coord_tokens):
    n = len(xi)
    if missing := [coord_tokens[i] for i in range(n) if not Q.has(xi[i])]:
        raise ParseError("no differential found for " + ", ".join(missing) +
                         "; write it as d" + missing[0] + " and check the coordinate list")
    try:
        P = sp.Poly(Q, *xi)
    except sp.PolynomialError as exc:
        raise ParseError("ds² must be a polynomial in the differentials") from exc
    if bad := [m for m in P.monoms() if sum(m) != 2]:
        raise ParseError("every term of ds² must contain exactly two differentials; found a term of degree "
                         f"{sum(bad[0])}")
    g = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            c = P.coeff_monomial(xi[i] * xi[j])
            g[i, j] = g[j, i] = simplify_component(c if i == j else c / 2)
    return g
