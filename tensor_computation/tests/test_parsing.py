"""Reading line elements written in different ways."""
import pytest
import sympy as sp

from grtensor import EXAMPLES, ParseError, parse_line_element

t, r, th, ph, M, a = sp.symbols("t r theta phi M a", positive=True)
SCHWARZSCHILD = sp.diag(-(1 - 2 * M / r), 1 / (1 - 2 * M / r), r**2, r**2 * sp.sin(th)**2)


def same(g, h):
    return all(sp.simplify(e) == 0 for e in g - h)


@pytest.mark.parametrize("text", [
    r"ds^2 = -\left(1-\frac{2M}{r}\right)dt^2 + \left(1-\frac{2M}{r}\right)^{-1}dr^2"
    r" + r^2\left(d\theta^2+\sin^2\theta\,d\phi^2\right)",
    r"-(1-2M/r)\,\mathrm{d}t^2 + \frac{\mathrm{d}r^2}{1-2M/r} + r^2\,\mathrm{d}\theta^2 + r^2\sin^{2}{\theta}\,\mathrm{d}\phi^{2}",
    r"ds^{2}=-\left(1-\frac{2M}{r}\right)dt\,dt+\frac{1}{1-\frac{2M}{r}}dr^2+r^2 d\Omega^2",
    r"-\left(1-\frac{2M}{r}\right)\dd{t}^2 + \frac{\dd{r}^2}{1-\frac{2M}{r}} + r^2\,\dd{\Omega}^2",
])
def test_latex_spellings_of_schwarzschild(text):
    x, g = parse_line_element(text, r"t, r, \theta, \phi")
    assert x == [t, r, th, ph]
    assert same(g, SCHWARZSCHILD)


def test_sympy_syntax():
    x, g = parse_line_element("-(1-2*M/r)*dt**2 + dr**2/(1-2*M/r) + r**2*dOmega**2", "t r theta phi", syntax="sympy")
    assert same(g, SCHWARZSCHILD)


def test_kerr_from_report_matches_standard_form():
    ex = EXAMPLES["Kerr (Boyer–Lindquist, as in the report)"]
    _, g = parse_line_element(ex["line_element"], ex["coords"], ex["functions"], ex["definitions"])
    S = r**2 + a**2 * sp.cos(th)**2
    kerr = sp.zeros(4)
    kerr[0, 0] = -(1 - 2 * M * r / S)
    kerr[0, 3] = kerr[3, 0] = -2 * M * a * r * sp.sin(th)**2 / S
    kerr[1, 1] = S / (r**2 - 2 * M * r + a**2)
    kerr[2, 2] = S
    kerr[3, 3] = (r**2 + a**2 + 2 * M * a**2 * r * sp.sin(th)**2 / S) * sp.sin(th)**2
    assert same(g, kerr)


def test_functions_with_and_without_arguments():
    for text in (r"-e^{2\Phi}dt^2 + dr^2", r"-e^{2\Phi(r)}dt^2 + dr^2", r"-\exp(2\Phi)dt^2 + dr^2"):
        _, g = parse_line_element(text, "t, r", functions=r"\Phi(r)")
        assert g[0, 0] == -sp.exp(2 * sp.Function("Phi")(r))


def test_off_diagonal_terms_are_split():
    _, g = parse_line_element(r"-dt^2 + 2\omega\,dt\,d\phi + d\phi^2", r"t, \phi")
    assert g[0, 1] == g[1, 0] == sp.Symbol("omega", positive=True)


def test_symbols_are_real_when_not_positive():
    x, g = parse_line_element(r"-dt^2 + dx^2", "t, x", positive=False)
    assert x[0].is_real and x[0].is_positive is None


@pytest.mark.parametrize("text, coords, message", [
    (r"-dt^2 + dr^2 + r^2 d\theta^2", r"t, r, \theta, \phi", "no differential"),
    (r"-dt^3 + dr^2", "t, r", "exactly two differentials"),
    (r"-dt^2 + dr^2 + r", "t, r", "exactly two differentials"),
    (r"-e^{2\Phi}dt^2 + dr^2", "t, r", None),
])
def test_errors_are_explained(text, coords, message):
    if message is None:
        parse_line_element(text, coords)  # an undeclared Phi is a constant symbol, not an error
        return
    with pytest.raises(ParseError, match=message):
        parse_line_element(text, coords)


def test_function_argument_must_be_a_coordinate():
    with pytest.raises(ParseError, match="not a coordinate"):
        parse_line_element(r"-e^{2\Phi}dt^2 + dr^2", "t, r", functions=r"\Phi(s)")
