# Curvature Calculator (`grtensor`)

Type a metric as a line element ds² in LaTeX and get its Christoffel symbols, Riemann, Ricci and Einstein tensors,
Ricci and Kretschmann scalars, geodesic equations and conserved momenta, all symbolically (SymPy).
It replaces the earlier `general_relativity_tensors.ipynb`, which had wrong signs and a wrong scalar curvature
(see "Tests" below).

## Setup

```bash
conda create -n apm426-tensor python=3.12
conda activate apm426-tensor
pip install -r requirements.txt
```

## Interactive app

```bash
streamlit run app.py
```

This opens the calculator in the browser (http://localhost:8501). Load an example from the sidebar or type your own
line element; the app first shows the metric as it was read, so you can check it before pressing **Compute**.
Every result can be copied as LaTeX or downloaded as one `.tex` file.

## From Python

```python
from grtensor import from_line_element, render

kerr = from_line_element(
    r"ds^2 = -\frac{\Delta - a^2\sin^2\theta}{\Sigma}dt^2 - 2a\sin^2\theta\,\frac{r^2+a^2-\Delta}{\Sigma}\,dt\,d\phi"
    r" + \frac{\Sigma}{\Delta}dr^2 + \Sigma\,d\theta^2"
    r" + \frac{(r^2+a^2)^2-\Delta a^2\sin^2\theta}{\Sigma}\sin^2\theta\,d\phi^2",
    coords=r"t, r, \theta, \phi",
    definitions="\\Sigma = r^2 + a^2\\cos^2\\theta\n\\Delta = r^2 + a^2 - 2Mr",
)
kerr.christoffel[1][0][0]   # Γ^r_tt
kerr.ricci                  # zero matrix: Kerr is a vacuum solution
render.kretschmann(kerr)    # LaTeX lines
```

A metric can also be given as a matrix: `Metric([t, r, theta, phi], g)`. `examples.ipynb` has worked examples
(Schwarzschild, a static star with unknown Φ(r) and Λ(r), Kerr as written in the report, FLRW).

## Writing the line element

- Differentials: `dt`, `d\theta`, `\mathrm{d}r` or `\dd{r}`; products such as `dt\,d\phi`. Each coordinate needs
  its differential somewhere in ds².
- `d\Omega^2` means `d\theta^2 + \sin^2\theta\, d\phi^2` in the last two coordinates.
- Abbreviations, one per line: `\Sigma = r^2 + a^2\cos^2\theta`. They may use each other.
- Unknown functions: `\Phi(r), \Lambda(r)` or `a(t)`. In ds² they can be written with or without the argument.
- `e^{...}` is the exponential.
- SymPy syntax also works: `-(1-2*M/r)*dt**2 + dr**2/(1-2*M/r) + r**2*dOmega**2` with coordinates `t r theta phi`.
- By default every symbol is treated as positive, which lets √(r²) simplify to r. Turn this off to assume only
  real symbols.

## Conventions

As in Schutz, Carroll and MTW:

- Γ^a_bc = ½ g^ad (∂_b g_dc + ∂_c g_db − ∂_d g_bc)
- R^a_bcd = ∂_c Γ^a_db − ∂_d Γ^a_cb + Γ^a_ce Γ^e_db − Γ^a_de Γ^e_cb
- R_bd = R^a_bad, R = g^bd R_bd, G_ab = R_ab − ½ g_ab R, K = R_abcd R^abcd

## How components are simplified

Each component is brought to a rational normal form in which sin θ appears at most linearly (using
sin² = 1 − cos², and likewise cosh² − sinh² = 1), then factored. In this form a component is zero exactly when it
simplifies to zero, and every component shares the same denominators, which keeps later sums fast: the full
curvature of Kerr takes about 10 seconds. SymPy's general `simplify` is not used because it can take minutes for a
single Kerr component. For identities outside trigonometry (for example with square roots), a component that does not
cancel symbolically is also checked numerically at random points; components containing unknown functions such as
Φ(r) are never declared zero that way. When printing, numerator and denominator each use the shorter of their cos-
and sin-forms.

## Tests

```bash
pytest
```

21 tests, about 10 seconds. They check closed-form results: the 2-sphere (also in coordinates with an off-diagonal
metric), Schwarzschild (9 Christoffel symbols, Kretschmann 48M²/r⁶), the static star against Schutz eqs.
(10.14)–(10.15), Reissner–Nordström, FLRW (R = 6(aä + ȧ² + k)/a²) and Kerr (Ricci flat, Kretschmann
48M²(r⁶ − 15a²r⁴cos²θ + 15a⁴r²cos⁴θ − a⁶cos⁶θ)/Σ⁶, exactly). Three tests cover the errors of the old notebook:
it contracted the wrong pair of Riemann indices (the unit sphere came out with R = −2), summed only the diagonal
terms g^ii R_ii (wrong for non-diagonal metrics such as Kerr), and used the float `1/2`. The parsing tests read the
same metric written in several LaTeX styles and check the error messages.

## Files

| Path | Contents |
| --- | --- |
| `grtensor/parsing.py` | line element (LaTeX or SymPy syntax) → metric matrix |
| `grtensor/metric.py` | `Metric`: Christoffel, Riemann, Ricci, Einstein, Kretschmann, simplification |
| `grtensor/render.py` | LaTeX output listing the nonzero independent components |
| `grtensor/examples.py` | example metrics used by the app, notebook and tests |
| `app.py` | Streamlit interface |
| `examples.ipynb` | worked examples |
| `tests/` | known results and parsing |
