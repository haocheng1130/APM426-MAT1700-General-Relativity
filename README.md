# Dynamical Systems of Particles Around a Kerr Black Hole

Final project for APM426 / MAT1700 General Relativity, University of Toronto, Winter 2025, revised in October 2026.
The report derives the geodesics of a rotating black hole with the Hamilton–Jacobi method (including the Carter
constant), integrates them numerically, and discusses charged particles in a black-hole magnetosphere. Two tools
written alongside it are included.

| Path | Contents |
| --- | --- |
| [`APM426_Report.pdf`](APM426_Report.pdf) | the report |
| [`APM426_Presentation.pdf`](APM426_Presentation.pdf) | the slides |
| [`LaTeX/`](LaTeX) | LaTeX sources of the report and slides, with their figures |
| [`Simulation/`](Simulation) | interactive Kerr geodesic simulator, the script for the report's orbit figure, and checks |
| [`tensor_computation/`](tensor_computation) | curvature calculator: from a line element in LaTeX to Christoffel symbols, Riemann, Ricci and Einstein tensors |

## Kerr geodesic simulator

Download [`Simulation/kerr_geodesics.html`](Simulation/kerr_geodesics.html) and open it in a browser; nothing needs to
be installed. Choose the spin, the turning points and the inclination of a bound orbit (or the constants E, L, Q
directly, also for light rays) and watch the orbit in 3D next to its radial function R(r). It integrates the equations
of the report in Mino time, d²r/dσ² = R′(r)/2 and d²θ/dσ² = Θ′(θ)/2, with an adaptive Dormand–Prince method.

## Report figure and checks

```bash
conda create -n apm426-tensor python=3.12
conda activate apm426-tensor
pip install -r Simulation/requirements.txt
python Simulation/kerr_figure.py LaTeX/kerr_orbit.pdf    # needs a LaTeX installation for the labels
cd Simulation/checks && python check_figure.py && python check_simulator.py
```

The checks compare the figure orbit and the simulator's example orbits with an independent integration of Hamilton's
equations derived directly from the metric (see [`Simulation/checks`](Simulation/checks)).

## Curvature calculator

```bash
pip install -r tensor_computation/requirements.txt
cd tensor_computation && streamlit run app.py
```

Type a metric such as the Kerr line element of the report and get its curvature, with the results as LaTeX. See
[`tensor_computation/README.md`](tensor_computation/README.md).

## References

The papers cited in the report are listed in its bibliography and are not included here: in particular B. Carter,
Phys. Rev. 174, 1559 (1968); M. Visser, arXiv:0706.0622; M. Takahashi and H. Koyama, arXiv:0807.0277; and Y. Mino,
Phys. Rev. D 67, 084027 (2003).
