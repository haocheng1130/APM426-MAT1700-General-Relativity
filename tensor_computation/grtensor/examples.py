"""Example metrics, written the way they appear in textbooks (and in the APM426 report for Kerr)."""

EXAMPLES = {
    "Kerr (Boyer–Lindquist, as in the report)": dict(
        line_element=r"ds^2 = -\frac{\Delta - a^2\sin^2\theta}{\Sigma}dt^2"
                     r" - 2a\sin^2\theta\,\frac{r^2+a^2-\Delta}{\Sigma}\,dt\,d\phi"
                     r" + \frac{\Sigma}{\Delta}dr^2 + \Sigma\,d\theta^2"
                     r" + \frac{(r^2+a^2)^2-\Delta a^2\sin^2\theta}{\Sigma}\sin^2\theta\,d\phi^2",
        coords=r"t, r, \theta, \phi",
        definitions="\\Sigma = r^2 + a^2\\cos^2\\theta\n\\Delta = r^2 + a^2 - 2Mr",
        functions="",
    ),
    "Schwarzschild": dict(
        line_element=r"ds^2 = -\left(1-\frac{2M}{r}\right)dt^2 + \frac{dr^2}{1-\frac{2M}{r}} + r^2\,d\Omega^2",
        coords=r"t, r, \theta, \phi", definitions="", functions="",
    ),
    "Static spherically symmetric (Schutz §10.1)": dict(
        line_element=r"ds^2 = -e^{2\Phi}dt^2 + e^{2\Lambda}dr^2 + r^2\,d\Omega^2",
        coords=r"t, r, \theta, \phi", definitions="", functions=r"\Phi(r), \Lambda(r)",
    ),
    "Reissner–Nordström": dict(
        line_element=r"ds^2 = -f\,dt^2 + \frac{dr^2}{f} + r^2\,d\Omega^2",
        coords=r"t, r, \theta, \phi", definitions=r"f = 1 - \frac{2M}{r} + \frac{Q^2}{r^2}", functions="",
    ),
    "FLRW (curvature k)": dict(
        line_element=r"ds^2 = -dt^2 + a(t)^2\left(\frac{dr^2}{1-kr^2} + r^2\,d\Omega^2\right)",
        coords=r"t, r, \theta, \phi", definitions="", functions=r"a(t)",
    ),
    "2-sphere of radius R_0": dict(
        line_element=r"ds^2 = R_0^2\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)",
        coords=r"\theta, \phi", definitions="", functions="",
    ),
}
