# Checks

Independent checks of the Kerr geodesic simulator and of the figure in Section 4 of the report. Both compare against
`hamilton.py`, which integrates Hamilton's equations for H = ½ g^μν p_μ p_ν in the affine parameter, with g^-1 and
∂H/∂x^μ derived symbolically from the Boyer–Lindquist metric. It does not use the separated functions R(r) and Θ(θ),
so agreement also confirms those formulas.

Run from this folder, in an environment with the packages of `../requirements.txt`:

```bash
python check_figure.py      # the report's figure orbit; the report quotes 3e-9
python check_simulator.py   # the eight example orbits of kerr_geodesics.html (needs Node.js)
```

`check_simulator.py` runs the page's own physics code through `dump_trajectories.js`. For each example it checks
how the orbit ends (bound, captured or escaped), the drift of the first integrals (dr/dσ)² = R and (dθ/dσ)² = Θ,
and the distance to the reference solution. The zoom-whirl orbit passes close to an unstable orbit, so it gets a
looser tolerance, and the unstable photon-sphere orbit is only checked for its fate. The example list in
`dump_trajectories.js` mirrors `PRESETS` in the page; update both together.
