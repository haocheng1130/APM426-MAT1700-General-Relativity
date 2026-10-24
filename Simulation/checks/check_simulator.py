"""Check kerr_geodesics.html: run its physics code on the example orbits and compare with Hamilton's equations.

Needs Node.js for dump_trajectories.js.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import hamilton  # noqa: E402

# Orbits near an unstable orbit amplify any difference, so they get a looser tolerance or are only checked for
# their fate: the zoom-whirl orbit passes close to the separatrix, and the photon sphere is unstable by design.
TOLERANCE = {"zoom-whirl": 1e-3, "photon sphere": None}

cases = json.loads(subprocess.run(["node", str(HERE / "dump_trajectories.js")], check=True, capture_output=True,
                                  text=True).stdout)
failed = False
for case in cases:
    samples = np.array(case["samples"])
    tol = TOLERANCE.get(case["name"], 1e-5)
    line = f"{case['name']:14s} {case['status']:9s} drift {case['drift']:.0e}"
    ok = case["status"] == case["expect"] and case["drift"] < 1e-7
    if tol is not None:
        ref, _ = hamilton.positions(case["a"], case["E"], case["L"], case["y0"], samples[:, 0])
        diff = hamilton.relative_difference(samples[:, 1:4], ref).max()
        line += f"   vs Hamilton {diff:.0e} (tolerance {tol:.0e})"
        ok = ok and diff < tol
    print(("ok    " if ok else "FAIL  ") + line)
    failed |= not ok
sys.exit(1 if failed else 0)
