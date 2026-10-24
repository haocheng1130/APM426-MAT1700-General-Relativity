"""Symbolic curvature of spacetime metrics: from a line element ds² to Christoffel symbols, Riemann, Ricci and Einstein."""
from .examples import EXAMPLES
from .metric import Metric
from .parsing import ParseError, parse_line_element


def from_line_element(text, coords, functions="", definitions="", syntax="latex", positive=True):
    """Build a Metric directly from a line element; see parse_line_element for the arguments."""
    x, g = parse_line_element(text, coords, functions, definitions, syntax, positive)
    return Metric(x, g)


__all__ = ["EXAMPLES", "Metric", "ParseError", "from_line_element", "parse_line_element"]
