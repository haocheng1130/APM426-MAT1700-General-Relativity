"""Curvature Calculator: type a metric as a line element ds² and get its curvature.

Run with:  streamlit run app.py
"""
import time

import streamlit as st

from grtensor import EXAMPLES, ParseError, from_line_element, render

st.set_page_config(page_title="Curvature Calculator", page_icon="∇", layout="wide")

# label -> (how the result is rendered, what to print when every component vanishes)
SECTIONS = {
    "Christoffel symbols": (render.christoffel, r"\Gamma^{a}_{bc} = 0"),
    "Riemann tensor": (None, r"R^{a}{}_{bcd} = 0"),
    "Ricci tensor": (render.ricci, r"R_{ab} = 0"),
    "Ricci scalar": (render.ricci_scalar, "R = 0"),
    "Einstein tensor": (render.einstein, r"G_{ab} = 0"),
    "Geodesic equations": (render.geodesics, None),
    "Conserved quantities": (render.conserved, None),
    "Kretschmann scalar": (render.kretschmann, r"R_{abcd}R^{abcd} = 0"),
}
INPUTS = ("line_element", "coords", "functions", "definitions")


def load_example():
    example = EXAMPLES[st.session_state.example]
    for key in INPUTS:
        st.session_state[key] = example[key]
    st.session_state.syntax = "LaTeX"


@st.cache_resource(show_spinner=False, max_entries=20)
def build_metric(line_element, coords, functions, definitions, syntax, positive):
    # cached, so tensors already computed for these inputs are reused
    return from_line_element(line_element, coords, functions, definitions, syntax.lower(), positive)


if "line_element" not in st.session_state:
    st.session_state.example = "Schwarzschild"
    load_example()

with st.sidebar:
    st.selectbox("Load an example", list(EXAMPLES), key="example", on_change=load_example)
    st.radio("Input syntax", ["LaTeX", "SymPy"], key="syntax", horizontal=True,
             help="SymPy syntax: -(1-2*M/r)*dt**2 + dr**2/(1-2*M/r) + r**2*dOmega**2, coordinates t r theta phi")
    positive = st.checkbox("Treat all symbols as positive", value=True,
                           help="Lets square roots such as √(r²) simplify to r. Turn off to assume only real symbols.")
    chosen = st.multiselect("Compute", list(SECTIONS), default=list(SECTIONS))
    lower = st.radio("Riemann components", ["R^a_bcd", "R_abcd (independent)"], horizontal=True) != "R^a_bcd"

st.title("Curvature Calculator")
st.caption(r"Type a metric as a line element $\mathrm{d}s^2$ in LaTeX. Conventions as in Schutz, Carroll and MTW: "
           r"$\Gamma^a_{bc} = \tfrac12 g^{ad}(\partial_b g_{dc} + \partial_c g_{db} - \partial_d g_{bc})$, "
           r"$R^a{}_{bcd} = \partial_c\Gamma^a_{db} - \partial_d\Gamma^a_{cb} + \Gamma^a_{ce}\Gamma^e_{db} - \Gamma^a_{de}\Gamma^e_{cb}$, "
           r"$R_{bd} = R^a{}_{bad}$, $G_{ab} = R_{ab} - \tfrac12 g_{ab}R$.")

st.text_area("Line element ds²", key="line_element", height=100)
left, right = st.columns(2)
left.text_input("Coordinates, in order", key="coords", placeholder=r"t, r, \theta, \phi")
right.text_input("Unknown functions", key="functions", placeholder=r"\Phi(r), \Lambda(r)")
st.text_area("Abbreviations used above, one per line", key="definitions", height=70,
             placeholder="\\Sigma = r^2 + a^2\\cos^2\\theta")

inputs = tuple(st.session_state[k] for k in INPUTS) + (st.session_state.syntax, positive)
try:
    with st.spinner("Reading the metric…"):
        metric = build_metric(*inputs)
except ParseError as exc:
    st.error(f"Could not read the line element: {exc}")
    st.stop()
except ValueError as exc:
    st.error(f"Not a valid metric: {exc}")
    st.stop()

st.subheader("Metric as read")
st.caption("Check that this is the metric you meant before computing.")
st.latex(render.line_element(metric))
for line in render.metric(metric):
    st.latex(line)
with st.expander("Inverse metric"):
    for line in render.inverse_metric(metric):
        st.latex(line)

request = (inputs, tuple(chosen), lower)
if st.button("Compute", type="primary", disabled=not chosen):
    st.session_state.request = request

if st.session_state.get("request") == request:
    results = []
    with st.status("Computing…", expanded=False) as status:
        for label in chosen:
            status.update(label=f"Computing {label.lower()}…")
            start = time.perf_counter()
            make, _ = SECTIONS[label]
            lines = render.riemann(metric, lower=lower) if label == "Riemann tensor" else make(metric)
            results.append((label, lines, time.perf_counter() - start))
        status.update(label="Done", state="complete")

    for tab, (label, lines, seconds) in zip(st.tabs(chosen), results):
        with tab:
            zero = SECTIONS[label][1]
            if label == "Riemann tensor" and lower:
                zero = "R_{abcd} = 0"
            count = f"{len(lines)} nonzero independent component{'s' if len(lines) != 1 else ''} · " if zero else ""
            st.caption(f"{count}{seconds:.1f} s")
            if label == "Conserved quantities":
                st.caption("Each coordinate the metric does not depend on gives a Killing vector and a conserved "
                           "momentum (as E = −p_t and L = p_φ in Lemma 3.1 of the report); the last line is the "
                           "normalization, with κ = −1 for massive particles and 0 for light.")
            for line in lines or [zero]:
                st.latex(line)
            with st.expander("LaTeX source"):
                st.code("\n".join(lines or [zero]), language="latex")

    tex = render.document(metric, [(label, lines) for label, lines, _ in results])
    st.download_button("Download all results (.tex)", tex, file_name="curvature.tex", mime="text/x-tex")
