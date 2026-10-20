# Approximate units: Theorem X's toy unit, and the counterexample to the factor 2

Section 10.

Outputs of the original runs: [data/exactification/](../../data/exactification/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `exact_toy.py` | Theorem X's toy composed exactly to depth 6 | `python3 exact_toy.py` | 1 s | under 100 MB | numpy | identical |
| `skeptic_exactify_count.py` | counterexample to the factor 2 | `python3 skeptic_exactify_count.py` | under 1 s | under 100 MB | numpy | re-run; output new |

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
