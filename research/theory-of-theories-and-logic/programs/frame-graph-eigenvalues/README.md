# Exact eigenvalues of the frame graphs from their intersection numbers, and the cut caps 2/(1 − θ1/k)

Section 3 of the [section README](../../README.md): 5/14 and 13/30 in the engine game at blocks 3 and 4, and the derived cap 28/11 for cuts in the chi game.

Outputs of the original runs: [data/frame-graph-eigenvalues/out_dpg_eig_check.txt](../../data/frame-graph-eigenvalues/out_dpg_eig_check.txt).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `dpg_eig_check.py` | exact dual-polar-graph eigenvalues and cut caps (5/14, 13/30, 28/11) | `python3 dpg_eig_check.py` | under 1 s | under 100 MB | Python 3 only | identical |

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
