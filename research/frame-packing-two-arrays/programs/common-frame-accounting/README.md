# common-frame-accounting: frame geometry and what a unit would need

Run every command from this folder with Python 3.9 or later and numpy. The outputs of the original runs are in `../../data/common-frame-accounting/` with the same names. Some programs write a JSON file next to themselves or into `out/`; run them from a copy if you want to keep the shipped files.

| program | what it does | command | time |
|---|---|---|---|
| `frames.py` | frame facts from actual matrices at blocks 1 to 3: free letters fix the raw frame, a paid move is a step to a neighbour, frame distances, the route operators, and the cut ratios (point cuts and the spectral ceiling) | `python3 frames.py` | 1 s |
| `needs.py` | what a unit with saving 0.1 or 0.01 needs at each block: moves saved per array, overlap per array in a common-frame design | `python3 needs.py` | 1 s |

`ctlword.out` and `ctlword.json` are the algebra track's check of the order-3 control word. Their program used the
moonshot rounds' word engine and is not included; the same word is checked independently by
`../skeptic-checks/slack/vctl.py`.
