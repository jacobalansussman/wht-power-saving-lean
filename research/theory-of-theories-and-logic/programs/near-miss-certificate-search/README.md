# The search for a sharper certificate at the block-4 near miss, and its skeptic

Section 7. `hookrel.py`, `mcuts.py`, `sp8lib.py` and `proj_counts.py` are the search; `vclib.py`, `vc_run.py`, `vc_ab.py` and `vc_1434.py` are the skeptic's independent library and drivers. `vc_1434.py` imports `vclib.py` from its own folder, so do not run it with `python3 -I`.

Outputs of the original runs: [data/near-miss-certificate-search/](../../data/near-miss-certificate-search/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `vc_1434.py` | orbit cut 1434 and its mirror: the counterexample to the first closing test | `python3 vc_1434.py 1434,854` | 14 s | about 3 GB | numpy | reproduced |
| `proj_counts.py` | projection counts behind the 2.849 threshold | `python3 proj_counts.py` | under 1 s | under 100 MB | Python 3 only | re-run; output new |
| `mcuts.py` | orbit cuts of the block-4 hook module | `python3 mcuts.py all OUT.json` | hours | not measured | numpy | not re-run |
| `hookrel.py` | hook-relative coordinate families | `python3 hookrel.py` | not measured | not measured | numpy | not re-run |

Other files here: `sp8lib.py`, `vc_ab.py`, `vc_run.py`, `vclib.py`.

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
