# Lifts: the block-3 gluing bound E ≥ 68, the block-4 and block-5 interval bounds, and the constants of Theorem SG

Section 6. `check_e3.py` reads `e3_Y.json` from the current folder, so run it from `data/lift-gluing/`. `dualcert.py` re-uses the row builder of `mycover.py` (it reads that file from the current folder) and needs scipy's HiGHS solver.

Outputs of the original runs: [data/lift-gluing/](../../data/lift-gluing/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `counts3.py` | block-3 frame graph counts behind E >= 68 | `python3 counts3.py` | under 1 s | under 100 MB | Python 3 only | identical |
| `kappa.py` | the constants kappa_m of Theorem SG | `python3 kappa.py` | under 1 s | under 100 MB | Python 3 only | identical |
| `check_e3.py` | the exact interval lift with excess 3 | `cd ../../data/lift-gluing && python3 ../../programs/lift-gluing/check_e3.py` | under 1 s | under 100 MB | Python 3 only | identical |
| `dualcert.py` | exact LP dual certificates for b_4 >= 18 and b_5 >= 141 (reads mycover.py) | `python3 dualcert.py 4 1` | 7 s at block 4 | 90 MB | numpy, scipy | block 4 identical apart from float digits |

Other files here: `mycover.py`.

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
