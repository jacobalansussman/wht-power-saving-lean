# The ceiling 169/57 for every certificate of the chi game, and the chi-game lift of excess 18

Section 8. `verify_q.py` and the skeptic's `skeptic_mycert.py` (with `geo.py`) each rebuild the 44-frame interval over F_4 and the rows of the certificate from scratch and check the inequality coefficient by coefficient. `skeptic_mycheck_lift.py` and `skeptic_mylift.py` recompute the excess of the lift in `data/chi-game-ceiling/lift_alpha_s2.txt`.

Outputs of the original runs: [data/chi-game-ceiling/](../../data/chi-game-ceiling/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `verify_q.py` | exact check of the 13-row 169/57 certificate over F_4 | `python3 verify_q.py ../../data/chi-game-ceiling/certificate_q4_m3_K5_u_1D.json` | 16 s | 240 MB | numpy | identical |
| `skeptic_mycert.py` | the skeptic's independent check of the same certificate (uses geo.py) | `python3 skeptic_mycert.py ../../data/chi-game-ceiling/certificate_q4_m3_K5_u_1D.json` | 13 s | 200 MB | numpy | identical |
| `skeptic_mycheck_lift.py` | recomputes the excess 18 of the chi-game interval lift | `python3 skeptic_mycheck_lift.py ../../data/chi-game-ceiling/lift_alpha_s2.txt` | 5 s | under 100 MB | Python 3 only | identical apart from the input path |

Other files here: `geo.py`, `skeptic_mylift.py`.

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
