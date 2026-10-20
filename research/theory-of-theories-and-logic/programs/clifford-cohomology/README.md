# Theorem E+: the hook bounds carry over to locality-preserving Clifford chain maps

Section 9. `indep.py`, `ctrl.py` and `ring.py` are the skeptic's independent checks (they import nothing from the author's `cqca.py`); `ctrl.py` is the control showing the test has power. `cqca.py` and `t1_bulk.py` are the author's library and first test.

Outputs of the original runs: [data/clifford-cohomology/](../../data/clifford-cohomology/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `indep.py` | Theorem E+: random Clifford automata, multiplicativity, defect zero | `python3 indep.py 1 2 3 4` | 1 s | under 100 MB | Python 3 only | identical |
| `ctrl.py` | control: wrong collapse orders fail | `python3 ctrl.py` | 1 s | under 100 MB | Python 3 only | identical |
| `ring.py` | odd-ring restriction depends on L mod 4 | `python3 ring.py` | under 1 s | under 100 MB | Python 3 only | identical |

Other files here: `cqca.py`, `t1_bulk.py`.

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
