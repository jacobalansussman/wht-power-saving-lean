# ZX rewriting: the sum-free tests through PyZX, and the PyZX 0.10.7 bug

Section 11. `zx1_sumfree.py` to `zx4_xw_sumfree.py` ran on a Linux machine of Google's Colab service with PyZX 0.10.7; their outputs are in `data/zx-sum-free/`. `skeptic_pyzx_bug_repro.py` imports `zx1_sumfree.py` from the current folder.

Outputs of the original runs: [data/zx-sum-free/](../../data/zx-sum-free/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `skeptic_pyzx_bug_repro.py` | the PyZX 0.10.7 basic_optimization bug | `python3 skeptic_pyzx_bug_repro.py` | 3 s | under 100 MB | numpy, pyzx==0.10.7 | reproduced |
| `zx1_sumfree.py` | sum-free circuits through PyZX (162 circuits, 378 runs) | `python3 zx1_sumfree.py` | not measured (ran on a Colab machine) | not measured | numpy, pyzx==0.10.7 | not re-run |

Other files here: `skeptic_small_checks.py`, `zx2_slot.py`, `zx3_ancilla.py`, `zx4_xw_sumfree.py`.

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
