# Theorem F, the angle law and Theorem H

Section 7.

Outputs of the original runs: [data/angle-law/](../../data/angle-law/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `omega_exact.py` | Theorem F exactly up to block 40 | `python3 omega_exact.py 40` | under 1 s | under 100 MB | Python 3 only | identical |
| `frob_ceiling.py` | the spectral ceiling table | `python3 frob_ceiling.py` | under 1 s | under 100 MB | Python 3 only | identical |
| `hook_cells.py` | Theorem H on hooks | `python3 hook_cells.py` | 1 s | 50 MB | numpy | identical |
| `frames_check.py` | frame-family checks of Theorem F | `python3 frames_check.py` | 14 s | 640 MB | numpy | identical |

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
