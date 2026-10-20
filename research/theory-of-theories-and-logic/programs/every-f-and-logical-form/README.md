# Logic, part 1: the finite engine groups, the every-f bound, and the transfer counterexample

Section 14. `engine_group_exact.py` and `transfer_and_cs.py` are the skeptic's; the rest are the author's. `transfer_and_cs.py`, `block3_traces.py` and `widened_check.py` import a sibling file from their own folder.

Outputs of the original runs: [data/every-f-and-logical-form/](../../data/every-f-and-logical-form/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `fbound.py` | engine groups at blocks 1 and 2, r', the block-1 ranks | `python3 fbound.py` | 3 s | 75 MB | numpy | identical |
| `engine_group_exact.py` | the skeptic's exact engine groups | `python3 engine_group_exact.py` | 90 s | 200 MB | Python 3 only | identical |
| `transfer_and_cs.py` | counterexample to the transfer criterion; infinite controlled-S element | `python3 transfer_and_cs.py` | under 1 s | under 100 MB | numpy | identical |
| `hankel.py` | why twiddle factors are not finite-state | `python3 hankel.py` | under 1 s | under 100 MB | numpy | identical |
| `block3_traces.py` | block-3 trace moduli (r'(3) >= 28) | `python3 block3_traces.py` | 1 s | under 100 MB | numpy | identical |

Other files here: `widened_check.py`.

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
