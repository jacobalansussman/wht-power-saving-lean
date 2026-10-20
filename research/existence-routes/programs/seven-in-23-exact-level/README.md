# The block-4 survivor 7 in 23 at the exact level

**What it shows.** `lift.py` lifts the symmetric survivor (order O O F O F) to actual 16 × 16 Clifford letters, using `w5lib/`
(the blueprint's library, copied); `f2blind.py` shows that the f = 2 widened test (`rep2.py`, from an earlier wave) is
blind when the move contents are unknowns: one array with two moves passes, where four are needed; `run1*.py`,
`d4test.py`, `d4gen.py`, `fullrow.py` and `flat.py` are the per-representation tests of the generic and the flat lift;
`ooffo*.py` and `foofo.py` the other move orders; `kmax.py`, `ubad.py`, `kals.py`, `kfree.py` and `crit.py` the rank
wall in the 14-dimensional representation. `skeptic/` is the AI skeptic's code (`vd4V.py`, `vd4addr.py`, `vf2blind.py`,
`voofo.py`; `critrun2/` the criterion on a second lift). The exclusion this folder first claimed was refuted by that
skeptic; the later closure of the construction is in `../../../moonshot-block4-block5/`.

**Commands.** `python3 f2blind.py`; in `skeptic/`: `python3 vf2blind.py`; the other programs as their headers say.

**Needs.** numpy.

**Time and memory.** `f2blind.py` about a second; the representation tests minutes, some on Colab machines.

**Expected output.** `out_f2blind.txt`, `skeptic/out_vf2blind.txt`.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 f2blind.py` (0.9 s): same verdict (residual 2.6e-30, i.e. it passes); the differences are floating-point noise and numpy warnings of the original run. Output: `rerun_2026-10-19_f2blind.txt`.
- `python3 vf2blind.py` (in `skeptic/`) (1.5 s): same verdict; floating-point noise in one residual (6.1e-16 against 6.5e-16). Output: `skeptic/rerun_2026-10-19_vf2blind.txt`.

The lift, the representation tests and the rank-wall programs were not re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
