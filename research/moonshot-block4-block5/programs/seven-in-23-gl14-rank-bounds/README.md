# 7 in 23 at block 4: the rank bounds in the 14-dimensional representation

**What it shows.** The 7-in-23 group-algebra construction (Z7 symmetry; orders O O F O F and O O F F O) must satisfy, at
distance 4 from the input, a rank condition in every character of the symmetry and every representation of the letters.
In gl14 (the 14-dimensional representation that comes from how address maps permute the 15 nonzero addresses) the
slot module splits into constituents (dimension, multiplicity), stored in `dec_s21_gl14.npz`.
- `pair.py` computes the joint-copy bound r <= dim Z' - max(delta - m n). For O O F O F in the first shared move's own
  character it gives r <= 11, where 12 is needed (the constituent (14,2)).
- `one65.py` checks the single-constituent obstruction for O O F F O when only the four reachable slots {8, 9, 19, 28}
  are free: on the constituent (6,5) a 4-dimensional part must fit into at most 3 dimensions, so the rank is at most 13
  where 14 is needed. With every slot off distance 4 free (the weaker test) it passes.
- `reach.py` recomputes the reachable slots; `decomp.py` (with `d4gen.py`, `d4test.py`, `lift.py` and the copies of the
  wave-5 library in `w5lib/`) rebuilds the decompositions; `gn.py` and `mrc.py` reach the ranks 11 and 13 numerically.

**Commands.** numpy (and scipy for `mrc.py`).
- `python3 pair.py dec_s21_gl14.npz p 0 d4 all` (under a second): last line `JOINT-COPY BOUND r <= ... = 11`;
- `python3 one65.py dec_s21_gl14.npz 1 all30 1 20` (under a second): dim A/P = 4 against rank 3 for every v;
- `python3 one65.py dec_s21_gl14.npz 1 d4 1 5`: 2 against 2 (the weak test passes, as it must);
- `python3 decomp.py 21 1 gl14 <ngens> out.npz` rebuilds a decomposition (see the docstring; seconds to minutes).

**Expected output.** `rerun_2026-10-19.txt`.

**Re-run for this revision:** `pair.py` and `one65.py` as above, on 19 October 2026; the numbers match the track's record
(11; 20, 16, 4 against 3; 2 against 2). The decompositions were not rebuilt; they are floating-point objects whose rank gaps
are 13 orders of magnitude or more, and the skeptic rebuilt them with its own code on two lifts.
