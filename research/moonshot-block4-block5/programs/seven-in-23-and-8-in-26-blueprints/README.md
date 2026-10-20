# The blueprint targets 7 in 23 and 8 in 26 at block 4

**What it holds.**
- `wave3/`: the first search (bp.py, bpgen.py and the character-split tests chidiag.py, chidiag2.py; junkdist.py measures
  how far junk sits from the inputs). It showed that the blueprint as written (both shared moves first) and its natural
  repair fail. Build the small C library that `bm.py` loads before the first run: `cc -O2 -shared -fPIC -o libmp_mac.so mp.c` on a Mac, or `cc -O2 -shared -fPIC -o libmp_lin.so mp.c` on Linux, in the folder that holds `mp.c`. (for the `bm.py` in this folder).
- `skeptic-wave3/`: the skeptic's own code. `symp2.py` found the 3,024 configurations of order O O F O F that escape
  the claimed theorem (out_symp2.txt); `cex_t_Z7.npy` is one escaping term. These escapes were later shown spurious.
- `wave5/`: the frame model (`sp8.py` with `vlib/`), the exact frame-level relaxation `relax.py`, the single-term
  exclusion (`t2_d3.py`), the pair runs (`t7c_pairchars.py`), and in `s2/` the group-algebra construction for 7 in 23
  (`t9_galg.py`, `t10_build.py`, `t10b_legit.py`, the witness `t13_witness.py` with its two `.npy` files, `t15_subalg.py`).
  `t14_z2cube_allsubs.py` is the 8-in-26 test.
- `skeptic-wave5/own.py`: the skeptic's independent rebuild from binary symplectic data only. Mode `z7` solves the
  reduced 7-in-23 problem (it passes, but generically: 180 equations, 208 unknowns, full row rank), with controls; mode
  `z2` runs the 8-in-26 test over all 11,391 elementary abelian subgroups of order 8 (every one fails in at least 4 of
  the 7 characters tested); `a0` shows that 7 in 22 fails structurally (distance 0.4564 from a fixed subspace).

**Commands** (numpy). In `skeptic-wave5/`: `python3 own.py 0 z7 0` (8 seconds on an Apple M4 Mac mini (16 GB)), `python3 own.py 0 z2 0`
(about 8 minutes), `python3 own.py 0 a0`. The wave-5 scripts take the arguments given in their docstrings.

**Expected output.** `skeptic-wave5/out_own_*.txt`; `wave5/out_*.txt`, `wave5/s2/out_*.txt`.

**Re-run for this revision:** `own.py 0 z7 0` and `own.py 0 z2 0` on 19 October 2026; every one of the 11,391 subgroups fails in at least 4 characters, as
recorded (`skeptic-wave5/rerun_own_*_2026-10-19.txt`). The wave-3 and wave-5 track scripts were not re-run.
