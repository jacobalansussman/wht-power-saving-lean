# Engine game, blocks 5 and 6: 2.99172, 2.97290, and 2.9447 from block 4 on

**What it shows.** The induction track carries the block-4 hook upward through projections between frame graphs,
controlled by a coupling argument, and gets at least 2.99172 moves per array at block 5, 2.97290 at block 6 and 2.9447
at every block from 4 on, in the engine game, with moves counted by rank. The inputs are the exact cap and mod-p
lower bounds on demands (the safe direction).

**Programs** (numpy where imported):
- `induct.py [M]`: the induction table up to block M (default 60).
- `nrs.py`, `nrs_closed.py`: the coupling numbers; `tail.py`: the tail of the bound; `checks.py`: consistency checks;
  `prices.py`: indicative caps at engine prices (they need upper bounds on grouped moves that nobody has);
  `cutceil.py`: ceilings of cut certificates.
- Run each as `python3 -I NAME.py`; each takes seconds. Original outputs in
  `../../data/engine-bounds-blocks-5-and-6/`.

The argument is written by hand (argued, not checked); two AI skeptics reproduced every number with their own code.

**Re-run for this revision:** yes, on 19 October 2026: `induct.py` (2.944709 at every block from 5 to 60) and
`checks.py` (block-5 and block-6 profiles ending at 2.9917 and 2.9729) give output identical to the original
(`../../data/engine-bounds-blocks-5-and-6/rerun_2026-10-19.txt`). The other programs were not re-run.
