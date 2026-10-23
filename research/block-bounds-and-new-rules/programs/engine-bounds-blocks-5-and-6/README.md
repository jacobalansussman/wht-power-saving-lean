# Engine game, blocks 5 and 6: 2.99172, 2.97290, and 2.9447 from block 4 on

**What it shows.** The induction track carries the block-4 hook upward through projections between frame graphs,
controlled by a coupling argument, and gets at least 2.99172 moves per array at block 5, 2.97290 at block 6 and 2.9447
at every block from 4 on, in the engine game, with moves counted by rank. The inputs are the block-4 hook's exact
cap 27,456 and lower bounds on what two, three and four moves in independent directions do to it: 48,576, 68,320 and
86,592 (the safe direction). Each lower bound is a nonsingular minor modulo a prime, stored with a checker that
shares no code with the program that found it, in `../hook-certificates-blocks-2-to-6/`: 48,576 and 68,320 on all 16
characters, 86,592 on one character together with the hand-argued lemma that all 16 have equal rank.

**Programs** (numpy where imported):
- `induct.py [M]`: the induction table up to block M (default 60).
- `nrs.py`, `nrs_closed.py`: the coupling numbers; `tail.py`: the tail of the bound; `checks.py`: consistency checks;
  `prices.py`: indicative caps at engine prices (they need upper bounds on grouped moves that nobody has);
  `cutceil.py`: ceilings of cut certificates.
- Run each as `python3 -I NAME.py`; each takes seconds. Original outputs in
  `../../data/engine-bounds-blocks-5-and-6/`. Part (a) of `nrs.py` compares the closed form with an exact table of
  projection counts at blocks 2 to 6 from my earlier notes, now copied into
  `../../data/engine-bounds-blocks-5-and-6/exact-projection-table/` (before 23 October 2026 `nrs.py` looked for it in
  a folder that is not in this repository and stopped there).
- `fact-a/engine_m.py M` (pure Python): Fact A of the proof by matrices at blocks 3, 4 and 5: the engine's letters
  act on frames as the proof says, and the letter group is transitive on frames.
- `skeptic-checks/`: the two skeptics' own programs (pure Python). `fcheck.py` recomputes 27,456, 48,576 and 68,320
  by another route (a block of j moves splits over the untouched bits into ranks of at most 280 × 280); `sk_nrs.py`
  brute-forces the projection counts up to block 7; `sk_claim1.py` tests the recursion on families with no symmetry;
  `sk_tail.py` recomputes the bounds up to block 200; `sk_kappa.py` checks κ_m up to block 40; `tensor_cx.py` is the
  counterexample to Proposition 7(c). Outputs, and their re-runs of the programs above, in
  `../../data/engine-bounds-blocks-5-and-6/skeptic-checks/`.

**The written proof** is [../../proofs/induction-across-block-sizes.md](../../proofs/induction-across-block-sizes.md)
(argued by hand, not checked by a person), and the two AI skeptics' reports, which upheld the bounds and asked for
four corrections to side claims and wording (listed at the top of the proof), are in [../../proofs/induction-across-block-sizes-skeptics.md](../../proofs/induction-across-block-sizes-skeptics.md).

**Re-run for this revision:** yes, on 19 October 2026: `induct.py` (2.944709 at every block from 5 to 60) and
`checks.py` (block-5 and block-6 profiles ending at 2.9917 and 2.9729) give output identical to the original
(`../../data/engine-bounds-blocks-5-and-6/rerun_2026-10-19.txt`). The other programs were not re-run. On
23 October 2026, from a copy of this section as published: `nrs.py` (now with its table), `fact-a/engine_m.py` at
blocks 3, 4 and 5 and all six skeptic programs give their stored output again
(`../../data/engine-bounds-blocks-5-and-6/rerun_2026-10-23.txt`; `sk_nrs.py` takes about 70 seconds, the others
under 10).
