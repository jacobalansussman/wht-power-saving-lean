# The hook at block 3 in the widened game

**What it shows.** The certificate behind 14W/5 when every address permutation (Toffoli included) and every phase
pattern is free: the subspace N = span{e_x ⊗ e_T : x in T} of V ⊗ Λ^4 V at f = 5 (dimension 280) is kept by every
monomial map, one paid move raises it by at most 80 (200 explicit vectors), and one finished array needs 224. At
block 2 (f = 3) the same construction gives cap 4 and demand 8, so 2 moves per array.

**Programs** (numpy):
- `verify_hook.py M F`: the author's check: exact Gaussian-integer vectors for the cap, a mod-p rank for the
  demand. `python3 -I verify_hook.py 3 5` (5 seconds) and `python3 -I verify_hook.py 2 3` (1 second).
- `skeptic-recomputation/hookcheck.py M F [exact] [allblocks]`: a skeptic's independent recomputation (no shared
  code): ranks modulo two primes, and with `exact = 1` exactly over Z[i] by fraction-free elimination (the exact
  kernel rank takes about 6 minutes). `python3 -I hookcheck.py 3 5 0` (4 seconds).
- `skeptic-keptness/kept_check.py`: a second skeptic types the letters from the Lean files and checks, exactly, that
  all shifts, all 168 additive bijections, all phase letters, the generators of all permutations, Toffoli, CCZ,
  controlled-S and other diagonals keep N; then the 200 cap vectors and the values f for every direction and every
  ordered pair of directions (5 to 8 minutes, 1.3 GB). `stress.py SEED STEPS` composes random words and re-checks
  keptness.

**Expected output.** `RESULT m=3 f=5: cap = 80 ..., dem >= 224: every unit needs >= 14/5 W`; the skeptic's values 80,
144, 224; keptness with 0 violations. Original outputs in `../../data/widened-game-hook-block3/`.

**Re-run for this revision:** yes, on 19 October 2026 on an Apple M4 Mac mini: `verify_hook.py 3 5` and `2 3`,
`hookcheck.py 3 5 0` (modulo two primes; the exact run was not repeated) and `kept_check.py` (464 seconds under load)
gave the same values and no violation (`../../data/widened-game-hook-block3/rerun_2026-10-19.txt`). `stress.py` was
not re-run.
