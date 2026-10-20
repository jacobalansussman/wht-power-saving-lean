# Block 4: at least 41/13 moves per array

**What it shows.** The k = 8 hook at f = 9 (dimension 102,960) has demand at least 86,592 per array against a cap of
27,456 per move, so every block-4 unit needs at least 41W/13 paid moves. `hookG.py` gets there in 14 seconds by
splitting the demand with one more symmetry (a group of bit swaps that commute with the transform and keep the hook):
each character's rank is computed modulo a prime (a certified lower bound), and the rank check of the column set is
exact.

**Programs** (numpy; `ent.c` is a small C helper for exact integer determinants):
- Build the helper once: `cc -O3 -shared -fPIC -o libent.so ent.c -lpthread`.
- `python3 -I hookG.py 4 8 01,23` (14 seconds, under 1 GB): prints the four character parts and
  `f(K) >= 86592 (certified)`. `b56.py` is the general module it builds on (any exterior-power product, any block).
- `exact-entry-demand/demcalc.py`: another agent's exact-entry demand calculator for coordinate certificates (used to
  screen other certificate families at block 4).
- `skeptic-rebuild/mydem.py`, `formulas.py`: the skeptic's independent rebuild from exact Hadamard minors.

**Expected output.** `RESULT hookG m=4 k=8: dim N 102960, rank check 102960 EXACT, f(K) >= 86592 (certified); delta
upper bound 91520`. Original outputs in `../../data/block4-hook/` (also the block-5 hook k = 4, the rank-2 block
values and the engine-price cap).

The separate certificate of the same demand (a 5,412 × 5,412 minor) and its checker are in
`../hook-certificates-blocks-2-to-6/`.

**Re-run for this revision:** yes, on 19 October 2026 on an Apple M4 Mac mini: the same four parts (1,382 + 1,366 +
1,366 + 1,298 = 5,412) and `f(K) >= 86592 (certified)` in 48 seconds (`../../data/block4-hook/rerun_2026-10-19.txt`).
The skeptic's rebuild and `demcalc.py` were not re-run.
