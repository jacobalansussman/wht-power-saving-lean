# heights: Theorem C as a filter, and the test units of the Proposition 6 round

`heightfilter.py` decides, for any concrete word, whether a height system exists: integers h_l(x) per input l and
frame x, 0 at the raw frame, m at the transformed frame, changing by at most one per step, such that after every
paid move all terms of the moved array sit at one height. If one exists, the word cannot save (Theorem C), at any
block and with any number of arrays. It is a system of difference constraints, solved exactly by Bellman–Ford in
milliseconds. `best_bound` gives the quantitative form (a lower bound on M as a sum of shortest-path distances).

These programs import `../ledger/ledger.py`.

Run every command from this folder with Python 3.9 or later and numpy. The outputs of the original runs are in `../../data/heights/` with the same names. Some programs write a JSON file next to themselves or into `out/`; run them from a copy if you want to keep the shipped files.

| program | what it does | command | time |
|---|---|---|---|
| `heightfilter.py` | the filter (`feasible`, `best_bound`, `word_outputs`) and its self-test on the plain schedule, HOLO3, the order-3 control C6O3 and SPLIT1 | `python3 heightfilter.py` | 1 s |
| `checks.py` | the named units in the true model (all four column-count blocks): plain, gauge exchange, SPLIT1 (one input, m + 1 moves, two frames in every basis: refutes Proposition S), HOLO3 (three arrays, holonomy with one term per input: refutes Proposition H3), C6O3; which hypotheses of Theorems A, B and C each satisfies | `python3 checks.py` | 1 s |
| `t_chain.py` | random exact test of the two mechanical claims inside Theorem C's proof | `python3 t_chain.py 1` (and seeds 2, 3) | 2 s each |
| `t_scalar.py` | the `ledger.py` warning: a three-line word that `ledger.py` accepts and the true model rejects | `python3 t_scalar.py` | 1 s |
| `frames.py` | the frame graph and Sp(2m, 2) acting on it, up to block 4 (a module) | imported | |
| `tokgame.py`, `tokgame2.py` | exhaustive single-term two-input games at blocks 1 and 2 (calibration only) | `python3 tokgame.py`, `python3 tokgame2.py` | seconds |
| `tokbeam.py` | beam searches for a violation of L1 + L2 ≤ M at blocks 3 and 4 (a search, not a proof) | `python3 tokbeam.py 3 300 150 1` (seeds 1 to 3); `python3 tokbeam.py 4 150 100 1` (seeds 1, 2) | minutes |
