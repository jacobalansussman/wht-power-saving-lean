# single-term-game: the single-term two-input class, solved exactly as a game

Take two inputs, any number of scratch arrays, and suppose every array a paid move acts on holds at most one term of
each input. Then what a program can do depends only on where its terms sit: single inputs at frames, and pairs (one
term of each input tied by a fixed twist). A pair releases its inputs only when a second object meets it. `stg.py`
plays every move of that game exactly, with symmetry reduction. It errs on one side only, by design: it treats every
cycle of pairs as a release, which can only make a negative result safer. `toylab.py` and `sweep.py` hold the
abstract frame games (any finite group in place of the engine's) and the small games used here.

Run every command from this folder with Python 3.9 or later and numpy. The outputs of the original runs are in `../../data/single-term-game/` with the same names. Some programs write a JSON file next to themselves or into `out/`; run them from a copy if you want to keep the shipped files.

| program | what it does | command | time |
|---|---|---|---|
| `t_stg_toys.py` | exhaustive game on a toy frame game, every paid class | `python3 t_stg_toys.py D5 0 1500000` (names in `sweep.py`) | seconds to hours |
| `t_stg_sub.py` | the engine sub-games: the cube B_m and the Hamming game H(k, 3), optionally split into parts | `python3 t_stg_sub.py B4:0 7 3000000 4/0/6 two` (part 0 of 6 of the block-4 cube 2 in 7) | about 10 min per part on a Colab machine |
| `t_stg_2m1.py` | decides exactly whether a single-term 2 in 5 exists in a game with the target three steps away | `python3 t_stg_2m1.py S5D5 S5S3` | seconds to minutes |
| `t_stgcal.py` | calibration of the game against actual modules in the group algebra | `python3 t_stgcal.py D7:0:30:10` (see its header for the spec) | minutes |
| `t_inv.py` | tests candidate invariants for a proof of Conjecture T on every reachable state | `python3 t_inv.py S4Z2:0:6 B3:0:6` | seconds to minutes |
| `t_phi.py` | the dimension potential Φ along random plays (needs scipy) | `python3 t_phi.py Sp4:0:single:30:6` | minutes |
| `t_h43num.py` | a numerical no-scratch check in the Hamming game; it could not be calibrated at block 3, so its results mean nothing at block 4 and are not used | see its header | |

`../../data/single-term-game/runN/` holds the outputs of the machine runs, numbered as in the round. Run 16 is the
block-4 cube 2 in 7, in six parts. Runs 12 to 14 are the cube and Hamming games to 5 and 6 moves. Run 19 is the 2 in 5
decisions, runs 15 and 17 the invariant tests, run 18 the Φ tests, and run 10 the toy games. Writing a compiled version
of `stg.py`, so that it can reach the block-4 Hamming game at 6 to 7 moves and the block-5 cube, is the first open
experiment of the README.
