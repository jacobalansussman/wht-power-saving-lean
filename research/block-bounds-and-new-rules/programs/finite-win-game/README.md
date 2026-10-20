# OpenAI's finite-win game

**What it shows.**
- `zeta-kernel-8-entries/`: for the subset-sum kernel T = [[1,1],[0,1]] on 8 entries, no word with 11 or fewer
  calls computes T ⊗ T ⊗ T, so Yates's 12 additions are optimal there. `cutip.py` gives a cut bound for every order of
  the outputs (an output-order class up to permutations of the 3 address bits); `lp_pis.py` keeps the 1,300 of 6,828
  classes whose bound is at most 11; `lp_cert.py` writes and checks exact integer certificates for the other 5,528;
  `lp_duals.py` adds dual bounds used for pruning; `tdfs.c` is the exhaustive depth-first search over words with
  those prunes; `cfree_bfs.py` shows that words without cancellation need 12 calls; `algebra.py`, `yates_check.py` and
  `lgvfilter.c` are the exact tests of candidate leaves. `run_all.sh` runs the whole chain:
  `PYTHON=/path/to/python3 sh run_all.sh` (a C compiler and numpy with scipy 1.9 or later; about 15 minutes, under
  200 MB). Expected: 5,528 classes excluded by certificate, all 1,300 searched classes complete with no survivor at
  11 calls (and none at 9 or 10), and both 12-call Yates words passing every prune.
- `hadamard-16-entries/`: tools for the Hadamard kernel on 16 entries, where 28 to 31 calls are open: `std_check.py`
  (checks a word), `solve_h.py words.txt [starts]` (least squares for the coefficients of a word shape; a failure is
  evidence only), `check_found.py`, and the annealing search `sa4.c`, which found nothing below 32 calls but also
  failed to re-find a 32-call word with its strongest objective, so its failures carry no weight. Outputs in
  `../../data/finite-win-game/`.

**Re-run for this revision:** yes for the zeta kernel: `run_all.sh` on 19 October 2026 on an Apple M4 Mac mini, in 13
minutes: 5,528 classes excluded by exact certificates (none failed), all 1,300 searched classes complete at 11 calls
with no survivor, none at 9 or 10, both 12-call Yates words pass every prune, and the identity order still has 1,451
leaves at 12 calls, as it must (`../../data/finite-win-game/rerun_zeta_2026-10-19.txt`). The 16-entry tools were not
re-run.
