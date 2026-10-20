# Shape R: the exact test against every fourth and fifth move

**What it shows.** `shapeR.py` reduces "do some 4th and 5th moves, with any contents, finish both arrays?" to a small
exact system (30 unknowns per generator, by averaging over the group of letters that pass through the paid move) and
decides it with coefficient sums on K's frame. `run_fam.py` runs ten families of shape-R prefixes; `theoremR_check.py`
types each prefix as in Theorem R and checks the theorem's prediction; `ctl.py` runs the controls (two private chains
must pass; Theorem A's prefixes must fail).

**Commands.** Build the small C library that `bm.py` loads before the first run: `cc -O2 -shared -fPIC -o libmp_mac.so mp.c` on a Mac, or `cc -O2 -shared -fPIC -o libmp_lin.so mp.c` on Linux, in the folder that holds `mp.c`. Then (numpy):
- `python3 ctl.py 1 3` (about a second): 9 of 9 target pairs pass for the private chains, 0 for Theorem A's prefixes;
- `python3 theoremR_check.py FAMILY 101 40` for FAMILY in A, RAWS, RAWR, SECR, MIX, COS, MIXED, MIXCOS, TYPEM, TYPEM1
  (5 to 10 seconds each): every line ends "AS THEOREM R" or "outside" (prefixes with two terms on one frame);
- `python3 run_fam.py FAMILY NPREFIX SEED OUT.jsonl` for the family sweeps behind `summary_final.txt`.

**Expected output.** `theoremR_check_s101.txt`, `theoremR_check_TYPEMR_s103.txt`, `summary_final.txt` (6,820 prefixes,
70,228 target pairs, 0 passes).

**Re-run for this revision:** `ctl.py 1 3` and `theoremR_check.py` for families A and TYPEM, on 19 October 2026; both
match the saved lines. The full family sweeps (minutes to hours, partly on Linux machines of Google's Colab service)
were not re-run.
