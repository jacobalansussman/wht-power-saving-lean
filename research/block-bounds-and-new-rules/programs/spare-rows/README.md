# Spare rows: certificates for rank-2 and rank-3 targets

**What it shows.** `comm_cert.py`: in the game OpenAI's network lives in, rank-2 targets need 2 moves per array and
rank-3 targets 3, by a small certificate from quadratic forms over the two-element field (the proof is by hand; the
program checks the certificate). `cut_lp.py`, `venn_mask.py M A B K J RANK`, `hook4odd_b3.py K`: certificates in the
full game with every phase and address map free (1.8485 and 1.7692 for rank 2 at 3 and 4 rows; 2.4883 for rank 3 at 4
rows). `spect_tiny.py KIND M W P MINUTES`: a small equal-frame search over the field with three elements (it cannot
rule out units over the rationals); it imports the search modules `frames.py` and `game.py` of an earlier round,
which are not published here, so it does not run from this folder alone. `test_rank.py`: rank helpers.

**Commands** (numpy): `python3 -I comm_cert.py`, `python3 -I cut_lp.py`, and as in each file's usage line. Original
outputs in `../../data/spare-rows/`; some contain floating-point warnings from the float screens.

**Re-run for this revision:** no.
