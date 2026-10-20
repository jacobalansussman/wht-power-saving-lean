# 5 in 14 at block 3: the single-term gain table

**What it shows.** In the engine game at block 3, a move whose content is a single free image of an earlier move's
output raises the hook by at most 72 (12 of the 14 classes of paid letter) or 0 (the other 2), below the cap of 80. So
a 5-in-14 unit, which is exactly tight on the hook, has no such move. The certificate is a set of exact vectors and
one nonsingular minor per class.

**Commands.** `python3 -I single-term-gain-certificate/check_r1.py ../../data/five-in-fourteen/cert_r1.json` (pure
Python; a few minutes). `gen_r1.py` (numpy) regenerates the certificate.

**Expected output.** 14 `[PASS]` lines and `RESULT: PASS; certified upper bounds on g2 over the 14 classes: [0, 0, 72,
..., 72]`.

The term-by-term filter and the greedy program that passes every test at f = 5 and fails at f = 6 are not published
as programs: they depend on a 105 MB pickle of depth-3 terms and on several helper modules that load pickles. Their
findings are reported in the section README as argued by hand (no skeptic) and computed by one program.

**Re-run for this revision:** yes, on 19 October 2026: `RESULT: PASS`, the same table, in 4 seconds
(`../../data/five-in-fourteen/rerun_2026-10-19.txt`).
