# Hook certificates for blocks 2 to 6, with separate checkers

**What it shows.** The numbers behind the lower bounds: the cap of every hook used (blocks 2 to 6: 4, 80, 27,456,
64,960, 438,480, 60,512) for every direction, by exact vectors over the Gaussian integers; the demands (block 2: 8;
block 3: 224; block 4: 86,592 = 16 × 5,412; blocks 5 and 6 in the widened game: 138,880 and 122,976) by nonsingular
minors modulo a prime; the cap of two moves in a row (144 at block 3, 48,576 at block 4); the exact gain of the second
of two plain moves at block 4 (24,288). Each generator is untrusted; each checker shares no code with its generator.
A nonzero minor modulo p proves the rank over Q(i) is at least its size, which is the safe direction for a demand.

**Programs** (C compiled with `cc -O2`, and pure Python):
- `gen_m3.py`, `check_m3.py`: block 3 (cap ≤ 80 for all 7 directions, demand ≥ 224, two moves ≤ 144). Run
  `python3 -I check_m3.py` in a folder that holds `cert_m3_dem.json`, `cert_m3_B2.json`, `cert_m3_M0.json` (1 second).
- `gen_m0.c`, `check_m0.c`: demand certificates (one character block M_0). `./check_m0 cert_m2_M0.txt` (block 2,
  instant); `./check_m0 cert_m4_M0.txt` (block 4: a 5,412 × 5,412 determinant modulo 2,147,483,629; about an hour,
  under 250 MB). `run_m4.sh` regenerates and checks the block-4 certificate.
- `check_capA.c`: caps for every direction, exact, no modulus: `./check_capA M K` with (2,2), (3,4), (4,8), (5,4),
  (5,5), (6,3) (seconds each).
- `check_yint.c`: the second step of a plain chain at block 4: `./check_yint 4 8 1 2` (1 second).
- `gen_b2.c`, `check_b2.c`: two moves in a row: `./check_b2 cert_b2_m3.txt`, `./check_b2 cert_b2_m4.txt` (under a
  second).
- `gen_big.c`, `check_big.c`, `run_b56.sh`, `run_full.sh`: blocks 5 and 6 in the widened game (full-row certificates
  in `../../data/hook-certificates/full-rows/`; about 6 minutes per block-6 class and 20 minutes for block 5 on a
  Linux machine). `run_m4_alpha.sh` checks a second character block at block 4 (same rank 5,412).

**Data.** `../../data/hook-certificates/` holds every certificate (`cert_*`), three corrupted certificates that the
checker must reject (`neg_m3_*`), and the original outputs (`outputs/`). Copy the certificates next to the programs
before running them.

**Expected output.** `RESULT: PASS` (or `[PASS]` lines and a `SUMMARY`) for every certificate, `RESULT: FAIL` for the
three negative controls; the block-4 determinant is 1797789153 modulo 2,147,483,629.

**Re-run for this revision:** yes, on 19 October 2026 on an Apple M4 Mac mini: `check_m3.py`, `check_m0` on the
block-2 certificate, `check_capA` for all six (block, k) pairs, `check_yint 4 8 1 2` and `check_b2` at blocks 3 and 4
all PASS, and the three negative controls FAIL as they should
(`../../data/hook-certificates/rerun_light_2026-10-19.txt`). The block-4 demand certificate was checked again the same
day: determinant 1797789153 modulo 2,147,483,629, `RESULT: PASS`, in 56 minutes at 235 MB
(`../../data/hook-certificates/rerun_check_m4_2026-10-19.txt`). The block-5 and block-6 full-row certificates were not
re-checked.

These are the proof-logging track's certificates of 13-14 October 2026; their written proofs (the cap lemma, the
splitting of the demand into equal character blocks) are argued by hand, and the Lean files of
`../lean-block2-block3/` prove the cap lemma and the block-3 demand independently.
