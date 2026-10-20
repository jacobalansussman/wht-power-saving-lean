# The commutative slice: exclusions for every W at once

**What it shows.** In the slice with only shifts, paid moves and mixing (no address maps, no phases), an earlier
reduction turns "a unit exists" into a test on the multiset of directions used. One probability distribution over
test functionals per block size then settles every W at once: at block 3 every multiset with fewer than 3W
directions fails (no saving), at block 4 every one with fewer than 15W/4 fails, at block 5 fewer than 31W/8. Separate
certificates settle 2 in 7, 3 in 10 and 4 in 13 at block 4 case by case. The reduction itself was not re-derived.

**Programs.**
- `check_lp.py CERT.json` (pure Python, exact fractions; seconds): checks a distribution certificate
  (`../../data/commutative-slice/cert_lp_m2.json` to `cert_lp_m5.json`). `gen_lp.py m OUT.json` (numpy and scipy)
  writes one.
- `check_matroid.py m W n CERTFILE` (pure Python): checks a case-by-case certificate, for example
  `python3 -I check_matroid.py 4 2 7 ../../data/commutative-slice/cert_m4_W2_n7.bin` (2 in 7: 0 of 116,280 multisets
  pass). `gen_matroid.py m W n OUTFILE` (numpy) writes one; `run_small.sh`, `run_3in10.sh`, `run_4in13.sh` regenerate
  and check the cases (set `PYTHON` to a Python with numpy; 4 in 13 takes about 10 minutes and writes a 40 MB file).

**Not published (size):** the certificates for 3 in 10 (3.9 MB) and 4 in 13 (40 MB); `run_3in10.sh` and `run_4in13.sh`
regenerate them. Their original check outputs are `out_3in10.txt` and `out_4in13.txt` in the data folder.

**Re-run for this revision:** yes, on 19 October 2026: `check_lp.py` on all four distribution certificates and
`check_matroid.py` on 2 in 7 (116,280 of 116,280 multisets excluded) all PASS
(`../../data/commutative-slice/rerun_2026-10-19.txt`). The 3-in-10 and 4-in-13 certificates were not regenerated.
