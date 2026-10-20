# The block-5 butterfly: the exact census behind Theorem W3

**What it shows.** In the butterfly 2 in 7 at block 5 (private starts, a shared trunk of three moves, private ends), the
top layer needs a special-position merge. `census_sp.py` lists them exactly over the two-element field in the dual polar
graph of Sp(10,2): per trunk path, 1,536 ordered anchor pairs with a special-position merge, 393,216 configurations, and
the leftover sits three or four layers out, each case equally often. `trans_sp.py` shows that all 3,360 trunk paths look
alike; `second_sp.py` looks at second branches. `skeptic/v_geo.py` is the skeptic's code, written from scratch, that
repeats the census over every path directly and checks the exact fact used in the proof (0 violations over 1,720,320
pairs).

**Commands.** `python3 census_sp.py 60 1` (7 seconds on an Apple M4 Mac mini (16 GB)); in `skeptic/`: `python3 v_geo.py 7` (6 seconds) and
`python3 v_geo2.py` (standard library).

**Expected output.** `out_census_sp.txt`, `skeptic/out_v_geo_M7.txt`, `skeptic/out_v_geo2.txt`.

**Re-run for this revision:** `census_sp.py 60 1` and `v_geo.py 7`, on 19 October 2026; identical apart from timings.
