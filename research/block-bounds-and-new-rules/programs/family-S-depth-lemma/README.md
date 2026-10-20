# Family S: single-term contents need 2m moves for two arrays

**What it shows.** If every paid move acts on a content that holds at most one term of each input, two arrays need at
least 2m moves at every block size, whatever the free group. b4-depth4's integer programmes (`integer-programmes/`)
encode the constraints of such schemes; `check_farkas.py` checks Farkas certificates showing that their minimum is at
least 2m (and 2m + 1 for the stricter programme) for every m, and model-checks the constraints by enumeration at
m = 2 and 3.

**Command.** `python3 -I check_farkas.py ../../data/family-S/farkas_familyS.json 3` (pure Python; about 3 minutes).
**Expected output:** `[PASS]` lines and `RESULT: PASS` (835,424 and 4,115,761 enumerated configurations, 0
violations).

The integer programmes need scipy (`familyS_ip.py`, `familyS_ip2.py`; `chain9.py` and `depthctl.py` are the depth
controls, and `chain9.py` imports the exact rank routine `fastrank.py` of the same folder); their outputs are in
`../../data/family-S/`. The constraints themselves were written by the b4-depth4 track and were not re-derived by the
certificate's author; one skeptic's own program agrees on the hardest case.

**Re-run for this revision:** yes, on 19 October 2026: `RESULT: PASS` in 190 seconds
(`../../data/family-S/rerun_2026-10-19.txt`). The integer programmes were not re-run.
