# Block 5: the three middle frames

**What it shows.** In the block-5 trunk with several branches, every target frame touches exactly three possible branch
frames, at distances 2, 4, 4 from the other array's anchor when the leftover is three layers out, and 3, 4, 4 when it is
four layers out; so each array has exactly one close branch. `v_census.py` is the skeptic's exhaustive check over all
anchor pairs of a walk (exact over the two-element field; no symmetry assumed).

**Command.** `python3 v_census.py 1` (or `2 3`), standard library, 8 seconds per walk on an Apple M4 Mac mini (16 GB).

**Expected output.** `out_v_census_1.txt`, `out_v_census_23.txt`. Only the three patterns (2,4,4), (3,4,4) and (4,4,4)
occur.

**Re-run for this revision:** `v_census.py 1`, on 19 October 2026 (`rerun_v_census_1_2026-10-19.txt`). Every per-row
count matches the saved file. The saved file's last total line reads 28,672 for the pattern (4,4,4), while its own rows
add up to 28,672 + 860,160 = 888,832; the re-run prints 888,832. The census statement is unaffected.
