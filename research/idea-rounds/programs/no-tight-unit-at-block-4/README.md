# No exactly tight unit at block 4

The recombiner's check of 16 October, with the earlier block-4 filter it reuses. Python with numpy; `nullsp.c` needs a
C compiler. Saved outputs: `../../data/no-tight-unit-at-block-4/`.

| Program | What it shows | Command |
|---|---|---|
| `b4_admissible.py` | enumerates the engine positions at block 4 by depth (orbits of the free group; calibrated on the counts 1, 30, 280) and runs the termwise filter on one word per position, exactly modulo 33,554,393; controls (identity, target, target then free letters) must pass | depth 1: `python3 b4_admissible.py > out_b4_depth1.txt` (about 2 minutes); depth 2 in two halves: `DEPTH_MAX=2 NE=1 ONLY_DEPTH=2 SLICE=0/2 python3 b4_admissible.py > out_b4_depth2_s0.txt`, and the same with `SLICE=1/2` (about 3 minutes each) |
| `enum_test.py` | the enumeration alone; at depth 3 it gives 1, 30, 280, 960 | `python3 enum_test.py 3` |
| `dep_supports.py` | the smallest sets of directions whose rank-one matrices z·zᵀ add up to zero: 7 elements, exactly the 3-flats (1, 15 and 155 of them at blocks 3, 4 and 5) | `python3 dep_supports.py > out_dep_supports.txt` |

**The filter's data file.** `b4_admissible.py` imports `make-E/filter9.py`, which reads `E_r1_p33554393.npy`, a basis of
the character-0 block of E = N ∩ K⁻¹N for the f = 9 hook (1,023 vectors of length 6,435, about 53 MB). It is not
published, because of its size; rebuild it in `make-E/`:

```
cd make-E
cc -O3 -o nullsp nullsp.c
python3 dumpK.py dump 1 E            # writes ME_r1_p33554393.bin, the 6,435 x 6,435 block of K mod p
./nullsp ME_r1_p33554393.bin nullE_r1_p33554393.bin
python3 dumpK.py load 1 E            # checks M x = 0 on every basis vector, writes E_r1_p33554393.npy
python3 filter9.py 3                 # the filter's self-test
```

The file used on 16 October has sha256 `483065e5b2cec24dd0730531d36e79c2ae1b2bdebc4bdc45646703616358dcc8`.
For this revision I rebuilt it with the first four commands on an Apple M4 Mac mini: the result is byte-identical
(same sha256), in about 4 minutes and 0.6 GB. Set `E_DIR` to another folder if the file lives elsewhere. `blk9.py` builds the hook blocks; it was
written for an earlier block-4 round (see also `../../../moonshot-block4-block5/`).

Re-run for this revision: `b4_admissible.py` at depth 1 and in both depth-2 halves, `enum_test.py 3` and `dep_supports.py` re-run: same results (0 of 30, 0 of 140 and 0 of 140; the depth-1 output has one extra progress line that the saved run, made with an earlier copy of the program, did not print); the data file E rebuilt from `make-E/`, byte-identical.
