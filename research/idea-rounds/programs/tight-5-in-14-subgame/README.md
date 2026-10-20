# The singular sub-game of a tight 5 in 14

Programs of the track of 16 October on 5 in 14 at block 3. Python with numpy. `lib514.py` builds Sp(6,2), the free
group P0, the Levi group, the moves and its own 56-dimensional irreducible; `libs8.py` holds the S8 dictionary and
Young's orthogonal form; `libcut.py` the isotypic decompositions and cuts. Run the programs from this folder in
numerical order: `c1_calibrate.py` writes `v56_basis.npz`, `c2_singular.py` writes `singular_lag.npy`, `c3_s8.py`
writes `s8_points_h_K.npy` and `c10_dgame.py` writes `dgame.npz`, which later programs read (copies are included so
that any program can also run on its own). Each run takes under a minute and well under 1 GB. Saved outputs:
`../../data/tight-5-in-14-subgame/` (`c2b_orbits.py` has none).

| Program | What it shows |
|---|---|
| `c1_calibrate.py` | the 56-dimensional irreducible, its character, the hook (cap 10, profile 10, 18, 28, ratio 14/5) and the cut profile of one plain chain |
| `c2_singular.py`, `c2b_orbits.py` | the 30 admissible positions are the 30 maximal totally singular subspaces of x·z; their graph is the point-plane incidence graph of PG(3,2); the symmetry group has order 40,320 |
| `c3_s8.py` | the dictionary with S8: a move is a transposition, the target a product of three disjoint ones, the free letters AGL(3,2) |
| `c4_cuts.py`, `c5_labels.py` | cuts in every S8 irreducible and the hook's pencil of labels: nothing beats 2.8 |
| `c6_symscreen.py`, `c7_symscreen2.py`, `c8_s4.py`, `c9_subcal.py`, `c14_105a.py`, `c15_cal105.py`, `c16_survivors.py` | the symmetry screens (per-irreducible hook balance, the 105-dimensional certificate, the orbit count), their calibrations on real symmetric programs, and the surviving shapes (one D8, six S3) |
| `c10_dgame.py`, `c11_staircase.py` | the sub-game's one-piece relaxation as a small linear-algebra game, and an explicit depth-1 "staircase" scheme that satisfies it |
| `c12_framegraded.py`, `c13_overlap.py` | frame-graded hooks (best 2.8) and the free-group types shared by the hook and the frame module (none) |
| `c17_s8labels.py`, `c18_s8sets.py` | generic labelled certificates on pairs and sets of S8 irreducibles: inconclusive (generic labels miss special subspaces) |

Command for each: `python3 cN_name.py > out_cN_name.txt`. Re-run for this revision: all 19 programs, in order: same output for the 18 with saved outputs (`c11_staircase.py` differs only in the source lines that numpy's overflow warnings print).
