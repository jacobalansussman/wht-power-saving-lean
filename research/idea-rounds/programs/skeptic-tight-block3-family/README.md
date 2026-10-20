# Skeptic 1: the tight block-3 family (lead D1)

The first skeptic's programs of 16 October. Python with numpy. They read the matrices of the small irreducible
representations of Sp(6,2) from `../../data/skeptic-tight-block3-family/irreps_small.npz` (3.7 MB, computed in the theory
round of 14 October; set `IRREPS_NPZ` to use another copy). Each run takes under two minutes and under 2.5 GB. Saved
outputs: `../../data/skeptic-tight-block3-family/`.

| Program | What it shows | Command |
|---|---|---|
| `s1_depth_cuts.py` | the depth each relaxation sees (calibrated on the frame module, which needs depth 3); all cuts of the 56- and 105-dimensional irreducibles, with the refusals of 5 in 14 and of the tight mixed shapes | `python3 s1_depth_cuts.py depth > out_s1_depth.txt`; `python3 s1_depth_cuts.py cuts > out_s1_cuts.txt` |
| `s2_labels.py` | 3,486 labelled certificates on the coupled pair 56 + 105a: none beats 2.8 | `python3 s2_labels.py > out_s2_labels.txt` |
| `s3_singular.py` | from scratch: the 30 singular frames, their step graph, the group S8 and its stabiliser AGL(3,2), the 8 sectors | `python3 s3_singular.py > out_s3_singular.txt` |
| `s4_typegraph.py` | the type graphs of the move and the target in the irreducibles behind the depth numbers | `python3 s4_typegraph.py > out_s4_typegraph.txt` |

Re-run for this revision: all four programs (`s1_depth_cuts.py` in both modes), same output.
