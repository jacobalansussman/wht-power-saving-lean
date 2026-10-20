# Root systems, index codes and reflection forms

Programs of the root-system track of 16 October. Python with numpy. Saved outputs:
`../../data/root-systems-and-index-codes/`. The probes (`s4`, `s5`, `s6`) are numerical searches from random starts:
evidence only, each calibrated on a known positive and negative case first.

| Program | What it shows | Command |
|---|---|---|
| `s1_e7_dictionary.py` | the block-3 frame game is W(E7) up to sign: E7's 126 roots onto the 63 nonzero vectors of F_2^6, its 135 frames onto the 135 Lagrangians, a paid move is a "J − 2I" block, the 36 A7 and 28 E6 quadrics | `python3 s1_e7_dictionary.py > out_s1_e7_dictionary.txt` (seconds) |
| `s2_index_code.py` | the crossing centre as a linear index code: copy ratios E8 0.600, E7 0.778, E6 1.167; E8 and E7 are perfect codes | `python3 s2_index_code.py > out_s2_index_code.txt` (under a second) |
| `s3_subspace.py` | the block-5 7-array problem is a three-qubit Pauli problem inside E8, with no symmetry in GL(5,2) | `python3 s3_subspace.py > out_s3_subspace.txt` (about 6.5 minutes: the scan of all of GL(5,2)) |
| `s4_refl_probe.py` | the no-scratch reflection form at block 4: finds plain 4 in 16, refuses an excluded 4 in 13, and stops at the trivial split on 4 in 15 | `python3 s4_refl_probe.py pos 1`, `... neg 1`, `... probe 1` (and seeds 2, 3) |
| `s5_sub7_probe.py` | the same at block 5: its positive control fails, so no block-5 claim is made from it | `python3 s5_sub7_probe.py pos 7 1 1` |
| `s6_core_probe.py` | block-4 cores with alternating or mixed targets: calibrated on 4 alternating targets in 20 gates; the cores stop at the trivial split | `python3 s6_core_probe.py altpos 1 3`, `... alt4 1`, `... mix 1` |

Re-run for this revision: `s1`, `s2`, and `s4` with `pos 1` and `neg 1` re-run: same output; not re-run: `s3_subspace.py` (6.5 minutes), the `s4` probe seeds, `s5`, `s6`.
