# Skeptic 2: targets on a subspace (lead D2, refuted)

The second skeptic's programs of 16 October. Python with numpy. `k1`, `k2`, `k2b` and `k4` import the Reed-Muller
round's modules (`t_free.py`, `t_small_alt.py`, `tile_alt.py`, `xbase.py`) and read its survivor list
`surv_fe_alt_415_all.txt`; by default they look for them in `../../../reed-muller-and-one-move/programs/phase-letters/ends-only-block4/`
and `../../../reed-muller-and-one-move/data/phase-letters/ends-only-block4/` (set `RM_PHASE_DIR` and `RM_PHASE_DATA` to
override). Each run takes seconds. Saved outputs: `../../data/skeptic-targets-on-a-subspace/`.

| Program | What it shows | Command |
|---|---|---|
| `k1_sumzero.py` | none of the 1,612 surviving 4-in-15 classes has four targets summing to zero (needed without scratch); the 448 targets rebuilt independently | `python3 k1_sumzero.py > out_k1_sumzero.txt` |
| `k2_sumzero_nonalt.py` | all 12,333 zero-sum sets of non-alternating targets fail the earlier exact tests | `python3 k2_sumzero_nonalt.py > out_k2_sumzero_nonalt.txt` |
| `k2b_calib.py` | the harness passes all 1,612 survivors and the plain 4 in 16, and fails 5 in 15 | `python3 k2b_calib.py > out_k2b_calib.txt` |
| `k3_core_hook.py` | the Theorem R table: every no-scratch core with outside gates spanning at most 4 directions is excluded (15 gates: 414,688 needed against 411,840) | `python3 k3_core_hook.py > out_k3_core_hook.txt` |
| `k4_sumzero_alt.py` | all 14,923 zero-sum sets with an alternating target fail too, without the survivor file | `python3 k4_sumzero_alt.py > out_k4_sumzero_alt.txt` |

Re-run for this revision: all five, same output.
