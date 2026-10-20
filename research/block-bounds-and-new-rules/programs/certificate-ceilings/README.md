# Ceilings of the certificate families

**What it shows.** How far each kind of certificate can ever go: the hook family's ceiling 4(n - 1)/(n + 2), the
coordinate ceiling 4(n - 1)/n, two-factor Venn sets (never better than the hook), and the arithmetic of the
block-size table of the section README (savings ln(mW/M)/ln m, the bound needed to rule out 0.1 at each block).

**Programs.**
- `numbers.py` (pure Python, exact fractions): the arithmetic of the table. `python3 -I numbers.py`.
- `table.py` (numpy): the ceiling table by block size. `venn2ceil.py M [exact]`: ceilings of two-factor Venn sets
  (float by default; `exact` uses Python integers). `mucheck.py`, `hook_kinv.py M K`: checks of the hook's demand
  upper bound (the copy of Λ^(k-1) inside N ∩ KN).
- Original outputs in `../../data/certificate-ceilings/`.

The ceilings are argued by hand, not checked (upheld by three AI skeptics, who corrected several labels); the programs
check the formulas and screen families.

**Re-run for this revision:** yes, on 19 October 2026: `numbers.py` and `table.py` give output identical to the
original (`../../data/certificate-ceilings/rerun_2026-10-19.txt`). `venn2ceil.py`, `mucheck.py` and `hook_kinv.py`
were not re-run.
