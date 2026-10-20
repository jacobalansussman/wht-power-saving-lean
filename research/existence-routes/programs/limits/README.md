# Limits in the number of digit columns

**What it shows.** `e1/xw2lib.py` holds free words of the everything-free game at block 2 (copied from an earlier wave);
`e1/e1_gram.py` computes exact Gaussian-integer Gram matrices of menus of free words at f columns, `e1/e1_paid.py` the
same with one paid letter, `e1/e1_analyze.py` their smallest eigenvalues (does an overlap decay as f grows?). `e1c/`,
`e1p/` and `e1an/` are later variants and the analyses (`e1an/out_*.txt`). `skeptic/` holds the AI skeptic's
independent rebuild (`vex2lim_rb*/`: `myxw.py`, `gram.py`, `paid.py`, `ctl.py`, `analyze*.py`) and `local/twist_check.py`
(the global 8th-root phase).

**Commands.** `python3 e1/e1_gram.py <menu> <f>` then `python3 e1/e1_analyze.py ...` (see the headers of the programs).

**Needs.** numpy, scipy.

**Time and memory.** The Gram runs took minutes to hours on Colab machines and wrote about 195 MB of matrices, which are not published.

**Expected output.** The analyses in `e1an/`.

**Re-run for this revision.** No.

Nothing in this folder was re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
