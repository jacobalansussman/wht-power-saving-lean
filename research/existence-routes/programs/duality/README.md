# Strassen's duality, repaired, and boosting pieces

**What it shows.** `calib.py` calibrates the job model (it imports `frontier_flow.py` from `../tau-theorem-and-laser/`);
`lp_combo*.py` and `colab_lp/` are fluid linear programs that combine relations (the made-up booster: 53.16[1] ≥ [243],
exponent 0.723 against 0.834). `skeptic/a5_exact.py` checks exactly, in Gaussian integers for blocks 1 to 6, that one
direction costs at most two big calls on a register with spare rows, and certificate A's prices; `a5_certA.py`,
`wd_cert.py` (certificate A's analogue in the widened game), `mylp.py` (an independent fluid program), `v2_roots.py` (the
sets of exponents where each relation's balance is non-negative) and `v2_pipeline.py` (the pipelining schedule).

**Commands.** `python3 calib.py`; `python3 lp_combo.py ...`; in `skeptic/`: `python3 a5_exact.py`, `python3 v2_roots.py`, `sh run1.sh`, `sh run2.sh`.

**Needs.** numpy and scipy.

**Time and memory.** Seconds to a few minutes each.

**Expected output.** Their `out_*.txt`.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 a5_exact.py` (in `skeptic/`) (0.5 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_a5_exact.txt`.
- `python3 v2_roots.py` (in `skeptic/`) (0.2 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_v2_roots.txt`.

The fluid programs were not re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
