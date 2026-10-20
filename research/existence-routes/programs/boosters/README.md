# Boosting pieces: the screen at blocks 4 to 8

**What it shows.** `certlib.py` is the certificate library (certificate A and B, hooks, projections, distance sets, isotropic
subspaces); `rank3.py` and `wronskian.py` are the exact identities behind Theorem B1 and Corollary B2; `lp*.py`,
`dual4*.py`, `dualpolar.py`, `lines6.py`, `projall.py`, `check_dist.py` and `certB.py` the linear programs and new
certificates; `intenum.py` enumerates relations and `check_survivor.py` checks the survivors S6, S7 and S8; `proj/`
holds exact projection counts. `calib.py` calibrates on two community words and needs the hub's price profile of T74
(`T74E_JSON`, not shipped; it is the hub's data). `skeptic/` is the AI skeptic's code (`vB1.py`, `vsurv.py` with about
65,000 tests per survivor, `vcert.py`, `vlp*.py`, `vlines.py`, `vbig.py`, `vcmp.py`; `cjob/` and `cjob2/` its larger
runs).

**Commands.** `python3 rank3.py`; `python3 wronskian.py`; `python3 check_survivor.py`; `python3 intenum.py 4:8 5:10`; in `skeptic/`: `python3 vB1.py`, `python3 vsurv.py`.

**Needs.** Python 3 for the identities; numpy and scipy for the linear programs.

**Time and memory.** Seconds to two minutes each (the linear programs up to about 2 minutes).

**Expected output.** Their `out_*.txt`.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 rank3.py` (0.1 s): matches the shipped output. Output: `rerun_2026-10-19_rank3.txt`.
- `python3 wronskian.py` (0.1 s): matches the shipped output. Output: `rerun_2026-10-19_wronskian.txt`.
- `python3 check_survivor.py` (4.3 s): matches the shipped output. Output: `rerun_2026-10-19_check_survivor.txt`.
- `python3 vB1.py` (in `skeptic/`) (0.1 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_vB1.txt`.

The linear programs, the enumerations and the skeptic's screens were not re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
