# Counting unknowns against equations

**What it shows.** `exact_counts.py`: the exact counts of Lemma C (in units of the free group's order |B|): unknowns, equations,
the naive dimension, the orbit and the dominance threshold of each shape; the plain 1 in 2 at block 2 has 6|B| = 18,432
unknowns against 15|B| = 46,080 equations. `dimcount.py` (with `batch.py`, the `specs_*.json` shape lists and `run*.sh`)
computes generic Jacobian ranks of the program map in small representations (stand-ins); `plainpoint.py` and
`specialpoint.py` find plain and special solutions there (`bsearch.py` supplies the engine's letters); `fullness.py`
tests whether the image fills the space; `mktable.py`, `summ.py` and `table_m2.md` summarise. `skeptic-compute/` is the
compute skeptic's independent code (exact ranks modulo two primes: `counts.py`, `jac.py`, `common.py`, `plain_modp.py`,
`orbit_modp.py`, `groups.py`, `widened*.py`, `block3_*.py`, `special_float.py`); `skeptic-logic/` is the logic skeptic's
(`indep.py`, `solve23.py`, `solve35.py`).

**Commands.** `python3 exact_counts.py`; `python3 batch.py specs_m2.json out.jsonl` (and the other spec files); `python3 skeptic-compute/counts.py`.

**Needs.** numpy; scipy for `specialpoint.py` and `special_float.py`.

**Time and memory.** The counts take under a second. The rank runs took minutes to an hour each on a Colab machine; one block-3 run of the skeptic ran out of memory and was stopped.

**Expected output.** Their `out_*` files.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 exact_counts.py` (0.1 s): matches the shipped output. Output: `rerun_2026-10-19_exact_counts.txt`.
- `python3 counts.py` (in `skeptic-compute/`) (0.2 s): matches the shipped output. Output: `skeptic-compute/rerun_2026-10-19_counts.txt`.

The rank runs (`dimcount.py`, `plainpoint.py`, `specialpoint.py`, `fullness.py` and the skeptics' rank programs) were not re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
