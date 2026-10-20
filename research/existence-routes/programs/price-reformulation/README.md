# OpenAI's price language, and the refutation of the 0.036 ceiling

**What it shows.** `calib.py` evaluates the price language on the recorded unit values; `sparse/` holds the sparsity searches.
`skeptic/check_prop3.py` is the refutation: the hypothesis of Proposition 3 (OpenAI's laws with the tensor law replaced
by homogeneity) gives one matrix two different prices at every exponent below 1, so the 0.036 ceiling follows from no
possible price; it also recomputes the 14 recorded values and shows that 6 belong to impossible shapes.
`skeptic/sparse_check.py` checks the sparsity lemma.

**Commands.** `python3 calib.py`; in `skeptic/`: `python3 check_prop3.py`, `python3 sparse_check.py 1` and `2`.

**Needs.** Python 3 and numpy.

**Time and memory.** Seconds.

**Expected output.** Their `out_*.txt`.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 check_prop3.py` (in `skeptic/`) (2.1 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_check_prop3.txt`.

`calib.py`, the sparsity searches and `sparse_check.py` were not re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
