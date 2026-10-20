# Restriction of scalars: a unit with cube roots of 2 becomes one with rational mixing

**What it shows.** An AI skeptic's toy for the refutation of "a unit's numbers must come from the supplied root's field": a unit
with coefficients in Q(2^(1/3)) (whose Galois closure is not abelian, so these numbers are never available to the
program) is rewritten with each value stored as three rational coordinates. 2 paid moves for 1 array become 6 for 3,
the same ratio, and the restricted word is checked exactly at f = 1, 2, 3.

**Commands.** `python3 weil.py`

**Needs.** Python 3 only.

**Time and memory.** Under a second.

**Expected output.** `out_weil.txt`.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 weil.py` (0.4 s): matches the shipped output. Output: `rerun_2026-10-19_weil.txt`.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
