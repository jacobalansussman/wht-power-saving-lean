# The τ-theorem, the Φ barrier and the packing loss

**What it shows.** `frontier_flow.py L` builds the balanced flow of Theorem 1 for a mixed unit at block 3 in exact fractions and
checks that every job size is produced exactly as often as it is consumed; `gamma_bound.py` bounds the constant of the
induction over sampled sizes up to 10^7 (for 1 + 4 in 10 its largest value is 1.2132, at n = 71); `lockstep_sim.py`
simulates the lockstep order of Lemma 2 literally and checks that it never deadlocks; `kron_check.py` shows that
Kronecker powers drift towards exponent 1 (Proposition 4); `macro_lp.py B...` prices a mixed unit packed into one big
word at block B with today's engine (a linear program); `grouped_lp.py` tests the transfer for grouped paid blocks in a
fluid model. `reprice_hub.py` is the what-if on a community profile whose premise (every idle block at the start of
its array) later failed; it is kept as the record of that failure.

**Commands.** `python3 frontier_flow.py 480`; `python3 gamma_bound.py`; `python3 kron_check.py`; `python3 macro_lp.py 3 6 12 24 48 96`; `python3 grouped_lp.py`; `python3 lockstep_sim.py`; `python3 reprice_hub.py`.

**Needs.** Python 3; numpy and scipy for `macro_lp.py` and `grouped_lp.py`.

**Time and memory.** Seconds, except `gamma_bound.py` (about 50 s), `macro_lp.py` (16 to 46 s), `grouped_lp.py` (about 76 s). `lockstep_sim.py`'s fifth case (7 + 38 in 88) grows past 1.2 GB; it was stopped in the original round, and its four finished cases are the shipped output.

**Expected output.** Their `out_*.txt`.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 frontier_flow.py 480` (1.0 s): matches the shipped output. Output: `rerun_2026-10-19_frontier_flow_480.txt`.
- `python3 gamma_bound.py` (48.9 s): matches the shipped output. Output: `rerun_2026-10-19_gamma_bound.txt`.
- `python3 kron_check.py` (0.1 s): matches the shipped output. Output: `rerun_2026-10-19_kron_check.txt`.
- `python3 reprice_hub.py` (0.1 s): matches the shipped output. Output: `rerun_2026-10-19_reprice_hub.txt`.
- `python3 macro_lp.py 3 6 12 24 48 96` (15.7 s): matches the shipped output. Output: `rerun_2026-10-19_macro_lp_3_6_12_24_48_96.txt`.

`lockstep_sim.py` and `grouped_lp.py` were not re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
