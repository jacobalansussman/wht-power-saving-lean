# The no-phase game as a group algebra, and the uniform family with address maps

**What it shows.** `phi_check.py` checks Theorem A's lift on real engine words at block 3 with the exact formal-content engine
`hb.py` (copied from the theory guide of an earlier round); `raw_check.py` checks lifted words on raw arrays at f = 1 to 6;
`t_lemmaR.py` checks Lemma R; `vclique.py`, `vrc.py` and `vtest5.py` are the exact tests of Lemma V (deletion, tight
cliques, rank cover) at blocks 4 and 5; `orbits_gl.py` and `tr_sym.py` the orbits of quadratic forms; `tf.py`, `lm.py`
and `sec.py` the numerical searches for 4 in 15 (none came near zero). `skeptic/` is the AI skeptic's code:
`my_conj_sweep.py` (every conjugate of the address groups), `my_cert_check.py` (the clique certificates, pair by pair),
`my_vtests.py`, `my_b5_conj.py`, `my_raw.py`, `my_orbits.py`, `my_gen_test.py`.

**Commands.** `python3 phi_check.py`; `python3 raw_check.py`; `python3 vclique.py ...`; in `skeptic/`: `python3 my_cert_check.py 6,15,8,12`, `python3 my_conj_sweep.py`.

**Needs.** numpy; torch for the numerical searches.

**Time and memory.** The exact checks take seconds.

**Expected output.** Their `out_*.txt`.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 phi_check.py` (0.7 s): matches the shipped output. Output: `rerun_2026-10-19_phi_check.txt`.
- `python3 raw_check.py` (4.6 s): T1 to T7 exact at f = 1 to 6 with errors of order 1e-16, as shipped; the two control lines at the end of the shipped file came from a separate one-line command that is not kept. Output: `rerun_2026-10-19_raw_check.txt`.
- `python3 my_cert_check.py 6,15,8,12` (in `skeptic/`) (0.1 s): reproduces the first certificate of the shipped file (|Y| = 27 > 12); the shipped file also holds runs with other address maps. Output: `skeptic/rerun_2026-10-19_my_cert_check_6-15-8-12.txt`.

The tests of Lemma V, the conjugate sweeps and the numerical searches were not re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
