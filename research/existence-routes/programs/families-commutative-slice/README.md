# Families of units and the commutative slice (Theorems 10, 17, 18)

**What it shows.** `refl.py` and `refl2.py` are the libraries (binary quadratic forms, the slice's point evaluations, a complex
Levenberg-Marquardt solver). `t_theorem10.py` tests Theorem 10 (the matrix E); `even_test.py`, `cover_test.py` and
`clique_test.py` are the exact tests; `rc_lp*.py` the linear programs; `exact_cert.py` turns an optimal dual into an
exact certificate; `verify_cert_b4.py` with `rc_cert_b4.json` and `verify_cert_b5.py` with `rc_cert_b5.json` and
`rc_vars_b5.json` check the certificates of Theorems 17 (M ≥ 4W at block 4) and 18 (M ≥ 155W/36 at block 5) in exact
fractions; `hand_b4.py` and `hand_b5.py` check the hand proofs; `box_enum*`, `cross_*`, `n1*`, `oddfamily.py`,
`piece_ilp.py` and `lemma11/` are the searches and small lemmas behind the other statements of the proof file.
`skeptic/` is the AI skeptic's own code (`my_rclp.py`, `my_verify_their_b5.py`, `my_hand_b4.py`, `my_thm10_check.py`,
...). Compiled binaries are not shipped; build them with `cc -O2 -o box_enum2 box_enum2.c` and the same in `lemma11/`.

**Commands.** `python3 verify_cert_b4.py`; `python3 verify_cert_b5.py`; `python3 hand_b4.py`; `python3 hand_b5.py`; in `skeptic/`: `python3 my_verify_their_b5.py`, `python3 my_hand_b4.py`.

**Needs.** numpy and scipy for the linear programs and searches; the certificate checks need Python 3 only.

**Time and memory.** Each check under 10 seconds.

**Expected output.** Their `out_*.txt`.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 verify_cert_b4.py` (0.2 s): matches the shipped output. Output: `rerun_2026-10-19_verify_cert_b4.txt`.
- `python3 verify_cert_b5.py` (7.2 s): matches the shipped output. Output: `rerun_2026-10-19_verify_cert_b5.txt`.
- `python3 hand_b4.py` (0.2 s): same results; the re-run prints four lines of intermediate output that the shipped file lacks. Output: `rerun_2026-10-19_hand_b4.txt`.
- `python3 hand_b5.py` (9.3 s): matches the shipped output. Output: `rerun_2026-10-19_hand_b5.txt`.
- `python3 my_hand_b4.py` (in `skeptic/`) (0.1 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_my_hand_b4.txt`.
- `python3 my_verify_their_b5.py` (in `skeptic/`) (0.2 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_my_verify_their_b5.txt`.

The linear programs and searches were not re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
