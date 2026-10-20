# The star lemma and the commutative slice at blocks 5 and 6 (Theorems A, C, D, E)

**What it shows.** `check_star_lemma.py` and `check_star_ineq.py` test the star lemma exactly; `check_geometry.py` the geometry of
the 16 directions with z1 = z2 = 1 at block 6; `recursion.py` computes Theorem D's thresholds (W ≥ 13 of 16 at block 6,
26 of 32 at block 7, 52 of 64 at block 8, 209 of 256 at block 10); `exact_checks.py` and `circ5.py` the floor of the
local argument (bent-function matrices); `allchars.py` and `allchars2.py` Theorem A (the corrected 2 in 7 at all 512
characters); `nearU6.py` and `nearU6b.py` the near-uniform multisets at block 6; `verify_star_cert.py` checks Theorem E's
exact certificate (`theorem-E-certificate/`, found by the linear program in `starlp2/`); `loc.py`, `lamsolve.py`,
`singer_*` and `lm4*/`, `lam*/` are numerical searches, which were unreliable (a known exact solution was found in 1 of
21 runs). `skeptic/` is the AI skeptic's own code (`v_star.py`, `v_recursion.py`, `v_kron_odd.py`, `v_cert5.py`,
`v_geometry.py`, `v_thmA*.py`, `v_forms6.py`, `v_nearU6.py`, `v_mixed.py`, ...).

**Commands.** `python3 check_star_lemma.py`; `python3 check_star_ineq.py`; `python3 check_geometry.py`; `python3 recursion.py`; `python3 exact_checks.py`; `python3 circ5.py`; `python3 verify_star_cert.py theorem-E-certificate/star_cert_b5.json theorem-E-certificate/star_vars_b5.json`; in `skeptic/`: `python3 v_star.py`, `v_recursion.py`, `v_kron_odd.py`, `v_cert5.py ../theorem-E-certificate`, `v_geometry.py`.

**Needs.** numpy; scipy for the numerical searches.

**Time and memory.** Each check under 10 seconds.

**Expected output.** Their `out_*.txt`.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 check_star_lemma.py` (7.6 s): matches the shipped output. Output: `rerun_2026-10-19_check_star_lemma.txt`.
- `python3 check_star_ineq.py` (6.2 s): matches the shipped output. Output: `rerun_2026-10-19_check_star_ineq.txt`.
- `python3 check_geometry.py` (3.8 s): matches the shipped output. Output: `rerun_2026-10-19_check_geometry.txt`.
- `python3 recursion.py` (0.1 s): matches the shipped output. Output: `rerun_2026-10-19_recursion.txt`.
- `python3 exact_checks.py` (4.8 s): matches the shipped output. Output: `rerun_2026-10-19_exact_checks.txt`.
- `python3 circ5.py` (3.4 s): same rank (8) with zero autocorrelation; the random search found a different truth table. Output: `rerun_2026-10-19_circ5.txt`.
- `python3 verify_star_cert.py theorem-E-certificate/star_cert_b5.json theorem-E-certificate/star_vars_b5.json` (9.8 s): matches the shipped output. Output: `rerun_2026-10-19_verify_star_cert_theorem-E-certificate-star-cert-b5-json_theorem-E-certificate-s.txt`.
- `python3 v_star.py` (in `skeptic/`) (2.8 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_v_star.txt`.
- `python3 v_recursion.py` (in `skeptic/`) (0.1 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_v_recursion.txt`.
- `python3 v_kron_odd.py` (in `skeptic/`) (0.2 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_v_kron_odd.txt`.
- `python3 v_geometry.py` (in `skeptic/`) (0.5 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_v_geometry.txt`.
- `python3 v_cert5.py ../theorem-E-certificate` (in `skeptic/`) (0.1 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_v_cert5_-theorem-E-certificate.txt`.

The numerical searches, `allchars*.py`, `nearU6*.py` and the skeptic's `v_thmA*.py`, `v_forms6.py`, `v_nearU6.py` and `v_mixed.py` were not re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
