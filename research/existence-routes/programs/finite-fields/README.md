# From solutions modulo a prime to complex units, and exact impossibility certificates

**What it shows.** `zsys.py` is the library (systems in the exported format, arithmetic modulo p, Hensel lifting). `modp_exist.py`,
`special_exist.py`, `padic_exist.py` and `padic_nf.py` find a solution modulo p (special position located from a
Gröbner basis modulo p), lift it p-adically, recognise its numbers as rationals or algebraic numbers of degree up to 6,
and check every equation exactly; only that check counts. `many_primes.py` is the many-primes route (Krick, Pardo and
Sombra) and `kps_table.py` the number of primes it needs. `lift_cert.py make` writes a Singular script that proposes
cofactors g_i, and `lift_cert.py check` verifies 1 = Σ g_i·f_i exactly in Python; `combine_cert.py`, `subst_cert.py`,
`extras_cert.py` and `verify_certs.sh` handle the derived systems; `crosscheck.py` re-reads raw exports; `spurious_scan.py`
refuses stray solutions. `out/` holds the calibration, the 17 relaxation hits and the certificate summaries.
`skeptic/` is the AI skeptic's independent code (`vlib.py`, `vcert.py`, `vall.py`, `vcover.py`, `vpoints.py`, ...);
`vrelaxQ.py`, `vcompare*.py`, `vscreens.py` and `fullsysQ.py` rebuild the systems from the skeleton recipes of another
round and need its builders (`NAG_RERUN_DIR`), which are not part of this section.

**Commands.** Unpack `../../data/finite-fields-systems-and-certificates.tar.gz` here first (see `../../data/README.md`). Then
`python3 kps_table.py`; in `skeptic/`: `python3 vall.py`, `python3 vcover.py`, `python3 vpoints.py`. To propose new
certificates: `python3 lift_cert.py make PREFIX OUT.sing`, run Singular on it, then `python3 lift_cert.py check PREFIX OUT.txt`.

**Needs.** Python 3 for every check. Singular to propose certificates; msolve 0.6.5 for Gröbner bases modulo p (write powers as `x^k`: it reads `x*x` as `x`).

**Time and memory.** `vall.py` about 3 minutes, the others seconds; all under 200 MB.

**Expected output.** `skeptic/out_vall.txt`: 240 certificates OK. `skeptic/out_vcover.txt`: 121 of the 122 exclusion steps of the 8 skeletons certified, the last one a constant system. `skeptic/out_vpoints.txt`: every claimed point exact.

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 kps_table.py` (0.1 s): same table; Python 3.12 adds a warning about a backslash in the docstring. Output: `rerun_2026-10-19_kps_table.txt`.
- `python3 vall.py` (in `skeptic/`) (183.8 s): 240 of 240 certificates OK, as shipped; the 13 FAIL lines are the documented failed attempts and differ only in the path they print. Output: `skeptic/rerun_2026-10-19_vall.txt`.
- `python3 vcover.py` (in `skeptic/`) (0.1 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_vcover.txt`.
- `python3 vpoints.py` (in `skeptic/`) (1.7 s): matches the shipped output. Output: `skeptic/rerun_2026-10-19_vpoints.txt`.

The searches (`modp_exist.py`, `padic_*.py`, `many_primes.py`, the Singular runs) were not re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
