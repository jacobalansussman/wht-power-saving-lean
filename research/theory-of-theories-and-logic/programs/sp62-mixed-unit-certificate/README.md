# Theorem M and the cut enumeration in the 30 irreducibles of Sp(6,2)

Section 4. `verify_exact.py` recomputes the exact ranks of the saved integer matrices of the certificate; `skeptic_e7_check.py` is a skeptic's rebuild from scratch, sharing no code. `exact_cert.py`, `fgcore.py`, `sp6.py`, `fg_build.py`, `fg_cuts.py` and the GAP scripts in `gap/` produced the matrices and the cut lists; `fgcore.py` and `fg_build.py` read `gap/ctbl.json` (the character table exported from GAP 4.12.1 with T. Breuer's library) relative to the current folder, and `fg_build.py` writes a 228 MB cache of the irreducibles that `fg_cuts.py` reads.

Outputs of the original runs: [data/sp62-mixed-unit-certificate/](../../data/sp62-mixed-unit-certificate/).

| program | what it shows | command (from this folder) | time | memory | needs | re-run on 19 Oct 2026 |
|---|---|---|---|---|---|---|
| `verify_exact.py` | Theorem M: exact rational ranks of the saved certificate matrices | `python3 verify_exact.py ../../data/sp62-mixed-unit-certificate/exact_mats_lam3.json` | 5 s | under 100 MB | Python 3 only | identical |
| `skeptic_e7_check.py` | Theorem M rebuilt independently through W(E7) | `python3 skeptic_e7_check.py 1` | 8 s | 50 MB | numpy | identical |
| `exact_cert.py` | producer of the certificate matrices (with fgcore.py and gap/ctbl.json) | `python3 exact_cert.py 3 8 28` | not measured | not measured | numpy | not re-run |
| `fg_build.py` | builds the irreducibles of Sp(6,2) (writes a 228 MB cache) | `python3 fg_build.py irreps.npz` | not measured | several GB | numpy | not re-run |
| `fg_cuts.py` | all cuts of the irreducibles (24,609) | `python3 fg_cuts.py irreps.npz out.json` | not measured | several GB | numpy | not re-run |

Other files here: `dpg_eig.py`, `fgcore.py`, `gap`, `sp6.py`.

Programs that need numpy (and scipy, mpmath or pyzx where listed) were run with
`uv run --no-project --with numpy python PROGRAM ...`; any Python 3.9+ with those packages works.
