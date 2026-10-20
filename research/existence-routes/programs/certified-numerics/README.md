# Certified numerics, and the exact block-2 stand-in schemes

**What it shows.** `certnum.py` is the library: polynomial systems with Gaussian-rational coefficients, a damped Gauss-Newton search
(not rigorous), Krawczyk's interval test on a square subsystem, rigorous box exclusion, and exact refusal certificates
(1 = Σ g_i·f_i of low degree). `selftest.py` checks it on systems with known answers. `cal1.py` and `calS.py` are the
calibration on exported relaxations (`sys/*.json`); `bulk.py`, `bulksum.py`, `runbulk.sh`, `bulk2.sh` ran the bulk
refusals on the exports of another round (not shipped; their results are `out/bulk*.jsonl`). `sx/` holds the stand-in
construction (`construct.py`), the exact schemes `sx/out/scheme_2in3.json` and `scheme_3in5.json`, and their second
checker `verify_scheme.py`. `v56/` is the 56-dimensional relaxation at block 3; `gct/` the capacity and stability
programs of the relaxation pieces. `skeptic/` is the AI skeptic's code (`vstandin.py` checks the stored schemes
exactly; `vgame.py`, `vconstruct.py`, `vbulk.py`, `vnulla*.py`, `vpoints.py`; `pieces/` rebuilds the piece
decompositions in two independent models).

**Commands.** `python3 selftest.py`; in `sx/`: `python3 verify_scheme.py out/scheme_2in3.json` and `... out/scheme_3in5.json`; in `skeptic/`: `python3 vstandin.py ../sx/out/scheme_2in3.json` and the same for 3 in 5.

**Needs.** python-flint 0.9 (`pip install python-flint==0.9.0`) and numpy; sympy for `sx/construct.py`.

**Time and memory.** Each check under a minute and under 1 GB. The constructions and the bulk runs took minutes to an hour on a Colab machine.

**Expected output.** `sx/out/verify_scheme_2in3.txt`, `sx/out/out_verify_3in5.txt` (INDEPENDENT EXACT CHECK: PASS) and `skeptic/out_vstandin_*.txt` (VERDICT PASS).

**Re-run for this revision.**
On 19 October 2026, on a Linux machine of Google's Colab service (8 processors, about 50 GB; Python 3.12, numpy 2.1) unless the line says otherwise:
- `python3 selftest.py`: every line as expected. Output: `rerun_2026-10-19_selftest.txt`.
- `python3 vstandin.py ../sx/out/scheme_2in3.json` (in `skeptic/`): VERDICT PASS, identical apart from a timing line. Output: `skeptic/rerun_2026-10-19_vstandin_2in3.txt`.
- `python3 vstandin.py ../sx/out/scheme_3in5.json` (in `skeptic/`): VERDICT PASS, identical apart from a timing line. Output: `skeptic/rerun_2026-10-19_vstandin_3in5.txt`.
- `python3 verify_scheme.py out/scheme_2in3.json` (in `sx/`): PASS, identical to the shipped output. Output: `sx/rerun_2026-10-19_verify_scheme_2in3.txt`.
- `python3 verify_scheme.py out/scheme_3in5.json` (in `sx/`): PASS, identical. Output: `sx/rerun_2026-10-19_verify_scheme_3in5.txt`.

The constructions, the calibration runs, the bulk runs and `v56/` were not re-run for this revision.

Every program here is the program as it ran in October 2026, with private paths removed; working names inside the files are mapped in `../../proofs/README.md`.
