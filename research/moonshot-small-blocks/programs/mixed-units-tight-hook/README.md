# Mixed units exactly tight on the block-3 hook (Theorem T): the table of first steps

Theorem T (no exact unit exactly tight on the block-3 hook can give any array a two-row target) is a written proof
(argued by hand, not checked; one AI skeptic re-derived it). Its finite input is a table: no single first step D·b
into a two-row array passes the tightness condition. These programs compute that table, twice, with unrelated code.

| File | What it shows |
|---|---|
| `v1_common.py`, `v1_entry.py` | The skeptic's own code (nothing shared with the author's). It builds the hook module V ⊗ Λ⁴V at f = 5 (dimension 560, its subspace N of dimension 280), enumerates the 1,680 classes of single first steps D·b of the widened game structurally (105 perfect matchings times signs), and checks the tightness condition for the two-row target R and for the kernel K. Result: 0 of 1,680 classes pass for R (smallest residual 2.8284 = √8), 105 pass for K. |
| `v2_contphase.py` | The skeptic's check that continuous phases do not help either: over all 105 matchings and all phases (an exact 32⁴ grid by FFT, then local minimisation), the smallest residual is again √8. Needs scipy. |
| `hooklib.py`, `mixlib.py` | The author's hook module (exact 8 × 8 Gaussian-dyadic group elements, arithmetic modulo a prime) and mixed-target programs. |
| `c10_entry_exhaustive.py` | The author's table: 14 engine classes and 1,680 widened classes (found by sampling until the count stopped growing), modulo two primes (argument 0 or 1). |
| `c24_defects.py` | The rank of each single-step failure (the "defect"): at least 8, attained by the plain first step; this is why a unit with hook room 16 can afford two plain-type defects. |
| `c23_open_mixed.py` | The list of mixed targets above the bar that Theorem T leaves open (hook room 16), with their savings. |

Commands (run from this folder; numpy needed, scipy for `v2_contphase.py`):

```
python3 v1_entry.py                  # 1 to 4 minutes, under 100 MB
python3 v2_contphase.py              # minutes
python3 c10_entry_exhaustive.py 0    # under a minute, about 0.4 GB; argument 1 for the second prime
python3 c24_defects.py               # about a minute
python3 c23_open_mixed.py            # seconds
```

Expected output: `../../data/mixed-units-tight-hook/` (`out_v1_entry.txt`: the class counts
`{(False, True, False, True): 105, (False, False, False, True): 1575}` and "smallest float residual ... 2.8284").
Re-runs of 19 October 2026: `rerun-2026-10-19/` there.

Comments inside the programs use the agents' track labels ("ms-mixed" for the author, "v-ms-mixed" for its skeptic)
and refer to their working notes; the section's Sources maps the labels.
