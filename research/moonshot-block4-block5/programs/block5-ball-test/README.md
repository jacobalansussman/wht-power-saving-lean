# Block 5: the radius-2 ball systems, and the open s1501 systems

**What it shows.** `vball.py` is the skeptic's own writer of the one-shot ball test: from a design's layered systems it
takes every layer-3 row within r steps of the seed column (no pruning), optionally every layer-4 row (`D4=1`), a
normalisation of the seed, and a_k b_k = 1 for every injection coefficient (so none can be 0). It writes a Singular
script and an msolve input. If 1 lies in the ideal, no choice of coefficients with every injection switched on works.
- `systems/m<design>_k2_d0|d1.sing|.ms` (modulo p = 1,073,741,789) and `q<design>_k2_d0|d1.sing` (over the rationals)
  for the designs s801, s1401 and s1402 at radius 2 (d0 = the ball alone, d1 = with the layer-4 rows);
- `verdicts_msolve.txt` and `verdicts_singular.txt`: every radius-2 system above is inconsistent ([-1] in msolve, "1 in
  ideal" in Singular); the radius-1 controls gave no verdict in 300 seconds, the slow pattern of consistent systems;
- `systems/m1501_k2_d1.*`, `q1501_k2_d1.sing` (radius 2) and `m1501_k3_d1.*`, `q1501_k3_d1.sing` (radius 3) for the one
  undecided design s1501, written for this revision so that anyone can try the exact verdict.

**Commands.** `Singular -q systems/<file>.sing` (Singular 4) and `msolve -t 1 -f systems/<file>.ms -o out.res`.
Writing the systems needs the layered systems as Python pickles, which are not published (pickle files can run code
when loaded); they are rebuilt by `../block5-sparse-designs/dumpsys.py`. Then
`D4=1 python3 vball.py sys_d4_<s>.pkl sys_d3_<s>.pkl <radius> <prime or 0> <out prefix>`.

**Re-run for this revision:** on 19 October 2026 I regenerated all 18 radius-2 systems for s801, s1401 and s1402 from
the skeptic's pickles with `vball.py`; every file is byte-for-byte identical to the one that was decided. Singular and
msolve were not re-run here (not installed on the Mac).
