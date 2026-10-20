#!/usr/bin/env python3
"""tables.py (fp-capacity): what the proved floors allow in packing terms, and what savings 0.1 and 0.01 need.
Floors (moves per array, engine game, rank-1 moves): block 2: 2 [Lean]; 3: 14/5 [Lean + certificate];
4: 41/13 [certificate]; 5: 2.99172; 6: 2.97290; m >= 7: 2.9447 (pf-induction, upheld).
  kappa_m = m / L_m - 1   : the most net saved moves per paid move any unit can have (Delta/M <= kappa_m)
  s_max   = ln(m/L_m)/ln m
  need(s) = m^s - 1       : live ascending packing per paid move a unit with saving s needs (Theorem A(a))
  per array: m - m^(1-s)  : saved moves per array, hence live ascending packing per array
Two-array targets: P* >= 2m - M in every gauge (Theorem A(b))."""
from fractions import Fraction as Fr
from math import log
L = {2: Fr(2), 3: Fr(14, 5), 4: Fr(41, 13), 5: 2.99172, 6: 2.97290}
rows = []
print("| block m | proved floor L_m | s_max | kappa_m = m/L_m - 1 (net packing per move, ceiling) | 0.1 needs per move | 0.1 per array | 0.01 per move | 0.01 per array |")
print("|---|---|---|---|---|---|---|---|")
for m in (2, 3, 4, 5, 6, 8, 16, 25, 41, 73):
    Lm = float(L.get(m, 2.9447))
    kap = m / Lm - 1; smax = log(m / Lm) / log(m)
    n1 = m ** 0.1 - 1; a1 = m - m ** 0.9; n2 = m ** 0.01 - 1; a2 = m - m ** 0.99
    ex1 = "" if n1 <= kap + 1e-12 else " (excluded)"; ex2 = "" if n2 <= kap + 1e-12 else " (excluded)"
    print("| %d | %.5g | %.4f | %.4f | %.4f%s | %.3f | %.4f%s | %.3f |" % (m, Lm, smax, kap, n1, ex1, a1, n2, ex2, a2))
print()
print("| two-array target | block | saving | Delta = 2m - M | irreducible live packing P* needed (every gauge) | live ascending packing per move (raw) |")
print("|---|---|---|---|---|---|")
for (m, M, note) in [(4, 7, ""), (4, "2,2,2,1", "grouped, total rank 7"), (5, 6, ""), (5, 7, ""), (5, 8, ""), (6, 8, ""), (6, 9, ""), (6, 10, "")]:
    if isinstance(M, str):
        Mi = 7; s = log(2 * m / 7) / log(m)   # rank counting
    else:
        Mi = M; s = log(2 * m / M) / log(m)
    print("| 2 in %s %s | %d | %.3f | %d | >= %d | >= %.3f |" % (M, note, m, s, 2 * m - Mi, 2 * m - Mi, (2 * m - Mi) / Mi))
print()
print("| W >= 3 target | block | saving | Delta | raw live ascending packing needed | irreducible packing needed |")
print("|---|---|---|---|---|---|")
for (W, M, m, note) in [(5, 14, 3, ""), (3, 10, 4, ""), (4, 13, 4, ""), (3, 7, 5, " (widened game only: the engine floor 2.99172 per array excludes it)")]:
    s = log(W * m / M) / log(m); D = W * m - M
    print("| %d in %d%s | %d | %.3f | %d | >= %d | >= %d - gamma_m(h) (coding gain of the layout; 0 for one-source/one-sink layouts) |" % (W, M, note, m, s, D, D, D))
