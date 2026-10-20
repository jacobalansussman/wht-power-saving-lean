#!/usr/bin/env python3
"""m=2, phi = id (old equal-frame model), timed convention: an explicit schedule in which ONE total serves the orthogonal
ports u_S=e0, u_T=e1 with all interference cancelled directly by the sources (pre-cancellation), at copy cost 3.
Frames as subspaces U of F_2^2 (subspace frames L_U); distance in G_sub = dim U + dim V - 2 dim(U&V).
twist.py's search only allows one read frame and ONE delivery frame for both targets; it reports phi=identity timed: None."""
import itertools
def span(vs):
    s = {0}
    for v in vs:
        if v not in s: s |= {x ^ v for x in s}
    return frozenset(s)
def dim(U): return len(U).bit_length() - 1
d = lambda U, V: dim(U) + dim(V) - 2 * dim(U & V)
Z, E0, E1, F = span([]), span([1]), span([2]), span([1, 2])
def interval(A, B): return {U for U in {Z, E0, E1, span([3]), F} if d(A, U) + d(U, B) == d(A, B)}
XS, XT = interval(E0, F), interval(E1, F)              # sources <u> -> F
YS, YT = interval(Z, E1), interval(Z, E0)              # targets 0 -> u^perp  (e0^perp = <e1>, e1^perp = <e0>)
# schedule (time order); every register walk is a monotone geodesic of its port interval
sched = [
  ("Y_T moves 0 -> <e0> (end), Y_S moves 0 -> <e1> (end); both empty", None),
  ("Y_T -= X_S  (Y_T at <e0>, X_S at <e0> = its start)", (E0 in YT and E0 in XS)),
  ("Y_S -= X_T  (Y_S at <e1>, X_T at <e1> = its start)", (E1 in YS and E1 in XT)),
  ("X_S moves <e0> -> F, X_T moves <e1> -> F", None),
  ("R := X_S + X_T  (both at F)", (F in XS and F in XT)),
  ("R walks F -> <e0>; Y_T += R (Y_T = x_T)", (d(F, E0) == 1 and E0 in YT)),
  ("R walks <e0> -> 0 -> <e1>; Y_S += R (Y_S = x_S)", (d(E0, Z) == 1 and d(Z, E1) == 1 and E1 in YS)),
]
for s, ok in sched: print(("ok  " if ok else "    ") + s)
assert all(ok is not False for _, ok in sched)
print("copy cost of R = 1 + 2 = 3 (two separate totals: 1 + 1 = 2)")
