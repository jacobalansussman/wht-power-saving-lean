#!/usr/bin/env python3
"""v_geometry.py (v-e3-slice6): in AG(n,2), the (d+1)-flats through a d-flat Psi number 2^(n-d)-1 and partition
AG(n,2) minus Psi into pieces of 2^d points. Checked for every d-flat, n = 2..5."""
import itertools
def flats(n, d):
    N = 1 << n; subs = {frozenset([0])}
    for _ in range(d):
        subs = {frozenset(S | {s ^ b for s in S}) for S in subs for b in range(1, N) if b not in S}
    return {frozenset(a ^ s for s in S) for S in subs for a in range(N)}
for n in range(2, 6):
    FL = {d: flats(n, d) for d in range(0, n + 1)}
    bad = 0; cnt = 0
    for d in range(0, n - 1):
        for Psi in FL[d]:
            ups = [G for G in FL[d + 1] if Psi <= G]
            pieces = [G - Psi for G in ups]
            union = set().union(*pieces)
            ok = (len(ups) == 2 ** (n - d) - 1 and all(len(p) == 2 ** d for p in pieces)
                  and sum(len(p) for p in pieces) == len(union) == (1 << n) - len(Psi) and not (union & Psi))
            bad += not ok; cnt += 1
    print("n=%d: d-flats checked (d <= n-2): %d, failures %d" % (n, cnt, bad))
