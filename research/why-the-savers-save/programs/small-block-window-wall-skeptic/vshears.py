#!/usr/bin/env python3
"""v-an2-smallnet / vshears.py (own code, stdlib; python3 -I).  Independent check of an2-smallnet NOTE 2.5:
banks hold vectors (rows) over Q; a call is 'bank i += c * bank j' (i != j), c in a coefficient set; a bank's
participation counts once per call it is in (source or target).  Goal: every bank ends as a nonzero multiple of a
DIFFERENT bank's original data (a fixed-point-free scaled permutation).  DFS over call sequences with a per-bank
participation cap; reports whether any exchange exists.  Coefficients {+-1, +-2, +-1/2} (wider than the claim's +-1)."""
import sys, json
from fractions import Fraction as Fr

def search(B, cap, coeffs, maxcalls):
    start = tuple(tuple(Fr(1) if i == j else Fr(0) for j in range(B)) for i in range(B))
    def goal(st):
        used = set()
        for i, row in enumerate(st):
            nz = [j for j, x in enumerate(row) if x != 0]
            if len(nz) != 1 or nz[0] == i or nz[0] in used: return False
            used.add(nz[0])
        return True
    seen = set()
    def dfs(st, part, ncalls):
        if goal(st): return True
        if ncalls == maxcalls: return False
        key = (st, part)
        if key in seen: return False
        seen.add(key)
        for i in range(B):
            if part[i] >= cap: continue
            for j in range(B):
                if j == i or part[j] >= cap: continue
                for c in coeffs:
                    row = tuple(a + c * b for a, b in zip(st[i], st[j]))
                    if all(x == 0 for x in row): continue
                    ns = st[:i] + (row,) + st[i+1:]
                    np_ = list(part); np_[i] += 1; np_[j] += 1
                    if dfs(ns, tuple(np_), ncalls + 1): return True
        return False
    return dfs(start, tuple([0]*B), 0)

if __name__ == '__main__':
    coeffs = [Fr(int(c)) if "/" not in c else Fr(c) for c in sys.argv[2].split(",")]
    out = {}
    for B, cap in [tuple(map(int, p.split("x"))) for p in sys.argv[3].split(",")]:
        r = search(B, cap, coeffs, (B * cap) // 2)
        print(f"banks {B}, at most {cap} calls per bank: exchange exists = {r}", flush=True)
        out[f"{B}|{cap}"] = r
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
