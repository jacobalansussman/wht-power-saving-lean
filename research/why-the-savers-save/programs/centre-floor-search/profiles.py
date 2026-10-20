#!/usr/bin/env python3
"""ctr-search / profiles.py: dimension profiles that Problem D (O1, one-sided charge) allows at cost <= 63.
A decoder with m directions has label dims d_1 <= ... <= d_m (cost sum d_j).  Necessary (NC_k): for every k-set K of
directions, the ports usable only with K have root rank <= k, so the other m - k directions' port sets cover the
complement of a rank-k flat of the 112 lines.  Upper bounds on what a dim-d direction covers:
  n(d) = 1,2,4,8,12,20,36,64,112 (max ports in a d-dim label subspace, nd.py);
  k = 7: co-flats Gamma_a (53) / N_i (56) with per-dim coverage bounds from coflat.py, other co-flats >= 70 (D7 bound);
  k < 7: rank-k flats hold <= 1,3,6,12,20,36 lines (A1,A2,A3,D4,D5,E6 in E8), so co-flats >= 112 - h(k).
The binding sets for NC_k are the m - k smallest dims.  Lists every profile (m <= MMAX) passing all bounds."""
import sys, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
n = {1: 1, 2: 2, 3: 4, 4: 8, 5: 12, 6: 20, 7: 36, 8: 64, 9: 112}
cf = json.load(open(os.path.join(HERE, 'coflat.json')))
mG = {int(k): v for k, v in cf['Gamma_a'].items()}; mG[9] = 53
mN = {int(k): v for k, v in cf['N_i'].items()}; mN[9] = 56
h = {0: 0, 1: 1, 2: 3, 3: 6, 4: 12, 5: 20, 6: 36}
INC = True
def ok_prefix(ds, m):
    """ds sorted ascending, complete profile of length m"""
    # incidence: every port is usable with >= 1 direction, and with exactly one only if that direction is parallel to
    # its root (one port per direction), so sum_k |S_W_k| >= 2*112 - m
    if INC and sum(n[d] for d in ds) < 224 - m: return False
    for k in range(0, 8):
        J = ds[:m - k]
        if not J: return False
        if k == 7:
            if not (sum(min(mG[d], 53) for d in J) >= 53 or sum(min(mN[d], 56) for d in J) >= 56 or sum(n[d] for d in J) >= 70):
                return False
        else:
            if sum(n[d] for d in J) < 112 - h[k]: return False
    return True
def gen(m, budget=63):
    out = []
    def rec(ds, lo, rem, left):
        if left == 0:
            if ok_prefix(ds, m): out.append(tuple(ds))
            return
        for d in range(lo, 10):
            if d * left > rem: break
            rec(ds + [d], d, rem - d, left - 1)
    rec([], 1, budget, m)
    return out
if __name__ == '__main__':
    MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    if len(sys.argv) > 2 and sys.argv[2] == 'noinc': INC = False
    res = {}
    for m in range(1, MMAX + 1):
        P = gen(m)
        res[m] = P
        print(f'm = {m}: {len(P)} profiles pass the counting bounds' + (f'; e.g. {P[:6]}' if P else ''), flush=True)
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(HERE, 'profiles.json'), 'w'))
