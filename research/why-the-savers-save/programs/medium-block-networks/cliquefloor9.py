#!/usr/bin/env python3
"""cx-medium / cliquefloor9.py (own code, stdlib; python3 -I).  The clique floor of cliquefloor.py for windows up to
k = 9 (the hub's 112-port design, E8), by an exhaustive enumeration of the subspaces W of F_2^k of ONE dimension d
(reduced row echelon forms), with cliques as big-int bitmasks.  For each d it reports Tmax(d) = the largest number of
maximum cliques that the ports with labels in a d-dim W meet, and the floor contribution
    omega * #cliques * d / Tmax(d);
the centre floor per call is the minimum over d.  usage: cliquefloor9.py DESIGN d [d ...]   (E7 reproduces 35)."""
import sys, os, json, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from designs_cx import design, labels

def cliques(n, adj):
    """all maximum cliques (bitmask adjacency), by Bron-Kerbosch with pivoting; returns (omega, list of frozensets)"""
    best = [0]; res = []
    nb = [0] * n
    for i in range(n):
        for j in adj[i]: nb[i] |= 1 << j
    def bk(R, P, X, size):
        if P == 0 and X == 0:
            if size > best[0]: best[0] = size; res.clear()
            if size == best[0]: res.append(R)
            return
        if size + bin(P).count('1') < best[0]: return
        u = (P | X).bit_length() - 1
        cand = P & ~nb[u]
        while cand:
            v = cand.bit_length() - 1; cand &= ~(1 << v)
            bk(R | (1 << v), P & nb[v], X & nb[v], size + 1)
            P &= ~(1 << v); X |= 1 << v
    bk(0, (1 << n) - 1, 0, 0)
    return best[0], res

def rref_subspaces(k, d):
    """yield bases (lists of ints) of every d-dim subspace of F_2^k, once each (RREF: pivots = highest bits)"""
    for piv in itertools.combinations(range(k - 1, -1, -1), d):        # pivot bit positions, descending
        pivset = set(piv)
        free = []                                                        # (row, bit) positions that are free
        for r, p in enumerate(piv):
            free += [(r, b) for b in range(p) if b not in pivset]
        for mask in range(1 << len(free)):
            rows = [1 << p for p in piv]
            for i, (r, b) in enumerate(free):
                if mask >> i & 1: rows[r] |= 1 << b
            yield rows

if __name__ == '__main__':
    name = sys.argv[1]; dims = [int(x) for x in sys.argv[2:]]
    L = design(name); labs, k, B = labels(L); n = len(L)
    adj = [[j for j in range(n) if j != i and B[i][j] == 1] for i in range(n)]
    t0 = time.time()
    om, cl = cliques(n, adj)
    ncl = len(cl)
    portmask = [0] * n
    for ci, c in enumerate(cl):
        m = c
        while m:
            p = m.bit_length() - 1; m &= ~(1 << p)
            portmask[p] |= 1 << ci
    labmask = {}
    for p, u in enumerate(labs): labmask[u] = labmask.get(u, 0) | portmask[p]
    print(f"{name}: n={n} k={k} omega={om} maximum cliques={ncl} ({time.time()-t0:.1f} s)", flush=True)
    out = dict(design=name, n=n, k=k, omega=om, cliques=ncl, dims={})
    for d in dims:
        t1 = time.time(); Tmax = 0; cnt = 0
        for rows in rref_subspaces(k, d):
            cnt += 1
            # all elements of the span
            els = [0]
            for r in rows: els += [x ^ r for x in els]
            m = 0
            for x in els:
                lm = labmask.get(x)
                if lm: m |= lm
            T = bin(m).count('1')
            if T > Tmax: Tmax = T
        fl = om * ncl * d / Tmax if Tmax else float('inf')
        out['dims'][d] = dict(subspaces=cnt, Tmax=Tmax, floor=fl)
        print(f"  d={d}: {cnt} subspaces, Tmax={Tmax} of {ncl} -> floor contribution {fl:.2f} per call ({time.time()-t1:.0f} s)", flush=True)
    out['floor'] = min(v['floor'] for v in out['dims'].values())
    print(f"{name}: centre floor over d in {dims}: {out['floor']:.2f} per call")
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'cliquefloor9_{name}_{"_".join(map(str, dims))}.json'), 'w'), indent=1)
