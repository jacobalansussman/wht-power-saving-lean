#!/usr/bin/env python3
"""v-cx-medium / vsub2.py (own code, stdlib + numpy; python3 -I).  The even-window sub-designs found by vsub.py that pass
the hub-type screen at one pair: E8's 56-line hyperplane section (k = 8, omega = 4: block 22 = 1 + 3*7) and E7's 32-line
hyperplane section (k = 6, omega = 4: block 16 = 1 + 3*5).  For each: two-angle check, junk per port, sigma (F_2 rank
of a port's junk partners' labels), the clique floor (cx-medium's argument, enumerated by closures as in vclique.py),
the coordinate decoder's cost (orthonormal basis of the span; totals gather the roots not orthogonal to f), and D per
replica at one pair = 2n - 3c."""
import sys, json, itertools
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from vclique import design, gram_B, rowint, add_basis, reduce_, max_cliques, ip

def section(L, want_n, want_om):
    n = len(L); B = gram_B(L); rows = [rowint(B[i]) for i in range(n)]
    basis = {}
    for r in rows: add_basis(basis, r)
    bv = list(basis.values()); k = len(bv)
    piv = {b.bit_length() - 1: (i, b) for i, b in enumerate(bv)}
    def coords(v):
        c = 0
        while v:
            p = v.bit_length() - 1; i, b = piv[p]; v ^= b; c ^= 1 << i
        return c
    C = [coords(r) for r in rows]
    for a in range(1, 1 << k):
        S = [i for i in range(n) if bin(C[i] & a).count('1') % 2 == 0]
        if len(S) != want_n: continue
        sub = [L[i] for i in S]; Bs = gram_B(sub)
        adj = [{j for j in range(len(S)) if j != i and Bs[i][j] == 1} for i in range(len(S))]
        om, cl = max_cliques(len(S), adj)
        if om == want_om: return sub
    raise SystemExit('not found')

def analyse(name, L, cap=3 * 10**6):
    n = len(L); B = gram_B(L); rows = [rowint(B[i]) for i in range(n)]
    bs = {}
    for r in rows: add_basis(bs, r)
    k = len(bs)
    junk = [sum(1 for j in range(n) if B[i][j] == 0) for i in range(n)]
    def rk(S):
        b = {}
        for p in S: add_basis(b, rows[p])
        return len(b)
    sig = [rk([j for j in range(n) if B[i][j] == 0]) for i in range(n)]
    adj = [{j for j in range(n) if j != i and B[i][j] == 1} for i in range(n)]
    om, cl = max_cliques(n, adj)
    # coordinate decoder: orthonormal basis of the span of the lines
    R = np.array([[float(x) for x in l] for l in L])
    U, s, Vt = np.linalg.svd(R, full_matrices=False)
    dimR = int((s > 1e-9).sum())
    best_coord = None
    rng = np.random.default_rng(0)
    # try the coordinate basis restricted to span when it lies in it, else SVD basis; and random rotations
    cands = []
    Q = Vt[:dimR]
    cands.append(('svd', Q))
    E = np.eye(R.shape[1])
    proj = Q.T @ Q
    if np.allclose(proj, np.eye(R.shape[1])): cands.append(('coord', E))
    for t in range(20):
        A, _ = np.linalg.qr(rng.normal(size=(dimR, dimR))); cands.append((f'rand{t}', A @ Q))
    for nm, Fm in cands:
        G = Fm @ R.T
        M = G.T @ G / 2
        assert all(abs(M[a, b] - (1.0 if a == b else 0.0)) < 1e-9 for a in range(n) for b in range(n) if B[a][b] == 1)
        cost = sum(rk([b for b in range(n) if abs(G[q, b]) > 1e-9]) for q in range(Fm.shape[0]))
        if best_coord is None or cost < best_coord[1]: best_coord = (nm, cost)
    # clique floor by closures
    of = [0] * n
    for ci, c in enumerate(cl):
        for p in c: of[p] |= 1 << ci
    def closure(S):
        b = {}
        for p in S: add_basis(b, rows[p])
        return frozenset(p for p in range(n) if reduce_(rows[p], b) == 0), len(b)
    seen = {}; frontier = [frozenset()]; best = None
    while frontier:
        nxt = []
        for S in frontier:
            for p in range(n):
                if p in S: continue
                Cc, d = closure(set(S) | {p})
                if Cc in seen: continue
                m_ = 0
                for q in Cc: m_ |= of[q]
                T = bin(m_).count('1'); seen[Cc] = 1
                if T and (best is None or d / T < best[0]): best = (d / T, d, T, len(Cc))
                nxt.append(Cc)
                if len(seen) > cap: raise SystemExit('cap')
        frontier = nxt
    floor = om * len(cl) * best[0]
    res = dict(n=n, k=k, sigma=sorted(set(sig)), junk_per_port=sorted(set(junk)), omega=om, cliques=len(cl), span_dim=dimR,
               clique_floor=floor, floor_best=best[1:], coord_decoder=best_coord, break_even_P1=2 * n / 3,
               block_P1=1 + 3 * (k - 1), D_at_floor=2 * n - 3 * floor, D_at_coord=2 * n - 3 * best_coord[1], closed_sets=len(seen))
    print(name, json.dumps(res), flush=True)
    return res

if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'E7-32': L = section(design('E7'), 32, 4)
    elif which == 'E8-56': L = section(design('E8'), 56, 4)
    r = analyse(which, L)
    json.dump(dict(r, lines=[[str(x) for x in l] for l in L]), open(__file__.replace('.py', f'_{which}.json'), 'w'), indent=1)
