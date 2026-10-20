#!/usr/bin/env python3
"""v-cx-medium / vclique.py (own code, stdlib; python3 -I).  Independent check of cx-medium's clique floor.
Differences from cx-medium: roots built in plain (undoubled, half-integer via Fractions) coordinates; E7 taken as the
E8 roots orthogonal to the root (1/2,...,1/2) (cx-medium used e1+e2); NO labels are constructed: the F_2 rank of a set
S of labels equals rank_F2(B[S,:]) (B = U U^T with U of full column rank, so U^T is onto), and the subspaces W that
matter are spans of labels, i.e. closed sets of ports {b : row B_b in rowspan B[S,:]}.  Enumerated by BFS from the
empty set (adding one port and closing).  floor = omega * #maxcliques * min_W dim(W)/T(W)."""
import itertools, sys, json
from fractions import Fraction as Fr

def lines(R):
    out, seen = [], set()
    for r in R:
        t = tuple(r); nt = tuple(-x for x in r)
        key = max(t, nt)
        if key not in seen: seen.add(key); out.append(key)
    return out

def E8():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for a in (1, -1):
            for b in (1, -1):
                v = [Fr(0)] * 8; v[i] = Fr(a); v[j] = Fr(b); R.append(tuple(v))
    for s in itertools.product((1, -1), repeat=8):
        if s.count(-1) % 2 == 0: R.append(tuple(Fr(x, 2) for x in s))
    return R

def ip(a, b): return sum(x * y for x, y in zip(a, b))

def design(name):
    if name == 'E8': return lines(E8())
    if name == 'E7':
        h = tuple([Fr(1, 2)] * 8)
        return lines([r for r in E8() if ip(r, h) == 0])
    if name[0] == 'D':
        n = int(name[1:]); R = []
        for i, j in itertools.combinations(range(n), 2):
            for a in (1, -1):
                for b in (1, -1):
                    v = [Fr(0)] * n; v[i] = Fr(a); v[j] = Fr(b); R.append(tuple(v))
        return lines(R)
    if name[0] == 'A':
        n = int(name[1:]); R = []
        for i, j in itertools.permutations(range(n + 1), 2):
            v = [Fr(0)] * (n + 1); v[i] = Fr(1); v[j] = Fr(-1); R.append(tuple(v))
        return lines(R)

def gram_B(L):
    n = len(L)
    for l in L: assert ip(l, l) == 2
    B = [[(1 + int(ip(L[i], L[j]))) % 2 for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n): assert ip(L[i], L[j]) in (0, 1, -1, 2, -2)
    return B

def rowint(r): return sum(b << j for j, b in enumerate(r))

def reduce_(v, basis):            # basis: dict pivot -> vector
    while v:
        p = v.bit_length() - 1
        if p in basis: v ^= basis[p]
        else: return v
    return 0

def add_basis(basis, v):
    v = reduce_(v, basis)
    if v: basis[v.bit_length() - 1] = v
    return v != 0

def max_cliques(n, adj):
    best = [0, []]
    def ext(R, P):
        if not P:
            if len(R) > best[0]: best[0], best[1] = len(R), [tuple(R)]
            elif len(R) == best[0]: best[1].append(tuple(R))
            return
        for v in sorted(P):
            ext(R + [v], {u for u in P if u > v and u in adj[v]})
    ext([], set(range(n)))
    return best[0], best[1]

def run(name, cap=10**7):
    L = design(name); n = len(L); B = gram_B(L)
    rows = [rowint(B[i]) for i in range(n)]
    basis = {}
    for r in rows: add_basis(basis, r)
    k = len(basis)
    adj = [{j for j in range(n) if j != i and B[i][j] == 1} for i in range(n)]
    om, cl = max_cliques(n, adj)
    of = [0] * n
    for ci, c in enumerate(cl):
        for p in c: of[p] |= 1 << ci
    def closure(S):
        bs = {}
        for p in S: add_basis(bs, rows[p])
        return frozenset(p for p in range(n) if reduce_(rows[p], bs) == 0), len(bs)
    seen = {}
    frontier = [frozenset()]
    best = None
    while frontier:
        nxt = []
        for S in frontier:
            for p in range(n):
                if p in S: continue
                C, d = closure(set(S) | {p})
                if C in seen: continue
                T = 0; m = 0
                for q in C: m |= of[q]
                T = bin(m).count('1')
                seen[C] = (d, T)
                if T and (best is None or Fr(d, T) < best[0]): best = (Fr(d, T), d, T, len(C))
                nxt.append(C)
                if len(seen) > cap: raise SystemExit('cap')
        frontier = nxt
    floor = om * len(cl) * best[0]
    res = dict(n=n, k=k, omega=om, cliques=len(cl), closed_sets=len(seen), best_dim=best[1], best_T=best[2],
               best_ports=best[3], floor_per_call=str(floor), floor_float=float(floor), per_port=float(floor / n),
               distinct_labels=len(set(rows)) == n)
    print(name, res, flush=True)
    return res

if __name__ == '__main__':
    out = {nm: run(nm) for nm in (sys.argv[1:] or ['A4', 'A5', 'D5', 'D6', 'E7'])}
    json.dump(out, open(__file__.replace('.py', '.json'), 'w'), indent=1)
