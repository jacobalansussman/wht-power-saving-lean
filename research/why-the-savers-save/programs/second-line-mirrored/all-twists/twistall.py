#!/usr/bin/env python3
"""cx-mirror / twistall.py (own code, stdlib).  EVERY twisted mirror, exhaustively.

qcheck.py: a twisted mirror v_b = s + j + A u_b (A orthogonal, u_b . A u_b = 1 on all ports) is a symmetric involution
A = I + M with M symmetric, zero diagonal, M^2 = 0, M j = 0; im M is a totally isotropic subspace of the even space E
(dim 0, 2 or 4).  So: M = 0 (the natural mirror); one M per isotropic plane P = <w, z> (M = w z^T + z w^T); and 28 per
Lagrangian L (M = B c B^T, c a nondegenerate alternating 4x4 matrix, B a basis of L).  Expected 1 + 5355 + 2295 * 28.
For each: N[a][b] = 1 + u_a . u_b + u_b . M u_a (a != b, symmetric); a drop-free x chain inside a support is an
independent set of {N = 1}; exact maximum by bitset branch and bound with a greedy colouring bound.
usage: python3 twistall.py PART NPARTS NPROC OUT.json"""
import json, sys, os, itertools
from multiprocessing import Pool
here = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(here, 'data.json'))); U = D['ports']; SK = D['supports']; v = len(U)
dot = lambda a, b: bin(a & b).count('1') & 1
E = [e for e in range(1, 512) if bin(e).count('1') % 2 == 0]
def span(vs):
    s = {0}
    for x in vs: s |= {y ^ x for y in s}
    return frozenset(s)
def planes():
    seen = set(); out = []
    for i, w in enumerate(E):
        for z in E[i + 1:]:
            if dot(w, z): continue
            P = span([w, z])
            if P in seen: continue
            seen.add(P); out.append((w, z))
    return out
def lagrangians():
    seen = set(); out = []
    iso3 = set()
    for (w, z) in planes():
        P = span([w, z])
        for x in E:
            if x in P or dot(x, w) or dot(x, z): continue
            T = span([w, z, x])
            if T in iso3: continue
            iso3.add(T)
            for y in E:
                if y in T or any(dot(y, t) for t in T if t): continue
                L = span([w, z, x, y])
                if L in seen: continue
                seen.add(L); out.append((w, z, x, y))
    return out
def Mvec(cols):
    """M as column images: M e_i, from a list of (p, q) pairs: M = sum p q^T + q p^T"""
    return cols
def build_M(pairs):
    col = [0] * 9
    for p, q in pairs:
        for i in range(9):
            c = 0
            if (q >> i) & 1: c ^= p
            if (p >> i) & 1: c ^= q
            col[i] ^= c
    return col
def applyM(col, x):
    y = 0
    for i in range(9):
        if (x >> i) & 1: y ^= col[i]
    return y
def alt4():
    out = []
    idx = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for bits in range(64):
        c = [[0] * 4 for _ in range(4)]
        for k, (i, j) in enumerate(idx):
            if (bits >> k) & 1: c[i][j] = c[j][i] = 1
        # Pfaffian over F_2: c01 c23 + c02 c13 + c03 c12
        if (c[0][1] & c[2][3]) ^ (c[0][2] & c[1][3]) ^ (c[0][3] & c[1][2]):
            out.append([(i, j) for (i, j) in idx if c[i][j]])
    return out
def all_M():
    Ms = [[0] * 9]
    for w, z in planes(): Ms.append(build_M([(w, z)]))
    A4 = alt4()
    for L in lagrangians():
        for c in A4: Ms.append(build_M([(L[i], L[j]) for i, j in c]))
    return Ms
def mis(nodes, adjb):
    """maximum independent set size in the subgraph induced by `nodes` (bitmask), adjacency bitmasks adjb"""
    best = 0
    def colour_bound(P):
        # greedy colouring of the complement graph's cliques = independent sets: number of independent sets
        # needed to cover P bounds an independent set from above? No: use the clique cover of P by cliques of G
        # (an independent set meets each clique once).
        k = 0; Q = P
        while Q:
            k += 1; x = (Q & -Q).bit_length() - 1; cl = 1 << x; cand = Q & adjb[x]
            while cand:
                y = (cand & -cand).bit_length() - 1; cl |= 1 << y; cand &= adjb[y]
            Q &= ~cl
        return k
    def bb(size, P):
        nonlocal best
        if not P:
            best = max(best, size); return
        if size + bin(P).count('1') <= best: return
        if size + colour_bound(P) <= best: return
        x = (P & -P).bit_length() - 1
        bb(size + 1, P & ~adjb[x] & ~(1 << x))
        bb(size, P & ~(1 << x))
    bb(0, nodes); return best
def chains(col):
    adjb = [0] * v
    for a in range(v):
        Ma = applyM(col, U[a])
        for b in range(v):
            if a != b and (1 ^ dot(U[a], U[b]) ^ dot(U[b], Ma)):
                adjb[a] |= 1 << b
    return [mis(sum(1 << x for x in S), adjb) for S in SK]
if __name__ == '__main__':
    part, nparts, nproc, outp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    Ms = all_M()
    mine = Ms[part::nparts]
    with Pool(nproc) as pool:
        res = pool.map(chains, mine, chunksize=64)
    hist = {}; best = None
    for col, c in zip(mine, res):
        hist[min(c)] = hist.get(min(c), 0) + 1
        if best is None or min(c) > min(best[1]): best = (col, c)
    maxany = max(max(c) for c in res)
    json.dump(dict(total_M=len(Ms), part=part, nparts=nparts, done=len(mine), worst_support_hist=hist,
                   best=dict(M=best[0], chains=best[1]), best_chain_any_support=maxany), open(outp, 'w'), indent=1)
    print('M total', len(Ms), '| this part', len(mine), '| worst-support chain histogram', dict(sorted(hist.items())),
          '| best', best[1], '| longest chain in any support', maxany)
