#!/usr/bin/env python3
# rebuild_and_field.py : rebuild the graph and local rule of arXiv:2610.10108v1 (Braverman, He) from the two lists in its appendix (see extract_lists.py),
# and measure (a) graph numbers, (b) the dimension of the edge-assignment space W over the field of 9 elements,
# (c) the dimension of the same integer system read over Gaussian numbers modulo a large prime p = 3 mod 4
#     (Z[i]/p is a field; nullity there is an upper bound for the nullity over the Gaussian rationals).
# This is a computation on MY reconstruction of the text. It is not a check of the preprint's theorem.
import re, sys, collections
import numpy as np
import json
Q = sys.argv[1]                          # the JSON file written by extract_lists.py
_d = json.load(open(Q))
R = _d["deleted_edges"]; PI = _d["vertex_order"]
print("parsed |R| =", len(R), "distinct", len(set(R)), "max", max(R), "; |pi| =", len(PI), "is permutation of 0..181:", sorted(PI) == list(range(182)))
R = set(R)
C = {102, 163, 164, 341, 386, 406, 476, 513, 862}; J = {103, 349, 470, 517, 603, 606, 674, 716, 861}
print("C,J subsets of R:", C <= R, J <= R, " |B| =", len(R - C - J))
B = R - C - J

# ---- arithmetic in Z[i]/p as pairs (a,b) = a + b*i ; p = 3 gives the field of 9 elements with theta = i
def mk(p):
    def add(x, y): return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)
    def sub(x, y): return ((x[0] - y[0]) % p, (x[1] - y[1]) % p)
    def mul(x, y): return ((x[0] * y[0] - x[1] * y[1]) % p, (x[0] * y[1] + x[1] * y[0]) % p)
    def inv(x):
        n = pow((x[0] * x[0] + x[1] * x[1]) % p, p - 2, p)
        return ((x[0] * n) % p, (-x[1] * n) % p)
    return add, sub, mul, inv
add3, sub3, mul3, inv3 = mk(3)
Z0 = (0, 0); ONE = (1, 0)
eta = [(r % 3, r // 3) for r in range(9)]
T = [(Z0, Z0, ONE)] + [(Z0, ONE, eta[r]) for r in range(9)] + [(ONE, eta[r], eta[s]) for r in range(9) for s in range(9)]
assert len(T) == 91 and len(set(T)) == 91
def dot3(x, y):
    s = Z0
    for a, b in zip(x, y): s = add3(s, mul3(a, b))
    return s
nbr = [[l for l in range(91) if dot3(T[p], T[l]) == Z0] for p in range(91)]   # symmetric form: same list for a line
assert all(len(n) == 10 for n in nbr)
edges = [(i, nbr[i][r]) for i in range(91) for r in range(10)]                 # e_{10 i + r} = {P_i, L_j}; vertex ids i and 91+j
eidx = {e: q for q, e in enumerate(edges)}
pos = {v: t for t, v in enumerate(PI)}
phys = [q for q in range(910) if q not in R]
print("candidate edges", len(edges), "physical", len(phys))
adjG = collections.defaultdict(list)
for q in phys:
    i, j = edges[q]; adjG[i].append(91 + j); adjG[91 + j].append(i)
deg = collections.Counter(len(adjG[v]) for v in range(182))
print("degrees of the physical graph G (degree: count):", sorted(deg.items()))
# orientation: physical earlier -> later, deleted later -> earlier
din = collections.Counter(); dout = collections.Counter(); cout = collections.Counter(); jin = collections.Counter()
for q, (i, j) in enumerate(edges):
    a, b = i, 91 + j
    early, late = (a, b) if pos[a] < pos[b] else (b, a)
    src, dst = (late, early) if q in R else (early, late)
    dout[src] += 1; din[dst] += 1
    if q in C: cout[src] += 1
    if q in J: jin[dst] += 1
print("in-degree = out-degree = 5 at every vertex:", all(din[v] == 5 and dout[v] == 5 for v in range(182)),
      "; d+_C = d-_J at every vertex:", all(cout[v] == jin[v] for v in range(182)))
def bfs(s):
    d = {s: 0}; fr = [s]
    while fr:
        nx = []
        for u in fr:
            for w in adjG[u]:
                if w not in d: d[w] = d[u] + 1; nx.append(w)
        fr = nx
    return d
dist = [bfs(v) for v in range(182)]
sd = collections.Counter(dist[edges[q][0]].get(91 + edges[q][1], -1) for q in B)
print("G-distance between the two ends of the 157 session edges (distance: count):", sorted(sd.items()), " sum =", sum(k * v for k, v in sd.items()))
print("G connected:", len(dist[0]) == 182, "; diameter:", max(max(d.values()) for d in dist))
allpairs = collections.Counter(x for d in dist for x in d.values())
print("ordered pairs of G by distance:", sorted(allpairs.items()))
print("moves 735 against 157*m: m=3:", 157 * 3, " m=4:", 157 * 4, " m=5:", 157 * 5)

# ---- the local rule, as integer data (Gaussian integers with digits 0,1,2), then read modulo p
def local_coords(z):                                   # over the field of 9 elements, as the text says
    nb = nbr[z]; u, w = T[nb[0]], T[nb[1]]
    out = []
    for y in nb:
        yv = T[y]; found = None
        for lam in eta:
            for mu in eta:
                if all(add3(mul3(lam, u[k]), mul3(mu, w[k])) == yv[k] for k in range(3)): found = (lam, mu)
        assert found is not None
        out.append(found)
    return out
def gpow(x, n):                                        # exact Gaussian integer power
    r = (1, 0)
    for _ in range(n): r = (r[0] * x[0] - r[1] * x[1], r[0] * x[1] + r[1] * x[0])
    return r
def gmul(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def evalrows(z, is_line, balanced):
    rows = []
    for lam, mu in local_coords(z):
        if balanced:
            f = lambda t: tuple(c - 3 if c == 2 else c for c in t)
            lam, mu = f(lam), f(mu)
        if is_line: lam, mu = gpow(lam, 3), gpow(mu, 3)
        rows.append([gmul(gpow(lam, 4 - r), gpow(mu, r)) for r in range(5)])
    return rows                                         # 10 x 5 exact Gaussian integers
def nullity(p, balanced, twist=True):
    add, sub, mul, inv = mk(p)
    Hr = np.zeros((910, 910), dtype=np.int64); Hi = np.zeros((910, 910), dtype=np.int64); row = 0
    for side in (0, 1):
        for z in range(91):
            M = [[(a % p, b % p) for a, b in r] for r in evalrows(z, side == 1 and twist, balanced)]
            A = [r[:] for r in M[:5]]                   # interpolate from the first five incident edges
            # invert A (5x5) by Gauss-Jordan over Z[i]/p
            aug = [A[k] + [ONE if k == c else Z0 for c in range(5)] for k in range(5)]
            for c in range(5):
                pr = next((k for k in range(c, 5) if aug[k][c] != Z0), None)
                if pr is None: return None
                aug[c], aug[pr] = aug[pr], aug[c]; iv = inv(aug[c][c]); aug[c] = [mul(iv, x) for x in aug[c]]
                for k in range(5):
                    if k != c and aug[k][c] != Z0:
                        f = aug[k][c]; aug[k] = [sub(x, mul(f, y)) for x, y in zip(aug[k], aug[c])]
            Ainv = [r[5:] for r in aug]
            qs = [eidx[(z, l)] for l in nbr[z]] if side == 0 else [eidx[(pp, z)] for pp in nbr[z]]
            for jrow in range(5, 10):                   # value on edge jrow = sum_k c_k * value on edge k
                c = [Z0] * 5
                for k in range(5):
                    for t in range(5): c[k] = add(c[k], mul(M[jrow][t], Ainv[t][k]))
                Hr[row, qs[jrow]] = 1
                for k in range(5): Hr[row, qs[k]] = (-c[k][0]) % p; Hi[row, qs[k]] = (-c[k][1]) % p
                row += 1
    # rank over the field Z[i]/p, vectorised elimination
    rank = 0; n = 910
    for col in range(n):
        nz = np.nonzero((Hr[rank:, col] != 0) | (Hi[rank:, col] != 0))[0]
        if len(nz) == 0: continue
        pr = rank + nz[0]
        if pr != rank: Hr[[rank, pr]] = Hr[[pr, rank]]; Hi[[rank, pr]] = Hi[[pr, rank]]
        iv = inv((int(Hr[rank, col]), int(Hi[rank, col])))
        a, b = Hr[rank].copy(), Hi[rank].copy()
        Hr[rank] = (a * iv[0] - b * iv[1]) % p; Hi[rank] = (a * iv[1] + b * iv[0]) % p
        fr, fi = Hr[rank + 1:, col].copy()[:, None], Hi[rank + 1:, col].copy()[:, None]
        pr_, pi_ = Hr[rank][None, :], Hi[rank][None, :]
        Hr[rank + 1:] = (Hr[rank + 1:] - (fr * pr_ - fi * pi_)) % p
        Hi[rank + 1:] = (Hi[rank + 1:] - (fr * pi_ + fi * pr_)) % p
        rank += 1
        if rank == n: break
    return 910 - rank
print("--- dimension of W (edge assignments obeying all 182 local rules); 910 unknowns, 910 conditions")
print("field of 9 elements (p=3), digits 0,1,2      : nullity", nullity(3, False), "  [the preprint states: at least 166]")
print("field of 9 elements (p=3), no cubing at lines: nullity", nullity(3, False, twist=False), "  [my variant, for comparison]")
for p in (10007, 1000003):
    assert p % 4 == 3
    for bal in (False, True):
        print(f"same integer rule read in Z[i]/{p} (a field; upper bound for the Gaussian rationals), digits {'-1,0,1' if bal else '0,1,2'}: nullity", nullity(p, bal))
    print(f"   ... and without the cubing at lines, Z[i]/{p}, digits 0,1,2: nullity", nullity(p, False, twist=False))
