#!/usr/bin/env python3
"""epat_check.py (w8-rmexact): independent re-derivation (no shared code) of the E-pattern exclusion of 4 in 15 with
the order-3 address group sigma = alpha^5 (alpha a root of x^4 + x + 1, polynomial basis), e3-rmbuild's U_4 pattern
(every direction orbit used 3 times; equivalently all 15 directions once).
Prints: the channel targets, the tight forms, the form F that fails, its odd count, the forced zeros of E on supp F,
the maximum matching, and the conclusion.  Pure Python; also brute-forces the matching over all assignments."""
import itertools
m = 4; POLY = 0b10011
idx = [(i, j) for i in range(m) for j in range(i, m)]
def gfmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a >> m: a ^= POLY
    return r
def Qval(Q, x):   # Q as dict of (i,j) -> 1
    return sum(((x >> i) & 1) * ((x >> j) & 1) for (i, j) in Q) & 1
def pairing(Q, S):   # S symmetric as dict (i,j)->bit for i<=j
    return sum(S.get(p, 0) for p in Q) & 1
g = 1
for _ in range(5): g = gfmul(g, 2)            # alpha^5, order 3
H = [1, g, gfmul(g, g)]
assert gfmul(H[2], g) == 1
def gram(h):      # matrix of x -> h*x, columns c_j = h * e_j ; (M M^T)_{ij} = sum_k M_ik M_jk = sum_k c_k[i] c_k[j]
    cols = [gfmul(h, 1 << j) for j in range(m)]
    return {(i, j): sum(((c >> i) & 1) * ((c >> j) & 1) for c in cols) & 1 for (i, j) in idx}
T = [gram(h) for h in H]
print("channel targets h h^T (as upper-triangular bits):", [[T[k][p] for p in idx] for k in range(3)])
forms = [tuple(p for b, p in enumerate(idx) if (F >> b) & 1) for F in range(1, 1 << len(idx))]
W = 4; C = 3          # virtual moves per direction (each direction lies in one orbit, used 3 times)
dirs = list(range(1, 16))
def supp(Q): return [x for x in dirs if Qval(Q, x)]
def nodd(Q): return sum(pairing(Q, S) for S in T)
tight = [Q for Q in forms if nodd(Q) > 0 and C * len(supp(Q)) == W * nodd(Q)]
print("tight forms: %d; support sizes %s; odd counts %s" % (len(tight), sorted({len(supp(Q)) for Q in tight}), sorted({nodd(Q) for Q in tight})))
zero = set()
for Q in tight:
    P = supp(Q)
    for x in P:
        for y in P: zero.add((x, y))
bad = []
for Q in forms:
    P = supp(Q); need = C * len(P) - W * nodd(Q)
    if need < 0: bad.append(('deletion', Q)); continue
    # structural rank of a C|P| x C|P| matrix with zero blocks on 'zero': max b-matching, brute force over the
    # direction-level bipartite graph (|P| <= 12; use Hall / augmenting on the expanded graph)
    rows = [(x, a) for x in P for a in range(C)]
    adj = {r: [(y, b) for y in P for b in range(C) if (r[0], y) not in zero] for r in rows}
    match = {}
    def aug(r, seen):
        for c in adj[r]:
            if c in seen: continue
            seen.add(c)
            if c not in match or aug(match[c], seen): match[c] = r; return True
        return False
    mm_ = sum(1 for r in rows if aug(r, set()))
    if mm_ < need: bad.append(('E-pattern', Q, P, need, mm_))
print("forms violating the tests: %d" % len(bad))
for b in bad[:6]: print("  ", b)
if bad:
    Q, P = bad[0][1], bad[0][2]
    print("first failing form F =", " + ".join("x%d" % (i + 1) if i == j else "x%dx%d" % (i + 1, j + 1) for (i, j) in Q),
          "support", P, "odd at", nodd(Q), "of 3 channel targets")
    print("forced zero blocks E[x,y] = 0 inside supp F (x row, y column):")
    for x in P: print("   ", x, ["0" if (x, y) in zero else "*" for y in P])
    print("tight forms through each pair of supp F:")
    for x, y in itertools.combinations_with_replacement(P, 2):
        ts = [Q2 for Q2 in tight if x in supp(Q2) and y in supp(Q2)]
        print("    (%d,%d): %d tight forms" % (x, y, len(ts)))
