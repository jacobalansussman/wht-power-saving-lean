# mycheck_lift.py LIFTFILE -- independent exact re-check of an alpha-interval lift (v-meta2-hermitian-logic).
# Shares no code with c2-hermitian: F_4 built here (0,1,2=a,3=a^2=a+1), N_s = (s | s^perp_herm), perp w.r.t. sum conj(x_i) z_i.
# Frames identified from the reductions themselves (no interval file used); edges = covering pairs of the subspace lattice.
# Ranks over Q(omega) computed by realification: a Q(omega)-row v gives two Q-rows (v, omega v) in Q^12; Q-rank = 2 * Q(omega)-rank.
import sys, itertools
from fractions import Fraction as Fr
def m4(a, b):
    # GF(4) multiply via log table: 1=a^0, 2=a^1, 3=a^2
    if a == 0 or b == 0: return 0
    lg = {1: 0, 2: 1, 3: 2}; ex = [1, 2, 3]
    return ex[(lg[a] + lg[b]) % 3]
def conj(a): return m4(a, a)
vecs = list(itertools.product(range(4), repeat=3))
def span(B):
    S = set()
    for co in itertools.product(range(4), repeat=len(B)):
        v = [0, 0, 0]
        for c, b in zip(co, B):
            for i in range(3): v[i] ^= m4(c, b[i])
        S.add(tuple(v))
    return frozenset(S)
def hperp(S): return frozenset(z for z in vecs if all(
    (m4(conj(x[0]), z[0]) ^ m4(conj(x[1]), z[1]) ^ m4(conj(x[2]), z[2])) == 0 for x in S))
# all subspaces of F_4^3
subs = set()
for k in range(4):
    for B in itertools.combinations([v for v in vecs if any(v)], k):
        S = span(list(B))
        subs.add(S)
subs = sorted(subs, key=lambda S: (len(S), sorted(S)))
assert len(subs) == 44, len(subs)
dimof = {1: 0, 4: 1, 16: 2, 64: 3}
def red(a, b): return (a % 2) + 2 * (b % 2)   # a + b w -> a + b alpha (digits: 1 -> 1, alpha -> 2, 1+alpha -> 3)
def rowspace6(rows):
    # F_4-span of 6-vectors
    S = {(0,) * 6}
    for r in rows:
        S = {tuple(x ^ m4(c, y) for x, y in zip(s, r)) for s in S for c in range(4)}
    return frozenset(S)
NS = {}
for S in subs:
    P = hperp(S)
    sp = frozenset(tuple(x) + tuple(z) for x in S for z in P)
    NS[sp] = S
L = open(sys.argv[1]).read().split("\n")
nf, ne, best = map(int, L[0].split())
B = []; which = []
for s in range(nf):
    v = list(map(int, L[1 + s].split())); assert len(v) == 36
    M = [[(v[2 * (6 * r + c)], v[2 * (6 * r + c) + 1]) for c in range(6)] for r in range(3)]
    R = [[red(a, b) for a, b in row] for row in M]
    sp = rowspace6(R)
    assert len(sp) == 64, "reduction not of rank 3"
    assert sp in NS, f"frame {s}: reduction is not any N_s"
    B.append(M); which.append(NS[sp])
assert len(set(which)) == 44, "frames not distinct"
def realify(M):
    out = []
    for row in M:
        r1 = []; r2 = []
        for a, b in row:
            r1 += [Fr(a), Fr(b)]               # v
            r2 += [Fr(-b), Fr(a - b)]          # omega * v : (a + b w) w = -b + (a - b) w
        out.append(r1); out.append(r2)
    return out
def qrank(rows):
    rows = [r[:] for r in rows]; rk = 0; n = len(rows[0])
    for c in range(n):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                f = rows[i][c] / rows[rk][c]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk
for s in range(nf):
    assert qrank(realify(B[s])) == 6
idx = {S: i for i, S in enumerate(which)}
edges = [(i, j) for i in range(44) for j in range(44) if which[i] < which[j] and dimof[len(which[j])] == dimof[len(which[i])] + 1]
assert len(edges) == 147, len(edges)
E = 0; bad = []
for i, j in edges:
    r = qrank(realify(B[i] + B[j])); assert r % 2 == 0
    e = r // 2 - 4; assert e >= 0
    E += e
    if e: bad.append((which[i], which[j], e))
# opposite pairs: rank must be 6 (automatic from reduction, re-check)
opp = 0
for i in range(44):
    for j in range(44):
        si, sj = which[i], which[j]
        if i < j and len(si & sj) == 1 and dimof[len(si)] + dimof[len(sj)] == 3:
            opp += 1; assert qrank(realify(B[i] + B[j])) == 12
print("file", sys.argv[1], "says", best, "; recomputed E =", E, "; bad edges", len(bad), "; opposite pairs checked", opp)
def dim(S): return dimof[len(S)]
for s, t, e in bad:
    lvl = (dim(s), dim(t))
    if lvl == (1, 2):
        p = next(v for v in s if any(v))
        absolute = (m4(conj(p[0]), p[0]) ^ m4(conj(p[1]), p[1]) ^ m4(conj(p[2]), p[2])) == 0
        tangent = (t == hperp(s))
        print("  flag", lvl, "excess", e, "p =", p, "absolute:", absolute, "line = p^perp:", tangent)
    else:
        print("  edge", lvl, "excess", e)
print("ratio 3/(1+E/147) =", Fr(3) / (1 + Fr(E, 147)), float(Fr(3) / (1 + Fr(E, 147))))
