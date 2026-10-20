# mycert.py CERT.json -- skeptic's own exact check of the 13-row alpha certificate (v-meta2-hermitian-compute).
# Own frames (geo.py: subspaces of F_4^3 as vector sets), own group (closure of generators: elementary matrices,
# diagonal alpha, Frobenius, the BILINEAR correlation s -> s^0, as permutations of the 44 frames; checked to preserve d),
# own row formulas (polymatroid; Ingleton in the form I(A;B) <= I(A;B|C) + I(A;B|D) + I(C;D); DFZ (65) re-typed from
# the local copy of arXiv:1311.4601, Theorem 8.6), each conditioned on T, own canonical forms.
import sys, os, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geo import *
from fractions import Fraction as Fr
import numpy as np

def mat_perm(M):          # M: 3x3 over F_4 (rows); v -> M v
    def app(v): return tuple(fmul(M[r][0], v[0]) ^ fmul(M[r][1], v[1]) ^ fmul(M[r][2], v[2]) for r in range(3))
    return tuple(IDX[frozenset(app(v) for v in SUBS[i])] for i in range(44))
gens = []
for i in range(3):
    for j in range(3):
        if i != j:
            for c in (1, 2):
                M = [[1 if r == k else 0 for k in range(3)] for r in range(3)]; M[i][j] = c
                gens.append(mat_perm(M))
for i in range(3):
    M = [[1 if r == k else 0 for k in range(3)] for r in range(3)]; M[i][i] = 2
    gens.append(mat_perm(M))
gens.append(tuple(IDX[frozenset(tuple(fconj(x) for x in v) for v in SUBS[i])] for i in range(44)))   # Frobenius
gens.append(tuple(bil_perp(i) for i in range(44)))                                                   # correlation
for g in gens:
    assert all(D[g[i]][g[j]] == D[i][j] for i in range(44) for j in range(44))
ident = tuple(range(44))
seen = {ident}; frontier = [ident]
while frontier:
    nf = []
    for h in frontier:
        for g in gens:
            k = tuple(g[h[i]] for i in range(44))
            if k not in seen: seen.add(k); nf.append(k)
    frontier = nf
G = np.array(sorted(seen), dtype=np.int16)
print("group order", len(G))
memo = {}
def canon(X):
    X = tuple(sorted(set(X)))
    if not X: return ()
    if X not in memo:
        im = np.sort(G[:, list(X)], axis=1).astype(np.int64)
        key = np.zeros(len(G), dtype=np.int64)
        for c in range(im.shape[1]): key = key * 44 + im[:, c]
        memo[X] = tuple(int(v) for v in im[int(np.argmin(key))])
    return memo[X]
def F(basis): return frame_of_basis(basis) if basis else IDX[frozenset({(0, 0, 0)})]
total = {}
def addterm(c, X, w):
    k = canon(X)
    if k == (): return
    total[k] = total.get(k, Fr(0)) + w * Fr(c)
def addexpr(expr, w):           # expr: list of (coef, set of frames); adds w * expr
    for c, X in expr: addterm(c, X, w)

def ingleton(A, B, C, D_, T):   # value >= 0
    r = lambda *ss: frozenset().union(T, *ss)
    # I(A;B|C)+I(A;B|D)+I(C;D)-I(A;B) with I(X;Y|Z) = r(XZ)+r(YZ)-r(XYZ)-r(Z)
    e = []
    def I(X, Y, Z=frozenset()):
        e.extend([(1, r(X, Z)), (1, r(Y, Z)), (-1, r(X, Y, Z)), (-1, r(Z))])
    def mI(X, Y, Z=frozenset()):
        e.extend([(-1, r(X, Z)), (-1, r(Y, Z)), (1, r(X, Y, Z)), (1, r(Z))])
    I(A, B, C); I(A, B, D_); I(C, D_); mI(A, B)
    return e
def dfz(Ls, T):                 # RHS - LHS of (65), letters A,B,C,W,X,Y,Z; H(P|Q) = r(PQT) - r(QT), H(P) = r(PT) - r(T)
    A, B, C, W, X, Y, Z = Ls
    r = lambda *ss: frozenset().union(T, *ss)
    e = []
    def H(c, *P): e.extend([(c, r(*P)), (-c, r())])
    def Hc(c, P, *Q): e.extend([(c, r(P, *Q)), (-c, r(*Q))])
    # LHS: 2H(A) + H(B) + 2H(C)
    H(-2, A); H(-1, B); H(-2, C)
    # RHS
    H(1, W); H(1, X); H(1, Y); H(1, Z)
    Hc(2, A, Z, Y); Hc(1, B, X, Z); Hc(2, C, A, X)
    Hc(3, X, W, Y); Hc(3, Z, W, C)
    Hc(5, W, A, B); Hc(5, Y, B, C)
    H(5, A); H(5, B); H(5, C); H(-5, A, B, C)
    return e

cert = json.load(open(sys.argv[1]))
types = {}
for row in cert["rows"]:
    lam = Fr(row["mult"]); assert lam >= 0
    t = row["type"]; types[t] = types.get(t, 0) + 1
    if t == "sub":
        S = frozenset(F(b) for b in row["S"]); a = F(row["a"]); b = F(row["b"])
        assert a not in S and b not in S and a != b
        e = [(1, S | {a}), (1, S | {b}), (-1, S | {a, b}), (-1, S)]
    elif t == "mono":
        S = frozenset(F(b) for b in row["S"]); a = F(row["a"]); assert a not in S
        e = [(1, S | {a}), (-1, S)]
    elif t in ("ing", "dfz"):
        Ls = [frozenset(F(b) for b in l) for l in row["letters"]]; T = frozenset(F(b) for b in row["T"])
        e = ingleton(*Ls, T) if t == "ing" else dfz(Ls, T)
    else:
        raise SystemExit("unknown row type " + t)
    addexpr(e, -lam)            # each row: e >= 0; we add -lam * e
bound = Fr(0)
P1 = [(s, t) for s in range(44) for t in range(44) if D[s][t] == 1]
P3 = [(s, t) for s in range(44) for t in range(44) if D[s][t] == 3]
for q in cert["eqs"]:
    mu = Fr(q["mult"]); dfn = q["def"]
    if dfn[0] == "single":
        a, b = F(dfn[1]), F(dfn[2]); addexpr([(1, {b}), (-1, {a})], mu)
    elif dfn[0] == "pair":
        (s, t), (s0, t0) = [F(x) for x in dfn[1]], [F(x) for x in dfn[2]]
        assert D[s][t] == D[s0][t0] and s != t and s0 != t0
        addexpr([(1, {s, t}), (-1, {t}), (-1, {s0, t0}), (1, {t0})], mu)
    elif dfn[0] == "norm":
        for s, t in P1: addexpr([(1, {s, t}), (-1, {t})], mu / len(P1))
        bound += mu
    else: raise SystemExit(dfn)
for s, t in P3: addexpr([(-1, {s, t}), (1, {t})], Fr(1, len(P3)))
neg = {k: v for k, v in total.items() if v < 0}
print("rows", types, "eqs", len(cert["eqs"]), "|P1|", len(P1), "|P3|", len(P3))
print("orbit variables touched", len(total), "; negative coefficients:", len(neg), list(neg.items())[:4])
print(("VERIFIED" if not neg else "NOT VERIFIED") + f": mean delta_3 <= {bound} * mean delta_1 = {float(bound)}")
# sanity: the F_4 family itself (Y_s = N_s, r = F_4-rank) must VIOLATE the certificate, i.e. the DFZ row must be < 0 on it
def rankF4(X):
    vecs = set()
    for i in X:
        for x in SUBS[i]:
            for z in SUBS[herm_perp(i)]: vecs.add(tuple(x) + tuple(z))
    n = len(vecs); d = 0
    while 4 ** d < n: d += 1
    # vecs is a union of subspaces, not a span: take the span
    S = {(0,) * 6}
    basis = []
    for v in vecs:
        if v not in S:
            S = {vadd(u, vscale(c, v)) for u in S for c in range(4)}
    n = len(S); d = 0
    while 4 ** d < n: d += 1
    return d
for row in cert["rows"]:
    if row["type"] == "dfz":
        Ls = [frozenset(F(b) for b in l) for l in row["letters"]]; T = frozenset(F(b) for b in row["T"])
        val = sum(c * rankF4(X) for c, X in dfz(Ls, T))
        print("DFZ row evaluated on the F_4 family:", val, "(must be < 0: the row is false in characteristic 2)")
