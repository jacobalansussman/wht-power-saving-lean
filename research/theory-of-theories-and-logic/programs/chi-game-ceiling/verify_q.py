# verify_q.py CERTIFICATE.json -- independent exact check of an F_4 (alpha) interval-LP dual certificate (meta2-hermitian).
# Shares no code with hq.py / cert_q.py / lpm.py.  Rebuilds PG(2,4) from scratch: points = normalised nonzero vectors
# of F_4^3 (F_4 = F_2[x]/(x^2+x+1), digits 0, 1, 2 = x, 3 = x+1), frames = subspaces as SETS OF POINTS (L0 = {}, M0 = all),
# d(s,t) = dim(join) - dim(meet); the symmetry group = all invertible 3x3 matrices over F_4, with or without the Frobenius,
# with or without the Hermitian polarity (every element preserves d; Lemma Sym makes averaging over it free); canonical
# form of a set of frames = lexicographically least sorted image.  Rows rebuilt from the explicit frames of the certificate:
#   sub, mono (polymatroid); ing (Ingleton); dfz (Dougherty-Freiling-Zeger arXiv:1311.4601 Thm 8.6, (65));
#   dfzd (its dual by Lemma D*: r*(S) = sum_{i in S} r(i) - r(E) + r(E \ S)), each conditioned on T.
# Equalities: r({s}) = r({s'}); r({s,t}) - r({t}) depends only on d(s,t); normalisation mean_{d=1} = 1.
# Checks  sum lam_i(-g_i) + sum mu_e eq_e - c >= 0 coefficientwise (r >= 0), lam >= 0, c = mean_{d=3}(r({s,t}) - r({t})).
import sys, json, itertools
from fractions import Fraction as Fr
import numpy as np

def gmul(a, b):                       # carry-less product mod x^2 + x + 1
    p = 0
    for i in range(2):
        if (b >> i) & 1: p ^= a << i
    if p & 4: p ^= 0b111
    return p
INV = {1: 1, 2: 3, 3: 2}
def norm(v):
    for c in v:
        if c: iv = INV[c]; return tuple(gmul(iv, x) for x in v)
    return None
vecs = [v for v in itertools.product(range(4), repeat=3) if any(v)]
POINTS = sorted({norm(v) for v in vecs}); assert len(POINTS) == 21
def span_pts(gens):
    out = set()
    for co in itertools.product(range(4), repeat=len(gens)):
        v = [0, 0, 0]
        for c, g in zip(co, gens):
            for i in range(3): v[i] ^= gmul(c, g[i])
        if any(v): out.add(norm(tuple(v)))
    return frozenset(out)
LINES = sorted({span_pts([p, r]) for p in POINTS for r in POINTS if p != r}, key=sorted); assert len(LINES) == 21
FRAMES = [frozenset()] + [frozenset({p}) for p in POINTS] + LINES + [frozenset(POINTS)]
where = {f: i for i, f in enumerate(FRAMES)}; NF = len(FRAMES)
DIM = {0: 0, 1: 1, 5: 2, 21: 3}
def join(a, b):
    u = a | b
    if len(u) <= 1: return frozenset(u)
    for l in LINES:
        if u <= l: return l
    return frozenset(POINTS)
def dist(i, j):
    a, b = FRAMES[i], FRAMES[j]
    return DIM[len(join(a, b))] - DIM[len(a & b)]
def frame_of(basis):
    return where[span_pts([tuple(v) for v in basis])] if basis else 0
conj = lambda c: gmul(c, c)
def dot(r, p):
    s = 0
    for i in range(3): s ^= gmul(conj(r[i]), p[i])
    return s
def polar(f):                          # Hermitian polarity on point sets: {} <-> all, point <-> line
    return frozenset(p for p in POINTS if all(dot(r, p) == 0 for r in f))
PP = [where[polar(f)] for f in FRAMES]
assert sorted(PP) == list(range(NF)) and all(PP[PP[i]] == i for i in range(NF))
pidx = {p: i for i, p in enumerate(POINTS)}
line_through = {}
for li, l in enumerate(LINES):
    for a in l:
        for b in l:
            if a != b: line_through[(pidx[a], pidx[b])] = 22 + li
two = [tuple(sorted(pidx[p] for p in l))[:2] for l in LINES]
G = set()
for cols in itertools.product(vecs, repeat=3):
    img = []
    for p in POINTS:
        v = tuple(gmul(cols[0][r], p[0]) ^ gmul(cols[1][r], p[1]) ^ gmul(cols[2][r], p[2]) for r in range(3))
        img.append(norm(v))
    if None in img or len(set(img)) != 21: continue
    for frob in (False, True):
        pi = [pidx[tuple(conj(c) for c in q)] if frob else pidx[q] for q in img]
        g = (0,) + tuple(1 + pi[i] for i in range(21)) + tuple(line_through[(pi[a], pi[b])] for a, b in two) + (NF - 1,)
        G.add(g); G.add(tuple(PP[g[i]] for i in range(NF)))
G = np.array(sorted(G), dtype=np.int64)
print("PG(2,4): frames", NF, "symmetry group order", len(G))
assert all(dist(int(G[k][i]), int(G[k][j])) == dist(i, j) for k in range(0, len(G), 997) for i in range(NF) for j in range(NF))
memo = {}
def canon(fs):
    fs = tuple(sorted(set(fs)))
    if not fs: return ()
    if fs not in memo:
        im = np.sort(G[:, list(fs)], axis=1)
        code = np.zeros(len(G), dtype=np.int64)
        for c in range(im.shape[1]): code = code * NF + im[:, c]
        memo[fs] = tuple(int(v) for v in im[int(np.argmin(code))])
    return memo[fs]
def lin(terms):
    d = {}
    for c, fs in terms:
        k = canon(fs)
        if k == (): continue
        d[k] = d.get(k, Fr(0)) + Fr(c)
    return d
def U(*sets):
    out = set()
    for s in sets: out |= set(s)
    return out
def ingleton(A, B, C, D, T):
    def r(*ls): return U(T, *ls)
    return lin([(1, r(A, B)), (1, r(A, C)), (1, r(A, D)), (1, r(B, C)), (1, r(B, D)),
                (-1, r(A)), (-1, r(B)), (-1, r(A, B, C)), (-1, r(A, B, D)), (-1, r(C, D))])
def dfz_terms():
    """(65) as RHS - LHS over letter indices 0..6 = A,B,C,W,X,Y,Z; H(P) = r(T+P) - r(T), H(P|Q) = r(T+P+Q) - r(T+Q)"""
    A, B, C, W, X, Y, Z = range(7); t = []
    def H(c, *ls): t.append((c, frozenset(ls))); t.append((-c, frozenset()))
    def Hc(c, x, *ys): t.append((c, frozenset((x,) + ys))); t.append((-c, frozenset(ys)))
    H(-2, A); H(-1, B); H(-2, C); H(1, W); H(1, X); H(1, Y); H(1, Z)
    Hc(2, A, Z, Y); Hc(1, B, X, Z); Hc(2, C, A, X); Hc(3, X, W, Y); Hc(3, Z, W, C); Hc(5, W, A, B); Hc(5, Y, B, C)
    H(5, A); H(5, B); H(5, C); H(-5, A, B, C)
    return t
def dual_terms(t):
    E7 = frozenset(range(7)); out = []
    for c, S in t:
        for i in S: out.append((c, frozenset({i})))
        out.append((-c, E7)); out.append((c, E7 - S))
    return out
DFZT = dfz_terms(); DFZDT = dual_terms(DFZT)
def rowform(terms, Ls, T):
    return lin([(c, U(T, *[Ls[i] for i in S])) for c, S in terms])
cert = json.load(open(sys.argv[1])); assert cert["q"] == 4 and cert["m"] == 3
total = {}
def add(d, w):
    for kk, v in d.items(): total[kk] = total.get(kk, Fr(0)) + w * v
counts = {}
for row in cert["rows"]:
    lam = Fr(row["mult"]); assert lam >= 0
    t = row["type"]; counts[t] = counts.get(t, 0) + 1
    if t == "sub":
        S = [frame_of(b) for b in row["S"]]; a = frame_of(row["a"]); b = frame_of(row["b"])
        assert a != b and a not in S and b not in S
        g = lin([(1, S + [a]), (1, S + [b]), (-1, S + [a, b]), (-1, S)])
    elif t == "mono":
        S = [frame_of(b) for b in row["S"]]; a = frame_of(row["a"]); assert a not in S
        g = lin([(1, S + [a]), (-1, S)])
    else:
        Ls = [set(frame_of(b) for b in l) for l in row["letters"]]; T = set(frame_of(b) for b in row["T"])
        if t == "ing": g = ingleton(*Ls, T)
        elif t == "dfz": g = rowform(DFZT, Ls, T)
        elif t == "dfzd": g = rowform(DFZDT, Ls, T)
        else: raise ValueError(t)
    add({kk: -v for kk, v in g.items()}, lam)
bound = Fr(0)
def pairs_at(d): return [(s, t) for s in range(NF) for t in range(NF) if s != t and dist(s, t) == d]
for e in cert["eqs"]:
    mu = Fr(e["mult"]); dfn = e["def"]
    if dfn[0] == "single": g = lin([(1, [frame_of(dfn[2])]), (-1, [frame_of(dfn[1])])]); b = 0
    elif dfn[0] == "pair":
        (s, t), (s0, t0) = [frame_of(x) for x in dfn[1]], [frame_of(x) for x in dfn[2]]
        assert dist(s, t) == dist(s0, t0) and s != t and s0 != t0
        g = lin([(1, [s, t]), (-1, [t]), (-1, [s0, t0]), (1, [t0])]); b = 0
    elif dfn[0] == "norm":
        P1 = pairs_at(1); g = {}
        for s, t in P1:
            for kk, v in lin([(1, [s, t]), (-1, [t])]).items(): g[kk] = g.get(kk, Fr(0)) + v / len(P1)
        b = 1
    else: raise ValueError(dfn)
    add(g, mu); bound += mu * b
P3 = pairs_at(3)
for s, t in P3:
    for kk, v in lin([(1, [s, t]), (-1, [t])]).items(): total[kk] = total.get(kk, Fr(0)) - v / len(P3)
bad = {kk: v for kk, v in total.items() if v < 0}
print("pairs at distance 1 / 3:", len(pairs_at(1)), len(P3))
print("rows by type:", counts, " equalities:", len(cert["eqs"]))
print("all coefficients >= 0:", not bad, list(bad.items())[:3] if bad else "")
print(("PROVED" if not bad else "FAILED") + f" (given the validity of each row type): mean delta_3 <= {bound} = {float(bound)} x mean delta_1"
      f"  (certificate claims {cert['bound']})")
