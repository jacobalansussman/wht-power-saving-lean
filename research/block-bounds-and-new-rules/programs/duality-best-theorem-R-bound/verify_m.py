# verify_m.py CERTIFICATE.json -- independent exact check of an interval-LP dual certificate at block m (pf-duality).
# Shares no code with lpm.py / cert_m.py.  Rebuilds: the frames (subspaces s of F_2^m, ordered by (size, sorted
# elements)), distances d(s,t) = dim(s+t) - dim(s & t), the group GL(m,2) x {s -> s^perp} as permutations, and the
# canonical form of a set of frames (lexicographically least sorted image).  Every row is rebuilt from its explicit
# frames as a linear form g(r) that is >= 0 for the rank function r(S) = dim sum_{s in S} Y_s of ANY family of
# subspaces over a field of characteristic != 2:
#   sub  : r(S+a) + r(S+b) - r(S+a+b) - r(S)        (submodularity)
#   mono : r(S+a) - r(S)                            (monotonicity)
#   ing  : Ingleton, letters A,B,C,D (unions of frames), every term joined with T (quotient by sum_{t in T} Y_t)
#   dfz  : Dougherty-Freiling-Zeger Thm 8.6 (arXiv:1311.4601, (65)), same conditioning
# Equalities (true for certificate families, PROOF.md Theorem I): r({s}) = r({s'}); r({s,t})-r({t}) depends only on
# d(s,t); normalisation mean_{d=1}(r({s,t})-r({t})) = 1.  Symmetric r (orbit-invariant) loses nothing (PROOF.md 6.1).
# Checks  sum lam_i(-g_i) + sum mu_e eq_e - c >= 0 coefficientwise (r >= 0), lam >= 0, c = mean_{d=m}(r({s,t})-r({t})).
import sys, json, itertools
from fractions import Fraction as Fr
import numpy as np

cert = json.load(open(sys.argv[1])); m = cert["m"]; n = 1 << m
frames = {frozenset({0})}; todo = [frozenset({0})]
while todo:
    nxt = []
    for S in todo:
        for v in range(1, n):
            T = frozenset(S | {x ^ v for x in S})
            if T not in frames: frames.add(T); nxt.append(T)
    todo = nxt
frames = sorted(frames, key=lambda S: (len(S), sorted(S))); NF = len(frames); where = {S: i for i, S in enumerate(frames)}
def dimS(S): return len(S).bit_length() - 1
def dist(i, j): return dimS(frames[i]) + dimS(frames[j]) - 2 * dimS(frames[i] & frames[j])
def perp(S): return frozenset(a for a in range(n) if all(bin(a & x).count('1') % 2 == 0 for x in S))
pp = [where[perp(S)] for S in frames]
G = []
for cols in itertools.product(range(1, n), repeat=m):
    img = [0] * n
    for x in range(1, n):
        low = x & -x; b = low.bit_length() - 1
        img[x] = img[x ^ low] ^ cols[b]
    if len(set(img)) != n: continue
    g = [where[frozenset(img[x] for x in S)] for S in frames]
    G.append(g); G.append([pp[g[i]] for i in range(NF)])
G = np.array(G, dtype=np.int64)
print("block", m, "frames", NF, "group order", len(G))
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
def dfz(Ls, T):
    A, B, C, W, X, Y, Z = Ls
    def H(*ls): return [(1, U(T, *ls)), (-1, set(T))]
    def Hc(x, *ys): return [(1, U(T, x, *ys)), (-1, U(T, *ys))]
    def k(c, t): return [(c * a, b) for a, b in t]
    lhs = k(2, H(A)) + H(B) + k(2, H(C))
    rhs = (H(W) + H(X) + H(Y) + H(Z) + k(2, Hc(A, Z, Y)) + Hc(B, X, Z) + k(2, Hc(C, A, X)) + k(3, Hc(X, W, Y))
           + k(3, Hc(Z, W, C)) + k(5, Hc(W, A, B)) + k(5, Hc(Y, B, C)) + k(5, H(A) + H(B) + H(C) + k(-1, H(A, B, C))))
    return lin(rhs + k(-1, lhs))
total = {}
def add(d, w):
    for kk, v in d.items(): total[kk] = total.get(kk, Fr(0)) + w * v
counts = {}
for row in cert["rows"]:
    lam = Fr(row["mult"]); assert lam >= 0
    t = row["type"]; counts[t] = counts.get(t, 0) + 1
    if t == "sub":
        S, a, b = row["witness"]; assert a != b and a not in S and b not in S
        g = lin([(1, S + [a]), (1, S + [b]), (-1, S + [a, b]), (-1, S)])
    elif t == "mono":
        S, a = row["witness"]; assert a not in S
        g = lin([(1, S + [a]), (-1, S)])
    elif t == "ing":
        g = ingleton(*[set(l) for l in row["letters"]], set(row["T"]))
    elif t == "dfz":
        g = dfz([set(l) for l in row["letters"]], set(row["T"]))
    else: raise ValueError(t)
    add({kk: -v for kk, v in g.items()}, lam)
bound = Fr(0)
for e in cert["eqs"]:
    mu = Fr(e["mult"]); dfn = e["def"]
    if dfn[0] == "single": g = lin([(1, [dfn[2]]), (-1, [dfn[1]])]); b = 0
    elif dfn[0] == "pair":
        (s, t), (s0, t0) = dfn[1], dfn[2]
        assert dist(s, t) == dist(s0, t0) and s != t and s0 != t0
        g = lin([(1, [s, t]), (-1, [t]), (-1, [s0, t0]), (1, [t0])]); b = 0
    elif dfn[0] == "norm":
        pairs = [(s, t) for s in range(NF) for t in range(NF) if s != t and dist(s, t) == 1]
        g = {}
        for s, t in pairs:
            for kk, v in lin([(1, [s, t]), (-1, [t])]).items(): g[kk] = g.get(kk, Fr(0)) + v / len(pairs)
        b = 1
    else: raise ValueError(dfn)
    add(g, mu); bound += mu * b
pairs = [(s, t) for s in range(NF) for t in range(NF) if s != t and dist(s, t) == m]
for s, t in pairs:
    for kk, v in lin([(1, [s, t]), (-1, [t])]).items(): total[kk] = total.get(kk, Fr(0)) - v / len(pairs)
bad = {kk: v for kk, v in total.items() if v < 0}
print("rows by type:", counts, " equalities:", len(cert["eqs"]))
print("all coefficients >= 0:", not bad, list(bad.items())[:3] if bad else "")
print(("PROVED" if not bad else "FAILED") + f" (given the validity of each row type): mean delta_{m} <= {bound} = {float(bound)} x mean delta_1"
      f"  (certificate claims {cert['bound']})")
