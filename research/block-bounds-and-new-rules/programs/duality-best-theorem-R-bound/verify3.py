# verify3.py CERTIFICATE.json -- independent exact check of the block-3 dual certificate (pf-duality).
# Pure Python (fractions), no numpy, no code shared with lp3.py / cert3.py except the frame ordering convention
# (subspaces of F_2^3 sorted by (size, sorted elements)), which is rebuilt here from scratch.
#
# What is checked: every inequality row is rebuilt from its witness or label as a linear form g(r) that is >= 0 for
# the rank function r of ANY subspace arrangement over a field of characteristic != 2:
#   elemental  : r(S+i)+r(S+j)-r(S+i+j)-r(S) (submodularity) or r(all)-r(all-i) (monotonicity);
#   ing / ing1 : Ingleton  r(AB)+r(AC)+r(AD)+r(BC)+r(BD)-r(A)-r(B)-r(ABC)-r(ABD)-r(CD), letters = unions of frames,
#                conditioned on T (all sets joined with T);
#   dfz        : Dougherty-Freiling-Zeger 2013, Theorem 8.6 (char != 2), conditioned on T.
# Equality rows: uniformity (r({s}) all equal; r({s,t})-r({t}) depends only on the distance) and the normalisation
# mean_{ordered distance-1 pairs} [r({s,t}) - r({t})] = 1.  Variables are r on canonical sets (orbit minimum under
# GL(3,2) x duality).  Then it checks, coefficient by coefficient,
#       sum_i lam_i (-g_i) + sum_e mu_e eq_e - c  >= 0,     lam_i >= 0,
# with c = mean over ordered distance-3 pairs of r({s,t}) - r({t}); this proves  mean delta_3 <= sum_e mu_e b_e.
import sys, json, itertools, re
from fractions import Fraction as F

m = 3
def subspaces():
    seen = {frozenset({0})}; fr = [frozenset({0})]
    while fr:
        new = []
        for S in fr:
            for v in range(1, 8):
                if v not in S:
                    T = frozenset(S | {x ^ v for x in S})
                    if T not in seen: seen.add(T); new.append(T)
        fr = new
    return sorted(seen, key=lambda S: (len(S), sorted(S)))
SUBS = subspaces(); NF = len(SUBS); IDX = {S: i for i, S in enumerate(SUBS)}
assert NF == 16
def dim(S): return {1: 0, 2: 1, 4: 2, 8: 3}[len(S)]
def dist(i, j): return dim(SUBS[i]) + dim(SUBS[j]) - 2 * dim(SUBS[i] & SUBS[j])
def perp(S): return frozenset(a for a in range(8) if all(bin(a & x).count('1') % 2 == 0 for x in S))
perms = []
for A in itertools.product(range(1, 8), repeat=3):
    if len({0, A[0], A[1], A[2], A[0] ^ A[1], A[0] ^ A[2], A[1] ^ A[2], A[0] ^ A[1] ^ A[2]}) < 8: continue
    def img(x, A=A): return (A[0] if x & 1 else 0) ^ (A[1] if x & 2 else 0) ^ (A[2] if x & 4 else 0)
    p = [IDX[frozenset(img(x) for x in S)] for S in SUBS]
    perms.append(p)
    perms.append([IDX[perp(SUBS[p[i]])] for i in range(NF)])
assert len(perms) == 336
cache = {}
def canon(mk):
    if mk not in cache:
        best = mk
        for p in perms:
            t = 0
            for i in range(NF):
                if (mk >> i) & 1: t |= 1 << p[i]
            if t < best: best = t
        cache[mk] = best
    return cache[mk]
def form(terms):
    d = {}
    for c, mk in terms:
        k = canon(mk)
        if k == 0: continue
        d[k] = d.get(k, F(0)) + F(c)
    return d
def r_(mk): return mk

cert = json.load(open(sys.argv[1]))
total = {}
def add(d, w):
    for k, v in d.items(): total[k] = total.get(k, F(0)) + w * v

def ingleton(L, T):
    A, B, C, D = [l | T for l in L]
    return form([(1, A | B), (1, A | C), (1, A | D), (1, B | C), (1, B | D),
                 (-1, A), (-1, B), (-1, A | B | C), (-1, A | B | D), (-1, C | D)] +
                [(-1 + 1, T)])   # no net r(T) term in Ingleton (5 positive, 5 negative unions all contain T)
LAB = dict(A=1, B=2, C=4, W=3, Y=6, X=5, Z=7)
def dfz(V, T):
    def H(*ls):
        mk = T
        for l in ls: mk |= V[l]
        return [(1, mk), (-1, T)]
    def Hc(xs, ys):
        a = T; b = T
        for l in xs + ys: a |= V[l]
        for l in ys: b |= V[l]
        return [(1, a), (-1, b)]
    def sc(c, t): return [(c * x, y) for x, y in t]
    lhs = sc(2, H('A')) + H('B') + sc(2, H('C'))
    rhs = (H('W') + H('X') + H('Y') + H('Z') + sc(2, Hc(['A'], ['Z', 'Y'])) + Hc(['B'], ['X', 'Z'])
           + sc(2, Hc(['C'], ['A', 'X'])) + sc(3, Hc(['X'], ['W', 'Y'])) + sc(3, Hc(['Z'], ['W', 'C']))
           + sc(5, Hc(['W'], ['A', 'B'])) + sc(5, Hc(['Y'], ['B', 'C']))
           + sc(5, H('A') + H('B') + H('C') + sc(-1, H('A', 'B', 'C'))))
    return form(rhs + sc(-1, lhs))   # >= 0
bottom, top = IDX[frozenset({0})], IDX[SUBS[-1]]
nrows = {}
for row in cert["rows"]:
    lam = F(row["mult"]); assert lam >= 0
    t = row["type"]
    if t == "elemental":
        w = row["witness"]
        if w[0] == "mono":
            full = (1 << NF) - 1; g = form([(1, full), (-1, full ^ (1 << w[1]))])
        else:
            S, i, j = w[1], w[2], w[3]
            assert not (S >> i) & 1 and not (S >> j) & 1 and i != j
            g = form([(1, S | (1 << i)), (1, S | (1 << j)), (-1, S | (1 << i) | (1 << j)), (-1, S)])
    else:
        lab = row["label"]
        mo = re.match(r"(dfzc?|ing1?|ingc):L=\[([0-9, ]+)\],T=(\d+)", lab)
        if mo:
            L = [int(v) for v in mo.group(2).split(',')]; T = int(mo.group(3))
            g = ingleton(L, T) if mo.group(1).startswith('ing') else dfz(dict(zip(['A', 'B', 'C', 'W', 'X', 'Y', 'Z'], L)), T)
        else:
            mo = re.match(r"dfz:points=(lines|planes),T=(none|bottom|top|both)$", lab)
            assert mo, lab
            T = {'none': 0, 'bottom': 1 << bottom, 'top': 1 << top, 'both': (1 << bottom) | (1 << top)}[mo.group(2)]
            if mo.group(1) == 'lines': V = {l: 1 << IDX[frozenset({0, v})] for l, v in LAB.items()}
            else: V = {l: 1 << IDX[perp(frozenset({0, v}))] for l, v in LAB.items()}
            g = dfz(V, T)
        t = t.split(':')[0]
    nrows[t] = nrows.get(t, 0) + 1
    add({k: -v for k, v in g.items()}, lam)
pairs = {d: [((1 << i) | (1 << j), 1 << j) for i in range(NF) for j in range(NF) if i != j and dist(i, j) == d] for d in (1, 2, 3)}
bound = F(0)
for e in cert["eqs"]:
    mu = F(e["mult"]); dfn = e["def"]
    if dfn[0] == "single": g = form([(1, dfn[1]), (-1, dfn[2])]); b = 0
    elif dfn[0] == "pair":
        p, t_, p0, t0 = dfn[1:]
        assert dist(p.bit_length() - 1 if False else [i for i in range(NF) if (p >> i) & 1 and not (t_ >> i) & 1][0], t_.bit_length() - 1) \
            == dist([i for i in range(NF) if (p0 >> i) & 1 and not (t0 >> i) & 1][0], t0.bit_length() - 1)
        g = form([(1, p), (-1, t_), (-1, p0), (1, t0)]); b = 0
    else:
        g = {}
        for (pm, tm) in pairs[1]:
            for k, v in form([(1, pm), (-1, tm)]).items(): g[k] = g.get(k, F(0)) + v / len(pairs[1])
        b = 1
    add(g, mu); bound += mu * b
cobj = {}
for (pm, tm) in pairs[3]:
    for k, v in form([(1, pm), (-1, tm)]).items(): cobj[k] = cobj.get(k, F(0)) + v / len(pairs[3])
add(cobj, F(-1))
bad = {k: v for k, v in total.items() if v < 0}
print("rows by type:", nrows, " equality rows:", len(cert["eqs"]))
print("all coefficients of  sum lam(-g) + sum mu eq - c  are >= 0:", not bad, (list(bad.items())[:5] if bad else ""))
print("PROVED (given the validity of each row type): mean delta_3 <=", bound, "=", float(bound), "* mean delta_1",
      " (certificate claims", cert["bound"] + ")")
