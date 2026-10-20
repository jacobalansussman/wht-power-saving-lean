# lp3.py -- the dual LP at block 3 (pf-duality, 13 Oct 2026).  Runs where numpy + scipy exist (Colab #1).
#
# Variables: r(S) = dim(sum of Y_s, s in S) for every set S of the 16 frames of the block-3 interval, taken up to the
# symmetry group of the interval (GL(3,2) acting on s, and duality s -> s-perp; order 336).  Every certificate gives a
# family that is invariant under this group (PROOF.md, Theorem I and section 6), so symmetrising loses nothing.
# Constraints (all valid for subspace arrangements over C):
#   * elemental polymatroid inequalities (monotone, submodular), r(empty) = 0;
#   * optional UNIFORM: r({s}) equal for all s and r({s,t}) - r({t}) = delta_{d(s,t)}  (true for certificates);
#   * Dougherty-Freiling-Zeger Theorem 8.6 (valid over every field of characteristic != 2), in conditional form
#     (quotient by the span of a set T of frames), for the instances listed in dfz_instances().
# Objective: maximise mean delta at distance 3 subject to mean delta at distance 1 = 1.  Without DFZ the optimum is 3
# (the GF(2) family is feasible).  Output: the optimum and the DFZ multipliers (the "fractional unit").
import sys, itertools, json
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack

m = 3
def subspaces(m):
    n = 1 << m; seen = {frozenset({0})}; fr = [frozenset({0})]
    while fr:
        new = []
        for S in fr:
            for v in range(1, n):
                if v not in S:
                    T = frozenset(S | {x ^ v for x in S})
                    if T not in seen: seen.add(T); new.append(T)
        fr = new
    return sorted(seen, key=lambda S: (len(S), sorted(S)))
def dimof(S): return len(S).bit_length() - 1
SUBS = subspaces(m); NF = len(SUBS); IDX = {S: i for i, S in enumerate(SUBS)}
def perp(S): return frozenset(a for a in range(1 << m) if all(bin(a & x).count('1') % 2 == 0 for x in S))
def dist(i, j):
    s, t = SUBS[i], SUBS[j]; return dimof(s) + dimof(t) - 2 * dimof(s & t)
# group: GL(3,2) as images of basis vectors, plus duality
def apply(A, x):  # A = images of e0,e1,e2
    y = 0
    for b in range(m):
        if (x >> b) & 1: y ^= A[b]
    return y
GL = []
for A in itertools.product(range(1, 1 << m), repeat=m):
    span = {0}
    for a in A: span |= {x ^ a for x in span}
    if len(span) == 1 << m: GL.append(A)
assert len(GL) == 168
PERMS = []
for A in GL:
    p = [IDX[frozenset(apply(A, x) for x in S)] for S in SUBS]
    PERMS.append(p)
    PERMS.append([IDX[perp(SUBS[p[i]])] for i in range(NF)])
PERMS = np.array(PERMS)
N = 1 << NF
masks = np.arange(N, dtype=np.int64)
canon = masks.copy()
for p in PERMS:
    img = np.zeros(N, dtype=np.int64)
    for i in range(NF):
        img |= ((masks >> i) & 1) << int(p[i])
    canon = np.minimum(canon, img)
reps, orb = np.unique(canon, return_inverse=True)   # orb[mask] = variable index; reps[0] = 0 (empty set)
NV = len(reps)
print("frames", NF, "group", len(PERMS), "orbits of subsets", NV); sys.stdout.flush()
EMPTY = orb[0]

def lin(terms):
    """terms: list of (coef, mask) -> dict var->coef (empty set dropped: r = 0)"""
    d = {}
    for c, mk in terms:
        v = int(orb[mk])
        if v == EMPTY: continue
        d[v] = d.get(v, 0) + c
    return {k: v for k, v in d.items() if v != 0}

# elemental inequalities, as rows  a.x >= 0  ->  -a.x <= 0
rows = set()
full = N - 1
for i in range(NF):
    rows.add(tuple(sorted(lin([(1, full), (-1, full ^ (1 << i))]).items())))
for i in range(NF):
    for j in range(i + 1, NF):
        S = masks[((masks >> i) & 1 == 0) & ((masks >> j) & 1 == 0)]
        a = orb[S | (1 << i)]; b = orb[S | (1 << j)]; c = orb[S | (1 << i) | (1 << j)]; d = orb[S]
        lo = np.minimum(a, b); hi = np.maximum(a, b)
        quad = np.unique(np.stack([lo, hi, c, d], axis=1), axis=0)
        for (x, y, z, w) in quad:
            dd = {}
            for v, cf in ((x, 1), (y, 1), (z, -1), (w, -1)):
                v = int(v)
                if v == EMPTY: continue
                dd[v] = dd.get(v, 0) + cf
            dd = tuple(sorted((k, v) for k, v in dd.items() if v != 0))
            if dd: rows.add(dd)
rows = [dict(r) for r in rows]
print("distinct elemental inequalities:", len(rows)); sys.stdout.flush()

def bit(i): return 1 << i
LINES1 = [i for i in range(NF) if dimof(SUBS[i]) == 1]
def frame_of_vec(v): return IDX[frozenset({0, v})]
# DFZ labelling of the Fano plane on F_2^3: A=e1,B=e2,C=e3,W=e1+e2,Y=e2+e3,X=e1+e3,Z=e1+e2+e3
LAB = dict(A=1, B=2, C=4, W=3, Y=6, X=5, Z=7)
def dfz_terms(V, T):
    """V: letter -> frame mask; T: conditioning mask. returns dict for  RHS - LHS >= 0"""
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
    return lin(rhs + sc(-1, lhs))
def dfz_instances(kind):
    out = []
    bottom, top = IDX[frozenset({0})], IDX[SUBS[-1]]
    Ts = {'none': 0, 'bottom': bit(bottom), 'top': bit(top), 'both': bit(bottom) | bit(top)}
    for tn, T in Ts.items():
        # Fano points = the 7 frames of 1-dim s (adjacent to the bottom), and = the 7 frames of 2-dim s via perp
        V1 = {l: bit(frame_of_vec(v)) for l, v in LAB.items()}
        V2 = {l: bit(IDX[perp(frozenset({0, v}))]) for l, v in LAB.items()}
        out.append((f"points=lines,T={tn}", dfz_terms(V1, T)))
        out.append((f"points=planes,T={tn}", dfz_terms(V2, T)))
        if kind == 'more':
            # points = line-frame joined with the bottom / top frame
            V3 = {l: bit(frame_of_vec(v)) | bit(top) for l, v in LAB.items()}
            out.append((f"points=lines+top,T={tn}", dfz_terms(V3, T)))
    return out

def solve(uniform, use_dfz, kind='basic'):
    A_ub = []; b_ub = []
    for r in rows:
        A_ub.append({k: -v for k, v in r.items()}); b_ub.append(0)
    names = []
    if use_dfz:
        for nm, r in dfz_instances(kind):
            A_ub.append({k: -v for k, v in r.items()}); b_ub.append(0); names.append(nm)
    A_eq = []; b_eq = []
    singles = [bit(i) for i in range(NF)]
    pairs = {d: [(bit(i) | bit(j), bit(j)) for i in range(NF) for j in range(NF) if i != j and dist(i, j) == d] for d in (1, 2, 3)}
    if uniform:
        s0 = singles[0]
        for s in singles[1:]:
            e = lin([(1, s), (-1, s0)])
            if e: A_eq.append(e); b_eq.append(0)
        for d in (1, 2, 3):
            p0 = pairs[d][0]
            for p in pairs[d][1:]:
                e = lin([(1, p[0]), (-1, p[1]), (-1, p0[0]), (1, p0[1])])
                if e: A_eq.append(e); b_eq.append(0)
    # normalisation: mean delta over ordered distance-1 pairs = 1; objective: mean over distance-3 pairs
    def mean(d):
        acc = {}
        for (pm, tm) in pairs[d]:
            for k, v in lin([(1, pm), (-1, tm)]).items(): acc[k] = acc.get(k, 0) + v / len(pairs[d])
        return acc
    A_eq.append(mean(1)); b_eq.append(1)
    obj = mean(3)
    c = np.zeros(NV)
    for k, v in obj.items(): c[k] = -v
    def mat(L):
        I = []; J = []; X = []
        for r_, row in enumerate(L):
            for k, v in row.items(): I.append(r_); J.append(k); X.append(v)
        return coo_matrix((X, (I, J)), shape=(len(L), NV)).tocsr()
    bounds = [(0, 0) if k == EMPTY else (0, None) for k in range(NV)]
    res = linprog(c, A_ub=mat(A_ub), b_ub=b_ub, A_eq=mat(A_eq), b_eq=b_eq, bounds=bounds, method='highs')
    out = {"uniform": uniform, "dfz": use_dfz, "kind": kind, "status": res.status, "message": res.message}
    if res.status == 0:
        out["optimum"] = -res.fun
        if use_dfz:
            marg = res.ineqlin.marginals[-len(names):]
            out["dfz_multipliers"] = {nm: float(-mu) for nm, mu in zip(names, marg)}
        x = res.x
        def val(terms): return sum(v * x[k] for k, v in lin(terms).items())
        out["delta"] = {d: float(sum(val([(1, p), (-1, t)]) for p, t in pairs[d]) / len(pairs[d])) for d in (1, 2, 3)}
        out["dim_single_mean"] = float(sum(val([(1, s)]) for s in singles) / NF)
    return out

if __name__ == "__main__":
    for uniform, dfz, kind in [(True, False, 'basic'), (True, True, 'basic'), (False, True, 'basic'), (True, True, 'more')]:
        o = solve(uniform, dfz, kind)
        print(json.dumps(o, default=str)); sys.stdout.flush()
