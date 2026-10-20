# lpm.py -- the interval LP at block m, restricted to sets of at most K frames (pf-duality, 14 Oct 2026).
#
# Frames: the 2^... subspaces s of F_2^m (the interval between L_0 and M_0, PROOF.md Lemma I1); d(s,t) = dim(s+t)-dim(s&t).
# Unknown: r(S) = dim(sum_{s in S} Y_s) for every set S of at most K frames, taken up to the symmetry group
# GL(m,2) x duality (s -> s^perp), which preserves distances; symmetrising a feasible point keeps it feasible.
# Rows (all valid for subspace arrangements over C, restricted to sets of size <= K):
#   submodularity r(S+a)+r(S+b) >= r(S+a+b)+r(S) (|S| <= K-2), monotonicity r(S+a) >= r(S) (|S| <= K-1), r(empty)=0;
#   UNIFORM: r({s}) all equal; r({s,t})-r({t}) depends only on d(s,t)  (true for certificates, Theorem I);
#   cuts: DFZ Theorem 8.6 (char != 2) and Ingleton, letters = frames, conditioned on a set T of frames.
# Objective: maximise mean delta_m subject to mean delta_1 = 1.
import sys, itertools, json, time
import numpy as np
from math import comb
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

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

class Interval:
    def __init__(self, m, K, log=print):
        t0 = time.time()
        self.m, self.K = m, K
        n = 1 << m
        SUBS = subspaces(m); self.SUBS = SUBS; NF = len(SUBS); self.NF = NF
        IDX = {S: i for i, S in enumerate(SUBS)}; self.IDX = IDX
        self.dim = np.array([len(S).bit_length() - 1 for S in SUBS])
        def perp(S): return frozenset(a for a in range(n) if all(bin(a & x).count('1') % 2 == 0 for x in S))
        self.perp = perp
        D = np.zeros((NF, NF), dtype=np.int64)
        for i in range(NF):
            for j in range(NF):
                D[i, j] = self.dim[i] + self.dim[j] - 2 * (len(SUBS[i] & SUBS[j]).bit_length() - 1)
        self.D = D
        PERP = [IDX[perp(S)] for S in SUBS]
        P = []
        for A in itertools.product(range(1, n), repeat=m):
            span = {0}
            for a in A: span |= {x ^ a for x in span}
            if len(span) != n: continue
            img = [0] * n
            for x in range(n):
                y = 0
                for b in range(m):
                    if (x >> b) & 1: y ^= A[b]
                img[x] = y
            p = [IDX[frozenset(img[x] for x in S)] for S in SUBS]
            P.append(p)
            P.append([PERP[p[i]] for i in range(NF)])
        self.P = np.array(P, dtype=np.int32)
        log(f"m={m} K={K} frames={NF} group={len(P)} ({time.time()-t0:.1f}s)")
        # binomial table and orbit ids per size
        self.Cb = np.array([[comb(a, j) for j in range(K + 2)] for a in range(NF + 1)], dtype=np.int64)
        self.oid, self.reps, self.off = [], [], []
        nv = 0
        for k in range(K + 1):
            N = comb(NF, k); oid = np.full(N, -1, dtype=np.int32); reps = []
            self.off.append(nv)
            if k == 0:
                oid[0] = 0; reps.append(()); nv += 1
            else:
                ptr = 0; CH = 1 << 20
                while True:
                    while ptr < N:
                        nz = np.flatnonzero(oid[ptr:ptr + CH] < 0)
                        if len(nz): ptr += int(nz[0]); break
                        ptr += CH
                    if ptr >= N: break
                    S = self.unrank(ptr, k)
                    imgs = np.sort(self.P[:, S], axis=1)
                    oid[self.rank(imgs)] = nv - self.off[k]
                    reps.append(tuple(S)); nv += 1
            self.oid.append(oid); self.reps.append(reps)
            log(f"  size {k}: {len(reps)} orbits ({time.time()-t0:.1f}s)")
        self.NV = nv
        self.EMPTY = 0

    def rank(self, S):  # S: (R, k) sorted ascending
        k = S.shape[1]; r = np.zeros(S.shape[0], dtype=np.int64)
        for i in range(k): r += self.Cb[S[:, i], i + 1]
        return r
    def unrank(self, r, k):
        S = []
        for i in range(k, 0, -1):
            c = i - 1
            while comb(c + 1, i) <= r: c += 1
            S.append(c); r -= comb(c, i)
        return sorted(S)
    def ids(self, S):
        """S: (R, k) array of frame indices, rows of DISTINCT entries (any order). returns global variable ids."""
        S = np.sort(np.asarray(S, dtype=np.int64), axis=1); k = S.shape[1]
        if k == 0: return np.zeros(S.shape[0], dtype=np.int64)
        return self.oid[k][self.rank(S)].astype(np.int64) + self.off[k]
    def id1(self, frames):
        fs = sorted(set(frames))
        if len(fs) > self.K: return None
        if not fs: return 0
        return int(self.ids(np.array([fs]))[0])

    def polymatroid_rows(self, log=print):
        """returns list of rows as (ids array (R,4), coefs) -- each row: sum coef*r >= 0."""
        rows = []
        K, NF = self.K, self.NF
        allf = np.arange(NF)
        for j in range(K - 1):
            for S in self.reps[j]:
                comp = np.setdiff1d(allf, np.array(S, dtype=np.int64))
                a, b = np.triu_indices(len(comp), 1); a = comp[a]; b = comp[b]
                base = np.tile(np.array(S, dtype=np.int64), (len(a), 1)) if S else np.zeros((len(a), 0), dtype=np.int64)
                iA = self.ids(np.column_stack([base, a])); iB = self.ids(np.column_stack([base, b]))
                iAB = self.ids(np.column_stack([base, a, b])); iS = self.ids(base)
                lo = np.minimum(iA, iB); hi = np.maximum(iA, iB)
                rows.append(np.stack([lo, hi, iAB, iS], axis=1))
        R = np.unique(np.concatenate(rows), axis=0)
        mono = []
        for j in range(K):
            for S in self.reps[j]:
                comp = np.setdiff1d(allf, np.array(S, dtype=np.int64))
                base = np.tile(np.array(S, dtype=np.int64), (len(comp), 1)) if S else np.zeros((len(comp), 0), dtype=np.int64)
                mono.append(np.stack([self.ids(np.column_stack([base, comp])), self.ids(base)], axis=1))
        Mo = np.unique(np.concatenate(mono), axis=0)
        out = []
        for lo, hi, ab, s in R:
            d = {}
            for v, c in ((lo, 1), (hi, 1), (ab, -1), (s, -1)):
                if v == 0: continue
                d[int(v)] = d.get(int(v), 0) + c
            d = {k: v for k, v in d.items() if v}
            if d: out.append(d)
        for a, s in Mo:
            d = {int(a): 1}
            if s != 0: d[int(s)] = d.get(int(s), 0) - 1
            d = {k: v for k, v in d.items() if v}
            if d: out.append(d)
        uniq = {tuple(sorted(d.items())) for d in out}
        log(f"polymatroid rows: {len(uniq)}")
        return [dict(t) for t in uniq]

    def pair_terms(self, d):
        """mean over ordered pairs (s,t) at distance d of r({s,t}) - r({t}), as dict id->coef"""
        I, J = np.nonzero(self.D == d)
        pst = self.ids(np.stack([I, J], axis=1)); pt = self.ids(J[:, None])
        acc = {}
        for a, b in zip(pst, pt):
            acc[int(a)] = acc.get(int(a), 0) + 1.0 / len(I)
            acc[int(b)] = acc.get(int(b), 0) - 1.0 / len(I)
        return {k: v for k, v in acc.items() if abs(v) > 1e-15 and k != 0}
    def uniform_eqs(self):
        eqs = []
        s = [self.id1([i]) for i in range(self.NF)]
        for v in sorted(set(s) - {s[0]}): eqs.append({v: 1, s[0]: -1})
        for d in range(1, self.m + 1):
            I, J = np.nonzero(self.D == d)
            pst = self.ids(np.stack([I, J], axis=1)); pt = self.ids(J[:, None])
            pairs = sorted(set(zip(pst.tolist(), pt.tolist())))
            p0 = pairs[0]
            for p in pairs[1:]:
                e = {}
                for v, c in ((p[0], 1), (p[1], -1), (p0[0], -1), (p0[1], 1)):
                    e[v] = e.get(v, 0) + c
                e = {k: v for k, v in e.items() if v}
                if e: eqs.append(e)
        return eqs

# ---- DFZ (Theorem 8.6) and Ingleton, as term lists over letters; each term (coef, frozenset of letters) = coef*r(T u letters)
A_, B_, C_, W_, X_, Y_, Z_ = range(7)
def _Hs(c, *ls): return [(c, frozenset(ls)), (-c, frozenset())]
def _Hc(c, xs, ys): return [(c, frozenset(xs) | frozenset(ys)), (-c, frozenset(ys))]
_D = _Hs(-2, A_) + _Hs(-1, B_) + _Hs(-2, C_) + _Hs(1, W_) + _Hs(1, X_) + _Hs(1, Y_) + _Hs(1, Z_)
_D += _Hc(2, [A_], [Z_, Y_]) + _Hc(1, [B_], [X_, Z_]) + _Hc(2, [C_], [A_, X_]) + _Hc(3, [X_], [W_, Y_]) + _Hc(3, [Z_], [W_, C_])
_D += _Hc(5, [W_], [A_, B_]) + _Hc(5, [Y_], [B_, C_]) + _Hs(5, A_) + _Hs(5, B_) + _Hs(5, C_) + _Hs(-5, A_, B_, C_)
_acc = {}
for _c, _s in _D: _acc[_s] = _acc.get(_s, 0) + _c
DFZ = [(c, s) for s, c in _acc.items() if c != 0]          # RHS - LHS >= 0
ING = [(1, frozenset({0, 1})), (1, frozenset({0, 2})), (1, frozenset({0, 3})), (1, frozenset({1, 2})), (1, frozenset({1, 3})),
       (-1, frozenset({0})), (-1, frozenset({1})), (-1, frozenset({0, 1, 2})), (-1, frozenset({0, 1, 3})), (-1, frozenset({2, 3}))]
FANO_LAB = [1, 2, 4, 3, 5, 6, 7]   # A,B,C,W,X,Y,Z as coefficient vectors over a basis (e1,e2,e3): W=A+B, X=A+C, Y=B+C, Z=A+B+C

def inst_row(I, terms, letters, T):
    """letters: list of frame-sets (tuples), T: tuple of frames. returns dict id->coef (value >= 0 valid) or None if too big"""
    d = {}
    for c, s in terms:
        fr = set(T)
        for l in s: fr |= set(letters[l])
        v = I.id1(fr)
        if v is None: return None
        if v == 0: continue
        d[v] = d.get(v, 0) + c
    return {k: v for k, v in d.items() if v}

def fano_letters(I, s_idx, basis):
    """frames s + span(coefficient combination of basis vectors), for the 7 DFZ letters (up-neighbour Fano at s)"""
    S = set(I.SUBS[s_idx]); out = []
    for lab in FANO_LAB:
        v = 0
        for b in range(3):
            if (lab >> b) & 1: v ^= basis[b]
        out.append((I.IDX[frozenset(S | {x ^ v for x in S})],))
    return out

def solve(I, rows, eqs, cuts, uniform=True, log=print):
    A_ub = [{k: -v for k, v in r.items()} for r in rows] + [{k: -v for k, v in r.items()} for _, r in cuts]
    A_eq = list(eqs) if uniform else []
    b_eq = [0.0] * len(A_eq)
    A_eq.append(I.pair_terms(1)); b_eq.append(1.0)
    c = np.zeros(I.NV)
    for k, v in I.pair_terms(I.m).items(): c[k] = -v
    def mat(L):
        Ii, J, X = [], [], []
        for r_, row in enumerate(L):
            for k, v in row.items(): Ii.append(r_); J.append(k); X.append(v)
        return coo_matrix((X, (Ii, J)), shape=(len(L), I.NV)).tocsr()
    bounds = [(0, 0) if k == 0 else (0, None) for k in range(I.NV)]
    t0 = time.time()
    res = linprog(c, A_ub=mat(A_ub), b_ub=np.zeros(len(A_ub)), A_eq=mat(A_eq), b_eq=b_eq, bounds=bounds, method='highs')
    if res.status != 0:
        log(f"LP status {res.status}: {res.message}"); return None
    marg = -res.ineqlin.marginals[len(rows):]
    active = sorted([(float(mu), nm) for (nm, _), mu in zip(cuts, marg) if mu > 1e-9], reverse=True)
    x = res.x
    delta = {}
    for d in range(1, I.m + 1):
        delta[d] = float(sum(v * x[k] for k, v in I.pair_terms(d).items()))
    return dict(opt=-res.fun, x=x, active=active, delta=delta, secs=time.time() - t0, res=res)

if __name__ == "__main__":
    m, K = int(sys.argv[1]), int(sys.argv[2])
    I = Interval(m, K)
    rows = I.polymatroid_rows()
    eqs = I.uniform_eqs()
    bottom = 0
    e = [1 << b for b in range(m)]
    cuts = []
    # up-neighbour Fano at the bottom: letters = 7 points of span(e1,e2,e3)
    for nmT, T in (("none", ()), ("bottom", (bottom,))):
        r = inst_row(I, DFZ, fano_letters(I, bottom, e[:3]), T)
        if r: cuts.append((f"dfz:points@bottom,T={nmT}", r))
    if m >= 4:
        x0 = I.IDX[frozenset({0, e[0]})]
        for nmT, T in (("none", ()), ("point", (x0,))):
            r = inst_row(I, DFZ, fano_letters(I, x0, e[1:4]), T)
            if r: cuts.append((f"dfz:lines@point,T={nmT}", r))
    for uni in (False, True):
        out = solve(I, rows, eqs, [], uniform=uni)
        print(json.dumps({"uniform": uni, "cuts": 0, "opt": out["opt"], "delta": out["delta"], "secs": round(out["secs"], 1)}))
        out = solve(I, rows, eqs, cuts, uniform=uni)
        print(json.dumps({"uniform": uni, "cuts": len(cuts), "opt": out["opt"], "delta": out["delta"], "active": out["active"], "secs": round(out["secs"], 1)}))
        sys.stdout.flush()
