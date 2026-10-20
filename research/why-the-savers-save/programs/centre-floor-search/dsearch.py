#!/usr/bin/env python3
"""ctr-search / dsearch.py: search in O1's one-sided normal form ("dictionary decoders").
Decoder: directions l_k in R^8 (a pool); scatter psi_k(a) = l_k . r_a / 4 (in col M); every port b picks a minimal
representation r_b = sum_{k in F_b} x_b(k) l_k, gather phi_k(b) = x_b(k).  Then sum_k psi_k(a) phi_k(b) = r_a.r_b/4 for
ALL a, b: an exact rank-8 realisation with today's junk (Gram/4), for every choice of representations.
Cost of term k (users A_k = {b : k in F_b}, readers R_k = {a : l_k . r_a != 0}):
   pair = rank_F2 B[R_k, A_k] (hub's Theorem C pairing rank),   star = dim span U(A_k) (one-sided, gather side).
The search picks representations (simulated annealing) to minimise the one-sided cost (default) or the pairing cost.
Pool families: F frame e_i (8), D D8-type e_i +- e_j (56), R the 112 roots, H half-spin (+-1)^8 (128 lines),
  T 2e_i + (+-1)^4 on a tetrad avoiding i (norm 8),  Q (+-1)^2 on pairs inside... (see pool()).
Representations enumerated: size 1-3 over the whole pool (exact integer Gram determinants), size 4 over the pool
directions inside span(e_T) of the port's tetrad (this contains the coordinate representation)."""
import sys, os, json, time, random, itertools, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *

def canon_dir(v):
    v = list(v)
    g = 0
    for x in v: g = np.gcd(g, abs(int(x)))
    v = [x // g for x in v]
    for x in v:
        if x != 0:
            return tuple(v) if x > 0 else tuple(-y for y in v)

def pool(fams):
    P = []
    if 'F' in fams: P += [tuple(1 if j == i else 0 for j in range(8)) for i in range(8)]
    if 'D' in fams:
        for i, j in itertools.combinations(range(8), 2):
            for s in (1, -1):
                v = [0] * 8; v[i] = 1; v[j] = s; P.append(tuple(v))
    if 'R' in fams: P += list(LINES)
    if 'H' in fams:
        for s in itertools.product((1, -1), repeat=7): P.append((1,) + s)
    if 'T' in fams:
        for i in range(8):
            for t in TET:
                if i in t: continue
                for s in itertools.product((1, -1), repeat=4):
                    v = [0] * 8; v[i] = 2
                    for a, x in zip(t, s): v[a] = x
                    P.append(tuple(v))
    if 'B' in fams:                                   # butterfly basis only (calibration 72)
        for a in range(4):
            for s in (1, -1):
                v = [0] * 8; v[2 * a] = 1; v[2 * a + 1] = s; P.append(tuple(v))
    out, seen = [], set()
    for v in P:
        c = canon_dir(v)
        if c not in seen: seen.add(c); out.append(c)
    return out

def det_int(M):
    """exact integer determinant (Bareiss)"""
    M = [list(map(int, r)) for r in M]
    n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k] != 0), None)
            if sw is None: return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[-1][-1]

def reps_for(P, maxsize=3, tet_size4=True, verbose=False, basis_rep=False):
    """for each port: list of minimal representations (tuples of pool indices)"""
    Pa = np.array(P, dtype=np.int64)          # p x 8
    p = len(P)
    GP = Pa @ Pa.T
    allreps = []
    t0 = time.time()
    for b in range(N):
        r = np.array(LINES[b], dtype=np.int64)
        g = Pa @ r; rr = int(r @ r)
        reps = []
        # size 1: l parallel to r  <=> Gram det 0
        d1 = GP.diagonal() * rr - g * g
        one = set(np.nonzero(d1 == 0)[0].tolist())
        reps += [(i,) for i in sorted(one)]
        # size 2: det Gram(li, lj) != 0, det Gram(li, lj, r) == 0, neither alone
        if maxsize >= 2:
            two = []
            for i in range(p):
                if i in one: continue
                j = np.arange(i + 1, p)
                gii, gjj, gij = GP[i, i], GP[j, j], GP[i, j]
                d2 = gii * gjj - gij * gij
                d3 = (gii * gjj * rr + 2 * gij * g[j] * g[i] - gii * g[j] ** 2 - gjj * g[i] ** 2 - rr * gij ** 2)
                ok = (d2 != 0) & (d3 == 0)
                for jj in j[ok].tolist():
                    if jj in one: continue
                    two.append((i, jj))
            reps += two
            twoset = set(two)
        if maxsize >= 3:
            # size 3: r in span(li, lj, lk), independent, minimal (no sub-pair / single spans r)
            cov2 = {}
            for (i, j) in two: cov2.setdefault(i, set()).add(j); cov2.setdefault(j, set()).add(i)
            three = []
            for i in range(p):
                if i in one: continue
                for j in range(i + 1, p):
                    if j in one or (i, j) in twoset: continue
                    gii, gjj, gij = GP[i, i], GP[j, j], GP[i, j]
                    if gii * gjj - gij * gij == 0: continue
                    k = np.arange(j + 1, p)
                    gkk, gik, gjk = GP[k, k], GP[i, k], GP[j, k]
                    # det Gram(li,lj,lk)
                    d3 = gii * (gjj * gkk - gjk * gjk) - gij * (gij * gkk - gjk * gik) + gik * (gij * gjk - gjj * gik)
                    # det Gram(li,lj,lk,r) via cofactor expansion along r row: use 4x4 formula
                    gi, gj, gk = g[i], g[j], g[k]
                    # 4x4 determinant [[gii,gij,gik,gi],[gij,gjj,gjk,gj],[gik,gjk,gkk,gk],[gi,gj,gk,rr]]
                    a11, a12, a13, a14 = gii, gij, gik, gi
                    a22, a23, a24 = gjj, gjk, gj
                    a33, a34 = gkk, gk
                    a44 = rr
                    # symmetric 4x4 det, expanded
                    d4 = (a11 * (a22 * (a33 * a44 - a34 * a34) - a23 * (a23 * a44 - a34 * a24) + a24 * (a23 * a34 - a33 * a24))
                          - a12 * (a12 * (a33 * a44 - a34 * a34) - a23 * (a13 * a44 - a34 * a14) + a24 * (a13 * a34 - a33 * a14))
                          + a13 * (a12 * (a23 * a44 - a34 * a24) - a22 * (a13 * a44 - a34 * a14) + a24 * (a13 * a24 - a23 * a14))
                          - a14 * (a12 * (a23 * a34 - a33 * a24) - a22 * (a13 * a34 - a33 * a14) + a23 * (a13 * a24 - a23 * a14)))
                    ok = (d3 != 0) & (d4 == 0)
                    for kk in k[ok].tolist():
                        if kk in one or (i, kk) in twoset or (j, kk) in twoset: continue
                        three.append((i, j, kk))
            reps += three
        if tet_size4:
            T = TET[TETOF[b]]
            inside = [i for i in range(p) if all(P[i][c] == 0 for c in range(8) if c not in T)]
            have = set(reps)
            small = [set(x) for x in reps]
            for F in itertools.combinations(inside, 4):
                if any(s <= set(F) for s in small): continue
                sub = Pa[list(F)]
                if np.linalg.matrix_rank(sub) == 4:                 # spans span(e_T) which contains r_b
                    reps.append(F)
        if basis_rep:
            Bm8 = [list(map(Fraction, [P[k][c] for k in range(8)])) for c in range(8)]
            x = solve_exact(Bm8, [Fraction(LINES[b][c]) for c in range(8)])
            F = tuple(k for k in range(8) if x[k] != 0)
            if F not in reps and not any(set(s) <= set(F) for s in reps): reps.append(F)
        allreps.append(reps)
        if verbose and b % 16 == 0:
            print(f'  port {b}: {len(reps)} reps ({time.time()-t0:.0f} s)', flush=True)
    return allreps

# ------------------------------------------------------------------ cost bookkeeping
def span_dim(labset):
    basis = []
    for v in labset:
        for x in basis: v = min(v, v ^ x)
        if v: basis.append(v)
    return len(basis), basis

class State:
    def __init__(self, P, reps, mode='star', rng=None):
        self.P, self.reps, self.mode = P, reps, mode
        self.p = len(P)
        Pa = np.array(P, dtype=np.int64)
        Ra = np.array(LINES, dtype=np.int64)
        self.read = [(Pa[k] @ Ra.T != 0) for k in range(self.p)]       # reader masks
        # Z_k = U(R_k)^perp : basis of label vectors orthogonal to every reader label
        self.Zperp = []
        for k in range(self.p):
            R = [a for a in range(N) if self.read[k][a]]
            d, basis = span_dim([LABS[a] for a in R])
            Z = [] if d == 9 else [c for c in range(1, 512) if all(dot(c, x) == 0 for x in basis)]
            self.Zperp.append(Z)
        self.users = [dict() for _ in range(self.p)]                   # label -> count
        self.choice = [0] * N
        self.costk = [0] * self.p
        self.rng = rng or random.Random(0)
    def term_cost(self, k):
        labs = list(self.users[k].keys())
        if not labs: return 0
        d, basis = span_dim(labs)
        if self.mode == 'star': return d
        # pairing: d - dim(U(A) cap Z), Z = U(R)^perp (dim 0 or 1 in P112)
        t = 0
        for z in self.Zperp[k]:
            v = z
            for x in basis: v = min(v, v ^ x)
            if v == 0: t += 1
        # dim of subspace = log2(t+1)
        return d - (t + 1).bit_length() + 1
    def assign(self, choice):
        self.choice = list(choice)
        self.users = [dict() for _ in range(self.p)]
        for b in range(N):
            for k in self.reps[b][self.choice[b]]:
                self.users[k][LABS[b]] = self.users[k].get(LABS[b], 0) + 1
        self.costk = [self.term_cost(k) for k in range(self.p)]
        return sum(self.costk)
    def move(self, b, newc):
        old = self.reps[b][self.choice[b]]; new = self.reps[b][newc]
        touched = set(old) | set(new)
        for k in old:
            self.users[k][LABS[b]] -= 1
            if self.users[k][LABS[b]] == 0: del self.users[k][LABS[b]]
        for k in new:
            self.users[k][LABS[b]] = self.users[k].get(LABS[b], 0) + 1
        delta = 0
        newcosts = {}
        for k in touched:
            c = self.term_cost(k); newcosts[k] = c; delta += c - self.costk[k]
        self.choice[b] = newc
        return delta, newcosts
    def commit(self, newcosts):
        for k, c in newcosts.items(): self.costk[k] = c
    def total(self): return sum(self.costk)

def anneal(st, iters=200000, T0=2.0, T1=0.05, init=None, seed=0, log=None):
    rng = random.Random(seed)
    if init is None: init = [rng.randrange(len(st.reps[b])) for b in range(N)]
    cur = st.assign(init)
    best, bestc = cur, list(st.choice)
    multi = [b for b in range(N) if len(st.reps[b]) > 1]
    if not multi: iters = 0
    for it in range(iters):
        T = T0 * (T1 / T0) ** (it / iters)
        b = rng.choice(multi)
        nc = rng.randrange(len(st.reps[b]))
        if nc == st.choice[b]: continue
        oc = st.choice[b]
        delta, newc = st.move(b, nc)
        if delta <= 0 or rng.random() < np.exp(-delta / T):
            st.commit(newc); cur += delta
            if cur < best: best, bestc = cur, list(st.choice)
        else:
            st.move(b, oc)          # revert (costs of touched return to the committed ones)
    st.assign(bestc)
    return best, bestc

def to_decoder(P, reps, choice):
    """explicit Phi (gather), Psi (scatter) as Fractions; terms = used directions"""
    used = sorted({k for b in range(N) for k in reps[b][choice[b]]})
    idx = {k: i for i, k in enumerate(used)}
    Phi = [[Fraction(0)] * N for _ in used]
    Psi = [[Fraction(ip(P[k], LINES[a]), 4) for a in range(N)] for k in used]
    for b in range(N):
        F = reps[b][choice[b]]
        A = np.array([P[k] for k in F], dtype=float).T            # 8 x |F|
        coef, *_ = np.linalg.lstsq(A, np.array(LINES[b], dtype=float), rcond=None)
        # exact rational coefficients: solve with Fractions on an independent row subset
        rows = [list(map(Fraction, [P[k][c] for k in F])) for c in range(8)]
        x = solve_exact(rows, [Fraction(LINES[b][c]) for c in range(8)])
        for k, v in zip(F, x): Phi[idx[k]][b] = v
    return used, Phi, Psi

def solve_exact(rows, rhs):
    """least-squares-free exact solve of rows x = rhs (consistent, full column rank)"""
    m = len(rows[0])
    A = [r[:] + [v] for r, v in zip(rows, rhs)]
    piv_rows = []
    r = 0
    for c in range(m):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: raise ValueError('rank deficient')
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [v / pv for v in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [a - f * bb for a, bb in zip(A[i], A[r])]
        r += 1
    for i in range(r, len(A)): assert A[i][m] == 0, 'inconsistent'
    return [A[i][m] for i in range(m)]

FLOOR_ONE_SIDED = 8 * 1260 * 6 / 1155      # cx-medium, proved (52.36)

def verify(P, reps, choice, tag):
    used, Phi, Psi = to_decoder(P, reps, choice)
    ok, bad = check_exact(Phi, Psi)
    c = costs(Phi, Psi)
    res = dict(tag=tag, exact=ok, bad=bad, n=c['n'], pair=c['pair'], star_gather=c['star_gather'],
               pair_terms=c['pair_terms'], star_terms=c['star_gather_terms'],
               directions=[P[k] for k in used])
    if ok and c['star_gather'] < FLOOR_ONE_SIDED:
        res['REFUSED'] = 'one-sided cost below the proved 52.36 floor: checker or model bug'
    return res

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--fams', default='F')
    ap.add_argument('--maxsize', type=int, default=3)
    ap.add_argument('--mode', default='star')
    ap.add_argument('--iters', type=int, default=200000)
    ap.add_argument('--restarts', type=int, default=4)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', default=None)
    ap.add_argument('--basisrep', action='store_true')
    ap.add_argument('--no4', action='store_true')
    a = ap.parse_args()
    P = pool(a.fams)
    print(f'pool {a.fams}: {len(P)} directions', flush=True)
    t0 = time.time()
    reps = reps_for(P, maxsize=a.maxsize, verbose=False, basis_rep=a.basisrep, tet_size4=not a.no4)
    nr = [len(r) for r in reps]
    print(f'reps per port: min {min(nr)} max {max(nr)} total {sum(nr)} ({time.time()-t0:.0f} s)', flush=True)
    st = State(P, reps, mode=a.mode)
    results = []
    for rs in range(a.restarts):
        best, bc = anneal(st, iters=a.iters, seed=a.seed + rs)
        v = verify(P, reps, bc, f'{a.fams}-{a.mode}-r{rs}')
        print(f'restart {rs}: annealed {a.mode} {best}; verified exact {v["exact"]} n {v["n"]} pair {v["pair"]} '
              f'star {v["star_gather"]} {v.get("REFUSED","")} ({time.time()-t0:.0f} s)', flush=True)
        print('   star terms', sorted(v['star_terms'], reverse=True)[:30], flush=True)
        results.append(v)
    if a.out:
        json.dump(results, open(a.out, 'w'), indent=1, default=str)
