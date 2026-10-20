# cut_lp.py (d-G3-SPECTATOR/more-spectators, wave 17).  Cut certificates in the Walsh-diagonal ("commutative") game.
# Group G = Sym_m(F_2); a rank-1 move is z z^T; the target B_k (x) I is M_K = sum_{i<k} e_i e_i^T.
# A linear functional lam on Sym_m is a quadratic form Q(z) = sum a_i z_i + sum_{i<j} b_ij z_i z_j (lam(zz^T) = Q(z)),
# lam(M_K) = sum_{i<k} a_i mod 2.  Cut certificate: V = (+) Q[G/ker lam] with half-space Y_0 per lam; f(g) = sum w_lam lam(g).
# (A) no address maps:  bound_A = max sum_lam w lam(M_K)  s.t.  sum_lam w Q_lam(z) <= 1 for all z != 0, w >= 0   (an LP)
# (B) all address maps free (GL(m,2) acts by congruence; averaging lemma): bound_B = max over single lam of
#     P[lam(aa^T + bb^T) = 1] / P[Q(z) = 1], a,b a uniform random independent pair, z uniform nonzero  [for k = 2]
#     and for k = 3: P[lam(aa^T+bb^T+cc^T) = 1] with a,b,c independent.
# (C) exhaustive check that k = 4 has no k forms with pairwise disjoint supports and odd sum_{i<k} a_i (perfect cube cert).
# Usage: python3 -I cut_lp.py
import itertools, sys
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction

def popc(x): return bin(x).count('1')
def forms(m):
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    nz = list(range(1, 2 ** m))
    for a in range(2 ** m):
        for b in range(2 ** len(pairs)):
            Qv = []
            for z in nz:
                v = popc(a & z) & 1
                for t, (i, j) in enumerate(pairs):
                    if (b >> t) & 1: v ^= ((z >> i) & 1) & ((z >> j) & 1)
                Qv.append(v)
            yield a, b, np.array(Qv, dtype=np.int8)

def lpA(m, k):
    rows = []; obj = []
    for a, b, Qv in forms(m):
        par = popc(a & ((1 << k) - 1)) & 1
        if par == 0: continue           # a form with lam(M_K) = 0 only adds to the constraints
        rows.append(Qv); obj.append(1.0)
    A = np.array(rows, dtype=float).T   # constraints x forms
    res = linprog(-np.array(obj), A_ub=A, b_ub=np.ones(A.shape[0]), bounds=(0, None), method='highs')
    return -res.fun, len(obj)

def boundB(m, k):
    n = 2 ** m; nz = list(range(1, n))
    # independent k-tuples (ordered) of vectors
    def indep(vs):
        span = {0}
        for v in vs:
            if v in span: return False
            span |= {s ^ v for s in span}
        return True
    tuples = [vs for vs in itertools.permutations(nz, k) if indep(vs)]
    best = (0, None)
    for a, b, Qv in forms(m):
        Q = {z: int(Qv[i]) for i, z in enumerate(nz)}
        p1 = Fraction(sum(Q.values()), len(nz))
        if p1 == 0: continue
        # lam(sum v v^T) = sum Q(v) mod 2  (lam linear)
        p2 = Fraction(sum((sum(Q[v] for v in vs) & 1) for vs in tuples), len(tuples))
        r = p2 / p1
        if r > best[0]: best = (r, (a, b, p1, p2))
    return best

def perfect_cube_search(k):
    cands = []
    for a, b, Qv in forms(k):
        if popc(a & ((1 << k) - 1)) & 1:
            cands.append(int(''.join(str(int(x)) for x in Qv), 2))
    # need k pairwise disjoint supports
    cands = sorted(set(cands))
    def rec(chosen, start):
        if len(chosen) == k: return chosen
        for i in range(start, len(cands)):
            c = cands[i]
            if all(c & d == 0 for d in chosen):
                r = rec(chosen + [c], i + 1)
                if r: return r
        return None
    return len(cands), rec([], 0)

if __name__ == '__main__':
    for k in (2, 3, 4):
        nc, sol = perfect_cube_search(k)
        print(f"(C) k={k}: {nc} forms with odd target value; k pairwise-disjoint ones: {'found ' + str([bin(s) for s in sol]) if sol else 'NONE'}", flush=True)
    for (m, k) in [(2, 2), (3, 2), (4, 2), (5, 2), (3, 3), (4, 3), (5, 3), (4, 4), (5, 4), (5, 5)]:
        v, nv = lpA(m, k)
        print(f"(A) no address maps: m={m} k={k}: best cut-certificate bound {v:.6f} moves per array ({nv} forms)", flush=True)
    for (m, k) in [(2, 2), (3, 2), (4, 2), (5, 2), (3, 3), (4, 3)]:
        r, info = boundB(m, k)
        print(f"(B) all address maps free: m={m} k={k}: best single-cut averaged bound {float(r):.6f} = {r}  (a={info[0]:b}, b={info[1]:b}, p1={info[2]}, p2={info[3]})", flush=True)
