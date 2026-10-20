"""r-riders / lb.py -- exact lower bound on the paid moves of ANY word (wider game, so also the equal-frame game)
that leaves the toy's live registers (and the delivered register Z) holding the target formal sums.

Argument (Lemma A of wht35/wider/w-size.md section 5, applied to an arbitrary target instead of kernel x identity):
on f digit columns, at the Walsh basis vector whose column labels have multiplicities m = (m_y), every letter of
this sub-game is a scalar on each array: dir(z) -> i^N_z(m), shift(w) -> (-1)^N_w(m), N_z(m) = sum_y m_y chi_z(y);
mixing is the same rational matrix at every point.  Take two points m, m' at which every shift has the same value
(N_w(m) = N_w(m') mod 2 for all w).  Then each direction z has the same value or the opposite one; let D be the set of
directions that flip.  The two products of letter matrices differ only in the factors 'dir z on one array' with z in D,
each of which differs by a rank-1 matrix, so   sum_{z in D} n_z  >=  rank( T(m) - T(m') )   for the live/output rows.
(Free per-register shifts at the end do not change the rank: they have equal values at m and m'.)
The best bound is the integer programme  min sum n_z  s.t.  sum_{z in D} n_z >= r_D  for every realisable D.
Ranks are computed in floating point for the search and every r_D used is RE-CHECKED EXACTLY over Z[i]
(Gaussian elimination over Q(i) with Fractions) on its witness pair.   usage: python3 -B lb.py A|B [extra]"""
import sys, itertools, json, os
pass  # numpy/scipy: install with pip (the original run used a local copy)
sys.dont_write_bytecode = True
import numpy as np
from fractions import Fraction as Fr
import toy

def target_matrix_fn(T):
    rows = sorted(T); cols = toy.INPUTS
    # for each (row, col) list of (coef, q)
    ent = [[[(c, q) for (inp, q), c in T[r].items() if inp == l] for l in cols] for r in rows]
    return rows, cols, ent

IPOW = np.array([1, 1j, -1, -1j])
def all_points():
    return np.array(list(itertools.product(range(4), repeat=7)), dtype=np.int64)   # m_1..m_7 (m_0 irrelevant)

def eval_T(ent, M):
    """M: (P,7) multiplicities; returns (P, R, C) complex."""
    P = M.shape[0]; R = len(ent); C = len(ent[0])
    out = np.zeros((P, R, C), dtype=complex)
    for r in range(R):
        for c in range(C):
            for coef, q in ent[r][c]:
                e = (M @ np.array(q[1:], dtype=np.int64)) % 4
                out[:, r, c] += coef * IPOW[e]
    return out

def N(m):  # N_z(m) for z = 1..7
    return [sum(m[y - 1] * toy.chi(z)[y] for y in range(1, 8)) for z in toy.DIRS]

def exact_rank(A):
    """A: list of rows of Gaussian integers given as (re, im) ints -> rank over Q(i)."""
    rows = [[(Fr(a), Fr(b)) for a, b in r] for r in A]
    def mul(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
    def subm(x, y): return (x[0] - y[0], x[1] - y[1])
    def inv(x):
        n = x[0] * x[0] + x[1] * x[1]; return (x[0] / n, -x[1] / n)
    rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] != (0, 0)), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        iv = inv(rows[rk][c]); rows[rk] = [mul(v, iv) for v in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != (0, 0):
                f = rows[i][c]; rows[i] = [subm(a, mul(f, b)) for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def exact_T(ent, m):
    out = []
    for r in range(len(ent)):
        row = []
        for c in range(len(ent[0])):
            re = im = 0
            for coef, q in ent[r][c]:
                e = sum(m[y - 1] * q[y] for y in range(1, 8)) % 4
                re += coef * [1, 0, -1, 0][e]; im += coef * [0, 1, 0, -1][e]
            row.append((re, im))
        out.append(row)
    return out

def bound(T, label):
    rows, cols, ent = target_matrix_fn(T)
    M = all_points(); TM = eval_T(ent, M)
    idx = {tuple(m): k for k, m in enumerate(M.tolist())}
    chiM = np.array([[toy.chi(z)[y] for y in range(1, 8)] for z in toy.DIRS], dtype=np.int64)  # (7 z, 7 y)
    best = {}
    deltas = [d for d in itertools.product(range(4), repeat=7)
              if all(v % 2 == 0 for v in (chiM @ np.array(d)) )]
    for d in deltas:
        nd = (chiM @ np.array(d)) % 4
        if not all(v in (0, 2) for v in nd): raise SystemExit('bad delta')
        D = frozenset(z for z, v in zip(toy.DIRS, nd) if v == 2)
        if not D: continue
        M2 = (M + np.array(d)) % 4
        k2 = np.array([idx[tuple(x)] for x in M2.tolist()])
        diff = TM - TM[k2]
        rk = np.linalg.matrix_rank(diff, tol=1e-7)
        j = int(np.argmax(rk)); r = int(rk[j])
        if r > best.get(D, (-1,))[0]:
            best[D] = (r, tuple(M[j].tolist()), tuple(M2[j].tolist()))
    # exact re-check of every witness
    for D, (r, m1, m2) in best.items():
        A = exact_T(ent, m1); B = exact_T(ent, m2)
        Dm = [[(a[0] - b[0], a[1] - b[1]) for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]
        er = exact_rank(Dm)
        assert er == r, ('exact rank mismatch', D, r, er)
        n1, n2 = N(m1), N(m2)
        assert all((a - b) % 2 == 0 for a, b in zip(n1, n2))
        assert frozenset(z for z, a, b in zip(toy.DIRS, n1, n2) if (a - b) % 4 == 2) == D
    # integer programme
    from scipy.optimize import milp, LinearConstraint, Bounds
    Ds = sorted(best, key=lambda D: (len(D), sorted(D)))
    A = np.array([[1 if z in D else 0 for z in toy.DIRS] for D in Ds]); lb = np.array([best[D][0] for D in Ds])
    res = milp(c=np.ones(7), constraints=LinearConstraint(A, lb, np.inf), integrality=np.ones(7), bounds=Bounds(0, np.inf))
    lp = milp(c=np.ones(7), constraints=LinearConstraint(A, lb, np.inf), bounds=Bounds(0, np.inf))
    print('[%s] realisable flip sets D: %d (sizes %s); max rank per size: %s' % (label, len(best),
          sorted({len(D) for D in best}), {s: max(best[D][0] for D in best if len(D) == s) for s in sorted({len(D) for D in best})}))
    print('[%s] LP bound %.4f ; INTEGER bound on total paid moves: %d   (n_z = %s)' % (label, lp.fun, round(res.fun),
          dict(zip(toy.DIRS, [round(v) for v in res.x]))))
    tight = [(sorted(D), best[D][0]) for D in Ds if abs(sum(res.x[z - 1] for z in D) - best[D][0]) < 1e-9]
    print('[%s] tight constraints at that optimum (first 12): %s' % (label, tight[:12]))
    return round(res.fun), best

if __name__ == '__main__':
    v = sys.argv[1] if len(sys.argv) > 1 else 'A'
    T = toy.targets(v)
    b, best = bound(T, v)
    with open('lb_%s.json' % v, 'w') as f:
        json.dump({'version': v, 'integer_bound': b, 'witnesses': [[sorted(D), r, m1, m2] for D, (r, m1, m2) in best.items()]}, f)
