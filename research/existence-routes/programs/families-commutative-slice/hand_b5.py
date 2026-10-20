# hand_b5.py: a symmetric certificate for Theorem 18? For every even 3-flat F at block 5 and each of its 21 ten-piece
# schemes, bound each piece by the AVERAGE slack over all odd 3-flats containing it; sum over the schemes of F (uniformly)
# and over F with one weight per S_5-orbit of F; ask (LP) for orbit weights making every direction's coefficient <= 1
# with the largest constant: sum_z c_z >= constant.
import itertools
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
pts = list(range(1, 32)); basis = {1, 2, 4, 8, 16}
def span_dirs(S):
    out = {0}
    for v in S: out |= {o ^ v for o in out}
    return out
flats3 = set()
for a in pts:
    for u, v, w in itertools.combinations(range(1, 32), 3):
        D = span_dirs([u, v, w])
        if len(D) != 8: continue
        F = frozenset(a ^ d for d in D)
        if 0 not in F: flats3.add(F)
flats3 = list(flats3)
odd = [F for F in flats3 if len(F & basis) % 2 == 1]; even = [F for F in flats3 if len(F & basis) % 2 == 0]
print("3-flats avoiding 0:", len(flats3), "odd", len(odd), "even", len(even))
through = {}
def avg_cover(U):
    U = frozenset(U)
    if U not in through:
        G = [Q for Q in odd if U <= Q]
        co = {z: Fr(0) for z in pts}
        for Q in G:
            for z in Q: co[z] += Fr(1, len(G))
        through[U] = co
    return through[U]
def schemes(F):
    P = sorted(F); a0 = P[0]; D = sorted({p ^ a0 for p in P} - {0}); out = []
    subs2 = {frozenset([0, u, w, u ^ w]) for u, w in itertools.combinations(D, 2)}
    for S2 in subs2:
        H = frozenset(p for p in P if (p ^ a0) in S2); Hp = frozenset(P) - H
        for t in sorted(S2 - {0}):
            def lines(X):
                o, seen = [], set()
                for p in sorted(X):
                    if p in seen: continue
                    l = frozenset([p, p ^ t]); seen |= l; o.append(l)
                return o
            K, L = lines(H); Kp, Lp = lines(Hp)
            out.append((H, Hp, [(X, Y) for X in (K, L) for Y in (Kp, Lp)]))
    return out
def contrib(F):
    co = {z: Fr(0) for z in pts}; sch = schemes(F)
    for (H, Hp, cross) in sch:
        for z in F: co[z] -= Fr(1, len(sch))
        for X in (H, Hp):
            for z, v in avg_cover(X).items(): co[z] += v / len(sch)
        for (X, Y) in cross:
            for z, v in avg_cover(X | Y).items(): co[z] += 2 * v / len(sch)
    return co
def orbit_key(F):
    best = None
    for p in itertools.permutations(range(5)):
        img = tuple(sorted(sum(((z >> j) & 1) << p[j] for j in range(5)) for z in F))
        best = img if best is None or img < best else best
    return best
orbits = {}
for F in even: orbits.setdefault(orbit_key(F), []).append(F)
print("S_5-orbits of even 3-flats:", len(orbits), [len(v) for v in orbits.values()])
cols = []
for key, Fs in orbits.items():
    tot = {z: Fr(0) for z in pts}
    for F in Fs:
        for z, v in contrib(F).items(): tot[z] += v
    cols.append(tot)
# LP: weights lam_o >= 0; constraints sum_o lam_o * coef_o(z) <= 1 for all z; maximise 10 * sum_o lam_o * |orbit o|
nO = len(cols)
A = np.array([[float(cols[o][z]) for o in range(nO)] for z in pts])
gain = np.array([10.0 * len(list(orbits.values())[o]) for o in range(nO)])
res = linprog(-gain, A_ub=A, b_ub=np.ones(len(pts)), bounds=[(0, None)] * nO, method="highs")
print("best bound from symmetric averaged schemes:", -res.fun, " (Theorem 18's LP: 155/36 =", 155 / 36, ")")
print("orbit weights:", [round(v, 6) for v in res.x])
