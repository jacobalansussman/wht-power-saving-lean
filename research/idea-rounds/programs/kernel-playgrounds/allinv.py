# All monomial-invariant subspaces of W = V (x) Lambda^k V (not only coordinate ones):
#   N  = span{e_x (x) e_T : x in T}                      (the hook)
#   A  = Lambda^(k+1) V  (antisymmetrisations of e_x (x) e_T, x not in T; GL-invariant)
#   C  = the complement of A inside span{x not in T} that S_(k+1) keeps (sum-zero combos per (k+1)-set)
# Every monomial-invariant subspace of W is a sum of some of N, A, C (weight-space argument in the NOTE).
# Theorem R for Y and the paid move / target, ranks mod p.
import sys, itertools
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from hookk import rank_mod, wedge, kron_list
P = 2147483647
def run(K, m, k):
    n = 2 ** m
    K = np.array(K, dtype=object); I = np.eye(2, dtype=object)
    K1 = kron_list([K] + [I] * (m - 1)); T = kron_list([K] * m)
    subsets = list(itertools.combinations(range(n), k)); sidx = {S: i for i, S in enumerate(subsets)}
    basis = [(x, S) for x in range(n) for S in subsets]; bidx = {b: i for i, b in enumerate(basis)}
    D = len(basis)
    vecsN, vecsA, vecsC = [], [], []
    for (x, S) in basis:
        if x in S:
            v = np.zeros(D, dtype=np.int64); v[bidx[(x, S)]] = 1; vecsN.append(v)
    for U in itertools.combinations(range(n), k + 1):
        # e_x (x) e_{U minus x} with the sign that makes sum_x sign * e_x (x) e_{U\x} the antisymmetrisation
        terms = []
        for pos, x in enumerate(U):
            rest = tuple(y for y in U if y != x)
            terms.append((bidx[(x, rest)], (-1) ** pos))
        a = np.zeros(D, dtype=np.int64)
        for (i, s) in terms: a[i] = s % P
        vecsA.append(a)
        for t in range(1, len(terms)):   # complement: differences sign_0 e_0 - sign_t e_t (sum-zero w.r.t. A's pairing)
            c = np.zeros(D, dtype=np.int64)
            c[terms[0][0]] = terms[0][1] % P; c[terms[t][0]] = (-terms[t][1]) % P
            vecsC.append(c)
    blocks = {'N': np.array(vecsN).T, 'A': np.array(vecsA).T, 'C': np.array(vecsC).T}
    mats = {}
    for name, g in (('cap', K1), ('dem', T)):
        gp = np.array([[int(v) % P for v in row] for row in g], dtype=np.int64)
        mats[name] = np.kron(gp, wedge(gp, k, P, subsets)) % P
    out = []
    for combo in (('N',), ('A',), ('C',), ('N', 'A'), ('N', 'C'), ('A', 'C')):
        Y = np.concatenate([blocks[c] for c in combo], axis=1) % P
        dY = rank_mod(Y.T.copy(), P)
        res = {}
        for name in ('cap', 'dem'):
            gY = mats[name].dot(Y) % P
            res[name] = rank_mod(np.concatenate([Y, gY], axis=1).T.copy(), P) - dY
        out.append(('+'.join(combo), dY, res['cap'], res['dem']))
    return out
for m, k in ((2, 2), (3, 4)):
    for name, K in (('H', [[1, 1], [1, -1]]), ('L', [[1, 0], [1, 1]])):
        for (combo, dY, cap, dem) in run(K, m, k):
            print(f"m={m} k={k} {name} Y={combo:5s} dimY={dY:4d} cap={cap:4d} dem={dem:4d} ratio={(dem/cap if cap else float('nan')):.4f}", flush=True)
