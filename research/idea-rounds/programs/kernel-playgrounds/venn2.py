# playgrounds/venn2.py -- monomial-invariant "Venn" certificates in Lambda^a V (x) Lambda^b V, V = C^(2^m).
# Kept subspace Y_S = span{ e_A (x) e_B : |A n B| in S } is kept by every monomial matrix (permutations of the
# column states and all diagonals), for every kernel. Theorem R: moves per array >= dem/cap.
# The hook of the project is (a, b, S) = (1, f-1, {1}).
# Usage: python3 -I venn2.py m   -> prints, for each kernel, the best (a,b,S) found and the hook for comparison.
import sys, itertools
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from hookk import rank_mod, wedge, kron_list

P = 2147483647

def run(K, m, a, b, p=P, maxdim=3000):
    n = 2 ** m
    K = np.array(K, dtype=object); I = np.eye(2, dtype=object)
    K1 = kron_list([K] + [I] * (m - 1)); T = kron_list([K] * m)
    SA = list(itertools.combinations(range(n), a)); SB = list(itertools.combinations(range(n), b))
    basis = [(A, B) for A in SA for B in SB]
    if len(basis) > maxdim: return None
    inter = np.array([len(set(A) & set(B)) for (A, B) in basis])
    vals = sorted(set(inter.tolist()))
    mats = {}
    for name, g in (('cap', K1), ('dem', T)):
        gp = np.array([[int(v) % p for v in row] for row in g], dtype=np.int64)
        mats[name] = np.kron(wedge(gp, a, p, SA), wedge(gp, b, p, SB)) % p
    out = []
    for r in range(1, len(vals)):
        for S in itertools.combinations(vals, r):
            inY = np.isin(inter, S)
            cap = rank_mod(mats['cap'][np.ix_(~inY, inY)], p)
            dem = rank_mod(mats['dem'][np.ix_(~inY, inY)], p)
            out.append((dem / cap if cap else 0.0, a, b, S, cap, dem, int(inY.sum())))
    return out

if __name__ == '__main__':
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    n = 2 ** m
    KS = {'H': [[1, 1], [1, -1]], 'L': [[1, 0], [1, 1]]}
    for name, K in KS.items():
        allres = []
        for a in range(1, n):
            for b in range(a, n):
                if a + b > n + 2: continue
                res = run(K, m, a, b)
                if res is None: continue
                allres += res
                best = max(res)
                print(f"{name} m={m} a={a} b={b}: best ratio {best[0]:.4f} S={best[3]} cap={best[4]} dem={best[5]} dimY={best[6]}", flush=True)
        allres.sort(reverse=True)
        print(f"== {name} m={m}: top 5:", [(round(r[0], 4),) + r[1:] for r in allres[:5]], flush=True)
