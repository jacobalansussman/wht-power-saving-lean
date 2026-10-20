# playgrounds/patterns.py -- exhaustive search over "pattern" certificates in V^(x)f at small n.
# A basis vector e_x, x in [n]^f, has the pattern = set partition of the f columns by equal values. Y_S = span{e_x :
# pattern(x) in S} is kept by every monomial matrix of V (all permutations of states, all diagonals). Theorem R:
# moves per array >= dem/cap. Exhaustive over all subsets S of patterns (small f), for kernels H and L.
# Usage: python3 -I patterns.py m f [Sym]   (Sym: also restrict to the S_f-antisymmetric-in-last-(f-1) part? no)
import sys, itertools
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from hookk import rank_mod, kron_list
P = 2147483647

def canon(x):
    lab = {}; out = []
    for v in x:
        if v not in lab: lab[v] = len(lab)
        out.append(lab[v])
    return tuple(out)

def run(K, m, f):
    n = 2 ** m
    K = np.array(K, dtype=object); I = np.eye(2, dtype=object)
    K1 = kron_list([K] + [I] * (m - 1)); T = kron_list([K] * m)
    xs = list(itertools.product(range(n), repeat=f))
    pats = sorted(set(canon(x) for x in xs)); pid = np.array([pats.index(canon(x)) for x in xs])
    mats = {}
    for name, g in (('cap', K1), ('dem', T)):
        gp = np.array([[int(v) % P for v in row] for row in g], dtype=np.int64)
        M = np.array([[1]], dtype=np.int64)
        for _ in range(f): M = np.kron(M, gp) % P
        mats[name] = M
    return mats, pid, pats

if __name__ == '__main__':
    m, f = int(sys.argv[1]), int(sys.argv[2])
    for name, K in (('H', [[1, 1], [1, -1]]), ('L', [[1, 0], [1, 1]])):
        mats, pid, pats = run(K, m, f)
        best = []
        npat = len(pats)
        for r in range(1, npat):
            for S in itertools.combinations(range(npat), r):
                inY = np.isin(pid, S)
                cap = rank_mod(mats['cap'][np.ix_(~inY, inY)], P)
                if cap == 0: continue
                dem = rank_mod(mats['dem'][np.ix_(~inY, inY)], P)
                best.append((dem / cap, cap, dem, [pats[i] for i in S]))
        best.sort(key=lambda t: -t[0])
        print(f"{name} m={m} f={f} #patterns={npat}: top ratios:", [(round(b[0], 4), b[1], b[2]) for b in best[:3]], flush=True)
        print("   best S:", best[0][3], flush=True)
