# twostep_quick.py -- non-Levi involutions with zero-divisor price 4 (invcost.py seed 4), second letters = pool + 200
# random widened letters (a subset: gives an UPPER bound on the best cost2)
import numpy as np
from vn2 import rho, OIDX, NIDX, E_of, rank, nullsp
from blk3 import D, mono, rand_widened, rand_engine
from fam import POOL
E = E_of('R'); RD = rho(D[1]); P = RD[np.ix_(OIDX, NIDX)]
def rhoN(g): R = rho(g); return R[np.ix_(NIDX, NIDX)]
def is_inv(u): return np.allclose(u @ u, np.eye(8)) and not np.allclose(u, np.eye(8)) and not np.allclose(u, -np.eye(8))
def inter(Ab, Bb):
    M = np.hstack([Ab, -Bb]); k = nullsp(M)
    return Ab @ k[:Ab.shape[1]] if k.shape[1] else np.zeros((Ab.shape[0], 0))
rng = np.random.default_rng(4)
invs = [u for u in POOL if is_inv(u)]
for _ in range(300):
    a = POOL[int(rng.integers(len(POOL)))]; b = POOL[int(rng.integers(len(POOL)))]
    u = a @ b @ np.linalg.inv(a)
    if is_inv(u): invs.append(u)
for _ in range(300):
    perm = np.arange(8); pairs = rng.permutation(8)
    for j in range(0, int(rng.integers(1, 5)) * 2, 2): perm[pairs[j]], perm[pairs[j + 1]] = pairs[j + 1], pairs[j]
    u = mono(list(perm), [0] * 8)
    for x in range(8):
        if perm[x] == x: u[x, x] = (-1) ** int(rng.integers(2))
        elif x < perm[x]:
            s = int(rng.integers(4)); u[perm[x], x] = 1j ** s; u[x, perm[x]] = 1j ** (-s)
    if is_inv(u): invs.append(u)
cands = []
for u in invs:
    Ru = rhoN(u); A1 = P @ E; Au = P @ Ru @ E
    for sgn in (1, -1):
        c1 = 0.5 * (rank(A1 + sgn * Au) + rank(A1 - sgn * Au))
        if c1 <= 4.0:
            K = inter(E, np.linalg.inv(Ru) @ E)
            if K.shape[1] < 136: cands.append((c1, sgn, u, K))
print('non-Levi price-4 candidates:', len(cands), [c[3].shape[1] for c in cands], flush=True)
r2 = np.random.default_rng(1)
V = list(POOL) + [rand_widened(r2) for _ in range(150)] + [rand_engine(r2) for _ in range(50)]
for (c1, sgn, u, K) in cands[:4]:
    Ru = rhoN(u); best = 99
    for v in V:
        Rv = rhoN(v)
        c2 = 0.5 * rank(np.vstack([P @ Rv @ K, P @ Rv @ Ru @ K]))
        best = min(best, c2)
    print('non-Levi u: cost1', c1, 'k', K.shape[1], 'best cost2 over', len(V), 'letters', best, 'TOTAL', c1 + best, flush=True)
