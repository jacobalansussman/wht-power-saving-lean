#!/usr/bin/env python3
"""ctr-search / xcheck.py: independent cross-checks.
(1) profile survivors by brute force (combinations_with_replacement, no pruning) for m <= 16, same bounds as
    profiles.py but re-implemented here;
(2) the one-sided reduction on random exact two-sided O1 decoders: build x_b = L c_b + n_b, y_a = L' c_a + n'_a
    (L'^T L = I, n_b in ker L'^T, n'_a in ker L^T and orthogonal to all n_b), check exactness, apply the reduction
    (C^T = A Phi, new scatter A^T C^T), check exactness again with identical gather supports."""
import sys, os, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *
n = {1: 1, 2: 2, 3: 4, 4: 8, 5: 12, 6: 20, 7: 36, 8: 64, 9: 112}
cf = json.load(open(os.path.join(HERE, 'coflat.json')))
mG = {int(k): v for k, v in cf['Gamma_a'].items()}; mG[9] = 53
mN = {int(k): v for k, v in cf['N_i'].items()}; mN[9] = 56
hk = [0, 1, 3, 6, 12, 20, 36]
def passes(p):
    m = len(p); s = sorted(p)
    if sum(n[d] for d in s) < 224 - m: return False
    for k in range(7):
        if sum(n[d] for d in s[:m - k]) < 112 - hk[k]: return False
    J = s[:m - 7]
    if not J: return False
    return sum(mG[d] for d in J) >= 53 or sum(mN[d] for d in J) >= 56 or sum(n[d] for d in J) >= 70
surv = {}
for m in range(1, 17):
    S = [p for p in itertools.combinations_with_replacement(range(1, 10), m) if sum(p) <= 63 and passes(p)]
    if S: surv[m] = S
print('(1) brute-force survivors:', {m: len(v) for m, v in surv.items()})
ref = json.load(open(os.path.join(HERE, 'profiles.json')))
same = all(sorted(map(tuple, ref.get(str(m), []))) == sorted(surv.get(m, [])) for m in range(1, 17))
print('    identical to profiles.py:', same)
# (2) reduction demo
rng = np.random.default_rng(5)
C = np.array(LINES, dtype=float).T / 2          # 8 x 112, c_b = r_b / 2
Ebool = E
for trial in range(5):
    m = 11
    L = rng.normal(size=(m, 8)); Lp = L @ np.linalg.inv(L.T @ L)        # Lp^T L = I
    KLp = np.linalg.svd(Lp.T)[2][8:].T          # basis of ker Lp^T (m x 3)
    KL = np.linalg.svd(L.T)[2][8:].T            # basis of ker L^T
    V = KLp[:, :1]                              # gather noise space (dim 1), inside ker Lp^T
    # scatter noise: inside ker L^T and orthogonal to V
    Q = KL - V @ np.linalg.lstsq(V, KL, rcond=None)[0] if False else KL
    Vp = KL @ np.linalg.svd(V.T @ KL)[2][1:].T  # vectors of ker L^T orthogonal to V
    X = L @ C + V @ rng.normal(size=(1, N))
    Y = Lp @ C + Vp[:, :1] @ rng.normal(size=(1, N))
    M = Y.T @ X
    ex1 = np.max(np.abs((M - C.T @ C)[Ebool]))
    A = C @ np.linalg.pinv(X)                   # C = A X  (8 x m)
    Y2 = A.T @ C
    ex2 = np.max(np.abs((Y2.T @ X - C.T @ C)[Ebool]))
    print(f'(2) trial {trial}: two-sided exact err {ex1:.1e}; after reduction err {ex2:.1e}; rank X {np.linalg.matrix_rank(X)}; '
          f'C = A X residual {np.max(np.abs(A @ X - C)):.1e}')
