# fulltouch4.py -- (fulltouch2.py in a NON-ABELIAN letter group of order 8) Lemma Q'' accounting for fresh-touch sequences on TWO two-row arrays (packing), exact in the regular
# representation of Gamma = <u1, u2> (|Gamma| = 4).  Question: does fully touching two packed arrays ever cost < 16?
import numpy as np, sys
from vn2 import rho, OIDX, NIDX, E_of, rank, nullsp
from blk3 import D
from fam import POOL
from fulltouch import is_inv, close, idx
E = E_of('R')
def run(els, seq, label):
    n = len(els)
    lam = []
    for g in els:
        L = np.zeros((n, n))
        for j, h in enumerate(els): L[idx(els, g @ h), j] = 1
        lam.append(L)
    RD = [rho(D[1] @ g)[:, NIDX] @ E for g in els]
    rows = np.array([r * n + h for r in OIDX for h in range(n)])
    Ub = np.eye(2 * 136 * n, dtype=complex); tot = 0; parts = []
    for (b1, b2) in seq:
        Z = np.zeros((560 * n, 136 * n), dtype=complex)
        M1 = sum((c * np.kron(RD[i], lam[i]) for i, c in b1.items()), Z.copy())
        M2 = sum((c * np.kron(RD[i], lam[i]) for i, c in b2.items()), Z.copy())
        M = np.hstack([M1, M2]); MU = M @ Ub
        r = rank(MU[rows]); tot += r; parts.append(round(r / n, 3))
        Ub = Ub @ nullsp(MU)
        if Ub.shape[1] == 0: break
    print(f'{label}: prices per point {parts}, total {tot / n:.3f} (two arrays), untouched left {Ub.shape[1] / n:.3f} of 272', flush=True)
invs = [u for u in POOL if is_inv(u)]
rng = np.random.default_rng(int(sys.argv[1]))
for trial in range(300):
    u1, u2 = invs[int(rng.integers(len(invs)))], invs[int(rng.integers(len(invs)))]
    if np.allclose(u1 @ u2, u2 @ u1): continue
    els = close([u1, u2], cap=8)
    if els is None or len(els) != 8: continue
    break
print('group order', len(els), flush=True)
i0, i1, i2 = 0, idx(els, u1), idx(els, u2); i12 = idx(els, u1 @ u2)
O = {}
seqs = {
  'plain: x1, x2': [({i0: 1}, O), (O, {i0: 1})],
  'packed x1 + x2, then x1': [({i0: 1}, {i0: 1}), ({i0: 1}, O)],
  'packed x1 + u1 x2, then x1 - u1 x2': [({i0: 1}, {i1: 1}), ({i0: 1}, {i1: -1})],
  'packed (1+u1)(x1 + x2), then x1, then x2': [({i0: 1, i1: 1}, {i0: 1, i1: 1}), ({i0: 1}, O), (O, {i0: 1})],
  'x1 + u1 x2, then u2 x1 + x2': [({i0: 1}, {i1: 1}), ({i2: 1}, {i0: 1})],
  '(1+u1) x1 + (1-u1) x2, then (1-u1) x1 + (1+u1) x2, then x1, x2': [({i0: 1, i1: 1}, {i0: 1, i1: -1}), ({i0: 1, i1: -1}, {i0: 1, i1: 1}), ({i0: 1}, O), (O, {i0: 1})],
  'random pair of rows': [({i: complex(rng.normal(), rng.normal()) for i in range(8)}, {i: complex(rng.normal(), rng.normal()) for i in range(8)}),
                          ({i: complex(rng.normal(), rng.normal()) for i in range(8)}, {i: complex(rng.normal(), rng.normal()) for i in range(8)})],
}
for lab, sq in seqs.items(): run(els, sq, lab)
