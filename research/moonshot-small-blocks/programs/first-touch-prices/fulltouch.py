# fulltouch.py -- Lemma Q'' accounting of a SEQUENCE of fresh touches of ONE two-row array, computed exactly in the regular
# representation of a small letter group Gamma = <u1, u2> (the f -> infinity model; averaged over residue classes of k
# when Gamma holds -1).  Hook, Y = N, D = dir(1).  Reports the total price per point and whether the array is fully touched.
import numpy as np, sys, itertools
from vn2 import rho, OIDX, NIDX, E_of, rank, nullsp
from blk3 import D
from fam import POOL
E = E_of('R')
def is_inv(u): return np.allclose(u @ u, np.eye(8)) and not np.allclose(u, np.eye(8)) and not np.allclose(u, -np.eye(8))
def close(gens, cap=32):
    els = [np.eye(8, dtype=complex)]
    keys = {np.round(els[0], 6).tobytes()}
    fr = list(els)
    while fr:
        nx = []
        for a in fr:
            for g in gens:
                b = g @ a; k = np.round(b, 6).tobytes()
                if k not in keys:
                    keys.add(k); els.append(b); nx.append(b)
                    if len(els) > cap: return None
        fr = nx
    return els
def idx(els, g):
    for i, e in enumerate(els):
        if np.allclose(e, g, atol=1e-8): return i
def run(els, seq, label):
    n = len(els)
    lam = []
    for g in els:
        L = np.zeros((n, n))
        for j, h in enumerate(els): L[idx(els, g @ h), j] = 1
        lam.append(L)
    RD = [rho(D[1] @ g)[:, NIDX] @ E for g in els]     # 560 x 136 each
    Ub = np.eye(136 * n, dtype=complex); tot = 0; parts = []
    for beta in seq:
        M = sum(c * np.kron(RD[i], lam[i]) for i, c in beta.items())          # (560 n) x (136 n), rows (x, h) ordering
        # out-of-N rows: rows of rho index r in OIDX, for every group coordinate
        rows = np.array([r * n + h for r in OIDX for h in range(n)])
        MU = M @ Ub
        r = rank(MU[rows]) if Ub.shape[1] else 0
        tot += r; parts.append(round(r / n, 3))
        k = nullsp(MU) if Ub.shape[1] else np.zeros((Ub.shape[1], 0))
        Ub = Ub @ k
        if Ub.shape[1] == 0: break
    print(f'{label}: |Gamma| {n}, prices per point {parts}, total {tot / n:.3f}, untouched left {Ub.shape[1] / n:.3f} (of 136)', flush=True)
if __name__ == '__main__':
    invs = [u for u in POOL if is_inv(u)]
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    done = 0
    for trial in range(200):
        u1, u2 = invs[int(rng.integers(len(invs)))], invs[int(rng.integers(len(invs)))]
        if np.allclose(u1, u2): continue
        els = close([u1, u2], cap=16)
        if els is None or len(els) < 4: continue
        i1, i2, i0 = idx(els, u1), idx(els, u2), 0
        i12 = idx(els, u1 @ u2)
        seqs = {
          '(1+u1), then 1': [{i0: 1, i1: 1}, {i0: 1}],
          '(1-u1), then 1': [{i0: 1, i1: -1}, {i0: 1}],
          '(1+u1), (1+u2), then 1': [{i0: 1, i1: 1}, {i0: 1, i2: 1}, {i0: 1}],
          '(1+u1)(1+u2), then (1-u1), then 1': [{i0: 1, i1: 1, i2: 1, i12: 1}, {i0: 1, i1: -1}, {i0: 1}],
          'projector of Gamma, then 1': [{i: 1.0 / len(els) for i in range(len(els))}, {i0: 1}],
          'random invertible beta': [{i: complex(rng.normal(), rng.normal()) for i in range(len(els))}],
        }
        for lab, sq in seqs.items():
            run(els, sq, lab)
        done += 1
        if done >= int(sys.argv[2] if len(sys.argv) > 2 else 2): break
