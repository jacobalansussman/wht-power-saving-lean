#!/usr/bin/env python3
"""vcore.py (v-t2-block4, skeptic): my own block-4 machinery, written without reading t2-block4's helpers' internals.
Conventions (mine): a direction is a 0/1 numpy vector of length m (index 1..15 <-> integer with bit i = coordinate i).
A binary quadratic form Q is an upper-triangular 0/1 m x m matrix; Q(x) = x^T Q x mod 2.
A symmetric matrix S is a full 0/1 m x m matrix; the pairing <Q, S> = sum_{i<=j} Q_ij S_ij mod 2, so <Q, x x^T> = Q(x).
Targets: S = h h^T over h in GL(m, 2). O_m = {g : g g^T = I} acts on targets by S -> g S g^T."""
import itertools
import numpy as np

def setup(m=4):
    dirs = np.array([[(x >> i) & 1 for i in range(m)] for x in range(1, 1 << m)], dtype=np.int64)   # (15, m)
    ut = [(i, j) for i in range(m) for j in range(i, m)]
    nq = len(ut)
    forms = []
    for bits in itertools.product([0, 1], repeat=nq):
        Q = np.zeros((m, m), dtype=np.int64)
        for b, (i, j) in zip(bits, ut): Q[i, j] = b
        forms.append(Q)
    forms = forms[1:]                     # drop Q = 0
    Qs = np.array(forms)                  # (1023, m, m)
    # support matrix: P[q, x] = Q_q(x)
    P = (np.einsum('xi,qij,xj->qx', dirs, Qs, dirs) % 2).astype(np.int64)
    # all invertible matrices over F_2 (as full 0/1 arrays)
    allm = np.array(list(itertools.product([0, 1], repeat=m * m)), dtype=np.int64).reshape(-1, m, m)
    def rank2(A):
        A = A.copy() % 2; r = 0; n, k = A.shape
        for c in range(k):
            piv = [i for i in range(r, n) if A[i, c]]
            if not piv: continue
            p = piv[0]; A[[r, p]] = A[[p, r]]
            for i in range(n):
                if i != r and A[i, c]: A[i] ^= A[r]
            r += 1
        return r
    # determinant over F_2 via float det parity is unsafe; use rank
    GL = np.array([g for g in allm if rank2(g) == m])
    return dict(m=m, dirs=dirs, Qs=Qs, P=P, GL=GL, rank2=rank2, ut=ut)

def pairing(Qs, S):
    m = S.shape[0]
    iu = np.triu_indices(m)
    return (Qs[:, iu[0], iu[1]] * S[iu[0], iu[1]]).sum(1) % 2

def key(S):
    return tuple(int(v) for v in S.ravel())

def targets_and_orth(D):
    GL = D['GL']; m = D['m']; I = np.eye(m, dtype=np.int64)
    T = {}
    for h in GL:
        S = (h @ h.T) % 2
        T.setdefault(key(S), S)
    O = [g for g in GL if np.array_equal((g @ g.T) % 2, I)]
    return T, O
