#!/usr/bin/env python3
"""vplain.py (v-rm-singer): exact (integer) check of the canonical Singer system on the plain schedule, own code.
System (Lemma C of e3-rmbuild, all paid moves in direction z0 = 1):
   X_Q = Atil Pi_Q (I + Gamma Pi_Q)^-1 Btil  must equal  I_W (x) diag(o_Q),
   Pi_Q = I_M (x) diag(pi), pi_k = Q(a^-k) (move channel k), o_Q[k] = <Q, s^-k s^-kT> (from the matrix).
Circulant convention: P|k> = |k+1> (left-regular action of the generator).  My derivation (RUNLOG): a move in
direction w = a^r is P^r [1 1^T] P^-r, so Atil_t = P^r_t, Btil_t = P^-r_t, Gamma_{s,t} = 2 P^(r_t - r_s) (s later).
Checks every nonzero form (not only orbit representatives).
"""
import sys, itertools
import numpy as np
sys.path.insert(0, '.')
from vfield import setup, channel_data


def circ_pow(n, r):
    P = np.zeros((n, n), dtype=np.int64)
    for k in range(n): P[(k + r) % n, k] = 1
    return P


def build(n, rs):
    M = len(rs)
    A = np.zeros((n, M * n), dtype=np.int64); B = np.zeros((M * n, n), dtype=np.int64)
    G = np.zeros((M * n, M * n), dtype=np.int64)
    for t, r in enumerate(rs):
        A[:, t * n:(t + 1) * n] = circ_pow(n, r)
        B[t * n:(t + 1) * n, :] = circ_pow(n, -r)
        for s in range(t + 1, M):
            G[s * n:(s + 1) * n, t * n:(t + 1) * n] = 2 * circ_pow(n, r - rs[s])
    return A, B, G


def unitri_solve(L, R):
    # L unit lower triangular (integer), solve L X = R exactly with Python ints
    X = R.astype(np.int64).copy()
    for i in range(L.shape[0]):
        X[i, :] -= L[i, :i] @ X[:i, :]
    assert np.abs(X).max() < 2 ** 40
    return X


def check(m, coords, rs, label):
    F, basis, C, sig, Sk, pairs = setup(m, coords)
    n = F.n; N = len(pairs); M = len(rs)
    A, B, G = build(n, rs)
    bad_forms = 0; good_forms = 0; maxerr = 0
    for q in itertools.product([0, 1], repeat=N):
        if not any(q): continue
        K, J = channel_data(F, C, Sk, pairs, q)
        pi = np.array(K * M, dtype=np.int64)
        Pi = np.diag(pi)
        L = np.eye(M * n, dtype=np.int64) + G @ Pi
        X = A @ (Pi @ unitri_solve(L, B))
        T = np.diag(np.array(J, dtype=np.int64))
        err = int(np.abs(X - T).max())
        maxerr = max(maxerr, err)
        if err: bad_forms += 1
        else: good_forms += 1
    print(f"m={m} {coords:4s} {label:38s} r={rs}: forms exact {good_forms}, violated {bad_forms}, max|err| {maxerr}",
          flush=True)


def main():
    for m in (3, 4):
        for coords in ('sd', 'poly'):
            F, basis, C, sig, Sk, pairs = setup(m, coords)
            # basis vectors that are orthonormal for the dot product in these coordinates:
            if coords == 'sd':
                ortho = basis                      # self-dual basis = coordinate basis
            else:
                ortho = [F.pw(i) for i in range(m)]  # 1, a, ..., a^(m-1) = coordinate basis
            rs = [F.log[w] for w in ortho]
            check(m, coords, rs, "plain (orthonormal basis)")
            check(m, coords, rs[:-1], "broken: last move dropped")
            other = [F.pw(i) for i in range(m)] if coords == 'sd' else selfdual_basis(F)
            check(m, coords, [F.log[w] for w in other], "broken: basis not orthonormal here")


def selfdual_basis(F):
    from vfield import selfdual
    return selfdual(F)


if __name__ == '__main__':
    main()
