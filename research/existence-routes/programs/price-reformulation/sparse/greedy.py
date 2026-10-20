#!/usr/bin/env python3
"""Is OpenAI's support for A (symmetric, ker A = E exactly) minimal?  Greedy removal of off-diagonal pairs from the
support of A0 (diagonal always allowed, it is free in their count), keeping a generic element of
{A symmetric, support in G, A K = 0} of full rank w/2.  Several random orders.  Usage: python3 greedy.py d restarts
"""
import sys, time
import numpy as np
from scipy.linalg import null_space
d = int(sys.argv[1]); restarts = int(sys.argv[2])
C = np.array([[0., 0.], [1., 0.]]); Z = np.diag([1., -1.]); X = np.array([[0., 1.], [1., 0.]])
def kron(*ms):
    out = np.array([[1.]])
    for m in ms: out = np.kron(out, m)
    return out
D = C.copy()
for m in range(1, d): D = kron(C, np.eye(2 ** m)) + kron(Z, D)
Gam = kron(X, np.eye(2 ** d)) + kron(Z, D)
N = 2 ** (d + 1); w = N * N; r = w // 2
J = np.array([[1.]])
for i in range(1, d + 2): J = np.kron(J, np.linalg.matrix_power(Z, i - 1) @ X)
L = np.kron(Gam, np.eye(N)) - np.kron(np.eye(N), Gam)
A0 = np.kron(J, J) @ L
K = null_space(L)
P0 = [(i, j) for i in range(w) for j in range(i + 1, w) if A0[i, j] != 0 or A0[j, i] != 0]
rng = np.random.default_rng(0)
def generic_rank(pairs):
    # unknowns: diagonal entries (w) + one per pair; constraints (A K)[i, c] = 0
    cols = [(i, i) for i in range(w)] + list(pairs)
    Mt = np.zeros((w * K.shape[1], len(cols)))
    for k, (i, j) in enumerate(cols):
        # A has entry at (i,j) and (j,i): contributes to row i: K[j,:], row j: K[i,:]
        Mt[i * K.shape[1]:(i + 1) * K.shape[1], k] += K[j, :]
        if i != j:
            Mt[j * K.shape[1]:(j + 1) * K.shape[1], k] += K[i, :]
    V = null_space(Mt)
    if V.shape[1] == 0: return 0
    y = V @ rng.standard_normal(V.shape[1])
    A = np.zeros((w, w))
    for k, (i, j) in enumerate(cols): A[i, j] = A[j, i] = y[k]
    return np.linalg.matrix_rank(A, tol=1e-9 * max(1e-300, np.abs(A).max()))
print(f"d={d} w={w}: |support of A0| = {len(P0)} pairs; generic rank on it = {generic_rank(P0)} (need {r})")
best = len(P0)
for rs in range(restarts):
    t0 = time.time()
    cur = list(P0); order = list(rng.permutation(len(P0)))
    for k in order:
        p = P0[k]
        trial = [q for q in cur if q != p]
        if generic_rank(trial) == r:
            cur = trial
    best = min(best, len(cur))
    print(f"  restart {rs}: removal-minimal support with {len(cur)} pairs ({len(cur)/w:.3f} w) [{time.time()-t0:.1f}s]")
print(f"RESULT d={d}: OpenAI {len(P0)} pairs ({len(P0)/w:.3f} w), smallest removal-minimal subset found {best} ({best/w:.3f} w)")
