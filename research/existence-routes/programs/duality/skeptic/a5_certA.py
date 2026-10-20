# v-ex2-duality a5_certA.py: independent checks (numpy complex; all entries are dyadic Gaussian rationals, so
# float arithmetic here is exact; we still compare with tolerance 1e-12).
# (1) A5: dir_e1 = S^1 K S^1 K^-1 S^3 with K = dir_e1...dir_em, K^-1 = K X_all, dir_z = ((1+i) I + (1-i) X_z)/2,
#     for m = 1..5, and K S_1 K^-1 = (dir S dir^-1) (x) I.
# (2) Certificate A (permutation module on nonzero Pauli labels): f(g) = |g L0 g^-1 \ L0| for g = rank-r block at
#     width w, computed from the matrices by conjugating Z-type Paulis; compare with 2^w - 2^(w-r); then the
#     spreading relations' balances at several widths.
import sys
import numpy as np
from itertools import product
from functools import reduce

def kron(*ms): return reduce(np.kron, ms)
I2 = np.eye(2, dtype=complex); X = np.array([[0, 1], [1, 0]], dtype=complex); Z = np.diag([1, -1]).astype(complex)
D = ((1 + 1j) * I2 + (1 - 1j) * X) / 2
S = np.diag([1, 1j])
def on(m, j, A):  # A on bit j (bit 0 = first tensor factor)
    return kron(*[A if t == j else I2 for t in range(m)])
for m in range(1, 6):
    K = reduce(np.matmul, [on(m, j, D) for j in range(m)])
    Xall = kron(*[X] * m)
    Kinv = K @ Xall
    ok_inv = np.allclose(K @ Kinv, np.eye(2 ** m), atol=1e-12)
    S1 = on(m, 0, S)
    W = S1 @ K @ S1 @ Kinv @ np.linalg.matrix_power(S1, 3)
    ok_a5 = np.allclose(W, on(m, 0, D), atol=1e-12)
    ok_fac = np.allclose(K @ S1 @ Kinv, on(m, 0, D @ S @ np.linalg.inv(D)), atol=1e-12)
    print(f"A5 m={m}: K*K^-1=I {ok_inv}; dir_e1 = S1 K S1 K^-1 S1^3 exactly {ok_a5}; K S1 K^-1 = (dir S dir^-1) x I {ok_fac}")

def pauli(m, x, z):
    return kron(*[np.linalg.matrix_power(X, (x >> t) & 1) @ np.linalg.matrix_power(Z, (z >> t) & 1) for t in range(m)])
def label_of(m, P):
    # find (x, z) and phase with P = c * pauli(x, z)
    for x in range(2 ** m):
        for z in range(2 ** m):
            Q = pauli(m, x, z)
            c = np.trace(Q.conj().T @ P) / 2 ** m
            if abs(abs(c) - 1) < 1e-9 and np.allclose(P, c * Q, atol=1e-9): return (x, z)
    return None
for w in (3, 4):
    L0 = [(0, z) for z in range(1, 2 ** w)]
    for r in range(1, w + 1):
        g = reduce(np.matmul, [on(w, j, D) for j in range(r)])  # rank-r block on rows 0..r-1
        gi = np.linalg.inv(g)
        img = set(label_of(w, g @ pauli(w, x, z) @ gi) for (x, z) in L0)
        moved = len(img - set(L0))
        print(f"cert A width {w}: rank-{r} block moves {moved} of the {2**w - 1} labels of L0 (formula {2**w - 2**(w - r)})")
def certA(w, r): return 2 ** w - 2 ** (w - r)
rels = {"Rs 2[2]<=[1]+[3]": ({2: 2}, {1: 1, 3: 1}), "R2 2[2]<=1.05[1]+[3]": ({2: 2}, {1: 1.05, 3: 1}),
        "R2b 2[4]<=[1]+[7]": ({4: 2}, {1: 1, 7: 1}), "R1 2[3]<=5[1]": ({3: 2}, {1: 5}),
        "A5 [1]<=2[m=4]": ({1: 1}, {4: 2}), "Kron [5]<=[2]+[3]": ({5: 1}, {2: 1, 3: 1})}
for name, (T, P) in rels.items():
    for w in (max(list(T) + list(P)), 8, 12):
        paid = sum(p * certA(w, r) for r, p in P.items()); tasks = sum(t * certA(w, k) for k, t in T.items())
        print(f"cert A, {name}, width {w}: paid {paid:.3f} tasks {tasks:.3f} -> {'allowed' if paid >= tasks - 1e-9 else 'EXCLUDED'}")
