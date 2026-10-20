"""v-a2-zx-compute small checks (numpy, < 50 MB)."""
import numpy as np, itertools
rng = np.random.default_rng(7)
H = np.array([[1, 1], [1, -1]]) / np.sqrt(2); S = np.diag([1, 1j]); X = np.array([[0, 1], [1, 0]]); Z = np.diag([1, -1])
I2 = np.eye(2)
a, b = (1 + 1j) / 2, (1 - 1j) / 2
D = np.array([[a, b], [b, a]])
def eqp(A, B):  # equal up to a scalar
    k = np.unravel_index(np.argmax(abs(B)), B.shape); s = A[k] / B[k]
    return abs(abs(s) - 1) < 1e-9 and np.allclose(A, s * B)
print("D == H S H exactly:", np.allclose(D, H @ S @ H))
print("S D S == e^{i pi/4} H:", np.allclose(S @ D @ S, np.exp(1j * np.pi / 4) * H))
print("D^2 == X:", np.allclose(D @ D, X))

# spider nest: prod_{v in U\0} dir(v) = I for dim U >= 3 (on m = 4 bits, f = 1)
def kronlist(L):
    out = np.eye(1)
    for M in L: out = np.kron(out, M)
    return out
def Xv(v): return kronlist([X if t else I2 for t in v])
def dirv(v): n = len(v); return a * np.eye(2 ** n) + b * Xv(v)
m = 4
def span(basis):
    S_ = set()
    for coeffs in itertools.product([0, 1], repeat=len(basis)):
        v = tuple(sum(c * bb[i] for c, bb in zip(coeffs, basis)) % 2 for i in range(m)); S_.add(v)
    return [v for v in S_ if any(v)]
for d, basis in [(2, [(1,0,0,0),(0,1,0,0)]), (3, [(1,0,0,0),(0,1,0,0),(0,0,1,0)]), (3, [(1,1,0,0),(0,1,1,0),(0,0,1,1)]), (4, [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)])]:
    P = np.eye(2 ** m)
    for v in span(basis): P = dirv(v) @ P
    print("spider nest dim %d: product == I: %s (== scalar*I: %s)" % (d, np.allclose(P, np.eye(2 ** m)), eqp(P, np.eye(2 ** m))))

# f = 1 degeneracy realised as a W = 2 scheme with 0 paid moves: slot (x) data, free = slot maps (any 2x2 invertible)
# and per-array monomial letters E_k (x) g.  CZ(s,d) = E0(x)I + E1(x)Z, CNOT(s->d) = E0(x)I + E1(x)X (per-array letters),
# CNOT(d->s) = (Hs(x)I) CZ (Hs(x)I) (slot map + per-array letter), SWAP = CNOT CNOT CNOT, I(x)H = SWAP (Hs(x)I) SWAP.
E0 = np.diag([1, 0]); E1 = np.diag([0, 1])
CZsd = np.kron(E0, I2) + np.kron(E1, Z); CXsd = np.kron(E0, I2) + np.kron(E1, X)
Hs = np.kron(H, I2); CXds = Hs @ CZsd @ Hs
SW = CXsd @ CXds @ CXsd
print("SWAP built from free pieces is SWAP:", np.allclose(SW, np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])))
print("I_2 (x) H == SWAP (H_slot (x) I) SWAP:", np.allclose(SW @ Hs @ SW, np.kron(I2, H)))
print("hence I_2 (x) D free at f = 1 (D = H S H, S per-array letter on both arrays):", np.allclose((SW @ Hs @ SW) @ np.kron(I2, S) @ (SW @ Hs @ SW), np.kron(I2, D)))

# Theorem 1 support count on random circuits: monomial gates + t one-qubit gates -> column support <= 2^t
def colsupp(A): return int((abs(A) > 1e-12).sum(axis=0).max())
n = 5; worst = 0
for trial in range(300):
    U = np.eye(2 ** n, dtype=complex); t = 0
    for step in range(25):
        if rng.random() < 0.3:
            q = rng.integers(n); g = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
            U = kronlist([g if k == q else I2 for k in range(n)]) @ U; t += 1
        else:
            perm = rng.permutation(2 ** n); ph = np.exp(2j * np.pi * rng.random(2 ** n))
            Pm = np.zeros((2 ** n, 2 ** n), dtype=complex); Pm[perm, np.arange(2 ** n)] = ph
            U = Pm @ U
    worst = max(worst, colsupp(U) / 2 ** min(t, n))
    assert colsupp(U) <= 2 ** t
print("Theorem 1 support bound s <= 2^t held on 300 random circuits; max s/2^min(t,n) =", worst)

# XW end-column letters are free and NOT monomial: at f = 1 and f = 2 every column is an end column, so H (and hence
# sqrtX = H S H) on every wire is a product of free boundary letters: K costs 0 paid moves there.
for f in (1, 2):
    for mm in (2, 3):
        n = mm * f
        Hall = kronlist([H] * n); Sall = kronlist([S] * n)
        print("XW, f=%d, m=%d: K = (H S H)^(x)n from end-column letters only: %s" % (f, mm, np.allclose(Hall @ Sall @ Hall, kronlist([D] * n))))
# Lemma 2 algebra: random commuting rank-1 idempotents are equal or annihilating; non-commuting ones break fusion at pi/4
from scipy.linalg import expm
W = 3
def r1(): 
    u = rng.normal(size=W) + 1j * rng.normal(size=W); v = rng.normal(size=W) + 1j * rng.normal(size=W); return np.outer(u, v.conj()) / (v.conj() @ u)
Q1, Q2 = r1(), r1()
Pd = Xv((1, 1)); one = np.eye(4)
A_ = np.kron(Q1, one - Pd); B_ = np.kron(Q2, one - Pd)
print("non-commuting rank-1 weights: fusion exp(i pi/4 A)exp(i pi/4 B) == exp(i pi/4 (A+B)):", np.allclose(expm(1j*np.pi/4*A_) @ expm(1j*np.pi/4*B_), expm(1j*np.pi/4*(A_+B_))))
# commuting orthogonal pair -> fused weight rank 2 idempotent; equal pair -> D^2 = P
T = rng.normal(size=(W, W)) + 1j * rng.normal(size=(W, W)); Ti = np.linalg.inv(T)
F1 = T @ np.diag([1, 0, 0]) @ Ti; F2 = T @ np.diag([0, 1, 0]) @ Ti
G = lambda Q: expm(1j * np.pi / 4 * np.kron(Q, one - Pd))
print("orthogonal pair fuses, fused weight idempotent rank 2:", np.allclose(G(F1) @ G(F2), G(F1 + F2)), np.linalg.matrix_rank(F1 + F2))
print("equal pair fuses to slot-weighted P (free):", np.allclose(G(F1) @ G(F1), np.kron(np.eye(W) - F1, one) + np.kron(F1, Pd)))
