"""clifford_rank.py -- in the engine game, is anything 'sparse but costly'?

Run: python3 clifford_rank.py   (numpy, seconds)

OpenAI's gap object Gamma has price (d/2) N with only (1 + d/2) N nonzeros: price linear in sparsity.
In the engine game a single array's operator is a Clifford U.  Claim [Bruhat decomposition w.r.t. the Siegel
parabolic, standard]: U = (monomial Clifford) H^(x)r (monomial Clifford), and r = log2(nonzeros per row).  So the
single-array cost equals log2 of the sparsity: sparse means exponentially cheap.

This script checks the two facts the lower half rests on, on random circuits at m = 2..5 qubits:
  (1) every free engine letter (CNOT, CZ, S, X, Z) leaves nonzeros-per-row unchanged;
  (2) every H changes log2(nonzeros per row) by exactly -1, 0 or +1;
and checks nonzeros per row = 2^rank_F2(X-part of U Z_j U^dagger) (the Bruhat index read off the tableau).
"""
import numpy as np

rng = np.random.default_rng(7)
H1 = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
S1 = np.diag([1, 1j])
X1 = np.array([[0, 1], [1, 0]])
Z1 = np.diag([1, -1])


def op1(g, j, m):
    out = np.array([[1]])
    for k in range(m):
        out = np.kron(out, g if k == j else np.eye(2))
    return out


def cnot(c, t, m):
    n = 2 ** m
    Mx = np.zeros((n, n))
    for x in range(n):
        bits = [(x >> (m - 1 - k)) & 1 for k in range(m)]
        if bits[c]:
            bits[t] ^= 1
        y = sum(b << (m - 1 - k) for k, b in enumerate(bits))
        Mx[y, x] = 1
    return Mx


def cz(a, b, m):
    n = 2 ** m
    dg = [(-1) ** (((x >> (m - 1 - a)) & 1) * ((x >> (m - 1 - b)) & 1)) for x in range(n)]
    return np.diag(dg)


def nnz_row(U):
    cnt = (np.abs(U) > 1e-9).sum(axis=1)
    assert len(set(cnt.tolist())) == 1
    return int(cnt[0])


def xrank_of_stabilizers(U, m):
    # X-part of U Z_j U^dag: read the Pauli by its action on basis states (support pattern of the column 0 image)
    rows = []
    for j in range(m):
        Pm = U @ op1(Z1, j, m) @ U.conj().T
        # Pauli string's X-part = the bit flip x -> x xor v, read from where column 0 is sent
        col = np.abs(Pm[:, 0])
        v = int(np.argmax(col))
        rows.append([(v >> (m - 1 - k)) & 1 for k in range(m)])
    A = np.array(rows, dtype=np.int64) % 2
    r = 0
    A = A.copy()
    for c in range(m):
        piv = [i for i in range(r, m) if A[i, c]]
        if not piv:
            continue
        A[[r, piv[0]]] = A[[piv[0], r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
    return r


tot_free = tot_h = viol = 0
mism = 0
for m in range(2, 6):
    for trial in range(40):
        U = np.eye(2 ** m)
        cur = 0
        for step in range(40):
            kind = rng.integers(0, 6)
            if kind == 0:
                g = op1(H1, rng.integers(m), m)
            elif kind == 1:
                g = op1(S1, rng.integers(m), m)
            elif kind == 2:
                g = op1(X1, rng.integers(m), m)
            elif kind == 3:
                g = op1(Z1, rng.integers(m), m)
            elif kind == 4:
                a, b = rng.choice(m, 2, replace=False)
                g = cnot(a, b, m)
            else:
                a, b = rng.choice(m, 2, replace=False)
                g = cz(a, b, m)
            U = g @ U
            r_new = int(round(np.log2(nnz_row(U))))
            if kind == 0:
                tot_h += 1
                if abs(r_new - cur) > 1:
                    viol += 1
            else:
                tot_free += 1
                if r_new != cur:
                    viol += 1
            cur = r_new
        if xrank_of_stabilizers(U, m) != cur:
            mism += 1
print(f"free letters applied: {tot_free}; H letters applied: {tot_h}; violations of (1)/(2): {viol}")
print(f"nonzeros per row = 2^(F2 rank of X-part of U Z U^dag) on all {4*40} final circuits: mismatches {mism}")
