"""feedback_state.py -- can the dense-pair feedback be lifted to f digit columns inside a unit?

Run: python3 feedback_state.py   (numpy; exact integer ranks mod a large prime; well under a minute)

The dense-pair lemma closes the 4x4 cell L = [[I, J],[I, 0]] (J = 4-cycle with weight lam) at z = 1 to get I + J.
A unit is right for every f, so the same trick would have to act on every digit column of a pair of directions,
i.e. produce (I + J)^(x)f.

Part A [proved by the Hankel argument, checked here]: the stream whose closure is (I+J)^(x)f per column, i.e. the
  transfer function (I + zJ)^(x)f, has McMillan degree exactly f * 4^f: its block Hankel matrix is block
  anti-triangular with the invertible block J^(x)f on the anti-diagonal.  So any state-space realisation carries
  f data-arrays of state.  A unit's array count may not grow with f.
Part B [computed]: the per-column cascade (OpenAI's 5.5) attains f * 4^f and gives (I+zJ)^(x)f; the single-delay
  closure of the simultaneous power L^(x)f (state 8^f - 4^f) does NOT give (I+J)^(x)f.
Part C [computed]: the reflection's sparse side is not stable under f-fold lifting: nnz per row of the generator of
  E_+(M^(x)f) grows like 2 (1 + d/2)^f, while the tensor price per coordinate of M^(x)f is f d.
"""
import numpy as np
from itertools import product

P = 2147483629
lam = 3


def Jcyc():
    J = np.zeros((4, 4), dtype=np.int64)
    for i in range(3):
        J[i, i + 1] = 1
    J[3, 0] = lam
    return J


def kron_list(ms):
    out = np.array([[1]], dtype=np.int64)
    for m in ms:
        out = np.kron(out, m)
    return out


def rank_modp(Mat):
    A = (Mat % P).astype(np.int64)
    r, c = A.shape
    rank = 0
    for col in range(c):
        nzr = np.nonzero(A[rank:, col] % P)[0]
        if len(nzr) == 0:
            continue
        piv = rank + nzr[0]
        A[[rank, piv]] = A[[piv, rank]]
        inv = pow(int(A[rank, col]), P - 2, P)
        A[rank] = (A[rank] * inv) % P
        f = A[:, col].copy()
        f[rank] = 0
        nz = np.nonzero(f)[0]
        if len(nz):
            # entries < 2^31, products < 2^62: safe in int64
            A[nz] = (A[nz] - (f[nz, None] * A[rank][None, :]) % P) % P
        rank += 1
        if rank == r:
            break
    return rank


def partA():
    print("Part A: McMillan degree of (I + zJ)^(x)f  (rank of the block Hankel matrix of its Markov parameters)")
    J = Jcyc()
    I4 = np.eye(4, dtype=np.int64)
    for f in range(1, 5):
        n = 4 ** f
        Pk = [np.zeros((n, n), dtype=np.int64) for _ in range(f + 1)]
        for S in product([0, 1], repeat=f):
            Pk[sum(S)] += kron_list([J if s else I4 for s in S])
        Hk = np.zeros((f * n, f * n), dtype=np.int64)
        for i in range(f):
            for j in range(f):
                k = i + j + 1
                if k <= f:
                    Hk[i * n:(i + 1) * n, j * n:(j + 1) * n] = Pk[k]
        rk = rank_modp(Hk)
        print(f"  f={f}: data width 4^f = {n:5d};  Hankel rank = {rk:5d};  f*4^f = {f*n:5d};  state arrays needed = {rk/n:.0f}")


def closure(K, ndata, z):
    M = K[:ndata, :ndata]; C = K[:ndata, ndata:]; B = K[ndata:, :ndata]; D = K[ndata:, ndata:]
    r = D.shape[0]
    return M + z * C @ np.linalg.solve(np.eye(r) - z * D, B)


def partB():
    print("\nPart B: simultaneous power versus per-column cascade (lam = 3, floats)")
    J = Jcyc().astype(float)
    I4 = np.eye(4)
    Z4 = np.zeros((4, 4))
    L = np.block([[I4, J], [I4, Z4]])          # data h (first 4), state l (last 4)
    print(f"  f=1: closure of L at z=1 equals I+J: {np.allclose(closure(L, 4, 1.0), I4 + J)}")
    for f in (2, 3):
        Lf = L
        for _ in range(f - 1):
            Lf = np.kron(Lf, L)
        # reorder: data = all factors in h (index < 4 in each factor)
        idx = list(product(range(8), repeat=f))
        data = [i for i, t in enumerate(idx) if all(x < 4 for x in t)]
        state = [i for i, t in enumerate(idx) if not all(x < 4 for x in t)]
        perm = data + state
        Lp = Lf[np.ix_(perm, perm)]
        H1 = closure(Lp, len(data), 1.0)
        target = I4 + J
        T = target
        for _ in range(f - 1):
            T = np.kron(T, target)
        print(f"  f={f}: simultaneous L^(x)f: state {len(state)} = 8^f - 4^f;  closure at z=1 equals (I+J)^(x)f: "
              f"{np.allclose(H1, T)};  max |difference| = {np.abs(H1 - T).max():.3g}")
    # per-column cascade for f = 2: state W1 (state in column 1), W2 (state in column 2); one shared delay
    f = 2
    n = 16
    # cell along column j acts on (data x, W_j): x -> x + J_j w_j ; w_j' = x
    J1 = np.kron(J, I4); J2 = np.kron(I4, J)
    # one tick: K1 on (x, W1), then K2 on (x', W2); assemble the full tick on (x, W1, W2)
    Id = np.eye(n); Zn = np.zeros((n, n))
    K1 = np.block([[Id, J1, Zn], [Id, Zn, Zn], [Zn, Zn, Id]])
    K2 = np.block([[Id, Zn, J2], [Zn, Id, Zn], [Id, Zn, Zn]])
    tick = K2 @ K1
    for z in (1.0, 0.37):
        H = closure(tick, n, z)
        T = np.kron(I4 + z * J, I4 + z * J)
        print(f"  f=2 per-column cascade (state 2*16 = 32): closure at z={z} equals (I+zJ)^(x)2: {np.allclose(H, T)}")


def partC():
    print("\nPart C: f-fold lift of the reflection.  nnz per row of the generator S(G^(x)f (x) I - I (x) G^(x)f) versus")
    print("        the tensor price per coordinate of (G(x)G)^(x)f (= f d).  G = OpenAI's Gamma, nnz per row 1 + d/2.")
    for d in (30, 60):
        cells = []
        for f in (1, 2, 3, 4):
            gen = 2 * (1 + d / 2) ** f
            price = f * d
            cells.append(f"f={f}: generator {gen:.3g}, price {price}, ratio h/price {gen/2/price:.3g}")
        print(f"  d={d}:  " + ";  ".join(cells))
    print("  The margin needs h/price < 2/3 (pairs at 3/2).  It holds only at f = 1.")


if __name__ == "__main__":
    partA()
    partB()
    partC()
