#!/usr/bin/env python3
"""lib3.py (t2-block3, own code): block-3 engine letters as 8x8 complex matrices, symplectic images, membership in the
free group B (monomial Clifford letters), the Levi group B cap K B K^-1, and the f = 5 hook module M = V (x) Lambda^4 V.

Engine conventions (as theory-guide/hb.py and t514-search/hooklib.py): dir z = (1+i)/2 I + (1-i)/2 X^z, K = dir1 dir2 dir4,
B = < Paulis, affine address maps, phase letters i^{(w.x+c) mod 2} > (all monomial with entries in <i> up to a scalar).
Entries of all group elements lie in Z[i][1/2]; float keys rounded at 1e-9 are exact in practice (finite group).
"""
import itertools
import numpy as np

N8 = 8
I8 = np.eye(N8, dtype=complex)


def X(a):
    M = np.zeros((N8, N8), dtype=complex)
    for x in range(N8):
        M[x ^ a, x] = 1
    return M


def Z(b):
    return np.diag([(-1) ** (bin(b & x).count("1") & 1) for x in range(N8)]).astype(complex)


CA, CB = (1 + 1j) / 2, (1 - 1j) / 2


def dirz(z):
    return CA * I8 + CB * X(z)


K = dirz(1) @ dirz(2) @ dirz(4)
Kinv = np.linalg.inv(K)


def perm_mat(f):
    """permutation matrix sending e_x to e_{f(x)}"""
    M = np.zeros((N8, N8), dtype=complex)
    for x in range(N8):
        M[f(x), x] = 1
    return M


def phase_letter(w, c=0):
    return np.diag([1j ** ((bin(w & x).count("1") + c) & 1) for x in range(N8)])


def apply_lin(A, x):
    """A: 3x3 0/1 matrix (list of rows), x: int 0..7 (bit j = coordinate j)"""
    v = [(x >> j) & 1 for j in range(3)]
    y = 0
    for i in range(3):
        if sum(A[i][j] * v[j] for j in range(3)) % 2:
            y |= 1 << i
    return y


def all_gl3():
    out = []
    for bits in itertools.product([0, 1], repeat=9):
        A = [list(bits[0:3]), list(bits[3:6]), list(bits[6:9])]
        # invertible over F2 iff the images of the basis span
        ims = {apply_lin(A, x) for x in range(8)}
        if len(ims) == 8:
            out.append(A)
    return out


def key(M, scale_free=True):
    """hash key; with scale_free the matrix is normalised by its first nonzero entry (projective class)"""
    if scale_free:
        flat = M.flatten()
        k = np.nonzero(np.abs(flat) > 1e-9)[0][0]
        M = M / flat[k]
    R = np.round(M.real * 2 ** 12).astype(np.int64)
    Iq = np.round(M.imag * 2 ** 12).astype(np.int64)
    return R.tobytes() + Iq.tobytes()


def is_monomial(M, tol=1e-9):
    nz = np.abs(M) > tol
    return bool((nz.sum(axis=0) == 1).all() and (nz.sum(axis=1) == 1).all())


def pauli(xv, zv):
    return X(xv) @ Z(zv)


def sympl(g):
    """symplectic image of a Clifford g as a 6x6 F2 matrix on labels (x|z) (column = image of a basis Pauli);
    returns None if g is not Clifford"""
    gi = np.linalg.inv(g)
    cols = []
    basis = [(1 << j, 0) for j in range(3)] + [(0, 1 << j) for j in range(3)]
    for (xv, zv) in basis:
        P = g @ pauli(xv, zv) @ gi
        found = None
        for a in range(8):
            for b in range(8):
                Q = pauli(a, b)
                # P = lambda Q  <=> tr(Q^dag P) = 8 lambda with |lambda| = 1
                t = np.trace(Q.conj().T @ P) / 8
                if abs(abs(t) - 1) < 1e-9:
                    found = (a, b)
                    break
            if found:
                break
        if found is None:
            return None
        a, b = found
        cols.append([(a >> j) & 1 for j in range(3)] + [(b >> j) & 1 for j in range(3)])
    return np.array(cols, dtype=np.int64).T % 2


def in_B(g):
    """free letter (up to a scalar): monomial, entries of equal modulus, ratios powers of i, and Clifford"""
    if not is_monomial(g):
        return False
    vals = g[np.abs(g) > 1e-9]
    r = vals / vals[0]
    if not np.allclose(np.abs(r), 1):
        return False
    ang = np.angle(r) / (np.pi / 2)
    if not np.allclose(ang, np.round(ang), atol=1e-9):
        return False
    return sympl(g) is not None


# ---------------------------------------------------------------- f = 5 hook module ------------------------------
S4 = list(itertools.combinations(range(N8), 4))
BASIS = [(x, T) for x in range(N8) for T in S4]
NIDX = np.array([i for i, (x, T) in enumerate(BASIS) if x in T])


def lambda4(g):
    U = np.array(S4)
    sub = g[U[:, None, :, None], U[None, :, None, :]]   # 70,70,4,4
    return np.linalg.det(sub)


def rho(g):
    return np.kron(g, lambda4(g))


def colspace(A, tol=1e-8):
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    r = int((s > tol * max(1.0, s[0])).sum())
    return U[:, :r]


def intersect(A, Bm, tol=1e-8):
    """basis of col(A) cap col(B) (A, B orthonormal columns)"""
    M = np.hstack([A, -Bm])
    U, s, Vh = np.linalg.svd(M)
    null = Vh[np.sum(s > tol):].conj().T if len(s) else None
    r = int(np.sum(s > tol))
    nullv = Vh[r:].conj().T
    v = A @ nullv[:A.shape[1]]
    return colspace(v) if v.shape[1] else v
