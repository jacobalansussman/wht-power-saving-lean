"""v-rm2-phase: own base code (written from the definitions in rm-theory 8.5-8.7, not from rm2-phase's files).
Block m. Directions z in F_2^m minus 0 (index = integer z-1). Forms Q: coefficient bits q_ij (i<=j),
Q(z) = sum_i q_ii z_i + sum_{i<j} q_ij z_i z_j. Pairing <Q,S> = sum_{i<=j} q_ij S_ij (linear in S), so
<Q, z z^T> = Q(z). Symmetric matrices S as m x m 0/1 numpy arrays."""
import itertools, numpy as np

def pairs(m):
    return [(i, j) for i in range(m) for j in range(i, m)]

def all_forms(m):
    P = pairs(m)
    return [dict(zip(P, bits)) for bits in itertools.product((0, 1), repeat=len(P))]

def vec(z, m):
    return np.array([(z >> i) & 1 for i in range(m)], dtype=np.int64)

def Qz_matrix(m):
    P = pairs(m); nf = 2 ** len(P)
    M = np.zeros((nf, 2 ** m - 1), dtype=np.int64)
    for f, bits in enumerate(itertools.product((0, 1), repeat=len(P))):
        for z in range(1, 2 ** m):
            v = vec(z, m)
            M[f, z - 1] = sum(b * v[i] * v[j] for b, (i, j) in zip(bits, P)) % 2
    return M

def pair_matrix(m, mats):
    """rows: forms, cols: matrices; <Q,S> = sum q_ij S_ij."""
    P = pairs(m)
    F = np.array(list(itertools.product((0, 1), repeat=len(P))), dtype=np.int64)  # nf x |P|
    Sv = np.array([[S[i, j] for (i, j) in P] for S in mats], dtype=np.int64)   # nt x |P|
    return (F @ Sv.T) % 2

def rank2(A):
    A = A.copy() % 2; r = 0; rows, cols = A.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]: piv = i; break
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(rows):
            if i != r and A[i, c]: A[i] ^= A[r]
        r += 1
    return r

def sym_nonsing(m):
    P = pairs(m); out = []
    for bits in itertools.product((0, 1), repeat=len(P)):
        S = np.zeros((m, m), dtype=np.int64)
        for b, (i, j) in zip(bits, P):
            S[i, j] = S[j, i] = b
        if rank2(S) == m: out.append(S)
    return out

def gl(m):
    out = []
    for cols in itertools.product(range(1, 2 ** m), repeat=m):
        A = np.array([vec(c, m) for c in cols], dtype=np.int64).T
        if rank2(A) == m: out.append(A)
    return out

def key(S):
    return tuple(S.flatten())
