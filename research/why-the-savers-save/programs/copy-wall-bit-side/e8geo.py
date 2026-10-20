"""e8geo.py (d-bitspan) -- E8 lines in the hub's s2-decoder ('geo') numbering, exact flats, cliques.
hamming_e8 is re-implemented here from the hub's e8units.py (same construction and order; checked against it by
check_numbering() when the hub copy is importable).  Everything exact (integers)."""
import itertools, math
import numpy as np

def hamming_e8():
    gen = np.array([[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 1, 1, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1],
                    [0, 1, 0, 1, 0, 1, 0, 1]])
    code = {tuple((np.array(c) @ gen) % 2) for c in itertools.product([0, 1], repeat=4)}
    w4 = sorted(c for c in code if sum(c) == 4)
    L = [tuple(2 * (np.arange(8) == a).astype(int)) for a in range(8)]
    for c in w4:
        idx = [i for i in range(8) if c[i]]
        for sg in itertools.product([1, -1], repeat=3):
            x = np.zeros(8, dtype=int); x[idx[0]] = 1
            for k in range(3): x[idx[k + 1]] = sg[k]
            L.append(tuple(x))
    L = np.array(L, dtype=np.int64)
    G = (L @ L.T) // 2
    assert len(L) == 120 and set(np.diag(G)) == {2}
    return L, G

L, G = hamming_e8()
NL = 120
PUNCT_S2 = [39, 52, 57, 71, 76, 78, 103, 109, 114]     # hub s2-decoder REPORT 'shipped bit puncture' (geo numbering)
PUNCT_WORD = [8, 17, 39, 51, 53, 58, 63, 72, 114]      # #589 word puncture [10,15,51,53,56,66,72,89,111] mapped (map_check.py)

def bits(mask):
    out = []; i = 0
    while mask:
        if mask & 1: out.append(i)
        mask >>= 1; i += 1
    return out

def _red(v):
    g = 0
    for x in v: g = math.gcd(g, abs(int(x)))
    if g > 1: v = v // g
    nz = np.nonzero(v)[0]
    if len(nz) and v[nz[0]] < 0: v = -v
    return v

def _normrows(C):
    g = np.gcd.reduce(np.abs(C), axis=1)
    g[g == 0] = 1
    C = C // g[:, None]
    first = np.argmax(C != 0, axis=1)
    s = np.sign(C[np.arange(len(C)), first]); s[s == 0] = 1
    return C * s[:, None]

def flats(maxrank=7, verbose=False):
    """All flats (sets of E8 lines = the lines inside a Q-subspace spanned by lines), rank 0..maxrank, exactly
    (integer normals; a line joins F+l iff its image in Q^8/F is parallel to l's).  Returns {rank: {mask: N}}."""
    import time
    out = {0: {0: np.eye(8, dtype=np.int64)}}
    LM = [1 << i for i in range(NL)]
    for r in range(1, maxrank + 1):
        t0 = time.time(); new = {}
        for mask, N in out[r - 1].items():
            C = _normrows(L @ N.T)
            nz = np.any(C != 0, axis=1)
            idx = np.nonzero(nz)[0]
            keys, inv = np.unique(C[idx], axis=0, return_inverse=True)
            inv = np.asarray(inv).reshape(-1)
            groups = [[] for _ in range(len(keys))]
            for a, gi in zip(idx.tolist(), inv.tolist()): groups[gi].append(a)
            for mem in groups:
                m2 = mask
                for i in mem: m2 |= LM[i]
                if m2 in new: continue
                l = L[mem[0]]
                t = N @ l
                nzt = [j for j in range(len(t)) if t[j] != 0]
                piv = min(nzt, key=lambda j: abs(t[j]))
                rows = [t[piv] * N[j] - t[j] * N[piv] for j in range(len(N)) if j != piv]
                N2 = _normrows(np.array(rows, dtype=np.int64).reshape(len(rows), 8)) if rows else np.zeros((0, 8), np.int64)
                new[m2] = N2
        out[r] = new
        if verbose: print('  flats rank %d: %d (%.1f s)' % (r, len(new), time.time() - t0), flush=True)
    return out

def check_flat_exact(mask, N):
    """exact: lines in mask are exactly the lines orthogonal to every row of N, and rank(lines) = 8 - len(N)."""
    inside = np.all((L @ N.T) == 0, axis=1) if len(N) else np.ones(NL, bool)
    ok1 = set(np.nonzero(inside)[0].tolist()) == set(bits(mask))
    from fractions import Fraction
    M = [[Fraction(int(x)) for x in L[i]] for i in bits(mask)]
    rk = 0
    if M:
        A = [row[:] for row in M]
        for c in range(8):
            p = next((k for k in range(rk, len(A)) if A[k][c] != 0), None)
            if p is None: continue
            A[rk], A[p] = A[p], A[rk]
            for k in range(len(A)):
                if k != rk and A[k][c] != 0:
                    f = A[k][c] / A[rk][c]
                    A[k] = [a - f * b for a, b in zip(A[k], A[rk])]
            rk += 1
    return ok1 and rk == 8 - len(N)

def cliques(live, size=8):
    """all sets of `size` pairwise NON-orthogonal live lines (bit interference cliques), as tuples."""
    adj = {i: sum(1 << j for j in live if j != i and G[i, j] != 0) for i in live}
    out = []
    def rec(cur, cand):
        if len(cur) == size:
            out.append(tuple(cur)); return
        c = cand
        while c:
            low = c & -c; i = low.bit_length() - 1; c ^= low
            if bin(cand >> i).count('1') + len(cur) < size: break
            rec(cur + [i], cand & adj[i] & ~((low << 1) - 1))
    rec([], sum(1 << i for i in live))
    return out

def f2_labels():
    """9-bit F_2 labels with U U^T = (1 + G) mod 2 (hub sym_factor_f2 construction, re-implemented)."""
    M = ((1 + G) % 2).astype(np.int64)
    # greedy independent rows
    basis, I = {}, []
    rows = [sum(int(M[i, j]) << j for j in range(NL)) for i in range(NL)]
    for i, x in enumerate(rows):
        y = x
        while y:
            h = y.bit_length() - 1
            if h in basis: y ^= basis[h]
            else: basis[h] = y; I.append(i); break
    N = M[np.ix_(I, I)]
    c = len(I)
    # orthonormal basis of N over F_2 by brute force search (c = 9: 511 vectors)
    vecs = [np.array([(x >> k) & 1 for k in range(c)]) for x in range(1, 1 << c)]
    def search(B):
        if len(B) == c: return B
        for v in vecs:
            if int(v @ N @ v) % 2 == 1 and all(int(v @ N @ b) % 2 == 0 for b in B):
                r = search(B + [v])
                if r: return r
        return None
    B = np.array(search([]))
    U = (M[:, I] @ B.T) % 2
    assert ((U @ U.T) % 2 == M).all()
    return U
