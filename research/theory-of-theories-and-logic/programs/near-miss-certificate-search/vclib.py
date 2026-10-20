# vclib.py -- v-cert-07 (skeptic of cert-07), written from scratch.
# Module M = Hom_Pauli(V, V (x) Lam^8 V), V = C^16 (block 4).  Basis f_T, T an 8-subset of F_2^4, with
# f_T(e_0) = e_sig(T) (x) e_T (sig = XOR of T).  For a Clifford G = g (x) Lam^8 g, G = C (x) rho(G) with C = lam*g,
# and comparing the coefficient of e_sig(T') (x) e_T' in G(e_sig(T) (x) e_T) (only the y = 0 Pauli translate can
# contribute, since sig(T'+y) = sig(T')) gives
#     rho(G)[T', T] = g[sig T', sig T] * det g[T', T] / (lam * g[0,0]).
# Depth element (DIFFERENT from cert-07's Hadamard): g_j = d^(x)j (x) I on address bits 0..j-1, with
# d = ((1+i)/2) I + ((1-i)/2) X (maps the Z-frame to the Y-frame on that bit).  Exact: d ~ [[1,-i],[-i,1]] (scalar
# dropped), entries in Z[i]; mod p (p = 1 mod 4) with i -> sqrt(-1) mod p.  g_j is block diagonal over the 2^(4-j)
# contiguous address blocks, so det g[T',T] = prod_z det(Dj[T'_z, T_z]) when |T'_z| = |T_z| for all z, else 0, and
# rho_j is block diagonal over these compositions.
import itertools
import numpy as np

SUBS = list(itertools.combinations(range(16), 8))
NT = len(SUBS)
IDX = {s: k for k, s in enumerate(SUBS)}
SA = np.array(SUBS)
SIG = np.array([int(np.bitwise_xor.reduce(np.array(s))) for s in SUBS])
HOOK = np.array([SIG[k] in SUBS[k] for k in range(NT)])

def orbits_gl42():
    """GL(4,2) orbits on 8-subsets via all transvections x_i += x_j; sorted by (size, min lexicographic index)."""
    gens = []
    for i in range(4):
        for j in range(4):
            if i != j:
                gens.append(np.array([x ^ ((((x >> j) & 1)) << i) for x in range(16)]))
    parent = list(range(NT))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for p in gens:
        im = np.sort(p[SA], axis=1)
        for k in range(NT):
            a, b = find(k), find(IDX[tuple(int(v) for v in im[k])])
            if a != b: parent[a] = b
    groups = {}
    for k in range(NT): groups.setdefault(find(k), []).append(k)
    orb = sorted(groups.values(), key=lambda v: (len(v), min(v)))
    lab = np.empty(NT, dtype=int)
    for q, v in enumerate(orb): lab[v] = q
    return orb, lab, gens

def sqrtm1(p):
    for g in range(2, 200):
        if pow(g, (p - 1) // 2, p) == p - 1:
            return pow(g, (p - 1) // 4, p)
    raise ValueError

def det_modp(A, p):
    A = [list(r) for r in A]; n = len(A); d = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] % p), None)
        if piv is None: return 0
        if piv != c: A[c], A[piv] = A[piv], A[c]; d = -d
        d = d * A[c][c] % p; inv = pow(A[c][c], p - 2, p)
        for r in range(c + 1, n):
            f = A[r][c] * inv % p
            if f:
                for k in range(c, n): A[r][k] = (A[r][k] - f * A[c][k]) % p
    return d % p

def element(j, mode, p=None):
    """Dj (2^j x 2^j): float unitary (mode 'f') or mod-p integer matrix (mode 'p', scalar dropped)."""
    if mode in ('f', 'r'):
        if mode == 'f':
            d = np.array([[(1 + 1j) / 2, (1 - 1j) / 2], [(1 - 1j) / 2, (1 + 1j) / 2]])
        else:   # real rotation r = S^-1 d S up to phase = [[1,1],[-1,1]]/sqrt2 (Z-frame -> X-frame), diagonal-conjugate of d
            d = np.array([[1.0, 1.0], [-1.0, 1.0]]) / np.sqrt(2)
        D = np.array([[1.0 + 0j]]) if mode == 'f' else np.array([[1.0]])
        for _ in range(j): D = np.kron(D, d)   # bit 0 is the last kron factor? see note below
        return D
    I = sqrtm1(p)
    d = np.array([[1, (-I) % p], [(-I) % p, 1]], dtype=object)
    D = np.array([[1]], dtype=object)
    for _ in range(j): D = np.kron(D, d) % p
    return D
# note: d^(x)j is symmetric under permuting its factors, so the bit order inside a block does not matter.

def blocks(j, mode, p=None):
    n = 1 << j; Z = 1 << (4 - j)
    D = element(j, mode, p)
    masks = range(1 << n); pc = [bin(c).count('1') for c in masks]
    if mode in ('f', 'r'):
        Mn = np.zeros((1 << n, 1 << n), dtype=complex if mode == 'f' else float)
    else:
        Mn = np.zeros((1 << n, 1 << n), dtype=np.int64)
    for c1 in masks:
        r = [i for i in range(n) if c1 >> i & 1]
        for c2 in masks:
            if pc[c2] != pc[c1]: continue
            cc = [i for i in range(n) if c2 >> i & 1]
            if not r:
                Mn[c1, c2] = 1
            elif mode in ('f', 'r'):
                Mn[c1, c2] = np.linalg.det(D[np.ix_(r, cc)])
            else:
                Mn[c1, c2] = det_modp([[int(D[a][b]) for b in cc] for a in r], p)
    codes = np.zeros((NT, Z), dtype=np.int64)
    for z in range(Z):
        inb = (SA >> j) == z
        codes[:, z] = (inb * (1 << (SA & (n - 1)))).sum(axis=1)
    comp = [tuple(bin(int(c)).count('1') for c in row) for row in codes]
    cls = {}
    for k, c in enumerate(comp): cls.setdefault(c, []).append(k)
    Df = D if mode in ('f', 'r') else np.array(D, dtype=np.int64)
    out = {}
    for c, rows in cls.items():
        rows = np.array(rows); cd = codes[rows]; sl = SIG[rows] & (n - 1)
        if mode in ('f', 'r'):
            G = Df[sl[:, None], sl[None, :]].copy()
            for z in range(Z): G = G * Mn[cd[:, z][:, None], cd[:, z][None, :]]
            G = G / Df[0, 0]     # lam = 1 up to a phase; rho must come out unitary (checked)
        else:
            G = Df[sl[:, None], sl[None, :]] % p
            for z in range(Z): G = (G * Mn[cd[:, z][:, None], cd[:, z][None, :]]) % p
        out[c] = (rows, G)
    return out

def rank_modp(A, p):
    A = np.array(A, dtype=np.int64) % p
    if A.shape[0] > A.shape[1]: A = A.T.copy()
    rws, cls_ = A.shape; r = 0
    for c in range(cls_):
        if r == rws: break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        k = r + nz[0]
        if k != r: A[[r, k]] = A[[k, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r, c:] = (A[r, c:] * inv) % p
        f = A[r + 1:, c]
        idx = np.nonzero(f)[0]
        if len(idx):
            rr = r + 1 + idx
            A[rr, c:] = (A[rr, c:] - (A[rr, c][:, None] * A[r, c:][None, :]) % p) % p
        r += 1
    return r

def psi_float(B, inS, tol=1e-7):
    psi = 0; phi = 0.0; lo = 9.0; hi = 0.0
    for c, (rows, G) in B.items():
        m = inS[rows]
        if m.all() or (~m).all(): continue
        s = np.linalg.svd(G[np.ix_(~m, m)], compute_uv=False)
        psi += int((s > tol).sum()); phi += float((s ** 2).sum())
        if (s > tol).any(): lo = min(lo, float(s[s > tol].min()))
        if (s <= tol).any(): hi = max(hi, float(s[s <= tol].max()))
    return psi, phi, lo, hi

def psi_modp(B, inS, p):
    psi = 0
    for c, (rows, G) in B.items():
        m = inS[rows]
        if m.all() or (~m).all(): continue
        psi += rank_modp(G[np.ix_(~m, m)], p)
    return psi

def ub3(B3, inS):
    return sum(min(int(inS[r].sum()), len(r) - int(inS[r].sum())) for r, G in B3.values())
