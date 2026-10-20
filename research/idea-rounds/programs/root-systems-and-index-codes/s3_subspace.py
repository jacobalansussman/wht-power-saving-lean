"""s3_subspace.py -- the 7-array subspace problem (rmfam N14) is a 3-qubit Pauli problem living in E8.

Setting (rmfam natural form, block 5): GF(32) with x^5+x^2+1, alpha = x, a self-dual basis f_1..f_5; directions
z in GF(32)^* in self-dual coordinates; array k is owed T_k = multiplication by alpha^(2k) (symmetric in this basis);
the program must deliver, for every quadratic form Q (Q(0) = 0), the sign (-1)^<Q,T_k> on array k.

Checks:
 S1  V = {0} u {alpha^k : k in {0,1,2,5,11,18,19}} is a 3-dim F_2-subspace, and T_(v+v') = T_v + T_v'.
 S2  the seven characters chi_k(Q) = <Q,T_k> are exactly the seven nonzero characters of a 3-dim quotient
     pi: Sym_5^* -> F_2^3 (so the target is D_Q = diag((-1)^{pi(Q).v}), v in V \ 0); kernel K0 has dimension 12.
 S3  with an 8th register (target 0, i.e. returned unchanged) and the free 8x8 Hadamard change of register basis,
     every target D_Q becomes the permutation matrix of the translation by pi(Q) on F_2^3: a Pauli X^pi(Q) on
     3 qubits.
 S4  in that Fourier basis the plain gates (reflections in the 7 non-trivial Hadamard rows) are E8 roots, and
     every target X^v (v != 0) is the product of 4 reflections in the E8 roots e_x - e_(x+v) (a D4 longest element);
     E8 root lines in this basis = real 3-qubit stabiliser states of support 2 or 8 (56 + 64 = 120).
 S5  the symmetry group of the problem inside GL(5,2) (h with h U_V h^T = U_V, U_V = span of the targets).
 S6  dependencies among the 31 directions: minimal supports (3-dim subspaces and 3-flats), for the outside-gate lemma.
"""
import itertools
import numpy as np

M_ = 5; Q_ = 32; N_ = 31; POLY = 0b100101
exp = [0] * (2 * N_); log = [None] * Q_
x = 1
for k in range(N_):
    exp[k] = x; log[x] = k; x <<= 1
    if x & Q_: x ^= POLY
for k in range(N_, 2 * N_): exp[k] = exp[k - N_]
def mul(a, b):
    return 0 if a == 0 or b == 0 else exp[log[a] + log[b]]
def tr(a):
    s = 0; y = a
    for _ in range(M_): s ^= y; y = mul(y, y)
    return s
# self-dual basis (deterministic search)
def self_dual():
    for cand in itertools.combinations(range(1, Q_), M_):
        if all(tr(mul(a, b)) == (1 if a == b else 0) for a in cand for b in cand):
            return list(cand)
F = self_dual()
def coords(z): return tuple(tr(mul(z, f)) for f in F)
def mult_matrix(g):
    return np.array([[tr(mul(g, mul(fi, fj))) for fj in F] for fi in F], dtype=np.int64)

S = [0, 1, 2, 5, 11, 18, 19]
pts = [exp[k] for k in S]
Vset = set([0] + pts)
assert all((a ^ b) in Vset for a in Vset for b in Vset), "not closed"
assert len(Vset) == 8
print("S1: {0} u {alpha^k : k in", S, "} is a 3-dim subspace V of GF(32): OK")
T = {exp[k]: mult_matrix(exp[(2 * k) % N_]) for k in S}
T[0] = np.zeros((5, 5), dtype=np.int64)
for a in Vset:
    assert (T[a] == T[a].T).all()
for a in Vset:
    for b in Vset:
        assert (((T[a] + T[b]) % 2) == T[a ^ b]).all()
print("    targets are symmetric and additive: T_(v+v') = T_v + T_v' (mod 2): OK")
for a in pts:
    M = T[a].copy() % 2
    # F_2 rank
    r = 0; Mm = M.copy()
    for c in range(5):
        p = next((i for i in range(r, 5) if Mm[i, c]), None)
        if p is None: continue
        Mm[[r, p]] = Mm[[p, r]]
        for i in range(5):
            if i != r and Mm[i, c]: Mm[i] ^= Mm[r]
        r += 1
    assert r == 5 and np.diag(M).any()
print("    every target invertible and non-alternating (congruent to I): OK")

# quadratic forms: Q(x) = sum_{i<=j} q_ij x_i x_j ; pairing <Q,T> = sum_i q_ii T_ii + sum_{i<j} q_ij T_ij
pairsij = [(i, j) for i in range(5) for j in range(i, 5)]
def pair(qbits, Tm):
    return sum(q * Tm[i, j] for q, (i, j) in zip(qbits, pairsij)) % 2
def Qval(qbits, z):
    c = coords(z)
    return sum(q * c[i] * c[j] for q, (i, j) in zip(qbits, pairsij)) % 2
# characters as 15-bit functionals
chi = {a: np.array([T[a][i, j] % 2 for (i, j) in pairsij], dtype=np.int64) for a in Vset}
# check <Q, zz^T> = Q(z)
for z in range(1, Q_):
    c = np.array(coords(z)); zz = np.outer(c, c) % 2
    for qbits in itertools.islice(itertools.product((0, 1), repeat=15), 0, 2 ** 15, 997):
        assert pair(qbits, zz) == Qval(qbits, z)
# the 7 characters span a 3-dim space and are its nonzero elements
basisV = [pts[0]]
for a in pts:
    span = {0}
    for b in basisV:
        span |= {s ^ b for s in span}
    if a not in span: basisV.append(a)
assert len(basisV) == 3
# pi(Q) = (chi_b1(Q), chi_b2(Q), chi_b3(Q)); chi_v(Q) = pi(Q) . v-coordinates
def vcoord(a):
    for bits in itertools.product((0, 1), repeat=3):
        s = 0
        for bit, b in zip(bits, basisV):
            if bit: s ^= b
        if s == a: return bits
allQ = np.array(list(itertools.product((0, 1), repeat=15)), dtype=np.int64)
PI = np.stack([(allQ @ chi[b]) % 2 for b in basisV], axis=1)   # 32768 x 3
for a in pts:
    vc = np.array(vcoord(a))
    assert (((allQ @ chi[a]) % 2) == ((PI @ vc) % 2)).all()
K0 = int((PI.sum(axis=1) == 0).sum())
print("S2: the 7 array characters are the 7 nonzero characters of pi: Sym_5^* -> F_2^3; |K0| =", K0, "= 2^12: OK")
assert K0 == 4096

# S3: Fourier picture
Vpts = list(itertools.product((0, 1), repeat=3))    # registers indexed by F_2^3 (v-coordinates)
H = np.array([[(-1) ** (sum(a * b for a, b in zip(u, v)) % 2) for v in Vpts] for u in Vpts], dtype=float) / np.sqrt(8)
seen = {}
for q in range(0, 2 ** 15, 1):
    p = tuple(PI[q])
    if p in seen: continue
    D = np.diag([(-1.0) ** (sum(a * b for a, b in zip(p, v)) % 2) for v in Vpts])
    X = H @ D @ H.T
    # translation by p
    Tp = np.zeros((8, 8))
    for i, v in enumerate(Vpts):
        w = tuple((a + b) % 2 for a, b in zip(v, p)); Tp[Vpts.index(w), i] = 1
    assert np.allclose(X, Tp)
    seen[p] = True
assert len(seen) == 8
print("S3: H D_Q H^T = translation by pi(Q) on F_2^3 for all 8 classes (a Pauli X^pi(Q) on 3 qubits): OK")

# S4: E8 in the Fourier basis: roots e_x +- e_y and (1/2) eps with an even number of minus signs
E8 = []
for i, j in itertools.combinations(range(8), 2):
    for si in (1, -1):
        for sj in (1, -1):
            v = np.zeros(8); v[i] = si; v[j] = sj; E8.append(v)
for s in itertools.product((0.5, -0.5), repeat=8):
    if sum(1 for t in s if t < 0) % 2 == 0: E8.append(np.array(s))
E8 = np.array(E8); assert len(E8) == 240
def is_root(v):
    v = np.asarray(v, float)
    return any(np.allclose(v, r) for r in E8)
had = [np.sqrt(8) * H[i] for i in range(8)]   # +-1 Hadamard rows, indexed by u
for i in range(1, 8):
    assert is_root(had[i] / 2), "Hadamard row %d not a root" % i
print("S4: the 7 non-trivial Hadamard rows (the plain gates) are E8 roots (h_u/2): OK")
def refl(r):
    r = np.asarray(r, float); return np.eye(8) - 2 * np.outer(r, r) / (r @ r)
for p in Vpts[1:]:
    P = np.eye(8)
    for i, v in enumerate(Vpts):
        w = tuple((a + b) % 2 for a, b in zip(v, p)); j = Vpts.index(w)
        if i < j:
            r = np.zeros(8); r[i] = 1; r[j] = -1; assert is_root(r); P = refl(r) @ P
    Tp = np.zeros((8, 8))
    for i, v in enumerate(Vpts):
        w = tuple((a + b) % 2 for a, b in zip(v, p)); Tp[Vpts.index(w), i] = 1
    assert np.allclose(P, Tp)
    # and also = product of the 4 Hadamard reflections with <u,p> = 1
    P2 = np.eye(8)
    for i, u in enumerate(Vpts):
        if sum(a * b for a, b in zip(u, p)) % 2 == 1: P2 = refl(had[i]) @ P2
    assert np.allclose(P2, Tp)
print("    each target X^v = product of 4 reflections in E8 roots e_x - e_(x+v) = product of the 4 Hadamard")
print("    reflections with <u,v> = 1 (two 4A1 frames of the D4 through v): OK")
# real 3-qubit stabiliser states by orbit of |000> under H_i, CZ_ij, X_i, Z_i
def kron(*ms):
    out = np.array([[1.0]])
    for m in ms: out = np.kron(out, m)
    return out
I2 = np.eye(2); Hh = np.array([[1, 1], [1, -1]]) / np.sqrt(2); Xx = np.array([[0, 1], [1, 0]]); Zz = np.diag([1, -1])
gens = []
for q in range(3):
    for g in (Hh, Xx, Zz):
        ops = [I2] * 3; ops[q] = g; gens.append(kron(*ops))
for a, b in itertools.combinations(range(3), 2):
    d = np.ones(8)
    for idx in range(8):
        bits = [(idx >> (2 - t)) & 1 for t in range(3)]
        if bits[a] and bits[b]: d[idx] = -1
    gens.append(np.diag(d))
def key(v):
    v = np.round(v, 9); k = np.nonzero(np.abs(v) > 1e-9)[0][0]
    if v[k] < 0: v = -v
    return tuple(v)
start = np.zeros(8); start[0] = 1
orbit = {key(start): start}; frontier = [start]
while frontier:
    nxt = []
    for v in frontier:
        for g in gens:
            w = g @ v; kk = key(w)
            if kk not in orbit: orbit[kk] = w; nxt.append(w)
    frontier = nxt
stab = list(orbit.values())
supp = {}
for v in stab:
    s = int((np.abs(v) > 1e-9).sum()); supp[s] = supp.get(s, 0) + 1
e8lines = sum(1 for v in stab if is_root(v * np.sqrt(2)) or is_root(-v * np.sqrt(2)))
print("    real 3-qubit stabiliser rays:", len(stab), "by support size:", dict(sorted(supp.items())),
      "; of these, E8 root lines (scaled):", e8lines)
assert e8lines == 120 and supp.get(2, 0) + supp.get(8, 0) == 120

# S5: symmetry group inside GL(5,2): h with h U_V h^T = U_V
UV = [T[a] % 2 for a in pts]
def key_m(Mm): return tuple((Mm % 2).flatten())
UVkeys = set(key_m(Mm) for Mm in UV)
T1 = UV[0]
# enumerate h in GL(5,2) with h T1 h^T in U_V: solve column by column (h maps the form T1 to a target)
cols = [np.array(c) for c in itertools.product((0, 1), repeat=5) if any(c)]
count = 0; syms = []
def gl_search(prefix):
    global count
    k = len(prefix)
    if k == 5:
        h = np.stack(prefix, axis=1)
        img = (h @ T1 @ h.T) % 2
        if key_m(img) in UVkeys:
            if all(key_m((h @ Mm @ h.T) % 2) in UVkeys for Mm in UV):
                syms.append(h)
        return
    # linear independence check
    span = {tuple([0] * 5)}
    for c in prefix:
        span |= {tuple((a + b) % 2 for a, b in zip(s, c)) for s in span}
    for c in cols:
        if tuple(c) in span: continue
        gl_search(prefix + [c])
gl_search([])
print("S5: |{h in GL(5,2) : h U_V h^T = U_V}| =", len(syms))

# S6: dependencies sum_{z in R} z z^T = 0 among the 31 directions; minimal supports
dirs = list(range(1, Q_))
vec = {z: np.outer(np.array(coords(z)), np.array(coords(z))) % 2 for z in dirs}
flat = {z: np.array([vec[z][i, j] for (i, j) in pairsij]) for z in dirs}
A = np.stack([flat[z] for z in dirs], axis=1)  # 15 x 31
def f2rank(Mm):
    Mm = Mm.copy() % 2; r = 0
    for c in range(Mm.shape[1]):
        p = next((i for i in range(r, Mm.shape[0]) if Mm[i, c]), None)
        if p is None: continue
        Mm[[r, p]] = Mm[[p, r]]
        for i in range(Mm.shape[0]):
            if i != r and Mm[i, c]: Mm[i] ^= Mm[r]
        r += 1
    return r
print("S6: rank of the 31 rank-one matrices z z^T:", f2rank(A), "(so the dependency code has dimension", 31 - f2rank(A), ")")
# weight-7 dependencies = 3-dim subspaces minus 0
sub3 = set()
for a, b, c in itertools.combinations(dirs, 3):
    sp = {0, a, b, c, a ^ b, a ^ c, b ^ c, a ^ b ^ c}
    if len(sp) == 8: sub3.add(frozenset(sp - {0}))
ok = all(((sum(vec[z] for z in s)) % 2 == 0).all() for s in sub3)
print("    the", len(sub3), "three-dim subspaces (minus 0) are dependencies:", ok)
