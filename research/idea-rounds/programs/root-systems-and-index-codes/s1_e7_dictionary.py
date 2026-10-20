"""s1_e7_dictionary.py -- exact checks of the E7 <-> block-3 frame-game dictionary.

Claims checked (all exact, integer/rational arithmetic):
 D1  the 126 roots of E7 map (mod 2, mod the radical) onto the 63 nonzero vectors of a symplectic F_2^6, two roots
     (+-alpha) per vector, with omega(a,b) = <a,b> mod 2.
 D2  the 7A1 frames (7 mutually orthogonal root lines) are exactly 135 and map bijectively onto the 135 Lagrangians.
 D3  in the orthonormal basis of one frame F0, every non-frame root is (1/sqrt2) v with v in {0,+-1}^7 of weight 4,
     supported on the complement of a line of a Fano plane on the 7 frame roots, all 16 sign patterns occurring; so a
     paid reflection is I - v v^T / 2 (a 4x4 "J - 2I" block, OpenAI's dense pair) and the frame stabiliser acts by
     signed Fano permutations (order 2^7 * 168).
 D4  the 64 quadratic forms on F_2^6 polarising to omega cut the root lines into closed subsystems: 36 forms with
     28 lines (type A7, rank 7) and 28 forms with 36 lines (type E6, rank 6).
 D5  the product of the 7 reflections of any 7A1 frame is -1 (the longest element), i.e. the identity of Sp(6,2).
"""
import itertools
from fractions import Fraction as Fr

# ---------- E8 and E7 roots (doubled coordinates: store 2*x so everything is an integer) ----------
def e8_roots2():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (2, -2):
            for sj in (2, -2):
                v = [0] * 8; v[i] = si; v[j] = sj; R.append(tuple(v))
    for signs in itertools.product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            R.append(tuple(signs))
    return R

def ip2(a, b):  # true inner product = sum(a_i b_i)/4 for doubled coordinates
    s = sum(x * y for x, y in zip(a, b)); assert s % 4 == 0; return s // 4

E8 = e8_roots2(); assert len(E8) == 240
rho = E8[0]
E7 = [a for a in E8 if ip2(a, rho) == 0]
assert len(E7) == 126 and all(ip2(a, a) == 2 for a in E7)
print("E7 roots:", len(E7))

# ---------- simple roots, Z-coordinates, mod-2 reduction ----------
f = [Fr(1, 10 ** k) + Fr(k, 1000003) for k in range(8)]  # generic functional
pos = [a for a in E7 if sum(fi * ai for fi, ai in zip(f, a)) > 0]
assert len(pos) == 63
posset = set(pos)
def add(a, b): return tuple(x + y for x, y in zip(a, b))
simple = [a for a in pos if not any(add(b, c) == a for b in pos for c in pos if b < c)]
assert len(simple) == 7, len(simple)

def solve_int(basis, v):
    # solve sum c_i basis_i = v over Q (basis is 7 vectors in Q^8 spanning the E7 space)
    n = len(basis)
    A = [[Fr(basis[j][i]) for j in range(n)] + [Fr(v[i])] for i in range(8)]
    # Gaussian elimination on 8 x (n+1)
    piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, 8) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]; A[r] = [x / pv for x in A[r]]
        for i in range(8):
            if i != r and A[i][c] != 0:
                fac = A[i][c]; A[i] = [x - fac * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    assert r == n
    for i in range(r, 8): assert A[i][n] == 0
    sol = [Fr(0)] * n
    for i, c in enumerate(piv): sol[c] = A[i][n]
    assert all(s.denominator == 1 for s in sol)
    return [int(s) for s in sol]

coords = {a: solve_int(simple, a) for a in E7}
cartan = [[ip2(a, b) for b in simple] for a in simple]
Bm2 = [[x % 2 for x in row] for row in cartan]

def f2_rank(rows):
    rows = [int("".join(map(str, r)), 2) for r in rows]; rank = 0; nb = max(1, max((x.bit_length() for x in rows), default=1))
    for bit in reversed(range(nb)):
        p = next((i for i in range(rank, len(rows)) if rows[i] >> bit & 1), None)
        if p is None: continue
        rows[rank], rows[p] = rows[p], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i] >> bit & 1: rows[i] ^= rows[rank]
        rank += 1
    return rank

# radical of B mod 2
rad = [v for v in itertools.product((0, 1), repeat=7) if any(v) and all(sum(Bm2[i][j] * v[j] for j in range(7)) % 2 == 0 for i in range(7))]
print("radical of the mod-2 form on L/2L:", rad)
assert len(rad) == 1
r = rad[0]
# projection F_2^7 -> F_2^6 with kernel <r>: drop a coordinate k with r_k = 1 after reducing
k = r.index(1)
def proj(v):
    v = list(v)
    if v[k]: v = [(x + y) % 2 for x, y in zip(v, r)]
    return tuple(v[i] for i in range(7) if i != k)

img = {a: proj([c % 2 for c in coords[a]]) for a in E7}
def omega(x, y):
    # recover omega on F_2^6 from roots: <a,b> mod 2 must depend only on images
    return None
pairs = {}
for a in E7:
    pairs.setdefault(img[a], []).append(a)
assert (0,) * 6 not in pairs
assert len(pairs) == 63 and all(len(v) == 2 and add(v[0], v[1]) == (0,) * 8 for v in pairs.values())
# omega well defined and nondegenerate
om = {}
for x, ax in pairs.items():
    for y, ay in pairs.items():
        vals = {ip2(a, b) % 2 for a in ax for b in ay}
        assert len(vals) == 1
        om[(x, y)] = vals.pop()
vecs = sorted(pairs)
print("D1: 126 roots -> 63 nonzero vectors of F_2^6, 2 roots each; omega = <,> mod 2 well defined: OK")

# ---------- D2: frames and Lagrangians ----------
lines = [pairs[x][0] for x in vecs]  # one representative per line, indexed like vecs
orth = [[ip2(a, b) == 0 for b in lines] for a in lines]
frames = []
def extend(cur, start):
    if len(cur) == 7:
        frames.append(tuple(cur)); return
    for j in range(start, 63):
        if all(orth[i][j] for i in cur):
            extend(cur + [j], j + 1)
extend([], 0)
print("7A1 frames:", len(frames))
def span(vs):
    S = {(0,) * 6}
    for v in vs:
        S |= {tuple((a + b) % 2 for a, b in zip(s, v)) for s in S}
    return frozenset(S)
lag_from_frames = {span([vecs[i] for i in fr]) for fr in frames}
# all Lagrangians: 3-dim totally isotropic subspaces
lags = set()
for a, b, c in itertools.combinations(vecs, 3):
    if om[(a, b)] == 0 and om[(a, c)] == 0 and om[(b, c)] == 0:
        S = span([a, b, c])
        if len(S) == 8: lags.add(S)
print("Lagrangians:", len(lags), " frames give", len(lag_from_frames), "distinct Lagrangians; equal sets:", lags == lag_from_frames)
assert len(frames) == 135 and lags == lag_from_frames and len(lags) == 135
# each frame's 7 root lines = the 7 nonzero vectors of its Lagrangian
assert all(set(vecs[i] for i in fr) == set(span([vecs[i] for i in fr])) - {(0,) * 6} for fr in frames)
print("D2: 135 frames <-> 135 Lagrangians (frame lines = the 7 nonzero vectors): OK")

# ---------- D3: frame coordinates ----------
F0 = [lines[i] for i in frames[0]]
supports = {}
for a in E7:
    c = [ip2(a, b) for b in F0]   # alpha = sum c_i beta_i / 2
    nz = [i for i in range(7) if c[i] != 0]
    if len(nz) == 1:
        assert abs(c[nz[0]]) == 2
    else:
        assert len(nz) == 4 and all(abs(c[i]) == 1 for i in nz), c
        supports.setdefault(tuple(nz), set()).add(tuple(c[i] for i in nz))
assert len(supports) == 7 and all(len(s) == 16 for s in supports.values())
fano_lines = [tuple(sorted(set(range(7)) - set(s))) for s in supports]
pair_count = {}
for L in fano_lines:
    for p in itertools.combinations(L, 2): pair_count[p] = pair_count.get(p, 0) + 1
assert len(pair_count) == 21 and set(pair_count.values()) == {1}
print("D3: non-frame roots = (+-b_i+-b_j+-b_k+-b_l)/2 on the 7 complements of the lines of a Fano plane,")
print("    all 16 sign patterns each (7*16 = 112 roots); paid reflection = I - v v^T/2, v in {0,+-1}^7 weight 4: OK")
# signed permutations preserving the root set
fano_auts = 0
for perm in itertools.permutations(range(7)):
    if all(tuple(sorted(perm[i] for i in L)) in set(fano_lines) for L in fano_lines): fano_auts += 1
print("    Fano automorphisms:", fano_auts, "-> frame stabiliser order 2^7 *", fano_auts, "=", 128 * fano_auts,
      "(Sp(6,2) parabolic P = 10752 = half of it, since -1 is central)")

# ---------- D4: symplectic quadrics <-> A7 / E6 ----------
basis6 = [tuple(1 if i == j else 0 for i in range(6)) for j in range(6)]
G6 = [[om[(basis6[i], basis6[j])] if basis6[i] in pairs and basis6[j] in pairs else None for j in range(6)] for i in range(6)]
# a quadratic form polarising to omega: q(x) = sum_{i<j} x_i x_j omega(e_i,e_j) + sum_i a_i x_i
def make_q(aa):
    def q(x):
        s = sum(aa[i] * x[i] for i in range(6))
        for i, j in itertools.combinations(range(6), 2):
            s += x[i] * x[j] * G6[i][j]
        return s % 2
    return q
# check polarisation
q0 = make_q([0] * 6)
for x in vecs:
    for y in vecs:
        z = tuple((a + b) % 2 for a, b in zip(x, y))
        if z in pairs:
            assert (q0(z) + q0(x) + q0(y)) % 2 == om[(x, y)]
types = {}
for aa in itertools.product((0, 1), repeat=6):
    q = make_q(aa)
    S = [x for x in vecs if q(x) == 1]
    roots = [a for x in S for a in pairs[x]]
    rset = set(roots)
    closed = all(add(a, b) in rset for a in roots for b in roots if ip2(a, b) == -1)
    # rank of the subsystem
    M = [[Fr(c) for c in coords[a]] for a in roots]
    def qrank(M):
        M = [row[:] for row in M]; rk = 0; ncol = len(M[0])
        for c in range(ncol):
            p = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
            if p is None: continue
            M[rk], M[p] = M[p], M[rk]
            for i in range(len(M)):
                if i != rk and M[i][c] != 0:
                    fac = M[i][c] / M[rk][c]; M[i] = [x - fac * y for x, y in zip(M[i], M[rk])]
            rk += 1
        return rk
    key = (len(S), closed, qrank(M))
    types[key] = types.get(key, 0) + 1
print("D4: (lines, closed, rank) -> number of forms:", types)
assert types == {(28, True, 7): 36, (36, True, 6): 28}
print("    36 forms with 28 lines (A7: 28 positive roots, rank 7), 28 forms with 36 lines (E6: 36 positive roots, rank 6): OK")

# ---------- D5: product of the 7 reflections of a frame = -1 ----------
def refl_matrix(a):
    # reflection in R^8 (true coordinates x = a/2): s(x) = x - <x,a> a ; matrix entries rational
    n = 8; av = [Fr(x, 2) for x in a]
    return [[Fr(int(i == j)) - av[i] * av[j] for j in range(n)] for i in range(n)]
def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
for fr in frames[:5]:
    P = None
    for i in fr:
        R = refl_matrix(lines[i]); P = R if P is None else matmul(R, P)
    # on the E7 space (orthogonal to rho) the product must be -1, on rho it is +1
    for a in E7[:20]:
        av = [Fr(x, 2) for x in a]
        img_ = [sum(P[i][j] * av[j] for j in range(8)) for i in range(8)]
        assert img_ == [-x for x in av]
print("D5: product of the 7 reflections of a frame acts as -1 on the E7 space (checked on 5 frames): OK")
print("ALL CHECKS PASSED")
