# v-ex2-certnum (skeptic): independent EXACT check of ex2-certnum's stand-in schemes and of the data behind Theorem 9.
# Own code, no sympy, no flint: Gaussian rationals as pairs of Fractions; the group B as monomial matrices (perm, phase
# exponents mod 4) by BFS; dim A bounded above by block-diagonality and below by a rank computed modulo p = 10^9+9
# (p = 1 mod 4, i -> sqrt(-1) mod p; a nonzero minor mod p is nonzero over Z[i], so rank_p <= rank over Q(i)).
# Usage: python3 vstandin.py path/to/scheme.json
import sys, json, itertools
from fractions import Fraction as Fr

n = 4
# ---------- Gaussian rationals ----------
Z = (Fr(0), Fr(0)); ONE = (Fr(1), Fr(0)); IU = (Fr(0), Fr(1))
def add(a, b): return (a[0] + b[0], a[1] + b[1])
def sub(a, b): return (a[0] - b[0], a[1] - b[1])
def mul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def inv(a):
    d = a[0] * a[0] + a[1] * a[1]; return (a[0] / d, -a[1] / d)
def isz(a): return a[0] == 0 and a[1] == 0
def mat(r, c): return [[Z] * c for _ in range(r)]
def mm(A, B):
    r, k, c = len(A), len(B), len(B[0]); C = mat(r, c)
    for i in range(r):
        Ai = A[i]; Ci = [Z] * c
        for l in range(k):
            a = Ai[l]
            if isz(a): continue
            Bl = B[l]
            for j in range(c):
                b = Bl[j]
                if not isz(b): Ci[j] = add(Ci[j], mul(a, b))
        C[i] = Ci
    return C
def madd(A, B): return [[add(x, y) for x, y in zip(ra, rb)] for ra, rb in zip(A, B)]
def kron(A, B):
    r1, c1, r2, c2 = len(A), len(A[0]), len(B), len(B[0])
    return [[mul(A[i // r2][j // c2], B[i % r2][j % c2]) for j in range(c1 * c2)] for i in range(r1 * r2)]
def rank(Mrows):
    M = [list(r) for r in Mrows]; rk = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if not isz(M[i][c])), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; iv = inv(M[rk][c]); M[rk] = [mul(iv, x) for x in M[rk]]
        for i in range(len(M)):
            if i != rk and not isz(M[i][c]):
                f = M[i][c]; M[i] = [sub(x, mul(f, y)) for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk
# ---------- the stand-in, typed from its definition (meta/ex-v-ex-dimension-logic/indep_core.py) ----------
def permmat(f):
    M = mat(n, n)
    for x in range(n): M[f(x)][x] = ONE
    return M
X0 = permmat(lambda x: x ^ 1); X1 = permmat(lambda x: x ^ 2)
I4 = [[ONE if i == j else Z for j in range(n)] for i in range(n)]
ca, cb = (Fr(1, 2), Fr(1, 2)), (Fr(1, 2), Fr(-1, 2))
D0 = [[add(mul(ca, I4[i][j]), mul(cb, X0[i][j])) for j in range(n)] for i in range(n)]
D1 = [[add(mul(ca, I4[i][j]), mul(cb, X1[i][j])) for j in range(n)] for i in range(n)]
R = kron(D0, D0); K = mm(D1, D0); Kr = kron(K, K); d = 16
# ---------- the group B as monomial matrices: g = (pi, e): g e_x = i^{e[x]} e_{pi[x]} ----------
gens = [((1, 0, 3, 2), (0,) * 4), ((2, 3, 0, 1), (0,) * 4),
        tuple([tuple(x ^ ((x & 1) << 1) for x in range(4)), (0,) * 4]),
        tuple([tuple(x ^ ((x >> 1) & 1) for x in range(4)), (0,) * 4]),
        ((0, 1, 2, 3), tuple(1 if x & 1 else 0 for x in range(4))),
        ((0, 1, 2, 3), tuple(1 if x & 2 else 0 for x in range(4))),
        ((0, 1, 2, 3), (1, 1, 1, 1))]
def gmul(h, g):  # (h g) e_x = h (i^{g.e[x]} e_{g.pi[x]}) = i^{g.e[x] + h.e[g.pi[x]]} e_{h.pi[g.pi[x]]}
    return (tuple(h[0][g[0][x]] for x in range(4)), tuple((g[1][x] + h[1][g[0][x]]) % 4 for x in range(4)))
idg = ((0, 1, 2, 3), (0, 0, 0, 0)); seen = {idg}; fr = [idg]
while fr:
    nf = []
    for g in fr:
        for h in gens:
            k = gmul(h, g)
            if k not in seen: seen.add(k); nf.append(k)
    fr = nf
print('|B| =', len(seen))
# rho(g) = g (x) g: rho(g) e_{a,b} = i^{e[a]+e[b]} e_{pi a, pi b}; index (a, b) -> 4a + b (same as kron)
def rho_entries(g):
    return {(4 * g[0][a] + g[0][b], 4 * a + b): (g[1][a] + g[1][b]) % 4 for a in range(4) for b in range(4)}
rhos = {}
for g in seen:
    e = rho_entries(g); rhos[tuple(sorted(e.items()))] = e
print('distinct rho(b):', len(rhos))
# adapted basis: U6 sym offdiag, U4 diag, L2 antisym. Membership in the block-diagonal algebra BD:
def piece(i, j):  # which piece/block of the coordinate change does a pair-index belong to
    return None
def blocks_of(M):
    # write M in the adapted basis: compute images of basis vectors and decompose. Basis vectors:
    B = []
    for a in range(4):
        for b in range(a + 1, 4):
            v = [Z] * 16; v[4 * a + b] = ONE; v[4 * b + a] = ONE; B.append(('U6', v))
    for a in range(4):
        v = [Z] * 16; v[5 * a] = ONE; B.append(('U4', v))
    for a in range(4):
        for b in range(a + 1, 4):
            v = [Z] * 16; v[4 * a + b] = ONE; v[4 * b + a] = (Fr(-1), Fr(0)); B.append(('L2', v))
    return B
BAS = blocks_of(None)
def coords(w):
    # coordinates of w in the adapted basis (exact): for a<b: w_ab = s + t, w_ba = s - t; diag: w_aa
    c = []
    for a in range(4):
        for b in range(a + 1, 4): c.append(('U6', mul((Fr(1, 2), Fr(0)), add(w[4 * a + b], w[4 * b + a]))))
    for a in range(4): c.append(('U4', w[5 * a]))
    for a in range(4):
        for b in range(a + 1, 4): c.append(('L2', mul((Fr(1, 2), Fr(0)), sub(w[4 * a + b], w[4 * b + a]))))
    return c
def adapted(M):  # 16x16 matrix in the adapted basis, labels
    cols = []
    for lab, v in BAS:
        w = [Z] * 16
        for i in range(16):
            s_ = Z
            for j in range(16):
                if not isz(v[j]) and not isz(M[i][j]): s_ = add(s_, mul(M[i][j], v[j]))
            w[i] = s_
        cols.append(coords(w))
    # Ma[r][c] = coordinate r of image of basis vector c
    return [[cols[c][r][1] for c in range(16)] for r in range(16)]
LAB = [lab for lab, _ in BAS]
def in_BD(M):
    Ma = adapted(M)
    return all(isz(Ma[r][c]) for r in range(16) for c in range(16) if LAB[r] != LAB[c])
# (1) every generator's rho is block diagonal => rho(B) and A inside BD (dim 36+16+36 = 88)
def rho_mat(g):
    M = mat(16, 16)
    for (r, c), e in rho_entries(g).items(): M[r][c] = [ONE, IU, (Fr(-1), Fr(0)), (Fr(0), Fr(-1))][e]
    return M
assert all(in_BD(rho_mat(g)) for g in gens); print('rho(generators) block diagonal in U6+U4+L2: yes (so A <= BD, dim 88)')
# (2) rank of span rho(B) modulo p >= 88 => dim A = 88 => A = BD
p = 1000000009; assert p % 4 == 1
ip = next(pow(a, (p - 1) // 4, p) for a in range(2, 100) if pow(pow(a, (p - 1) // 4, p), 2, p) == p - 1)
assert ip * ip % p == p - 1
pw = [1, ip, p - 1, p - ip]
rowsp = []
for e in rhos.values():
    v = [0] * 256
    for (r, c), k in e.items(): v[16 * r + c] = pw[k]
    rowsp.append(v)
def rank_mod(rows, p):
    M = [r[:] for r in rows]; rk = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(rk, len(M)) if M[i][c] % p), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; iv = pow(M[rk][c], p - 2, p); M[rk] = [x * iv % p for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % p:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk
rA = rank_mod(rowsp, p); print('rank of span rho(B) mod p =', rA, '(lower bound for dim_C A; need 88)'); assert rA == 88
# exact data of Theorem 9
Ra = adapted(R); Ka = adapted(Kr)
idx = {l: [r for r in range(16) if LAB[r] == l] for l in ('U6', 'U4', 'L2')}
def blk(Ma, rl, cl): return [[Ma[r][c] for c in idx[cl]] for r in idx[rl]]
print('R, K preserve Sym2 and L2:', all(isz(Ra[r][c]) and isz(Ka[r][c]) for r in range(16) for c in range(16)
      if (LAB[r] == 'L2') != (LAB[c] == 'L2')))
print('rank R[U6,U4] =', rank(blk(Ra, 'U6', 'U4')), ' rank R[U4,U6] =', rank(blk(Ra, 'U4', 'U6')))
print('need per array: rank K[U6,U4] =', rank(blk(Ka, 'U6', 'U4')), ' rank K[U4,U6] =', rank(blk(Ka, 'U4', 'U6')))
print('also rank K[U6,U6], K[U4,U4], R blocks:', rank(blk(Ka, 'U6', 'U6')), rank(blk(Ka, 'U4', 'U4')),
      rank(blk(Ra, 'U6', 'U6')), rank(blk(Ra, 'U4', 'U4')))
# ---------- the stored scheme ----------
if len(sys.argv) > 1:
    S = json.load(open(sys.argv[1])); W, s = S['W'], S['s']
    def rd(Mj): return [[(Fr(x[0]), Fr(x[1])) for x in row] for row in Mj]
    E = []
    for l in range(W):
        M = mat(16, 16 * W)
        for i in range(16): M[i][16 * l + i] = ONE
        E.append(M)
    srcs = list(E); nmix = 0; nz = 0
    assert len(S['moves']) == s and len(S['outputs']) == W
    for j in range(s):
        assert len(S['moves'][j]) == W + j
        c = mat(16, 16 * W)
        for i, Mj in enumerate(S['moves'][j]):
            M = rd(Mj); assert len(M) == 16 and all(len(r) == 16 for r in M)
            assert in_BD(M), 'move %d coefficient %d not in A' % (j, i); nmix += 1
            if any(not isz(x) for r in M for x in r): nz += 1; c = madd(c, mm(M, srcs[i]))
        srcs.append(mm(R, c))
    ok = True
    for t in range(W):
        assert len(S['outputs'][t]) == W + s
        o = mat(16, 16 * W)
        for i, Mj in enumerate(S['outputs'][t]):
            M = rd(Mj); assert in_BD(M), 'output %d coefficient %d not in A' % (t, i); nmix += 1
            if any(not isz(x) for r in M for x in r): nz += 1; o = madd(o, mm(M, srcs[i]))
        tgt = mat(16, 16 * W)
        for i in range(16):
            for j in range(16): tgt[i][16 * t + j] = Kr[i][j]
        bad = sum(1 for i in range(16) for j in range(16 * W) if o[i][j] != tgt[i][j])
        print('output %d: %d entries differ from rho(D1 D0) in block %d (exact)' % (t, bad, t)); ok &= bad == 0
    print('scheme %d in %d: %d mixing matrices, all in A (block diagonal, A = BD); %d nonzero; VERDICT %s' %
          (W, s, nmix, nz, 'PASS' if ok else 'FAIL'))
