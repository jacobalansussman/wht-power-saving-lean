# Independent skeptic rebuild of Theorem M (c2-finite-group): 105a of Sp(6,2) = W(E7)^+ inside 7 (x) Lambda^2 7 of the
# E7 reflection representation (integer matrices in simple-root coordinates). Shares no code with exact_cert.py.
import numpy as np, itertools, random, sys
from fractions import Fraction
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
n = 7
C = 2 * np.eye(n, dtype=np.int64)
for a, b in [(0, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1, 3)]:   # Bourbaki E7: 1-3-4-5-6-7, 2-4 (0-based)
    C[a, b] = C[b, a] = -1
I7 = np.eye(n, dtype=np.int64)
def refl(b):  # s_beta(v) = v - (beta, v) beta, (beta, v) = b^T C v
    b = np.asarray(b, dtype=np.int64); return I7 - np.outer(b, C @ b)
S = [refl(I7[i]) for i in range(n)]
# positive roots
roots = {tuple(I7[i]) for i in range(n)}; fr = list(roots)
while fr:
    nx = []
    for r in fr:
        for s in S:
            y = tuple(int(x) for x in s @ np.array(r))
            if y not in roots and tuple(-x for x in y) not in roots: roots.add(y); nx.append(y)
    fr = nx
pos = sorted(r if any(x > 0 for x in r) and all(x >= 0 for x in r) else tuple(-x for x in r) for r in roots)
pos = sorted(set(pos)); assert len(pos) == 63, len(pos)
PI = {r: i for i, r in enumerate(pos)}
def pair_index(v):
    t = tuple(int(x) for x in v)
    return PI[t] if t in PI else PI[tuple(-x for x in t)]
ip = lambda a, b: int(np.array(a) @ C @ np.array(b))
u = [tuple(int(x) % 2 for x in (C @ np.array(r))) for r in pos]   # mod-2 functional = image in the symplectic F_2^6
assert len(set(u)) == 63
# frame R0: 7 mutually orthogonal positive roots
def find_frame():
    def bt(chosen, start):
        if len(chosen) == 7: return chosen
        for j in range(start, 63):
            if all(ip(pos[j], pos[c]) == 0 for c in chosen):
                r = bt(chosen + [j], j + 1)
                if r: return r
        return None
    return bt([], 0)
R0 = find_frame(); R0set = frozenset(R0)
print("frame R0 (positive-root indices):", R0)
def xor(a, b): return tuple(x ^ y for x, y in zip(a, b))
lines = [t for t in itertools.combinations(range(7), 3) if xor(xor(u[R0[t[0]]], u[R0[t[1]]]), u[R0[t[2]]]) == (0,) * 7]
print("lines among the 7 frame roots:", len(lines))
assert len(lines) == 7
# Sp(6,2) element <-> det-1 element of W(E7)
def plus(g):
    d = round(np.linalg.det(g.astype(float))); assert d in (1, -1); return g if d == 1 else -g
def frame_img(g): return frozenset(pair_index(g @ np.array(pos[i])) for i in R0)
def dist(g): return {7: 0, 3: 1, 1: 2, 0: 3}[len(frame_img(g) & R0set)]
# representation 7 (x) Lambda^2 7
pairs = list(itertools.combinations(range(n), 2))
def lam2(g):
    M = np.zeros((21, 21), dtype=np.int64)
    for a, (b, c) in enumerate(pairs):
        for e, (d, f) in enumerate(pairs):
            M[a, e] = g[b, d] * g[c, f] - g[b, f] * g[c, d]
    return M
def rho(g): return np.kron(g, lam2(g))
D = 147
# class sums: 2a = 63 transvections t_u <-> -s_beta ; 2b = t_u t_v t_{u+v} over 315 totally isotropic lines
Z2a = np.zeros((D, D), dtype=np.int64)
for r in pos: Z2a += rho(-refl(r))
U = {u[i]: i for i in range(63)}
tl = set()
for i in range(63):
    for j in range(i + 1, 63):
        if ip(pos[i], pos[j]) % 2 == 0:
            k = U[xor(u[i], u[j])]; tl.add(frozenset((i, j, k)))
assert len(tl) == 315, len(tl)
Z2b = np.zeros((D, D), dtype=np.int64)
for t in tl:
    i, j, k = sorted(t)
    assert ip(pos[i], pos[j]) == 0 and ip(pos[i], pos[k]) == 0 and ip(pos[j], pos[k]) == 0
    Z2b += rho(-(refl(pos[i]) @ refl(pos[j]) @ refl(pos[k])))
assert (Z2a @ Z2b == Z2b @ Z2a).all()
Id = np.eye(D, dtype=np.int64)
assert ((Z2a + 21 * Id) @ (Z2a + 45 * Id) @ (Z2a + 9 * Id) == 0).all(), "Z2a not killed by (x+21)(x+45)(x+9)"
P = (Z2a + 45 * Id) @ (Z2a + 9 * Id)            # integer; columns span the -21 eigenspace
print("max |P| entry", int(np.abs(P).max()))
assert ((Z2a + 21 * Id) @ P == 0).all() and ((Z2b - 3 * Id) @ P == 0).all(), "P columns not in the (-21, 3) eigenspace"
p = (1 << 61) - 1
def rank_mod(M, want=False):
    A = [[int(x) % p for x in row] for row in M]; r = 0; rows = len(A); cols = len(A[0]) if rows else 0; piv = []
    for c in range(cols):
        pr = next((i for i in range(r, rows) if A[i][c]), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]; iv = pow(A[r][c], p - 2, p); A[r] = [x * iv % p for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c]:
                f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == rows: break
    return (r, piv) if want else r
r1 = rank_mod(P.T.tolist())
r2 = rank_mod(((Z2a + 21 * Id) @ (Z2a + 9 * Id)).T.tolist()); r3 = rank_mod(((Z2a + 21 * Id) @ (Z2a + 45 * Id)).T.tolist())
print("eigenspace ranks mod p (-21, -45, -9):", r1, r2, r3, " sum", r1 + r2 + r3)
assert r1 + r2 + r3 == D and r1 == 105
# integer basis of X: 105 independent columns of P
_, piv = rank_mod(P.tolist(), want=True)   # pivots = independent columns
Xb = P[:, piv]; assert Xb.shape[1] == 105
# A = <-s_alpha, alpha in R0>; 6 independent generators
Agen = [rho(-refl(pos[R0[i]])) for i in range(7)]
prod = Id.copy()
for a in Agen: prod = prod @ a
assert (prod == Id).all(), "product of the 7 A-generators is not 1"
def weight_proj(eps, Y):   # 64 * projector onto psi with psi(a_i) = eps_i (i < 6), applied to Y
    Y = Y.astype(object)
    for i in range(6):
        Y = Y + eps[i] * (Agen[i].astype(object) @ Y)
    return Y
chars = {}
for bits in range(64):
    eps = [(-1) ** ((bits >> i) & 1) for i in range(6)]
    e7 = 1
    for e in eps: e7 *= e
    Sset = frozenset(i for i in range(7) if (eps + [e7])[i] == -1)
    chars[Sset] = eps
assert len(chars) == 64 and all(len(s) % 2 == 0 for s in chars)
def has_line(s): return any(set(l) <= s for l in lines)
O = {'0': [], '2(21)': [], '4noline(7)': [], '4line(28)': [], '6(7)': []}
for s in chars:
    if len(s) == 0: O['0'].append(s)
    elif len(s) == 2: O['2(21)'].append(s)
    elif len(s) == 6: O['6(7)'].append(s)
    elif has_line(s): O['4line(28)'].append(s)
    else: O['4noline(7)'].append(s)
print({k: len(v) for k, v in O.items()})
wd = {}; WY = {}
for s, eps in chars.items():
    Y = weight_proj(eps, Xb); WY[s] = Y; wd[s] = rank_mod(Y.T.tolist())
print("weight-space dims per orbit (mod p):", {k: sum(wd[s] for s in v) for k, v in O.items()}, " total", sum(wd.values()))
assert sum(wd.values()) == 105
Ncols = []
for s in O['4line(28)']:
    Y = WY[s]; rk, pv = rank_mod(Y.tolist(), want=True)
    Ncols += [Y[:, c] for c in pv]
Nm = np.array(Ncols, dtype=object).T
dN, pvN = rank_mod(Nm.tolist(), want=True); Nm = Nm[:, pvN]
print("dim N =", dN, " (", Nm.shape, ")")
def bareiss_rank(cols):   # exact rank over Q of the integer matrix with these columns
    A = [list(map(int, r)) for r in np.array(cols, dtype=object).T.tolist()]   # rows = vectors
    rows = len(A); m = len(A[0]); r = 0; prev = 1
    for c in range(m):
        pr = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]
        for i in range(r + 1, rows):
            num = [A[r][c] * A[i][k] - A[i][c] * A[r][k] for k in range(m)]
            assert all(x % prev == 0 for x in num)
            A[i] = [x // prev for x in num]
        prev = A[r][c]; r += 1
        if r == rows: break
    return r
dN_exact = bareiss_rank(Nm.T)
print("dim N exact:", dN_exact)
# P_0-invariance: random elements of W^+ stabilising R0
def rand_w(L=40):
    g = I7.copy()
    for _ in range(L): g = S[random.randrange(7)] @ g
    return plus(g)
stab = []; tries = 0
while len(stab) < 6 and tries < 20000:
    tries += 1; g = rand_w()
    if frame_img(g) == R0set: stab.append(g)
print("found", len(stab), "stabiliser elements in", tries, "tries")
for g in stab:
    gN = rho(g).astype(object) @ Nm
    assert rank_mod(np.concatenate([Nm, gN], axis=1).T.tolist()) == dN, "N not kept by a P_0 element"
print("N kept by all sampled P_0 elements (mod p)")
# profile at distances 1, 2, 3, several elements each
found = {1: [], 2: [], 3: []}; tries = 0
while min(len(v) for v in found.values()) < 3 and tries < 5000:
    tries += 1; g = rand_w(); d = dist(g)
    if d in found and len(found[d]) < 3: found[d].append(g)
# plus the "named" elements: one transvection -s_beta with beta not orthogonal mod 2 to the frame
res = {}
for d in (1, 2, 3):
    vals = []
    for g in found[d]:
        gN = rho(g).astype(object) @ Nm
        M = np.concatenate([Nm, gN], axis=1)
        rp = rank_mod(M.T.tolist())
        vals.append(rp - dN)
    res[d] = vals
    print(f"distance {d}: f_N mod p over {len(vals)} random elements: {vals}", flush=True)
# exact ranks for one element per distance
ex = []
for d in (1, 2, 3):
    g = found[d][0]; gN = rho(g).astype(object) @ Nm
    rx = bareiss_rank(np.concatenate([Nm, gN], axis=1).T)
    ex.append(rx - dN_exact); print(f"distance {d}: exact f_N = {rx - dN_exact}", flush=True)
print("EXACT PROFILE", tuple(ex), " psi3+4psi2-10psi1 =", ex[2] + 4 * ex[1] - 10 * ex[0], " 7psi3+38psi2-88psi1 =", 7 * ex[2] + 38 * ex[1] - 88 * ex[0])
