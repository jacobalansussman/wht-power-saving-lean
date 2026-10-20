# v1_entry.py -- v-ms-mixed (skeptic of ms-mixed's Theorem T), 14 Oct 2026.  OWN CODE: nothing imported from
# ms-mixed / hooklib / mixlib.  Recomputes the key finite fact behind the corollary of Theorem T:
#   for every free monomial b (widened: S8 x Z4^8; engine a subset), the single term g = D b (D = dir(1))
#   FAILS T2 for the two-row target R = dir(1) dir(2):  pi_O rho(g) E_R != 0,  E_R = N cap rho(R)^-1 N,
#   in the hook M = V (x) Lambda^4 V (f = 5), N = span{e_x (x) e_T : x in T}.
# Classes: T2(D b) is invariant under b -> s b with D s D^-1 monomial; I enumerate the 1680 classes STRUCTURALLY
# (perfect matching of columns induced by b^-1 on the D-pairs {2k,2k+1}, times q^2 = +-1 per pair), not by sampling.
# Two arithmetics: complex floats (SVD) and exact mod-p at two primes not used by ms-mixed.
# Lambda^4 is built differently from hooklib: Lambda^4(D b) = Lambda^4(D) Lambda^4(b) with Lambda^4(b) combinatorial
# (monomial: permutation of 4-sets, sign of sorting, product of phases) and Lambda^4(D) by Gaussian-integer minors
# computed exactly with Python ints (D = A/2 with A in Z[i]).
import sys, itertools, random
import numpy as np

S4 = list(itertools.combinations(range(8), 4)); SI = {T: i for i, T in enumerate(S4)}
BAS = [(x, T) for x in range(8) for T in S4]
NI = [i for i, (x, T) in enumerate(BAS) if x in T]
OI = [i for i, (x, T) in enumerate(BAS) if x not in T]
assert len(NI) == 280 and len(OI) == 280

# ---------- exact Gaussian integers as (re, im) python ints -------------------------------------------------------
def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gadd(a, b): return (a[0]+b[0], a[1]+b[1])
def gneg(a): return (-a[0], -a[1])
def det4(m):  # exact 4x4 determinant over Z[i] by Laplace
    def det(mm):
        n = len(mm)
        if n == 1: return mm[0][0]
        tot = (0, 0)
        for j in range(n):
            if mm[0][j] == (0, 0): continue
            sub = [row[:j] + row[j+1:] for row in mm[1:]]
            t = gmul(mm[0][j], det(sub))
            tot = gadd(tot, t if j % 2 == 0 else gneg(t))
        return tot
    return det(m)

def dir2(z):   # 2*dir(z) = (1+i) I + (1-i) X^z, exact Z[i]
    A = [[(0, 0)] * 8 for _ in range(8)]
    for x in range(8):
        A[x][x] = gadd(A[x][x], (1, 1)); A[x ^ z][x] = gadd(A[x ^ z][x], (1, -1))
    return A
def gmatmul(A, B):
    return [[ (sum(A[r][k][0]*B[k][c][0] - A[r][k][1]*B[k][c][1] for k in range(8)),
               sum(A[r][k][0]*B[k][c][1] + A[r][k][1]*B[k][c][0] for k in range(8))) for c in range(8)] for r in range(8)]
def lam4_exact(A):   # 70x70 over Z[i]
    return [[det4([[A[r][c] for c in J] for r in I]) for J in S4] for I in S4]

# scaled exact matrices: D2 = 2 D, R4 = 4 R; Lambda^4 of scaled: 2^4 Lambda^4(D), 4^4 Lambda^4(R); rho scales by 2^5, 4^5
A_D = dir2(1); A_R = gmatmul(dir2(1), dir2(2)); A_K = gmatmul(A_R, dir2(4))
L_D = lam4_exact(A_D); L_R = lam4_exact(A_R); L_K = lam4_exact(A_K)

def to_c(A, scale): return np.array([[complex(a[0], a[1]) for a in row] for row in A]) / scale
def to_p(A, p, I):  return np.array([[(a[0] + a[1] * I) % p for a in row] for row in A], dtype=np.int64)

# monomial b: b e_x = i^ph[x] e_perm[x];  Lambda^4 b e_T = (prod phases) * sign * e_{perm T sorted}
def lam4_mono_entries(perm, ph):
    out = []   # list of (row, col, sign, phase_exp)
    for j, T in enumerate(S4):
        img = [perm[t] for t in T]
        srt = sorted(img)
        # sign of permutation sorting img
        inv = sum(1 for a in range(4) for b2 in range(a+1, 4) if img[a] > img[b2])
        out.append((SI[tuple(srt)], j, -1 if inv % 2 else 1, sum(ph[t] for t in T) % 4))
    return out

def rho_mono_sc(perm, ph):
    """rho(b) e_j = sgn_j * i^ex_j e_{sig_j}: returns (sig, sgn, ex) arrays of length 560"""
    sig = np.zeros(560, int); sgn = np.zeros(560, int); ex = np.zeros(560, int)
    L = lam4_mono_entries(perm, ph)
    for x in range(8):
        for (r, c, s, e) in L:
            j = x * 70 + c; sig[j] = perm[x] * 70 + r; sgn[j] = s; ex[j] = (e + ph[x]) % 4
    return sig, sgn, ex
def times_mono_c(G, perm, ph):      # G @ rho(b)
    sig, sgn, ex = rho_mono_sc(perm, ph)
    return G[:, sig] * (sgn * (1j ** ex))[None, :]
def times_mono_p(G, perm, ph, p, I):
    sig, sgn, ex = rho_mono_sc(perm, ph)
    ipow = np.array([1, I, p - 1, (p - I) % p], np.int64)
    return (G[:, sig] * ((sgn * ipow[ex]) % p)[None, :]) % p
# check the shortcut against a dense rho(b) once (float)
def rho_mono_dense(perm, ph):
    Rm = np.zeros((560, 560), complex)
    sig, sgn, ex = rho_mono_sc(perm, ph)
    Rm[sig, np.arange(560)] = sgn * (1j ** ex)
    return Rm

# ---------- floats -------------------------------------------------------------------------------------------------
def nullsp_c(A, tol=1e-9):
    u, s, vh = np.linalg.svd(A)
    r = int((s > tol * s[0]).sum())
    return vh[r:].conj().T, s
Dc = to_c(A_D, 2); Rc = to_c(A_R, 4); Kc = to_c(A_K, 8)
rhoD_c = np.kron(Dc, to_c(L_D, 16)); rhoR_c = np.kron(Rc, to_c(L_R, 256)); rhoK_c = np.kron(Kc, to_c(L_K, 4096))
# sanity: unitary
for nm, Rm in (('D', rhoD_c), ('R', rhoR_c), ('K', rhoK_c)):
    assert np.allclose(Rm @ Rm.conj().T, np.eye(560)), nm
ER_c, sR = nullsp_c(rhoR_c[np.ix_(OI, NI)]); EK_c, sK = nullsp_c(rhoK_c[np.ix_(OI, NI)])
ED_c, sD = nullsp_c(rhoD_c[np.ix_(OI, NI)])
def gap(s, k):  # smallest nonzero singular value
    nz = s[s > 1e-9 * s[0]]
    return nz.min()
print(f"floats: dim E_R = {ER_c.shape[1]} (gap {gap(sR,0):.3f}), dim E_K = {EK_c.shape[1]} (gap {gap(sK,0):.3f}), "
      f"dim N_D = {ED_c.shape[1]} -> cap = {280 - ED_c.shape[1]}, dem_R = {280 - ER_c.shape[1]}, dem_K = {280 - EK_c.shape[1]}")

# ---------- mod p --------------------------------------------------------------------------------------------------
def isprime(q):
    if q < 2: return False
    d = 2
    while d * d <= q:
        if q % d == 0: return False
        d += 1
    return True
def sqrtm1(p):
    for g in range(2, 1000):
        if pow(g, (p - 1) // 2, p) == p - 1: return pow(g, (p - 1) // 4, p)
def rank_p(A, p):
    A = A.copy() % p; r = 0; rows, cols = A.shape
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        pv = r + nz[0]
        if pv != r: A[[r, pv]] = A[[pv, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        col = A[:, c].copy(); col[r] = 0; nzr = np.nonzero(col)[0]
        if len(nzr): A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
        r += 1
    return r
def null_p(A, p):
    A = A.copy() % p; r = 0; rows, cols = A.shape; piv = []
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        pv = r + nz[0]
        if pv != r: A[[r, pv]] = A[[pv, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        col = A[:, c].copy(); col[r] = 0; nzr = np.nonzero(col)[0]
        if len(nzr): A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
        piv.append(c); r += 1
    free = [c for c in range(cols) if c not in set(piv)]
    Z = np.zeros((cols, len(free)), np.int64)
    for k, fc in enumerate(free):
        Z[fc, k] = 1
        for i, pc in enumerate(piv): Z[pc, k] = (-A[i, fc]) % p
    return Z
def mm_p(A, B, p):   # exact for p < 2^21 and inner dim <= 560
    return (A % p) @ (B % p) % p

primes = []
q = (1 << 21) - 1
while len(primes) < 2:
    if q % 4 == 1 and isprime(q): primes.append(q)
    q -= 2
PDATA = {}
for p in primes:
    I = sqrtm1(p); inv2 = pow(2, p - 2, p)
    Dp = to_p(A_D, p, I) * inv2 % p
    rhoD = np.kron(Dp, to_p(L_D, p, I) * pow(inv2, 4, p) % p) % p
    Rp = to_p(A_R, p, I) * pow(inv2, 2, p) % p
    rhoR = np.kron(Rp, to_p(L_R, p, I) * pow(inv2, 8, p) % p) % p
    Kp = to_p(A_K, p, I) * pow(inv2, 3, p) % p
    rhoK = np.kron(Kp, to_p(L_K, p, I) * pow(inv2, 12, p) % p) % p
    ER = null_p(rhoR[np.ix_(OI, NI)], p); EK = null_p(rhoK[np.ix_(OI, NI)], p); ND = null_p(rhoD[np.ix_(OI, NI)], p)
    # also Nperp-type: Y = O (complement coordinates)
    ERo = null_p(rhoR[np.ix_(NI, OI)], p); EKo = null_p(rhoK[np.ix_(NI, OI)], p)
    print(f"p = {p}: dim E_R = {ER.shape[1]}, dim E_K = {EK.shape[1]}, cap(N) = {280 - ND.shape[1]}; "
          f"Nperp: dim E_R = {ERo.shape[1]}, dim E_K = {EKo.shape[1]}", flush=True)
    PDATA[p] = dict(I=I, rhoD=rhoD, ER=ER, EK=EK, ERo=ERo, EKo=EKo)

# ---------- the 1680 classes, structurally ---------------------------------------------------------------------------
def matchings(s):
    if not s: yield []; return
    a = s[0]
    for j in range(1, len(s)):
        rest = s[1:j] + s[j+1:]
        for m in matchings(rest): yield [(a, s[j])] + m
MAT = list(matchings(list(range(8))))
assert len(MAT) == 105
reps = []
for m in MAT:
    for sg in itertools.product((0, 1), repeat=4):
        perm = [0] * 8; ph = [0] * 8
        for k, (x, y) in enumerate(m):
            perm[x] = 2 * k; perm[y] = 2 * k + 1; ph[x] = 0; ph[y] = sg[k]
        reps.append((perm, ph))
print("class representatives:", len(reps))

# cross-check the monomial shortcut against a dense float rho via 4x4 determinants, 3 random letters
def lam4_det(g):
    M = np.zeros((70, 70), complex)
    for i, A in enumerate(S4):
        for j, B in enumerate(S4): M[i, j] = np.linalg.det(g[np.ix_(A, B)])
    return M
_r = random.Random(1)
for _ in range(3):
    perm = list(range(8)); _r.shuffle(perm); ph = [_r.randrange(4) for _ in range(8)]
    b = np.zeros((8, 8), complex)
    for x in range(8): b[perm[x], x] = 1j ** ph[x]
    g = Dc @ b
    assert np.allclose(np.kron(g, lam4_det(g)), times_mono_c(rhoD_c, perm, ph)), "shortcut mismatch"
print("monomial shortcut agrees with dense determinant rho on 3 random letters")
import collections
tab = collections.Counter(); minres = 1e9; tabdef = collections.Counter()
for (perm, ph) in reps:
    G = times_mono_c(rhoD_c, perm, ph)
    A = G[np.ix_(OI, NI)]
    rR = np.linalg.norm(A @ ER_c); rK = np.linalg.norm(A @ EK_c)
    Ao = G[np.ix_(NI, OI)]
    row = [rR > 1e-8, rK > 1e-8]
    # mod p (both primes): defect ranks for R:N, K:N, R:Nperp, K:Nperp
    dfs = []
    for p, d in PDATA.items():
        Gp = times_mono_p(d['rhoD'], perm, ph, p, d['I'])
        Ap = Gp[np.ix_(OI, NI)]; Aop = Gp[np.ix_(NI, OI)]
        dfs.append((rank_p(mm_p(Ap, d['ER'], p), p), rank_p(mm_p(Ap, d['EK'], p), p),
                    rank_p(mm_p(Aop, d['ERo'], p), p), rank_p(mm_p(Aop, d['EKo'], p), p)))
    assert dfs[0] == dfs[1], (perm, ph, dfs)
    # float defect rank for R:N must agree
    sv = np.linalg.svd(A @ ER_c, compute_uv=False); frk = int((sv > 1e-8).sum())
    assert frk == dfs[0][0], (frk, dfs[0])
    assert (rR > 1e-8) == (dfs[0][0] > 0) and (rK > 1e-8) == (dfs[0][1] > 0)
    minres = min(minres, rR)
    tab[(dfs[0][0] == 0, dfs[0][1] == 0, dfs[0][2] == 0, dfs[0][3] == 0)] += 1
    tabdef[dfs[0]] += 1
print("(R:N pass, K:N pass, R:Nperp pass, K:Nperp pass) -> classes:", dict(tab))
print("smallest float residual ||pi_O rho(Db) E_R|| over classes:", round(minres, 4))
print("defect ranks (R:N, K:N, R:Nperp, K:Nperp) -> classes:")
for k, v in sorted(tabdef.items()): print("   ", k, v)
