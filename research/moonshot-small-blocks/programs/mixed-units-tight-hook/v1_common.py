import itertools, random
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


def matchings(s):
    if not s: yield []; return
    a = s[0]
    for j in range(1, len(s)):
        rest = s[1:j] + s[j+1:]
        for m in matchings(rest): yield [(a, s[j])] + m
MAT = list(matchings(list(range(8))))
assert len(MAT) == 105
