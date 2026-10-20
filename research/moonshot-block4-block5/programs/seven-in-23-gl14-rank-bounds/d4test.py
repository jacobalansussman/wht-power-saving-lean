"""d4test.py -- e3-bpexact: the EXACT distance-4 condition of the lifted survivor, irrep by irrep of the letter group B.

Derivation (RUNLOG 05:15; NOTE section 3).  Only the top orbit's F1-part (and F2) put terms at input distance 4.  With
tau_g(F1) = psi(g) c_g p_g^-1 F1 and a Fourier transform over the arrays, the distance-4 part of output 0 is right iff,
for every character a of Z7 (F2 aside),
      r_a . D3 . A . P_a . kappa . f0 |_d4  =  K        for some r_a in C[B],  ONE A in C[B] (the top's F1 coefficient),
P_a = sum_j a(j) p_j, kappa = m (sum_b y_b b) m^-1, f0 = m g2.  The elements D3 . (letter) . (letter) . f0 live in the left
C[B]-module  C[B D3 B] f0 = sum_d C[B] g_d f0  (30 cosets g_d = D3 b_d; slot d = cell of g_d f0, a neighbour of F).  Right
multiplication by a letter k:  g_d k = beta(d,k) g_(k.d).  In an irrep sigma of B this is the block-monomial matrix
R_sigma(k) (30 d_sigma square), and the condition becomes, with d1 the coset of D3 itself:
      X_a . E_d1 R_sigma(A) . R_sigma(P_a) R_sigma(kappa) . Pi_S  =  sigma(b0) E_d0       (S = the 16 slots at distance 4)
i.e. some w in Wsp = {E_d1 R_sigma(A)} has  w M_a Pi_(S - d0) = 0  and  w M_a Pi_d0 invertible, for every a.
sigma = trivial is the distance-4 part of the frame-sum relaxation.  sigma = V (the natural 16-dim representation of the
letters) and its complex conjugate are faithful.  Wsp is computed from the letters (no assumption that it is full).
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lift import *

class D4:
    def __init__(self, L, D3=None, verbose=True):
        self.L = L; P = L.P
        self.D3 = G.D[1] if D3 is None else D3
        N3 = L.cell(self.D3)
        self.nbZ = list(P.nb[L.iZ])
        assert N3 in self.nbZ
        self.di = {c: i for i, c in enumerate(self.nbZ)}
        # coset representatives g_d = D3 b_d
        self.b = [None] * 30
        self.b[self.di[N3]] = I16
        pool = [b for (b, s) in L.letters]
        for b in pool:
            c = L.cell(mm(self.D3, b)); i = self.di[c]
            if self.b[i] is None: self.b[i] = b
            if all(x is not None for x in self.b): break
        assert all(x is not None for x in self.b)
        self.g = [mm(self.D3, bb) for bb in self.b]
        self.gi = [minv(x) for x in self.g]
        self.d1 = self.di[N3]
        self.slot = [L.cell(mm(gd, L.f0)) for gd in self.g]
        assert sorted(self.slot) == sorted(P.nb[L.F])
        self.d0 = [d for d in range(30) if self.slot[d] == L.iK][0]
        self.S = [d for d in range(30) if P.dZ[self.slot[d]] == 4]
        self.b0 = mm(L.K, minv(mm(self.g[self.d0], L.f0)))
        assert is_letter(self.b0)
        assert len(self.S) == 16 and self.d0 in self.S
        self.Sx = [d for d in self.S if d != self.d0]
        self._cache = {}
        if verbose:
            print(f'[d4] d1 {self.d1}, d0 {self.d0}, |S| {len(self.S)}', flush=True)
    def act(self, k):
        """right multiplication by the letter k: list of (k.d, beta(d,k))"""
        kk = key(k)
        if kk in self._cache: return self._cache[kk]
        out = []
        for d in range(30):
            w = mm(self.g[d], k)
            dd = self.di[self.L.cell(w)]
            beta = mm(w, self.gi[dd])
            assert is_letter(beta), 'beta not a letter'
            out.append((dd, beta))
        self._cache[kk] = out
        return out
    def R(self, k, sig, ds):
        M = np.zeros((30 * ds, 30 * ds), dtype=complex)
        for d, (dd, beta) in enumerate(self.act(k)):
            M[d * ds:(d + 1) * ds, dd * ds:(dd + 1) * ds] = sig(beta)
        return M
    def Wspace(self, sig, ds, letters, tol=1e-9):
        """basis of {E_d1 R_sigma(A)}: rows block d1 of R_sigma(b) over the given letters, as vectors (ds*30ds)"""
        vecs = []
        for b in letters:
            dd, beta = self.act(b)[self.d1]
            M = np.zeros((ds, 30 * ds), dtype=complex); M[:, dd * ds:(dd + 1) * ds] = sig(beta)
            vecs.append(M.reshape(-1))
        A = np.array(vecs)
        U, s, Vh = np.linalg.svd(A, full_matrices=False)
        r = int((s > tol * s[0]).sum())
        return Vh[:r]          # orthonormal rows
    def test(self, sig, ds, ys, ks, chars_req, letters_for_W, rng, tol=1e-8, full_rows=False):
        """ks: letters k_b = m b m^-1 with coefficients ys (kappa).  Returns per-character (rank of the d0 block, needed ds)."""
        L = self.L
        Rk = sum(y * self.R(k, sig, ds) for y, k in zip(ys, ks))
        Rp = [self.R(p, sig, ds) for p in L.p]
        Ma = {}
        for a in range(7):
            Pa = sum(L.chars[a][j] * Rp[j] for j in range(7))
            Ma[a] = Pa @ Rk
        colx = np.concatenate([np.arange(d * ds, (d + 1) * ds) for d in self.Sx])
        col0 = np.arange(self.d0 * ds, (self.d0 + 1) * ds)
        res = {}
        if full_rows:
            # Wsp = all (ds x 30ds) matrices: rows independent; alpha in C^(30 ds)
            C = np.concatenate([Ma[a][:, colx] for a in chars_req], axis=1)     # 30ds x (|req| 15 ds)
            U, s, Vh = np.linalg.svd(C, full_matrices=True)
            r = int((s > tol * max(1, s[0])).sum())
            Nul = U[:, r:].conj().T                                                # rows alpha with alpha C = 0 (left null space)
            res['dim_null'] = Nul.shape[0]
            res['junk'] = float(np.abs(Nul @ C).max()) if Nul.shape[0] else 0.0     # must be ~0
            for a in range(7):
                B0 = Nul @ Ma[a][:, col0] if Nul.shape[0] else np.zeros((0, ds))
                rk = np.linalg.matrix_rank(B0, tol=1e-7 * max(1, np.abs(B0).max() if B0.size else 1)) if B0.size else 0
                res[a] = rk
            return res
        Wb = self.Wspace(sig, ds, letters_for_W)          # N x (ds * 30ds)
        N = Wb.shape[0]; res['dim_W'] = N
        # conditions: w Ma Pi_x = 0 for a in req;  w = sum x_i W_i
        rows = []
        for a in chars_req:
            T = Ma[a][:, colx]
            rows.append(np.array([(Wb[i].reshape(ds, 30 * ds) @ T).reshape(-1) for i in range(N)]))
        C = np.concatenate(rows, axis=1)               # N x (...)
        U, s, Vh = np.linalg.svd(C, full_matrices=True)
        r = int((s > tol * max(1, s[0])).sum())
        Nul = U[:, r:].conj().T
        res['dim_null'] = Nul.shape[0]
        if Nul.shape[0] == 0:
            for a in range(7): res[a] = 0
            return res
        x = Nul.T @ (rng.standard_normal(Nul.shape[0]) + 1j * rng.standard_normal(Nul.shape[0]))
        w = (x @ Wb).reshape(ds, 30 * ds)
        for a in range(7):
            B0 = w @ Ma[a][:, col0]
            res[a] = np.linalg.matrix_rank(B0, tol=1e-7 * max(1, np.abs(B0).max()))
        return res

# ---------------------------------------------------------------- representations of the letters
def sig_triv(b): return np.ones((1, 1), dtype=complex)
def sig_nat(b): return b
def sig_conj(b): return np.conj(b)
def affine_of(b):
    """b = monomial: e_y -> phase e_(A y + c).  Returns (A as list of column images, c)."""
    img = [int(np.argmax(np.abs(b[:, y]))) for y in range(16)]
    c = img[0]
    A = [img[1 << j] ^ c for j in range(4)]
    return A, c
def sig_gl15(b):
    A, c = affine_of(b)
    M = np.zeros((15, 15), dtype=complex)
    for v in range(1, 16):
        w = 0
        for j in range(4):
            if (v >> j) & 1: w ^= A[j]
        M[w - 1, v - 1] = 1
    return M
def sig_aff16(b):       # affine permutation action on the 16 addresses (factors through B -> AGL(4,2))
    img = [int(np.argmax(np.abs(b[:, y]))) for y in range(16)]
    M = np.zeros((16, 16), dtype=complex)
    for y in range(16): M[img[y], y] = 1
    return M
def sig_absnat(b):      # |b|: the permutation part of the natural representation (same as aff16)
    return np.abs(b).astype(complex)
