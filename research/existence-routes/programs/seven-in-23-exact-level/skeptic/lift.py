"""lift.py -- e3-bpexact: the w5-blueprint survivor (7 in 23, Z7, O O F O F) built with ACTUAL 16 x 16 Clifford matrices.

Every symplectic object of w5-blueprint/s2/t10_build.py is replaced by a matrix whose symplectic image (sp8.sigma) is
checked.  Conventions (GAME.md s5): a term [w] x_l is an element w of the matrix group Gamma applied to input l; a letter b
acts by [w] -> [b w]; tau_h([w] x_l) = [c_h w c'_h] x_(l+h), c'_h = K^-1 c_h^-1 K.  The cell of w is sigma(w)^-1 L_Z.

Construction (same configuration rule as t10_build.py):
  F    = first neighbour of the kernel cell at input distance 3;
  t1   = first neighbour of F at input distance 2 and distance 2 from the kernel;  lam1 = a neighbour of t1 at distance 1;
  O1   = [m1] x_0 (cell lam1),  O2 = [g2] x_0, g2 = m2 m1 (cell t1),  F1's move m (cell of m g2 = F), f0 = m g2;
  B_m  = letters b whose symplectic image fixes Lam_m = sigma(m)^-1 L_Z (so every term m b g2 sits in cell F);
  q    = an order-7 letter of B_m (exact order 7 as a matrix), p_j = m q^j m^-1 (letters), beta_h = q^h c_h^-1 (letters);
  F1   = m . sum_h psi(h)^-1 beta_h tau_h(sum_b y_b b g2 x_0)  =  sum_(h,b) psi(h)^-1 y_b [m q^h b g2 c'_h] x_h.
  Absorption: tau_g(F1) = psi(g) c_g p_g^-1 F1 (checked as formal sums).
"""
import sys, os, random, itertools, time
import numpy as np
W5 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'w5lib')   # unchanged copies of w5-blueprint/sp8.py, relax.py, vlib/
sys.path.insert(0, W5)
from sp8 import *          # G (= vlib Grp(4)), z7, sigma, Polar, LZ, table, img_mask, inv, mul, IDM, ...

I16 = np.eye(16, dtype=complex)
def key(w):
    a = np.round(w * 4096)
    return np.concatenate([a.real, a.imag]).astype(np.int64).tobytes()
def mm(*Ms):
    out = Ms[0]
    for M in Ms[1:]: out = snap(out @ M)
    return out
from grp import snap
def minv(M): return np.conj(M).T
def is_letter(M):
    return G.is_mono(M)

class Lift:
    def __init__(self, seed=21, psi=1, P=None, nletters=4000, verbose=True):
        t0 = time.time()
        self.rng = rng = random.Random(seed); self.psi = psi
        self.P = P or Polar()
        P = self.P
        self.K = G.K; self.Ki = minv(G.K)
        self.iZ = P.idx[LZ]
        self.iK = self.cell(self.K)
        # letters (matrices) with symplectic images
        self.letters = []
        for _ in range(nletters):
            b = G.rand_engine(rng); self.letters.append((b, sigma(b)))
        # configuration (t10_build's rule)
        D = P.D
        self.F = [j for j in P.nb[self.iK] if P.dZ[j] == 3][0]
        T = [Y for Y in P.nb[self.F] if P.dZ[Y] == 2 and D[Y, self.iK] == 2]
        self.t1 = T[0]
        L1c = [Y for Y in P.nb[self.t1] if P.dZ[Y] == 1]
        self.lam1 = L1c[0]
        self.m1 = self.find_move(I16, self.lam1)
        self.m2 = self.find_move(self.m1, self.t1)
        self.g2 = mm(self.m2, self.m1)
        self.m = self.find_move(self.g2, self.F)
        self.f0 = mm(self.m, self.g2)
        assert self.cell(self.m1) == self.lam1 and self.cell(self.g2) == self.t1 and self.cell(self.f0) == self.F
        # symmetry Z7 (exact order 7 as matrices)
        c = z7()[0][1]
        c7 = np.linalg.matrix_power(c, 7)
        lam = c7[0, 0]; assert np.allclose(c7, lam * I16)
        k = [kk for kk in range(4) if abs(1j ** kk - lam) < 1e-9][0]   # lam = i^k; (i^k c)^7 = i^(7k) lam = i^(8k) = 1
        c = snap((1j ** k) * c)
        assert np.allclose(np.linalg.matrix_power(c, 7), I16)
        self.c = [I16]
        for h in range(1, 7): self.c.append(mm(self.c[-1], c))
        self.cp = [mm(self.Ki, minv(ch), self.K) for ch in self.c]
        assert all(is_letter(x) for x in self.c + self.cp)
        # B_m: letters fixing Lam_m
        self.Lam_m = self.cell(self.m)          # = sigma(m)^-1 L_Z
        self.Bm = [(b, s) for (b, s) in self.letters if self.fixes(s, self.Lam_m)]
        # order-7 letter q in B_m
        self.q = None
        for (b, s) in self.Bm:
            o = self.sp_order(s)
            if o % 7 == 0:
                e = o // 7
                qs = b; qq = np.linalg.matrix_power(b, e) if e > 1 else b
                qq = snap(qq)
                # exact order 7 as a matrix: qq^7 is in the kernel (Paulis x Z4), a 2-group; raise to 2^a(2^-a mod 7)
                M7 = np.linalg.matrix_power(qq, 7); a2 = 0; Mp = M7
                while not np.allclose(Mp, I16):
                    Mp = snap(Mp @ Mp); a2 += 1
                    assert a2 < 10
                ee = (2 ** a2) * pow(2 ** a2, -1, 7)
                qf = snap(np.linalg.matrix_power(qq, ee))
                if np.allclose(np.linalg.matrix_power(qf, 7), I16) and not np.allclose(qf, I16):
                    self.q = qf; break
        assert self.q is not None, 'no order-7 letter in B_m'
        assert is_letter(self.q) and self.fixes(sigma(self.q), self.Lam_m)
        self.qp = [I16]
        for j in range(1, 7): self.qp.append(mm(self.qp[-1], self.q))
        mi = minv(self.m)
        self.p = [mm(self.m, qj, mi) for qj in self.qp]
        assert all(is_letter(x) for x in self.p), 'p_j not letters'
        self.beta = [mm(self.qp[h], minv(self.c[h])) for h in range(7)]
        assert all(is_letter(x) for x in self.beta)
        self.chars = [[np.exp(2j * np.pi * kk * h / 7) for h in range(7)] for kk in range(7)]
        if verbose:
            print(f'[lift] seed {seed}: cells lam1 {self.lam1} t1 {self.t1} F {self.F} K {self.iK}; |B_m sample| {len(self.Bm)} '
                  f'of {nletters}; q order 7 found; {time.time() - t0:.0f}s', flush=True)
    # ------------------------------------------------------------------
    def cell(self, w):
        return self.P.idx[img_mask(table(inv(sigma(w))), LZ)]
    def fixes(self, s, i):
        return self.P.idx[img_mask(table(s), self.P.L[i])] == i
    def sp_order(self, s):
        k = 1; h = s
        while h != IDM:
            h = mul(s, h); k += 1
            if k > 400: return -1
        return k
    def find_move(self, g, target, tries=20000):
        for _ in range(tries):
            z = self.rng.randrange(1, 16); b = self.rng.choice(self.letters)[0]
            mv = mm(G.D[z], b)
            if self.cell(mm(mv, g)) == target: return mv
        raise RuntimeError('no move found')
    # ------------------------------------------------------------------ formal sums: dict (l, key) -> [coef, matrix]
    @staticmethod
    def fs_add(u, l, w, c):
        k = (l, key(w))
        if k in u: u[k][0] += c
        else: u[k] = [c, w]
    def tau(self, h, u):
        out = {}
        for (l, kk), (c, w) in u.items():
            self.fs_add(out, (l + h) % 7, mm(self.c[h], w, self.cp[h]), c)
        return out
    def lmul(self, b, u):
        out = {}
        for (l, kk), (c, w) in u.items(): self.fs_add(out, l, mm(b, w), c)
        return out
    @staticmethod
    def clean(u, tol=1e-12):
        return {k: v for k, v in u.items() if abs(v[0]) > tol}
    @staticmethod
    def fs_diff(u, v):
        keys = set(u) | set(v)
        return max([abs((u[k][0] if k in u else 0) - (v[k][0] if k in v else 0)) for k in keys] + [0])
    # ------------------------------------------------------------------ the shared move F1
    def F1(self, ys, bs):
        """ys: coefficients, bs: letters of B_m (matrices).  Returns F1 as a formal sum."""
        u = {}
        ps = np.conj(self.chars[self.psi])
        for h in range(7):
            for y, b in zip(ys, bs):
                if abs(y) < 1e-14: continue
                self.fs_add(u, h, mm(self.m, self.qp[h], b, self.g2, self.cp[h]), ps[h] * y)
        return self.clean(u)
    def F1_from_def(self, ys, bs):
        """m . sum_h psi(h)^-1 beta_h tau_h(y), y = sum_b y_b [b g2] x_0  (independent rebuild)."""
        y = {}
        for yy, b in zip(ys, bs): self.fs_add(y, 0, mm(b, self.g2), yy)
        acc = {}
        ps = np.conj(self.chars[self.psi])
        for h in range(7):
            th = self.lmul(self.beta[h], self.tau(h, y))
            for k, (c, w) in th.items(): self.fs_add(acc, k[0], w, ps[h] * c)
        return self.clean(self.lmul(self.m, self.clean(acc)))
    def check_absorbed(self, F1):
        worst = 0
        for g in range(7):
            lhs = self.tau(g, F1)
            rhs = self.lmul(mm(self.c[g], minv(self.p[g])), F1)
            rhs = {k: [self.chars[self.psi][g] * v[0], v[1]] for k, v in rhs.items()}
            worst = max(worst, self.fs_diff(lhs, rhs))
        return worst

if __name__ == '__main__':
    L = Lift()
    bs = [b for (b, s) in L.Bm[:60]]
    ys = np.random.default_rng(1).standard_normal(len(bs)) + 0j
    F1 = L.F1(ys, bs); F1b = L.F1_from_def(ys, bs)
    print('F1 terms', len(F1), ' def-rebuild diff %.1e' % L.fs_diff(F1, F1b), ' absorption defect %.1e' % L.check_absorbed(F1))
    print('cells of array-0 terms:', sorted(set(L.cell(v[1]) for k, v in F1.items() if k[0] == 0)), ' F =', L.F)
    print('letters c, cp, p, beta, q all Z4-monomial: OK')
