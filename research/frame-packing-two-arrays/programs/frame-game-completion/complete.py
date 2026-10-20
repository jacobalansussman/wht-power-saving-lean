#!/usr/bin/env python3
"""complete.py (pk2-complete, packing round 2, 15 Oct 2026): the lift of a FRAME-GAME program to all four blocks.

Claim tested (PROOF.md Theorem 1): there is ONE algebra homomorphism

    Lambda : C[Sp(2m,2)] -> C[G],   Lambda(lambda) = [lambda^] . D(lambda^),   D(lambda^) in C[N]  (Pauli algebra),

with pi o Lambda = id (pi = the frame-game quotient), Lambda(C[P]) in C[B] and Lambda(T_z) in [dir z] C[B].  Here
D = sum_k D_k e_k over the four blocks k = f mod 4:
  k = 1 : D_1(g) = Pexp(g^-1)                    (fp-build's Lambda, Prop 1.3)
  k = 3 : D_3(g) = Pexp(conj(g)^-1) = Pexp(g^T)   (v-fp-build's Lambda_3), real Pauli basis X^a Z^b
  k = 0 : D_0(g) = sum_v eps(g, v)^-1 e0_v, g P_v g^-1 = eps(g, v) P_(gbar v)   (End(V) by conjugation; the 16:02
          Colab runs used the Hermitian basis H_v = i^|a&b| P_v instead, a conjugate lift: basis0='H')
  k = 2 : D_2(g) = sum_w del(g, w)^-1 e2_w, g P_w g^T = del(g, w) P_(g.w),  P_w = X^a Z^b  (Bell vectors vec(P_w))
where e{0,2}_w = 4^-m sum_u (-1)^omega(w,u) [X^a_u Z^b_u] e_k are the character idempotents of the (commutative) Pauli
algebra in the even blocks, and e_k = (1/4) sum_j i^(-jk) [i^j 1].  Everything in D is a FREE element.

Conventions follow fp-build/fpa.py: a label v = (a|b) is an int, a = bits 0..m-1 (x-part), b = bits m..2m-1;
a symplectic matrix is a tuple of 2m row bitmasks, (M v)_i = parity(M[i] & v); gbar is the action of g on labels by
conjugation g P g^-1.  Engine letters: dir z = (1+i)/2 1 + (1-i)/2 X^z (paid), monomial lifts of P (free).
"""
import sys, random, itertools
import numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__))     # fpa.py here is an unchanged copy of fp-build/fpa.py
if HERE not in sys.path: sys.path.insert(0, HERE)
from fpa import Sp, par, cadd, lmul, inp, end_solve, verify_combo, Fr

def pc(x): return bin(x).count("1")

class Eng:
    def __init__(s, m, basis0='P'):
        s.basis0 = basis0                 # 'P': block-0 characters realised on End(V) with basis X^a Z^b (default);
                                          # 'H': with the Hermitian basis i^|a&b| X^a Z^b (used by the 16:02 Colab runs)
        s.m = m; n = s.n = 1 << m; s.NL = 4 ** m; s.S = Sp(m)
        s.mask = n - 1
        s.Xm = [s._X(a) for a in range(n)]; s.Zm = [s._Z(b) for b in range(n)]
        s.P = [s.Xm[v & s.mask] @ s.Zm[v >> m] for v in range(s.NL)]          # scalar-free P_v = X^a Z^b
        s.hp = [pc((v & s.mask) & (v >> m)) % 4 for v in range(s.NL)]          # H_v = i^hp P_v Hermitian
        s.H = [(1j ** s.hp[v]) * s.P[v] for v in range(s.NL)]
        # symplectic form table omega(v, u) = a.b' + b.a'
        s.om = np.zeros((s.NL, s.NL), dtype=np.int8)
        for v in range(s.NL):
            for u in range(s.NL):
                s.om[v, u] = (pc((v & s.mask) & (u >> m)) + pc((v >> m) & (u & s.mask))) & 1
        s.chi = (1 - 2 * s.om.astype(np.int64)).astype(float)                  # chi_w(u) = (-1)^omega(w,u)
        s.K = np.eye(n, dtype=complex)
        for j in range(m): s.K = s.dirm(1 << j) @ s.K
    def _X(s, a):
        M = np.zeros((s.n, s.n), complex)
        for x in range(s.n): M[x ^ a, x] = 1
        return M
    def _Z(s, b): return np.diag([(-1.0) ** (pc(b & x) & 1) for x in range(s.n)]).astype(complex)
    def dirm(s, z): return (1 + 1j) / 2 * np.eye(s.n) + (1 - 1j) / 2 * s._X(z)
    # ------------------------------------------------------------------ identify c * X^a Z^b
    def ident(s, M):
        col = M[:, 0]; r = int(np.argmax(np.abs(col))); a = r          # X^a Z^b |0> = |a>
        Y = s.Xm[a] @ M                                                 # = c Z^b (diagonal)
        d = np.diag(Y); c = d[0]
        assert abs(abs(c) - 1) < 1e-9, "not a Pauli multiple"
        b = 0
        for i in range(s.m):
            if abs(d[1 << i] / c + 1) < 1e-9: b |= 1 << i
        v = a | (b << s.m)
        assert np.abs(M - c * s.P[v]).max() < 1e-9, "not a Pauli multiple"
        return c, v
    def symp(s, g):
        gi = np.linalg.inv(g)
        imgs = [s.ident(g @ s.P[1 << j] @ gi)[1] for j in range(2 * s.m)]   # column j = image of e_j
        return tuple(sum(((imgs[j] >> i) & 1) << j for j in range(2 * s.m)) for i in range(2 * s.m))
    # ------------------------------------------------------------------ the even-block cocycles
    def eps0(s, g):
        """block 0: g B_v g^-1 = eps(g, v) B_(gbar v) on End(V), B = P (default) or H; eps in mu4 (in +-1 for H)"""
        gi = np.linalg.inv(g); e = np.zeros(s.NL, complex); img = np.zeros(s.NL, dtype=np.int64)
        B = s.P if s.basis0 == 'P' else s.H
        for v in range(s.NL):
            c, w = s.ident(g @ B[v] @ gi)
            if s.basis0 == "H": c = c / (1j ** s.hp[w])
            assert abs(abs(c) - 1) < 1e-9 and min(abs(c - 1j ** k) for k in range(4)) < 1e-9
            e[v] = c; img[v] = w
        return e, img
    def del2(s, g):
        d = np.zeros(s.NL, complex); img = np.zeros(s.NL, dtype=np.int64)
        for w in range(s.NL):
            c, w2 = s.ident(g @ s.P[w] @ g.T)
            d[w] = c; img[w] = w2
        return d, img
    def pexp(s, M):
        """coefficients d_u with M = sum_u d_u P_u"""
        return np.array([np.trace(s.P[u].conj().T @ M) / s.n for u in range(s.NL)])
    def Dk(s, g):
        """D_k(g) for k = 0..3 as coefficient vectors over the scalar-free basis [P_u] (to be multiplied by e_k)"""
        e, _ = s.eps0(g); d, _ = s.del2(g)
        D0 = (s.chi @ (1 / e)) / s.NL                              # sum_v eps(v)^-1 e0_v ; e0_v = 4^-m sum_u chi_v(u) [P_u]
        D2 = (s.chi @ (1 / d)) / s.NL
        D1 = s.pexp(np.linalg.inv(g))
        D3 = s.pexp(g.T)
        return [D0, D1, D2, D3]
    def combine(s, Ds):
        """sum_k D_k e_k as a dict {(u, j): coef} over group elements i^j P_u"""
        out = {}
        for k, D in enumerate(Ds):
            for u in range(s.NL):
                if abs(D[u]) < 1e-13: continue
                for j in range(4):
                    out[(u, j)] = out.get((u, j), 0) + D[u] * (1j ** (-j * k)) / 4
        return {key: c for key, c in out.items() if abs(c) > 1e-13}
    def D(s, g): return s.combine(s.Dk(g))
    def Qout(s):
        """free left correction Q with Q Lambda(Kbar) = [K]: Q_k = [K] D_k(K)^-1 [K]^-1"""
        K = s.K; Ki = np.linalg.inv(K)
        e, _ = s.eps0(K); d, _ = s.del2(K)
        Dinv = [(s.chi @ e) / s.NL, s.pexp(K), (s.chi @ d) / s.NL, s.pexp(K.conj())]     # inverses of D_k(K)
        out = {}
        for k, Dv in enumerate(Dinv):
            for u in range(s.NL):
                if abs(Dv[u]) < 1e-13: continue
                c, w = s.ident(K @ s.P[u] @ Ki)            # K P_u K^-1 = c P_w, c in mu4
                jc = int(round(np.angle(c) / (np.pi / 2))) % 4
                assert abs(c - 1j ** jc) < 1e-9
                for j in range(4):
                    key = (w, (j + jc) % 4)
                    out[key] = out.get(key, 0) + Dv[u] * (1j ** (-j * k)) / 4
        return {key: c for key, c in out.items() if abs(c) > 1e-13}

    # ------------------------------------------------------------------ monomial lift of p in P (as t_lift.py)
    def lift_P(s, p):
        m, n = s.m, s.n
        A = [p[i] & s.mask for i in range(m)]
        C = [p[m + i] & s.mask for i in range(m)]
        Ainv = s.S._inv_m(A)
        Sm = []
        for i in range(m):
            acc = 0
            for k in range(m):
                if (C[i] >> k) & 1: acc ^= Ainv[k]
            Sm.append(acc)
        Mph = s._phase(Sm); Mad = np.zeros((n, n), complex)
        for x in range(n):
            y = sum(par(A[i] & x) << i for i in range(m)); Mad[y, x] = 1
        return Mph @ Mad
    def _phase(s, Srows):
        d = []
        for x in range(s.n):
            q = 0
            for i in range(s.m):
                if (x >> i) & 1:
                    if (Srows[i] >> i) & 1: q += 1
                    for j in range(i + 1, s.m):
                        if (x >> j) & 1 and (Srows[i] >> j) & 1: q += 2
            d.append(1j ** (q % 4))
        return np.diag(d)

# ====================================================================== raw-array evaluation at f columns
class Raw:
    def __init__(s, E, f):
        s.E = E; s.f = f; n = E.n; s.N = n ** f
        idx = np.arange(s.N, dtype=np.int64)
        s.idx = idx
        s.Arep = {}; s.Brep = {}
        s._cache = {}
    def rep(s, a):
        """a repeated in every base-n digit"""
        return sum(a << (s.E.m * k) for k in range(s.f))
    def pauli_idx(s, u):
        """(src, sign) with (X^a Z^b)^(x)f v = sign * v[src]; P|x> = (-1)^{b.x}|x^a>, so out[J] = (-1)^{b.(J^a)} v[J^a]"""
        if u in s._cache: return s._cache[u]
        a = u & s.E.mask; b = u >> s.E.m
        A = s.rep(a); B = s.rep(b)
        src = s.idx ^ A; t = src & B
        if hasattr(np, "bitwise_count"):
            par_ = (np.bitwise_count(t) & 1).astype(np.int8)
        else:
            par_ = np.zeros(s.N, dtype=np.int8)
            while t.any():
                par_ ^= (t & 1).astype(np.int8); t = t >> 1
        out = (src, (1 - 2 * par_).astype(np.int8))
        if s.N <= (1 << 20): s._cache[u] = out
        return out
    def pauli(s, u, v):
        """(X^a Z^b)^(x)f v, fast"""
        src, sg = s.pauli_idx(u)
        return sg * v[src]
    def mat(s, M, v):
        t = v.reshape((s.E.n,) * s.f)
        for ax in range(s.f):
            t = np.tensordot(M, t, axes=([1], [ax])); t = np.moveaxis(t, 0, ax)
        return t.reshape(-1)
    def nfree(s, Dd, v):
        """apply a C[N] element {(u, j): c}: [i^j P_u] acts as i^(j f) P_u^(x)f"""
        out = np.zeros_like(v); byu = {}
        for (u, j), c in Dd.items(): byu[u] = byu.get(u, 0) + c * (1j ** (j * s.f % 4))
        for u, c in byu.items():
            if abs(c) > 1e-14: out += c * s.pauli(u, v)
        return out
