"""t_lift.py (fp-build): the lift of a FRAME-GAME program to the engine game, tested on raw complex arrays.

Lift (PROOF.md, Proposition 2): Phi: C[G]_mu4 -> M_n(C) (x) C[Sp], g -> g (x) gbar, is an isomorphism; a frame-game
program lifts by sending each free p in P to Lambda(p) = Phi^-1(1 (x) p) = Pexp(phat^-1) . phat (phat a monomial lift,
Pexp = the Pauli expansion, a FREE combination) and each paid T_z to Lambda(T_z) = c_z . dir z with
c_z = ((1-i)[1] + (1+i)[X^z])/2 (e3-rmbuild Lemma R's corrected move).  The output 1 (x) Kbar is turned into K by the
free element Pexp(K).  Claim: the lifted word is right at every f = 1 mod 4 (the block where i acts as i); the other
blocks (f even, f = 3 mod 4) are NOT guaranteed.  Here: the order-3 control at block 3, its exact frame-game readout
from fpa.end_solve, lifted, replayed on random complex arrays at f = 1..5, and the same for the plain 2 in 6 with
random free letters inserted (a frame-game program that is not written with engine letters in mind)."""
import sys, random, itertools, warnings
warnings.filterwarnings("ignore")
import numpy as np
from fpa import *
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "ledger"))
import ledger as L

m = 3; n = 1 << m
S = Sp(m)
E = L.Engine(m)          # for symplectic images of matrices (convention check)

def X(a):
    M = np.zeros((n, n), complex)
    for x in range(n): M[x ^ a, x] = 1
    return M
def Z(b): return np.diag([(-1) ** (bin(b & x).count("1") & 1) for x in range(n)]).astype(complex)
PAULIS = [(a, b, X(a) @ Z(b)) for a in range(n) for b in range(n)]
def pexp(M):
    """Pauli expansion: list of (coef, Pauli matrix)"""
    out = []
    for a, b, P in PAULIS:
        c = np.trace(P.conj().T @ M) / n
        if abs(c) > 1e-12: out.append((c, P))
    return out
def dirm(z): return (1 + 1j) / 2 * np.eye(n) + (1 - 1j) / 2 * X(z)
def addr_mat(Arows):
    """permutation x -> A x, A given by row bitmasks (row i of A)"""
    M = np.zeros((n, n), complex)
    for x in range(n):
        y = sum(par(Arows[i] & x) << i for i in range(m))
        M[y, x] = 1
    return M
def phase_mat(Srows):
    d = []
    for x in range(n):
        q = 0
        for i in range(m):
            if (x >> i) & 1:
                if (Srows[i] >> i) & 1: q += 1
                for j in range(i + 1, m):
                    if (x >> j) & 1 and (Srows[i] >> j) & 1: q += 2
        d.append(1j ** (q % 4))
    return np.diag(d)
def lift_P(p):
    """monomial lift phat of p in P (p = ph(S) addr(A), S = C A^-1)"""
    A = [p[i] & (n - 1) for i in range(m)]                  # top-left block rows
    C = [p[m + i] & (n - 1) for i in range(m)]              # bottom-left block rows
    Ainv = S._inv_m(A)
    # S_ = C A^-1 (row i of product = xor of rows of Ainv selected by bits of C[i])
    Sm = []
    for i in range(m):
        acc = 0
        for k in range(m):
            if (C[i] >> k) & 1: acc ^= Ainv[k]
        Sm.append(acc)
    return phase_mat(Sm) @ addr_mat(A)
def symp_of(M):
    """my-convention row bitmasks of the symplectic image of a monomial/Clifford matrix (via ledger)"""
    SC = E.SC
    re = np.rint(M.real * SC).astype(np.int64); im = np.rint(M.imag * SC).astype(np.int64)
    Sl = E._symp((re, im))
    return tuple(sum(int(Sl[k, j]) << j for j in range(2 * m)) for k in range(2 * m))

def tens(M, v, f):
    t = v.reshape((n,) * f)
    for ax in range(f):
        t = np.tensordot(M, t, axes=([1], [ax])); t = np.moveaxis(t, 0, ax)
    return t.reshape(-1)
def apply_free_combo(combo, v, f):
    return sum(c * tens(P, v, f) for c, P in combo)
def Lam_P(p, v, f):
    ph = lift_P(p)
    return apply_free_combo(pexp(np.linalg.inv(ph)), tens(ph, v, f), f)
def Lam_T(z, v, f):
    return apply_free_combo(pexp(np.linalg.inv(dirm(z))), tens(dirm(z), v, f), f)

# convention checks
rng = random.Random(7)
for _ in range(20):
    p = S.rand_P(rng)
    assert symp_of(lift_P(p)) == p, "lift convention"
for z in range(1, n): assert symp_of(dirm(z)) == S.T(z), "dir convention"
print("conventions: monomial lifts and dir z have exactly the frame-game images (20 random p, all z)")

def run_lifted(prog_states_fn, sols, f, seed):
    """prog_states_fn(apply) -> list of state vectors; sols[l] = list of (state index, p, x)"""
    rngn = np.random.default_rng(seed)
    xs = [rngn.normal(size=n ** f) + 1j * rngn.normal(size=n ** f) for _ in range(2)]
    st = prog_states_fn(xs, f)
    KP = pexp(np.linalg.multi_dot([dirm(1 << j) for j in range(m)][::-1]) if m > 1 else dirm(1))
    Km = np.eye(n, dtype=complex)
    for j in range(m): Km = dirm(1 << j) @ Km
    KP = pexp(Km)
    errs = []
    for l, sol in enumerate(sols):
        out = np.zeros(n ** f, complex)
        for (r, p, x) in sol:
            out += float(x) * Lam_P(p, st[r], f)
        out = apply_free_combo(KP, out, f)
        tgt = tens(Km, xs[l], f)
        errs.append(float(np.abs(out - tgt).max() / np.abs(tgt).max()))
    return max(errs)

# ---- case 1: the order-3 control
cols = [1 << ((k + 1) % 3) if k < 3 else 1 << k for k in range(m)]
s3 = S.addr(cols)
y1, y2 = inp(S, 0), inp(S, 1)
states = [y1, y2]
for r in (cadd((1, y1), (1, y2)), cadd((1, y1), (-1, lmul(S, s3, y2)))):
    c = r
    for j in range(m):
        c = lmul(S, S.T(1 << j), c); states.append(c)
sols = [end_solve(S, states, {(l, S.K): Fr(1)})[1]["solution"] for l in (0, 1)]
def prog1(xs, f):
    st = [xs[0], xs[1]]
    for c0 in (xs[0] + xs[1], xs[0] - Lam_P(s3, xs[1], f)):
        c = c0
        for j in range(m):
            c = Lam_T(1 << j, c, f); st.append(c)
    return st
print("case 1 (order-3 control, frame-game readout, lifted): max relative error by f:",
      {f: float("%.1e" % run_lifted(prog1, sols, f, 11 + f)) for f in (1, 2, 3, 4, 5)})

