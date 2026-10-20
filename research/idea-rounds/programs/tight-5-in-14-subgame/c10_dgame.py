# c10_dgame.py -- the singular subgame's V_56 relaxation as a two-copy "D-game".
# V56|S8 = A + B + C with A=(2,1^6) ~ U7, B=(3,1^5) ~ U21, C=(2,2,2,1,1) ~ U7+U21 as B'-modules. On U := A+B the whole
# group algebra acts through the free algebra Alg = End(U7)+End(U21); transport C to U by the B'-isomorphism Phi.
# rho1(g) = action on A+B, rho2(g) = Phi^-1 rho_C(g) Phi. Gauge: E_g = rho1(g)^-1 rho2(g). Then (NOTE.md):
#   N_0 = 0;  move: N += Alg (E_h - I) Y  for any Y in N + Alg^W;  target: (E_K - I) e_l in N for all l.
# Here: E_h, E_K, the blocks of G = E_h - I and F = E_K - I, and the hook family as graphs (checks cap 10 / demand 28).
import numpy as np
from lib514 import *
import c8_s4 as C8
R = C8.R
Ab, Bb, C7a, C21b = C8.Ab, C8.Bb, C8.C7a, C8.C21b
U1 = np.hstack([Ab, Bb])           # basis of A+B (orthonormal), types: first 7 = U7, last 21 = U21
U2 = np.hstack([C7a, C21b])        # B'-equivariant image in C
def rho1(g): return U1.T @ R.r56(g) @ U1
def rho2(g): return U2.T @ R.r56(g) @ U2
# checks: rho1(g) block diagonal for all g (A, B are S8-irreducible); rho1(b) = rho2(b) for b in B'
for g in [h1, K, comp(h1, h2)]:
    M = rho1(g); print('rho1 block-off-diagonal norm:', np.abs(M[:7, 7:]).max(), np.abs(M[7:, :7]).max())
for b in C8.Bp_gens: print("B' generator: |rho1 - rho2| =", np.abs(rho1(b) - rho2(b)).max())
Eh = np.linalg.solve(rho1(h1), rho2(h1)); EK = np.linalg.solve(rho1(K), rho2(K))
G = Eh - np.eye(28); F = EK - np.eye(28)
r = lambda M: np.linalg.matrix_rank(M, tol=1e-9)
print('rank G =', r(G), ' rank F =', r(F))
print('blocks of G: rank G77 %d, G7,21 %d, G21,7 %d, G21,21 %d; rank G_7rows %d, G_21rows %d; rank G_7cols %d, G_21cols %d' %
      (r(G[:7, :7]), r(G[:7, 7:]), r(G[7:, :7]), r(G[7:, 7:]), r(G[:7]), r(G[7:]), r(G[:, :7]), r(G[:, 7:])))
print('blocks of F: rank F77 %d, F7,21 %d, F21,7 %d, F21,21 %d' % (r(F[:7, :7]), r(F[:7, 7:]), r(F[7:, :7]), r(F[7:, 7:])))
# eigen-structure: E_h is a product of two involutions
w = np.linalg.eigvals(Eh); print('eigenvalues of E_h (rounded):', sorted(np.round(w, 4), key=lambda z: (round(z.real, 4), round(z.imag, 4))))
# hook family as graphs Y(t) = {(y, T y)} with T = diag(t7 I7, t21 I21): cap = rank(G T) , demand = rank(F T)
for t7, t21 in [(-1/3, 1.0), (0.7, -2.0), (1.0, 1.0)]:
    T = np.diag([t7] * 7 + [t21] * 21)
    print(f'graph certificate t=({t7:.3f},{t21:.3f}): cap rank(G T) = {r(G @ T)}, demand rank(F T) = {r(F @ T)}')
# the (S, Zhat) types of the flag argument (cap, demand) for all 8 non-trivial one-dimensional types
P7 = np.zeros((28, 28)); P7[:7, :7] = np.eye(7); P21 = np.eye(28) - P7
def capdem(S, Z):
    """S, Z: subsets of {7,21}; cap = dim((G S + Z)/Z), demand = dim((F S + Z)/Z)"""
    def proj(sel): return sum((P7 if s == 7 else P21) for s in sel) if sel else np.zeros((28, 28))
    Sm = proj(S); Zm = proj(Z)
    def q(M):
        if not Z: return r(M @ Sm)
        return r(np.hstack([M @ Sm, Zm])) - r(Zm)
    return q(G), q(F)
for S in [(7,), (21,), (7, 21)]:
    for Z in [(), (7,), (21,)]:
        if not set(Z) <= set(S): continue
        c, dm = capdem(S, Z)
        print(f'type S={S} Z={Z}: cap {c}, demand {dm}', f'ratio {dm/c:.3f}' if c else '')
np.savez('dgame.npz', G=G, F=F, Eh=Eh, EK=EK)
