"""vf2blind.py -- v-e3-bpexact: independent check of e3-bpexact's f = 2 blindness construction.
Runs their f2blind.py for its matrices, then checks with MY OWN projectors (P_D = diagonal, P_S = (1+SWAP)/2 - P_D,
P_A = (1-SWAP)/2) that every content and readout is block-diagonal (so lies in the image of C[Z4 wr S16] at f = 2, whose
image is M16 + M120 + M120), that the identity K(x)K = R_x + sum_i R_i (D_i(x)D_i) A_i holds, and that each move uses ONE
D(z) on ONE array.  Also: does a single move with free content suffice (it must not, or the 2-move claim is moot)."""
import runpy, sys, io, contextlib, numpy as np
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    g = runpy.run_path(__import__('os').path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), '..', 'f2blind.py'), run_name='f2b')
K2, Rx, Rm, Am, Ds, G = g['K2'], g['Rx'], g['Rm'], g['Am'], g['Ds'], g['G']
n = 16
SW = np.zeros((256, 256))
for x in range(16):
    for y in range(16): SW[y * 16 + x, x * 16 + y] = 1
PD = np.zeros((256, 256))
for x in range(16): PD[x * 17, x * 17] = 1
PS = (np.eye(256) + SW) / 2 - PD; PA = (np.eye(256) - SW) / 2
Ps = {'D': PD, 'S': PS, 'A': PA}
def offdiag(M):
    return max(np.abs(Ps[a] @ M @ Ps[b]).max() for a in Ps for b in Ps if a != b)
print('projector sanity: ranks', [int(round(np.trace(P))) for P in Ps.values()], ' sum = I', np.allclose(sum(Ps.values()), np.eye(256)))
print('off-diagonal parts: Rx %.1e, R1 %.1e, R2 %.1e, A1 %.1e, A2 %.1e' % (offdiag(Rx), offdiag(Rm[0]), offdiag(Rm[1]), offdiag(Am[0]), offdiag(Am[1])))
print('moves are D(z)(x)D(z) for z = 1, 2:', np.allclose(Ds[0], np.kron(G.D[1], G.D[1])), np.allclose(Ds[1], np.kron(G.D[2], G.D[2])))
lhs = Rx + Rm[0] @ Ds[0] @ Am[0] + Rm[1] @ Ds[1] @ Am[1]
print('identity |K(x)K - (Rx + sum R_i M_i A_i)| = %.1e  (|K(x)K| max %.2f)' % (np.abs(lhs - K2).max(), np.abs(K2).max()))
# one move only, contents free: best possible = can block-diag R and A make K2 - R1 M1 A1 block-diagonal?  The off-diagonal
# D->S block of R M A is  R_SS (M A)_SD  with rank <= rank(M_SD-part)... check directly rank of M1's off-diagonal blocks
M1 = Ds[0]
print('rank of the D(1)(x)D(1) block S<-D: %d, D<-S: %d ; of K(x)K: %d, %d' % (
    np.linalg.matrix_rank(PS @ M1 @ PD, 1e-9), np.linalg.matrix_rank(PD @ M1 @ PS, 1e-9),
    np.linalg.matrix_rank(PS @ K2 @ PD, 1e-9), np.linalg.matrix_rank(PD @ K2 @ PS, 1e-9)))
