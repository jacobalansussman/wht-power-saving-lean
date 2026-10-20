"""s2_als.py -- masked alternating least squares on complex P112 (pattern 'cx': B[S,T] = 0 for orthogonal S != T,
B[S,S] = 1).  B = X Y^T, X[S,j] = weight of total j into target S (serve mask), Y[T,j] = coefficient of x_T in total j
(read mask).  usage: VENV -I s2_als.py <design> [restarts] [iters] [seed]
designs: coord (control: 8 coordinate totals, read = serve = C_{o_k}), free8 (8 unmasked totals),
         k8special (7 free + 1 confined to C_h, h one of the 8 omega-7 classes), k8c60 (7 free + 1 confined, |C_h| = 60),
         k8c48 (7 free + 1 confined, |C_h| = 48), k8two (6 free + 2 confined special classes)"""
import sys
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from base import *
import numpy as np
design = sys.argv[1]; R = int(sys.argv[2]) if len(sys.argv) > 2 else 4
IT = int(sys.argv[3]) if len(sys.argv) > 3 else 300; seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
Z = np.zeros((N, N), bool)
for i in range(N):
    for k in orth[i]: Z[i, k] = True
coord = [UI[o] for o in range(8)]
iso = [h for h in range(1, 512) if not dot9(h, h)]
def Cmask(h): return np.array([not dot9(lab[i], h) for i in range(N)])
special = [147, 169, 202, 240, 260, 318, 349, 359]
c60 = [h for h in iso if len(C_of([h])) == 60]
c48 = [h for h in iso if len(C_of([h])) == 48]
allm = np.ones(N, bool)
if design == 'coord':
    rm = [Cmask(o) for o in coord]; sm = rm
elif design == 'free8':
    rm = [allm] * 8; sm = rm
elif design == 'k8special':
    rm = [allm] * 7 + [Cmask(special[0])]; sm = rm
elif design == 'k8c60':
    rm = [allm] * 7 + [Cmask(c60[0])]; sm = rm
elif design == 'k8c48':
    rm = [allm] * 7 + [Cmask(c48[0])]; sm = rm
elif design == 'k8two':
    rm = [allm] * 6 + [Cmask(special[0]), Cmask(special[1])]; sm = rm
else:
    raise SystemExit('unknown design')
K = len(rm)
RM = np.array(rm).T; SM = np.array(sm).T        # N x K
rows = [np.where(Z[S])[0] for S in range(N)]
rng = np.random.default_rng(seed)
best = (1e9, None)
W = 3.0
for rs in range(R):
    X = rng.standard_normal((N, K)) * SM; Y = rng.standard_normal((N, K)) * RM
    for it in range(IT):
        for S in range(N):
            a = np.where(SM[S])[0]
            A = np.vstack([Y[rows[S]][:, a], W * Y[S, a]]); b = np.zeros(len(rows[S]) + 1); b[-1] = W
            X[S] = 0; X[S, a] = np.linalg.lstsq(A, b, rcond=None)[0]
        for T in range(N):
            a = np.where(RM[T])[0]
            A = np.vstack([X[rows[T]][:, a], W * X[T, a]]); b = np.zeros(len(rows[T]) + 1); b[-1] = W
            Y[T] = 0; Y[T, a] = np.linalg.lstsq(A, b, rcond=None)[0]
        B = X @ Y.T
        res = np.sqrt((B[Z] ** 2).sum() + ((np.diag(B) - 1) ** 2).sum())
        if res < 1e-11: break
    if res < best[0]: best = (res, (X.copy(), Y.copy()))
    # natural?  B = D Gram/2 restricted to P112 (row scaling): test rank of [B ; Gram] row spaces
    # natural up to rescaling (B = D1 Gram D2) forces B[S,T]B[T,S] = B[S,S]B[T,T] Gram[S,T]^2 / Gram[S,S]Gram[T,T]
    Gm = np.array([[G[P[i]][P[k]] for k in range(N)] for i in range(N)], float)
    dg = np.diag(B)
    nat = np.abs(B * B.T - np.outer(dg, dg) * Gm ** 2 / np.outer(np.diag(Gm), np.diag(Gm))).max()
    print('%s restart %d iters %d residual %.3e  rank(B) %d  naturality defect %.2e (0 = B = D1 Gram D2 test passes)'
          % (design, rs, it + 1, res, np.linalg.matrix_rank(B, tol=1e-6), nat), flush=True)
print('BEST', design, '%.3e' % best[0])
if best[0] < 1e-9:
    np.save(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'out/als_%s_X.npy') % design, best[1][0])
    np.save(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'out/als_%s_Y.npy') % design, best[1][1])
