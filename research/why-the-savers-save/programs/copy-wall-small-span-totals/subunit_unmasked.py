"""subunit_unmasked.py (d-R7): on the 35-port sub-unit U' (labels in V = ann(1,2), isolated port dropped; root rank 7)
and on P112, random-restart unmasked ALS at rank K; is each converged fitting B natural (rowspace contains Lambda_U D)?
usage: subunit_unmasked.py K restarts"""
import sys, time, collections; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
from als import als
K, NR = int(sys.argv[1]), int(sys.argv[2])
U = np.array(lines_of_dual([1, 2])); _O = ORTH[np.ix_(U, U)]; U = U[_O.sum(1) < len(U) - 1]
OU = ORTH[np.ix_(U, U)]; GU = GP[np.ix_(U, U)].astype(float)
u_, s_, vt_ = np.linalg.svd(GU); rU = int((s_ > 1e-9 * s_[0]).sum()); LAMU = u_[:, :rU] * s_[:rU]
def nat(B, tol=1e-6):
    U2, s, Vt = np.linalg.svd(B); r = int((s > tol * s[0]).sum()); Q = Vt[r:]
    A = np.einsum('it,tc->ict', Q, LAMU).reshape(-1, len(U))
    u, s2, vt = np.linalg.svd(A, full_matrices=False); d = vt[-1]
    return r, s2[-1] / s2[0], s2[-2] / s2[0], np.abs(d).min() / np.abs(d).max()
rng = np.random.default_rng(K); cls = collections.Counter(); t0 = time.time()
for rs in range(NR):
    res, X, R, h = als(np.ones((K, len(U)), bool), OU, iters=3000, rng=rng, ridge=1e-14)
    if res > 1e-8: cls['unconverged'] += 1; continue
    r, a, b, dm = nat(X @ R)
    c = 'natural' if a < 1e-5 and dm > 1e-3 else ('partial' if a < 1e-5 else 'NON-NATURAL')
    cls[(c, r)] += 1
    if c == 'NON-NATURAL': print('  rs %d rank %d smin %.2e snext %.2e dmin %.2e' % (rs, r, a, b, dm), flush=True)
print('U\' (%d ports, root rank %d) K=%d restarts %d: %s  %.0fs' % (len(U), rU, K, NR, dict(cls), time.time() - t0))
