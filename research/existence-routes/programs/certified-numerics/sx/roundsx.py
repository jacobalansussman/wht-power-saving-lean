# roundsx.py W s in.npz out.npz -- second phase: while the solution set (with the current pins) still has positive local
# dimension, pin a free coordinate to the nearest simple Gaussian dyadic rational (grid 2^-k, k = 0..4, nearest first),
# re-solve by damped Gauss-Newton, keep the pin if the residual returns to round-off.  Pinned values are kept in 'fix'.
import sys, time, numpy as np, warnings
warnings.filterwarnings('ignore')
exec(open('sxcore.py').read())
W, s = int(sys.argv[1]), int(sys.argv[2]); z = np.load(sys.argv[3]); p = z['p']; free = z['free'].copy()
S = Sys(W, s); rng = np.random.default_rng(5); t0 = time.time()
Js = np.concatenate([S.jac(rng.normal(size=S.P) + 1j * rng.normal(size=S.P)) for _ in range(2)], 1)
u, sv, vh = np.linalg.svd(Js, full_matrices=False); rT = int((sv > 1e-10 * sv[0]).sum()); Qt = u[:, :rT]
res0, jac0 = S.res, S.jac
S.res = lambda p: Qt.conj().T @ res0(p); S.jac = lambda p: Qt.conj().T @ jac0(p)
var = free.copy()            # coordinates still variable (nonzero, not pinned to a value)
def grid(x):
    out = []
    for k in range(0, 5):
        g = 2.0 ** -k
        c = (np.round(x.real / g) + 1j * np.round(x.imag / g)) * g
        out.append((abs(c - x), k, c))
    return sorted(out, key=lambda t: (t[0] > 1e-9, t[1], t[0]))
tried = set()
while True:
    fi = np.where(var)[0]
    J = S.jac(p)[:, fi]; u_, sv_, vh_ = np.linalg.svd(J); rk = int((sv_ > 1e-9 * sv_[0]).sum())
    N = vh_[rk:].conj().T                     # null directions in the variable coordinates
    ldim = len(fi) - rk
    if ldim == 0: break
    lev = np.linalg.norm(N, axis=1)           # how much each variable coordinate moves along the solution set
    order = sorted([j for j in range(len(fi)) if lev[j] > 1e-6 and fi[j] not in tried],
                   key=lambda j: min(abs(c - p[fi[j]]) for _, _, c in grid(p[fi[j]])[:2]))
    if not order: print('no pinnable coordinate left', flush=True); break
    done = False
    for j in order[:40]:
        k = fi[j]; tried.add(k)
        for dist, kk, c in grid(p[k])[:3]:
            q = p.copy(); q[k] = c; v2 = var.copy(); v2[k] = False
            q, nr = gn(S, q, v2, iters=40, tol=1e-13)
            if nr < 1e-12:
                p, var = q, v2; done = True
                print('pinned x%d = %s (grid 2^-%d, moved %.2e); local dim was %d; residual %.1e (%.0fs)' %
                      (k, c, kk, dist, ldim, nr, time.time() - t0), flush=True)
                break
        if done: break
    if not done: print('no pin accepted among 40 candidates', flush=True); break
np.savez(sys.argv[4], p=p, free=free, var=var)
print('FINAL local dim %d, variable coordinates %d, residual %.1e' % (ldim, var.sum(), np.linalg.norm(S.res(p))), flush=True)
