# sparsify.py W s seed [maxbatch] -- find a stand-in solution by damped Gauss-Newton, then move along the solution set
# pinning coordinates to 0 in batches (smallest |x| first) while the residual stays at round-off.  Saves sx_<W>in<s>_<seed>.npz
import sys, time, numpy as np, warnings
warnings.filterwarnings('ignore')
exec(open('sym2core.py').read()); Sys = Sys2
W, s, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); maxb = int(sys.argv[4]) if len(sys.argv) > 4 else 256
S = Sys(W, s); rng = np.random.default_rng(seed); t0 = time.time()
print('W %d s %d: P = %d unknowns, %d raw equations, coords per mixing matrix %d' % (W, s, S.P, len(S.Tv), nA), flush=True)
# compress the equations to the span of all outputs (the affine hull T of the image, through 0)
Js = np.concatenate([S.jac(rng.normal(size=S.P) + 1j * rng.normal(size=S.P)) for _ in range(2)], 1)
u, sv, vh = np.linalg.svd(Js, full_matrices=False); rT = int((sv > 1e-10 * sv[0]).sum()); Qt = u[:, :rT]
print('dim T = %d; target in T: %.1e' % (rT, np.linalg.norm(S.Tv - Qt @ (Qt.conj().T @ S.Tv))), flush=True)
res0 = S.res
S.res = lambda p: Qt.conj().T @ res0(p)
jac0 = S.jac
S.jac = lambda p: Qt.conj().T @ jac0(p)
best = None
for start in range(8):
    p = (rng.normal(size=S.P) + 1j * rng.normal(size=S.P)) * 0.3
    p, nr = gn(S, p, np.ones(S.P, bool), iters=200, tol=1e-13)
    print('start %d: residual %.2e max|p| %.2f (%.0fs)' % (start, nr, np.abs(p).max(), time.time() - t0), flush=True)
    if nr < 1e-12: best = p; break
p = best; free = np.ones(S.P, bool)
J = S.jac(p); svj = np.linalg.svd(J, compute_uv=False); rk = int((svj > 1e-9 * svj[0]).sum())
print('rank J at solution %d of %d (gap %.1e | %.1e)' % (rk, rT, svj[rk - 1] / svj[0], svj[rk] / svj[0] if rk < len(svj) else 0), flush=True)
b = maxb; stuck = set()
while True:
    cand = [k for k in np.argsort(np.abs(p)) if free[k] and k not in stuck]
    if not cand: break
    batch = cand[:b]
    q = p.copy(); q[batch] = 0; f2 = free.copy(); f2[batch] = False
    q, nr = gn(S, q, f2, iters=40, tol=1e-13)
    if nr < 1e-12:
        p, free = q, f2
        print('pinned %d to 0 -> %d free, residual %.1e, max|p| %.2f (%.0fs)' % (len(batch), free.sum(), nr, np.abs(p).max(), time.time() - t0), flush=True)
        b = min(maxb, b * 2)
    else:
        if b == 1:
            stuck.add(batch[0])
        b = max(1, b // 2)
    if len(stuck) > 0 and len(stuck) % 50 == 0 and b == 1: print('  stuck', len(stuck), flush=True)
np.savez('sym2_%din%d_%d.npz' % (W, s, seed), p=p, free=free)
fi = np.where(free)[0]
J = S.jac(p)[:, fi]; svj = np.linalg.svd(J, compute_uv=False); rk = int((svj > 1e-9 * svj[0]).sum())
print('FINAL: %d nonzero coordinates, rank of J on them %d (local solution dim %d), residual %.1e (%.0fs)' %
      (len(fi), rk, len(fi) - rk, np.linalg.norm(S.res(p)), time.time() - t0), flush=True)
for ks in range(len(S.slots)):
    nz = [(i, a) for (k2, i, a), k in zip(S.idx, range(S.P)) if k2 == ks and free[k]]
    print('slot %d (%s): %d nonzero' % (ks, 'move %d' % (ks + 1) if ks < s else 'output %d' % (ks - s), len(nz)), flush=True)
