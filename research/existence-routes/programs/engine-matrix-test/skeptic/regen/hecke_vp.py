#!/usr/bin/env python3
"""hecke_vp.py: realised fakes of the two-sided trivially pinned component by variable projection.
For links gamma (combinations of a family of elements of P) and input couplings b, the move outputs w_i are frame
vectors (x) inputs; the component holds for SOME outer couplings (a, tau) iff each target sigma_m (x) e_l lies in
span{ S(w_i) } + span{ delta_0 (x) e_l' }, S = distance-class sums.  Minimise that distance over (b, gamma).
usage: hecke_vp.py m W M [k_sub] [trials] [seed] [save.npz]"""
import sys, time
import numpy as np
from scipy.optimize import least_squares
from framefake import build, rand_P

def main(m, W, M, k_sub=12, trials=4, seed=0, save=None):
    D = build(m, seed + 1)
    n, i0, iK, dist = D["n"], D["i0"], D["iK"], D["dist"]
    tinv = np.argsort(D["t"])
    plist = rand_P(D, k_sub); pinvs = [np.argsort(p) for p in plist]; nb = len(plist)
    pairs = [(ip, i) for i in range(M) for ip in range(i)]
    Sd = np.array([(dist == d).astype(float) for d in range(m + 1)])
    npar = M * W + len(pairs) * nb
    T = np.zeros(((m + 1) * W, W))
    for l in range(W): T[m * W + l, l] = 1          # sigma_m (x) e_l, layout (d, l')
    E0 = np.zeros(((m + 1) * W, W))
    for l in range(W): E0[0 * W + l, l] = 1
    def ws(z):
        b = z[:M * W].reshape(M, W); g = z[M * W:].reshape(len(pairs), nb)
        pidx = {pr: q for q, pr in enumerate(pairs)}
        w = []
        for i in range(M):
            u = np.zeros((n, W)); u[i0] = b[i]
            for ip in range(i):
                c = g[pidx[(ip, i)]]
                u = u + sum(c[k] * w[ip][pinvs[k]] for k in range(nb))
            w.append(u[tinv] - u)
        return w
    def res(z):
        w = ws(z)
        Bm = np.hstack([E0] + [(Sd @ wi).reshape(-1, 1) for wi in w])   # (m+1)W x (W + M)
        coef, *_ = np.linalg.lstsq(Bm, T, rcond=None)
        return (T - Bm @ coef).ravel()
    rng = np.random.default_rng(seed); best = None
    for tr in range(trials):
        z0 = rng.normal(size=npar)
        r = least_squares(res, z0, method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=400)
        nr = np.linalg.norm(r.fun)
        if best is None or nr < best[0]: best = (nr, r.x)
        if nr < 1e-11: break
    nr = best[0]
    # condition number of the final span (to guard against a fake that only works by cancellation of huge terms)
    w = ws(best[1]); Bm = np.hstack([E0] + [(Sd @ wi).reshape(-1, 1) for wi in w])
    sv = np.linalg.svd(Bm, compute_uv=False)
    if save:
        np.savez(save, m=m, W=W, M=M, z=best[1], plist=np.array(plist), pairs=np.array(pairs), resid=nr)
    return nr, npar, n, nb, sv

if __name__ == "__main__":
    m, W, M = map(int, sys.argv[1:4])
    k = int(sys.argv[4]) if len(sys.argv) > 4 else 12; tr = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    sd = int(sys.argv[6]) if len(sys.argv) > 6 else 0; sv_ = sys.argv[7] if len(sys.argv) > 7 else None
    t0 = time.time(); nr, npar, n, nb, sv = main(m, W, M, k, tr, sd, sv_)
    print("hecke-vp m=%d W=%d M=%d frames=%d family=%d unknowns=%d: distance of targets %.2e %s (span sing. values %.1e..%.1e) [%.1fs]"
          % (m, W, M, n, nb, npar, nr, "REALISED FAKE" if nr < 1e-9 else "none found", sv.min(), sv.max(), time.time() - t0), flush=True)
