#!/usr/bin/env python3
"""hecke_abs_t.py: hecke_abs.py with exact Jacobians (torch.func.jacrev).  usage: hecke_abs_t.py m W M [trials]"""
import sys, time
import numpy as np, torch
from scipy.optimize import least_squares
from hecke_abs import setup
torch.set_default_dtype(torch.float64)
def solve(m, W, M, trials=6, seed=0):
    om, mu, kap, pairs, pb = setup(m, W, M)
    nA, nB = W * M, M * W; nh = sum(len(b) for b in pb); n = nA + nB + nh
    # basis tensor: for each pair, its spectral basis vectors (j = 1..m)
    rows, cols, js, coef, pidx = [], [], [], [], []
    k = 0
    Bmat = torch.zeros(nh, m, M, M)
    for (i, ip), b in zip(pairs, pb):
        for bb in b:
            Bmat[k, :, i, ip] = torch.tensor(bb); k += 1
    mu_t = torch.tensor(mu); kap_t = torch.tensor(kap)
    def res(z):
        A = z[:nA].reshape(W, M); B = z[nA:nA + nB].reshape(M, W); x = z[nA + nB:]
        L = torch.einsum("k,kjab->jab", x, Bmat)
        R = mu_t[:, None, None] * torch.eye(M) + L
        out = torch.einsum("wa,jab,bv->jwv", A, R, B) - (kap_t - 1)[:, None, None] * torch.eye(W)
        return out.reshape(-1)
    jr = torch.func.jacrev(res)
    f = lambda z: res(torch.tensor(z)).numpy(); J = lambda z: jr(torch.tensor(z)).numpy()
    rng = np.random.default_rng(seed); best = None
    for t in range(trials):
        z0 = rng.normal(size=n)
        meth = "lm" if m * W * W >= n else "trf"
        r = least_squares(f, z0, jac=J, method=meth, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
        nr = np.linalg.norm(r.fun)
        if best is None or nr < best[0]: best = (nr, r.x)
        if nr < 1e-12: break
    return best[0], float(np.linalg.norm(f(best[1])))
if __name__ == "__main__":
    m, W, M = map(int, sys.argv[1:4]); trials = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    t0 = time.time(); nr, nr2 = solve(m, W, M, trials)
    print("m=%d W=%d M=%d (bound max(2W, W+m-1) = %d): best residual %.2e (recheck %.2e) %s  [%.1fs]"
          % (m, W, M, max(2 * W, W + m - 1), nr, nr2, "ABSTRACT FAKE" if nr2 < 1e-10 else "none found", time.time() - t0), flush=True)
