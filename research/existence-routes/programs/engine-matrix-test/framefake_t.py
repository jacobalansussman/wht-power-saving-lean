#!/usr/bin/env python3
"""framefake_t.py: framefake.py with exact Jacobians (torch.func.jacrev) for the Colab; same model and output format.
usage: framefake_t.py mode m W M [k_sub] [trials] [seed]"""
import sys, time
import numpy as np, torch
from scipy.optimize import least_squares
from framefake import build, rand_P
torch.set_default_dtype(torch.float64)

def solve(mode, m, W, M, k_sub=12, trials=6, seed=0):
    D = build(m, seed + 1)
    n, i0, iK, dist = D["n"], D["i0"], D["iK"], D["dist"]
    tinv = torch.tensor(np.argsort(D["t"]))
    plist = rand_P(D, k_sub); pinvs = [torch.tensor(np.argsort(p)) for p in plist]; nb = len(plist)
    pairs = [(ip, i) for i in range(M) for ip in range(i)]
    sizes = [W * W, M * W, len(pairs) * nb, W * M if mode == "hecke" else W * M * nb]
    npar = sum(sizes)
    Sd = torch.tensor(np.array([(dist == d).astype(float) for d in range(m + 1)]))
    tgt = []
    for l in range(W):
        if mode == "hecke":
            T = torch.zeros(m + 1, W); T[m, l] = 1
        else:
            T = torch.zeros(n, W); T[iK, l] = 1
        tgt.append(T)
    def comb(c, v): return sum(c[k] * v[pinvs[k]] for k in range(nb))
    def res(z):
        k = 0
        tau = z[k:k + sizes[0]].reshape(W, W); k += sizes[0]
        b = z[k:k + sizes[1]].reshape(M, W); k += sizes[1]
        gam = z[k:k + sizes[2]].reshape(len(pairs), nb); k += sizes[2]
        a = z[k:k + sizes[3]]
        a = a.reshape(W, M) if mode == "hecke" else a.reshape(W, M, nb)
        pidx = {pr: q for q, pr in enumerate(pairs)}
        w = []
        e0 = torch.zeros(n, 1); e0[i0, 0] = 1
        for i in range(M):
            u = e0 * b[i].reshape(1, W)
            for ip in range(i):
                u = u + comb(gam[pidx[(ip, i)]], w[ip])
            w.append(u[tinv] - u)
        out = []
        for l in range(W):
            o = e0 * tau[l].reshape(1, W)
            for i in range(M):
                o = o + (a[l, i] * w[i] if mode == "hecke" else comb(a[l, i], w[i]))
            out.append(((Sd @ o) if mode == "hecke" else o) - tgt[l])
        return torch.cat([x.reshape(-1) for x in out])
    jr = torch.func.jacrev(res)
    f = lambda z: res(torch.tensor(z)).numpy()
    J = lambda z: jr(torch.tensor(z)).numpy()
    rng = np.random.default_rng(seed); best = None
    for tr in range(trials):
        z0 = rng.normal(size=npar) * 0.5
        meth = "lm" if len(f(z0)) >= npar else "trf"
        r = least_squares(f, z0, jac=J, method=meth, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        nr = np.linalg.norm(r.fun)
        if best is None or nr < best[0]: best = (nr, r.x)
        if nr < 1e-11: break
    return best[0], float(np.linalg.norm(f(best[1]))), npar, n, nb

if __name__ == "__main__":
    mode = sys.argv[1]; m, W, M = map(int, sys.argv[2:5])
    k_sub = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    trials = int(sys.argv[6]) if len(sys.argv) > 6 else 6
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 0
    t0 = time.time()
    nr, nr2, npar, n, nb = solve(mode, m, W, M, k_sub, trials, seed)
    print("%s m=%d W=%d M=%d frames=%d family=%d unknowns=%d: best residual %.2e (recheck %.2e) %s [%.1fs]"
          % (mode, m, W, M, n, nb, npar, nr, nr2, "REALISED FAKE" if nr2 < 1e-9 else "none found", time.time() - t0), flush=True)
