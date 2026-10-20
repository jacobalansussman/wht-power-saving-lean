#!/usr/bin/env python3
"""framefake.py: realisable fakes for the trivially pinned components of the engine identity, in the frame module.
Pauli-blind level: G -> Sp(2m,2), free group B -> P = Stab(L0), frames = Lagrangians.  A word with M moves (each = t up
to free letters, absorbed) acting on the base frame delta_{L0} (the P-fixed vector e):
   u_i = delta_{L0} b_i + sum_{i' < i} Gam_{i'i} w_{i'},   w_i = (M(t) - I) u_i,        (b_i in C^{1 x W})
   output row l:  o_l = tau_l delta_{L0} + sum_i Alpha_{l,i} w_i,
with Gam_{i'i}, Alpha_{l,i} in A_P = span{M(p) : p in P} (here: random sub-families, a valid restriction for fakes).
mode 'hecke' (two-sided e..e compression): Alpha_{l,i} = a_{l,i} I (only eps(alpha) survives); the conditions are the
   distance-class sums S_d(o_l) = [d = m] e_l^T  (d = 0..m).
mode 'column' (one-sided compression M(Z) delta_{L0} = delta_{K L0} (x) I_W): o_l = delta_{K L0} e_l^T exactly.
usage: framefake.py mode m W M [k_sub] [trials] [seed]"""
import sys, time
import numpy as np
from scipy.optimize import least_squares
from sp import Sp

def build(m, seed=1):
    S = Sp(m)
    L0, F = S.frames()
    idx = {L: i for i, L in enumerate(F)}; n = len(F)
    def perm_of(g): return np.array([idx[frozenset(S.apply(g, x) for x in L)] for L in F])
    Pg = [perm_of(g) for g in S.P_gens()]
    # random elements of P as products of generators (random walks), plus identity
    rng = np.random.default_rng(seed)
    t = perm_of(S.transvection(1))
    Kb = [((1 << j) if j < m else ((1 << j) ^ (1 << (j - m)))) for j in range(S.n)]
    iK = idx[frozenset(S.apply(Kb, x) for x in L0)]
    dist = np.array([S.dist(L, L0) for L in F])
    return dict(S=S, n=n, i0=idx[L0], iK=iK, t=t, Pg=Pg, dist=dist, rng=rng, perm_of=perm_of)

def rand_P(D, k):
    """identity, the coordinate permutations (x,b) -> (pi x, pi b) (so that plain words are representable), then
    k - 1 random elements of P (random walks in the generators)."""
    import itertools
    n = D["n"]; rng = D["rng"]; out = [np.arange(n)]
    S = D["S"]; m = S.m
    import os
    for pi in ([] if os.environ.get("NOPERM") else itertools.permutations(range(m))):
        if list(pi) == list(range(m)): continue
        g = []
        for j in range(S.n):
            g.append(1 << pi[j] if j < m else 1 << (m + pi[j - m]))
        out.append(D["perm_of"](g))
    for _ in range(k - 1):
        p = np.arange(n)
        for _ in range(30):
            g = D["Pg"][rng.integers(len(D["Pg"]))]
            p = g[p]
        out.append(p)
    return out

def pmat(p):
    n = len(p); M = np.zeros((n, n)); M[p, np.arange(n)] = 1; return M

def solve(mode, m, W, M, k_sub=12, trials=6, seed=0):
    D = build(m, seed + 1)
    n, i0, iK, dist = D["n"], D["i0"], D["iK"], D["dist"]
    t = D["t"]; tinv = np.argsort(t)
    def Tm(v): return v[tinv] - v                                # (M(t) - I) v, M(p) v = v[p^-1]
    plist = rand_P(D, k_sub)
    pinvs = [np.argsort(p) for p in plist]
    nb = len(plist)
    def comb(c, v): return sum(ck * v[pi] for ck, pi in zip(c, pinvs))
    pairs = [(ip, i) for i in range(M) for ip in range(i)]
    npar = W * W + M * W + len(pairs) * nb + (W * M if mode == "hecke" else W * M * nb)
    Sd = np.array([(dist == d).astype(float) for d in range(m + 1)])  # (m+1) x n
    def unpack(z):
        k = 0
        tau = z[k:k + W * W].reshape(W, W); k += W * W
        b = z[k:k + M * W].reshape(M, W); k += M * W
        gam = {}
        for pr in pairs: gam[pr] = z[k:k + nb]; k += nb
        if mode == "hecke": a = z[k:k + W * M].reshape(W, M); k += W * M
        else: a = z[k:k + W * M * nb].reshape(W, M, nb); k += W * M * nb
        return tau, b, gam, a
    def outputs(z):
        tau, b, gam, a = unpack(z)
        w = []
        for i in range(M):
            u = np.zeros((n, W)); u[i0] = b[i]
            for ip in range(i):
                u = u + comb(gam[(ip, i)], w[ip])
            w.append(Tm(u))
        outs = []
        for l in range(W):
            o = np.zeros((n, W)); o[i0] = tau[l]
            for i in range(M):
                if mode == "hecke": o = o + a[l, i] * w[i]
                else: o = o + comb(a[l, i], w[i])
            outs.append(o)
        return outs
    def res(z):
        outs = outputs(z); r = []
        for l, o in enumerate(outs):
            if mode == "hecke":
                tgt = np.zeros((m + 1, W)); tgt[m, l] = 1
                r.append((Sd @ o - tgt).ravel())
            else:
                tgt = np.zeros((n, W)); tgt[iK, l] = 1
                r.append((o - tgt).ravel())
        return np.concatenate(r)
    rng = np.random.default_rng(seed)
    best = None
    for tr in range(trials):
        z0 = rng.normal(size=npar) * 0.5
        meth = "lm" if len(res(z0)) >= npar else "trf"
        r = least_squares(res, z0, method=meth, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=1500 if meth == "trf" else 200 * npar)
        nr = np.linalg.norm(r.fun)
        if best is None or nr < best[0]: best = (nr, r.x)
        if nr < 1e-11: break
    if len(sys.argv) > 8:
        tau, b, gam, a = unpack(best[1])
        np.savez(sys.argv[8], mode=mode, m=m, W=W, M=M, tau=tau, b=b, a=a, plist=np.array(plist),
                 gam=np.array([gam[pr] for pr in pairs]), pairs=np.array(pairs), resid=best[0])
    return best[0], np.linalg.norm(res(best[1])), npar, n

if __name__ == "__main__":
    mode = sys.argv[1]; m, W, M = map(int, sys.argv[2:5])
    k_sub = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    trials = int(sys.argv[6]) if len(sys.argv) > 6 else 6
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 0
    t0 = time.time()
    nr, nr2, npar, n = solve(mode, m, W, M, k_sub, trials, seed)
    print("%s m=%d W=%d M=%d frames=%d k_sub=%d unknowns=%d: best residual %.2e (recheck %.2e) %s [%.1fs]"
          % (mode, m, W, M, n, k_sub, npar, nr, nr2, "REALISED FAKE" if nr2 < 1e-9 else "none found", time.time() - t0), flush=True)
