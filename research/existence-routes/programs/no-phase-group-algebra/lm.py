#!/usr/bin/env python3
"""lm.py (e3-rmbuild): Levenberg-Marquardt for the transfer-form system of tf.py (same model, same characters),
with the normal equations accumulated over chunks of characters (memory stays at (#unknowns)^2).
Usage as tf.py, plus --chunk.  Prints one JSON line per restart; saves the best point with --save.
"""
import argparse, json, time, sys
import numpy as np
import torch
from torch.func import jacrev, jacfwd
import tf


def unpack(theta, W, M, n):
    k = 0
    def take(*shape):
        nonlocal k
        sz = int(np.prod(shape))
        re = theta[k:k + sz]; im = theta[k + sz:k + 2 * sz]; k += 2 * sz
        return torch.complex(re, im).reshape(*shape)
    A = take(W, M, n); B = take(M, W, n)
    ii, jj = torch.tril_indices(M, M, -1)
    g = take(len(ii), n)
    G = torch.zeros(M, M, n, dtype=torch.complex128).index_put((ii, jj), g)
    return A, B, G


def nparams(W, M, n):
    return 2 * (W * M * n + M * W * n + (M * (M - 1) // 2) * n)


class Sys:
    def __init__(s, W, M, n, Cc, TAU, WT):
        s.W, s.M, s.n = W, M, n
        s.Cc, s.TAU, s.sw = Cc, TAU, torch.sqrt(WT)
        s.mask = torch.tril(torch.ones(M, M, dtype=torch.float64), -1).to(torch.complex128)[:, :, None]
        a = torch.arange(n)
        s.idx = (a[:, None] - a[None, :]) % n
    def big(s, X, r, c):
        return X[:, :, s.idx].permute(0, 2, 1, 3).reshape(r * s.n, c * s.n)
    def resid(s, theta, sl):
        W, M, n = s.W, s.M, s.n
        A, B, G = unpack(theta, W, M, n)
        Ab = s.big(A, W, M); Bb = s.big(B, M, W); Gb = s.big(G * s.mask, M, M)
        Cc = s.Cc[sl]; TAU = s.TAU[sl]
        I = torch.eye(M * n, dtype=torch.complex128)
        L = I[None] - Gb[None] * Cc[:, None, :]
        Y = torch.linalg.solve(L, Bb[None].expand(Cc.shape[0], -1, -1))
        X = Ab[None] @ (Cc[:, :, None] * Y)
        D = torch.diag_embed(TAU - 1)
        targ = torch.einsum("ab,qij->qaibj", torch.eye(W, dtype=torch.complex128), D).reshape(-1, W * n, W * n)
        R = (X - targ) * s.sw[sl][:, None, None]
        R = R.reshape(-1)
        return torch.cat([R.real, R.imag])


def lm_solve(S, theta, iters=200, chunk=64, mu=1e-2, verbose=False, maskfree=None):
    nQ = S.Cc.shape[0]
    P = theta.numel()
    def full_norm(th):
        tot = 0.0
        for a in range(0, nQ, chunk):
            r = S.resid(th, slice(a, a + chunk))
            tot += float((r * r).sum())
        return tot
    f = full_norm(theta)
    for it in range(iters):
        JTJ = torch.zeros(P, P, dtype=torch.float64); JTr = torch.zeros(P, dtype=torch.float64)
        for a in range(0, nQ, chunk):
            sl = slice(a, a + chunk)
            r = S.resid(theta, sl)
            J = jacfwd(lambda th: S.resid(th, sl))(theta)
            JTJ += J.T @ J; JTr += J.T @ r
        dg = torch.diagonal(JTJ).clone()
        while True:
            Am = JTJ + mu * (torch.diag(dg) + 1e-9 * torch.eye(P, dtype=torch.float64))
            try:
                step = torch.linalg.solve(Am, -JTr)
            except Exception:
                mu *= 10; continue
            th2 = theta + step
            f2 = full_norm(th2)
            if f2 < f:
                theta, f = th2, f2; mu = max(mu / 3, 1e-12); break
            mu *= 4
            if mu > 1e12:
                return theta, f
        if verbose:
            print("  it", it, "f", f, "mu", mu, flush=True)
        if f < 1e-24:
            break
    return theta, f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--n", type=int, default=1)
    ap.add_argument("--power", type=int, default=1)
    ap.add_argument("--W", type=int, required=True)
    ap.add_argument("--dirs", required=True)
    ap.add_argument("--chars", type=int, default=2)
    ap.add_argument("--restarts", type=int, default=4)
    ap.add_argument("--iters", type=int, default=150)
    ap.add_argument("--scale", type=float, default=0.7)
    ap.add_argument("--gscale", type=float, default=-1.0, help="init scale of the G block (default: --scale)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--chunk", type=int, default=64)
    ap.add_argument("--save", default="")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--orbits", type=int, default=0, help="use a random subset of this many character orbits (0 = all)")
    a = ap.parse_args()
    torch.set_num_threads(int(__import__("os").environ.get("OMP_NUM_THREADS", "1")))
    m, n = a.m, a.n
    sig = tf.singer_table(m, a.power) if n > 1 else list(range(1 << m))
    if n > 1:
        assert tf.order_of(sig) == n, (tf.order_of(sig), n)
    dirs = [int(d) for d in a.dirs.split(",")]
    M = len(dirs)
    Cc, TAU, WT = tf.build_char_data(m, n, sig, dirs, a.chars)
    Sfull = Sys(a.W, M, n, torch.tensor(Cc), torch.tensor(TAU), torch.tensor(WT))
    if a.orbits:
        rng = np.random.default_rng(a.seed + 1000)
        keep = np.sort(rng.choice(len(WT), size=min(a.orbits, len(WT)), replace=False))
        keep = np.union1d(keep, [0])          # always keep the trivial character's orbit
        Cc, TAU, WT = Cc[keep], TAU[keep], WT[keep]
    S = Sys(a.W, M, n, torch.tensor(Cc), torch.tensor(TAU), torch.tensor(WT))
    print(json.dumps(dict(m=m, n=n, W=a.W, M=M, chars=a.chars, orbits=len(WT), unknowns=nparams(a.W, M, n))), flush=True)
    best = None
    for r in range(a.restarts):
        g = torch.Generator().manual_seed(a.seed + r)
        theta = torch.randn(nparams(a.W, M, n), generator=g, dtype=torch.float64) * a.scale
        if a.gscale >= 0:
            ng = 2 * (2 * a.W * M * n)
            theta[ng:] *= a.gscale / a.scale
        t0 = time.time()
        theta, f = lm_solve(S, theta, iters=a.iters, chunk=a.chunk, verbose=a.verbose)
        ffull = sum(float((Sfull.resid(theta, slice(q, q + a.chunk)) ** 2).sum()) for q in range(0, Sfull.Cc.shape[0], a.chunk))
        print(json.dumps(dict(restart=r, loss=f, loss_all_chars=ffull, secs=round(time.time() - t0, 1))), flush=True)
        if best is None or f < best[0]:
            best = (f, theta)
    print("BEST", best[0], flush=True)
    if a.save:
        torch.save(dict(theta=best[1], dirs=dirs, m=m, n=n, power=a.power, W=a.W, chars=a.chars, loss=best[0]), a.save)


if __name__ == "__main__":
    main()
