#!/usr/bin/env python3
"""sec.py (e3-rmbuild): one SECTOR block of the transfer-form system, Levenberg-Marquardt.

Shifts X^a are free and commute with every slice letter, so the central idempotents
    e_0 (average of all shifts)  and  1 - e_0
of C[A x| H] (H any group of address maps) split every program into a sector-0 part and a nonzero-sector part with
independent coefficients (PROOF.md, Lemma S).  In the commutative slice (H = 1) every single sector lambda is its own
block, with SCALAR coefficients, and its characters are chi_{lambda,Q}(dir w) = i^{[lambda.w]} (-1)^{Q(w)}, Q a binary
quadratic form.  With a cyclic address group Z_n acting FREELY on the orbit of lambda, the block's coefficient algebra
is the full matrix algebra M_n(C) (shift diagonals + the cyclic shift generate it), the representation is
Ind chi_{lambda,Q} on C^n with channel k seeing chi(dir(sigma^-k w)).
In a nonzero sector no character is trivial on every move, so the free part A0 of the transfer form is an unknown too:
    T = A0 + Atil C (I - G C)^{-1} Btil  =  diag_k(chi(K at channel k)) (x) I_W.
--lam 0 gives sector 0 with scalar coefficients (= tf.py/lm.py at n = 1 with chars 2; A0 then also free, but the
trivial character forces A0 = I).
"""
import argparse, json, time
import numpy as np
import torch
from torch.func import jacfwd
import tf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--n", type=int, default=1)
    ap.add_argument("--power", type=int, default=1)
    ap.add_argument("--lam", type=int, required=True)
    ap.add_argument("--W", type=int, required=True)
    ap.add_argument("--dirs", required=True)
    ap.add_argument("--restarts", type=int, default=3)
    ap.add_argument("--iters", type=int, default=80)
    ap.add_argument("--scale", type=float, default=0.7)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--chunk", type=int, default=128)
    ap.add_argument("--save", default="")
    a = ap.parse_args()
    torch.set_num_threads(1)
    m, n, W = a.m, a.n, a.W
    N = 1 << m
    sig = tf.singer_table(m, a.power) if n > 1 else list(range(N))
    inv = [0] * N
    for v in range(N):
        inv[sig[v]] = v
    sk = [list(range(N))]
    for k in range(1, n):
        sk.append([inv[v] for v in sk[-1]])
    dirs = [int(d) for d in a.dirs.split(",")]
    M = len(dirs)
    forms = tf.binary_forms(m)          # entries 0/2 (exponent of i)
    I4 = np.array([1, 1j, -1, -1j])
    lam = a.lam
    Cl = []; TL = []
    for t in forms:
        def chi_exp(w):
            return ((tf.popc(lam & w) & 1) + t[w]) % 4
        c = np.array([I4[chi_exp(sk[k][z])] - 1 for z in dirs for k in range(n)])
        tau = np.array([I4[sum(chi_exp(sk[k][1 << j]) for j in range(m)) % 4] for k in range(n)])
        Cl.append(c); TL.append(tau)
    Cc = torch.tensor(np.array(Cl)); TAU = torch.tensor(np.array(TL))
    nQ = Cc.shape[0]
    Wn, Mn = W * n, M * n
    ii, jj = torch.tril_indices(M, M, -1)
    shapes = [("A", (Wn, Mn)), ("B", (Mn, Wn)), ("G", (len(ii), n, n)), ("A0", (Wn, Wn))]
    sizes = [int(np.prod(s)) for _, s in shapes]
    P = 2 * sum(sizes)
    print(json.dumps(dict(m=m, n=n, lam=lam, W=W, M=M, unknowns=P, chars=nQ)), flush=True)
    targ_all = torch.einsum("ab,qij->qaibj", torch.eye(W, dtype=torch.complex128),
                            torch.diag_embed(TAU)).reshape(nQ, Wn, Wn)

    def unpack(th):
        out = {}; k = 0
        for (nm, sh), sz in zip(shapes, sizes):
            out[nm] = torch.complex(th[k:k + sz], th[k + sz:k + 2 * sz]).reshape(sh); k += 2 * sz
        Gb = torch.zeros(M, M, n, n, dtype=torch.complex128).index_put((ii, jj), out["G"])
        out["Gbig"] = Gb.permute(0, 2, 1, 3).reshape(Mn, Mn)
        return out

    def resid(th, sl):
        p = unpack(th)
        C = Cc[sl]
        L = torch.eye(Mn, dtype=torch.complex128)[None] - p["Gbig"][None] * C[:, None, :]
        Y = torch.linalg.solve(L, p["B"][None].expand(C.shape[0], -1, -1))
        X = p["A0"][None] + p["A"][None] @ (C[:, :, None] * Y)
        R = (X - targ_all[sl]).reshape(-1)
        return torch.cat([R.real, R.imag])

    def fnorm(th):
        return sum(float((resid(th, slice(s, s + a.chunk)) ** 2).sum()) for s in range(0, nQ, a.chunk))

    best = None
    for r in range(a.restarts):
        g = torch.Generator().manual_seed(a.seed + r)
        th = torch.randn(P, generator=g, dtype=torch.float64) * a.scale
        f = fnorm(th); mu = 1e-2; t0 = time.time()
        for it in range(a.iters):
            JTJ = torch.zeros(P, P, dtype=torch.float64); JTr = torch.zeros(P, dtype=torch.float64)
            for s in range(0, nQ, a.chunk):
                sl = slice(s, s + a.chunk)
                rr = resid(th, sl); J = jacfwd(lambda x: resid(x, sl))(th)
                JTJ += J.T @ J; JTr += J.T @ rr
            dg = torch.diagonal(JTJ).clone()
            ok = False
            while mu < 1e12:
                step = torch.linalg.solve(JTJ + mu * (torch.diag(dg) + 1e-9 * torch.eye(P, dtype=torch.float64)), -JTr)
                f2 = fnorm(th + step)
                if f2 < f:
                    th, f, mu, ok = th + step, f2, max(mu / 3, 1e-12), True
                    break
                mu *= 4
            if not ok or f < 1e-24:
                break
        print(json.dumps(dict(restart=r, loss=f, iters=it + 1, secs=round(time.time() - t0, 1))), flush=True)
        if best is None or f < best[0]:
            best = (f, th)
    print("BEST", best[0], flush=True)
    if a.save:
        torch.save(dict(theta=best[1], dirs=dirs, m=m, n=n, power=a.power, lam=lam, W=W, loss=best[0]), a.save)


if __name__ == "__main__":
    main()
