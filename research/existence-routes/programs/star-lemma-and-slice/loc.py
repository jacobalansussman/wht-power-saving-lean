#!/usr/bin/env python3
"""loc.py -- the local rank problem on AG(n,2) (e3-slice6, 15 Oct 2026).

Question: is there a NONSINGULAR N x N matrix E (N = 2^n, rows/cols = points of AG(n,2)) such that for every
affine hyperplane H of AG(n,2) the principal block E[H,H] has rank <= s ?

Why: for U_m (one move per direction) with W = 2^(m-2) - s arrays, Theorem 10 makes E[F,F] nonsingular on the even
(m-2)-flat F = {z1 = z2 = 1} (= AG(m-2,2)), and every affine hyperplane of F lies in an odd (m-2)-flat avoiding 0
(Lemma 12'), whose E-block has rank = slack = s.  So "no such E" would EXCLUDE U_m at that W.  (m = 6: n = 4,
s = 3 is 13 in 63.)

Method: Levenberg-Marquardt on the polynomial system
    E[H,H] - U_H V_H^T = 0  (every H; U_H, V_H are 2^(n-1) x s),   E G - I = 0  (G is E's inverse),
real or complex (complex = real pairs).  Residual -> 0 with bounded G = numerical evidence of feasibility.
Usage: loc.py n s seed iters [complex] [extra]
"""
import sys, time
import numpy as np

def hyperplanes(n):
    N = 1 << n
    pts = np.arange(N)
    out = []
    for h in range(1, N):
        par = np.array([bin(h & x).count("1") & 1 for x in pts])
        for b in (0, 1):
            out.append(np.nonzero(par == b)[0])
    return out

def run(n, s, seed, iters, cplx=False, verbose=True):
    rng = np.random.default_rng(seed)
    N = 1 << n
    Hs = hyperplanes(n)
    h = N // 2
    nH = len(Hs)
    dt = complex if cplx else float
    def rnd(*sh):
        x = rng.standard_normal(sh)
        if cplx:
            x = x + 1j * rng.standard_normal(sh)
        return x
    E = rnd(N, N)
    G = np.linalg.inv(E)
    U = [rnd(h, s) for _ in Hs]
    V = [rnd(h, s) for _ in Hs]
    # variable layout: E (N*N), G (N*N), then per H: U (h*s), V (h*s)
    nE = N * N
    offU = 2 * nE
    per = 2 * h * s
    nv = offU + nH * per
    neq = nH * h * h + nE + 1

    def pack():
        return np.concatenate([E.ravel(), G.ravel()] + [np.concatenate([U[k].ravel(), V[k].ravel()]) for k in range(nH)])

    def unpack(x):
        E_ = x[:nE].reshape(N, N)
        G_ = x[nE:2 * nE].reshape(N, N)
        U_, V_ = [], []
        for k in range(nH):
            o = offU + k * per
            U_.append(x[o:o + h * s].reshape(h, s))
            V_.append(x[o + h * s:o + per].reshape(h, s))
        return E_, G_, U_, V_

    def resid(x):
        E_, G_, U_, V_ = unpack(x)
        r = [ (E_[np.ix_(H, H)] - U_[k] @ V_[k].T).ravel() for k, H in enumerate(Hs) ]
        r.append((E_ @ G_ - np.eye(N)).ravel())
        r.append(np.array([np.sum(E_ * E_) - N]))   # fixes the scale of E (rank conditions are homogeneous)
        return np.concatenate(r)

    def jac(x):
        E_, G_, U_, V_ = unpack(x)
        J = np.zeros((neq, nv), dtype=dt)
        row = 0
        for k, H in enumerate(Hs):
            o = offU + k * per
            for a in range(h):
                for b in range(h):
                    J[row, H[a] * N + H[b]] = 1.0
                    # - sum_t U[a,t] V[b,t]
                    J[row, o + a * s: o + a * s + s] = -V_[k][b, :]
                    J[row, o + h * s + b * s: o + h * s + b * s + s] = -U_[k][a, :]
                    row += 1
        # E G - I : d/dE_{i,l} -> G_{l,j}; d/dG_{l,j} -> E_{i,l}
        for i in range(N):
            for j in range(N):
                J[row, i * N:(i + 1) * N] = G_[:, j]
                J[row, nE + np.arange(N) * N + j] = E_[i, :]
                row += 1
        J[row, :nE] = 2 * E_.ravel()
        return J

    x = pack()
    lam = 1e-2
    r = resid(x)
    f = np.vdot(r, r).real
    t0 = time.time()
    for it in range(iters):
        J = jac(x)
        if cplx:
            # treat as real problem in (Re, Im): use complex normal equations (holomorphic residual)
            JH = J.conj().T
        else:
            JH = J.T
        A = JH @ J
        g = JH @ r
        while True:
            step = np.linalg.solve(A + lam * np.diag(np.diag(A).real + 1e-9), -g)
            xn = x + step
            rn = resid(xn)
            fn = np.vdot(rn, rn).real
            if fn < f:
                x, r, f = xn, rn, fn
                lam = max(lam / 3, 1e-12)
                break
            lam *= 4
            if lam > 1e12:
                break
        if verbose and (it % 10 == 0 or f < 1e-24):
            E_, G_, _, _ = unpack(x)
            sv = np.linalg.svd(E_, compute_uv=False)
            print(f"it {it} res2 {f:.3e} lam {lam:.1e} sminE {sv[-1]:.3e} smaxE {sv[0]:.3e} t {time.time()-t0:.1f}s", flush=True)
        if f < 1e-26 or lam > 1e12:
            break
    E_, G_, U_, V_ = unpack(x)
    return f, E_, Hs

if __name__ == "__main__":
    n, s, seed, iters = map(int, sys.argv[1:5])
    cplx = len(sys.argv) > 5 and sys.argv[5] == "complex"
    f, E, Hs = run(n, s, seed, iters, cplx)
    sv = np.linalg.svd(E, compute_uv=False)
    worst = max(np.linalg.svd(E[np.ix_(H, H)], compute_uv=False)[s] if s < len(H) else 0.0 for H in Hs)
    print(f"FINAL n={n} s={s} seed={seed} cplx={cplx} res2={f:.3e} sigma_min(E)={sv[-1]:.3e} cond={sv[0]/sv[-1]:.3e} "
          f"max_H sigma_(s+1)(E_H)={worst:.3e}")
    np.save(f"loc_E_n{n}_s{s}_seed{seed}{'_c' if cplx else ''}.npy", E)
