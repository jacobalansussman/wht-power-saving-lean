#!/usr/bin/env python3
"""v_kant.py (skeptic v-e3-engineE): is there an EXACT solution of the two-sided component's equations near a saved
numerical fake?  The equations F(z) = 0 (z = tau, b, gamma, a real; F = distance-class sums of the outputs minus the
targets, (m+1) W^2 real equations) are polynomial.  At the saved point: residual eta0 = |F|, Jacobian J (complex-step,
exact to rounding) of full row rank with smallest singular value s_min, and a local Lipschitz estimate Lj of J.
Newton-Kantorovich for underdetermined systems (Newton with the pseudo-inverse): if h = Lj * |J^+| * |J^+ F| < 1/2, an
exact real zero exists within 2 |J^+ F| of z.  The Lipschitz constant is an estimate from sampled points (not interval
arithmetic), so this is certified numerics in spirit, not a proof.  Own frame module (v_frames.py).
usage: v_kant.py fake.npz"""
import sys
import numpy as np
from v_frames import main as frames

d = np.load(sys.argv[1], allow_pickle=True)
m, W, M = int(d["m"]), int(d["W"]), int(d["M"])
R = frames(m)
n, i0, pt, pK, dist = R["n"], R["i0"], R["pt"], R["pK"], R["dist"]
plist = [np.asarray(p, int) for p in d["plist"]]; nb = len(plist)
pairs = [tuple(int(v) for v in x) for x in d["pairs"]]
Ind = np.array([(dist == q).astype(float) for q in range(m + 1)])
def pack(tau, b, gam, a): return np.concatenate([tau.ravel(), b.ravel(), gam.ravel(), a.ravel()])
z0 = pack(d["tau"], d["b"], d["gam"], d["a"])
def unpack(z):
    k = 0
    tau = z[k:k + W * W].reshape(W, W); k += W * W
    b = z[k:k + M * W].reshape(M, W); k += M * W
    gam = z[k:k + len(pairs) * nb].reshape(len(pairs), nb); k += len(pairs) * nb
    a = z[k:k + W * M].reshape(W, M)
    return tau, b, gam, a
def Mp(p, v): out = np.zeros_like(v); out[p] = v; return out
def F(z):
    tau, b, gam, a = unpack(z)
    gi = {pr: gam[q] for q, pr in enumerate(pairs)}
    w = []
    for i in range(M):
        u = np.zeros((n, W), dtype=z.dtype); u[i0] = b[i]
        for ip in range(i): u = u + sum(c * Mp(p, w[ip]) for c, p in zip(gi[(ip, i)], plist))
        w.append(Mp(pt, u) - u)
    r = []
    for l in range(W):
        o = np.zeros((n, W), dtype=z.dtype); o[i0] = tau[l]
        for i in range(M): o = o + a[l, i] * w[i]
        tgt = np.zeros((m + 1, W)); tgt[m, l] = 1
        r.append((Ind @ o - tgt).ravel())
    return np.concatenate(r)
def J(z):
    h = 1e-30; cols = []
    for k in range(len(z)):
        zc = z.astype(complex); zc[k] += 1j * h
        cols.append(F(zc).imag / h)
    return np.array(cols).T
F0 = F(z0); J0 = J(z0)
s = np.linalg.svd(J0, compute_uv=False)
Jp = np.linalg.pinv(J0)
step = np.linalg.norm(Jp @ F0)
rng = np.random.default_rng(5)
Lj = 0.0
for r in (1e-3, 1e-4):
    for _ in range(4):
        dz = rng.normal(size=len(z0)); dz *= r / np.linalg.norm(dz)
        Lj = max(Lj, np.linalg.norm(J(z0 + dz) - J0, 2) / r)
hk = Lj * np.linalg.norm(Jp, 2) * step
print("%s: m=%d W=%d M=%d: %d equations, %d unknowns; |F| = %.2e; Jacobian rank %d (s_min %.3e, s_max %.3e); "
      "Newton step %.2e; sampled Lipschitz of J %.2e; Kantorovich h = %.2e -> %s"
      % (sys.argv[1], m, W, M, len(F0), len(z0), np.linalg.norm(F0), int((s > 1e-10 * s[0]).sum()), s[-1], s[0], step, Lj, hk,
         "exact real zero within %.1e (h < 1/2)" % (2 * step) if hk < 0.5 else "no conclusion"))
