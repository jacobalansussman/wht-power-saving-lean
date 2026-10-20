#!/usr/bin/env python3
"""check_fake.py: independent re-check of a saved frame-module fake.  It does NOT reuse framefake.py's evaluation:
it rebuilds the frames and the free elements from sp.py, writes the fake out as an explicit register word
(registers hold frame-module vectors (x) inputs; steps: 'lin' = new register := sum_k c_k M(p_k) (register),
'move' = apply M(t), 'sub', 'add'), runs the word on W FRESH random input vectors y (complex) per input slot, and
compares the distance-class sums (hecke) or the full vectors (column) of every output with those of delta_{K L0} y_l.
usage: check_fake.py file.npz"""
import sys
import numpy as np
from sp import Sp
d = np.load(sys.argv[1], allow_pickle=True)
mode, m, W, M = str(d["mode"]), int(d["m"]), int(d["W"]), int(d["M"])
S = Sp(m); L0, F = S.frames(); idx = {L: i for i, L in enumerate(F)}; n = len(F)
dist = np.array([S.dist(L, L0) for L in F])
t = np.array([idx[frozenset(S.apply(S.transvection(1), x) for x in L)] for L in F])
Kb = [((1 << j) if j < m else ((1 << j) ^ (1 << (j - m)))) for j in range(S.n)]
iK = idx[frozenset(S.apply(Kb, x) for x in L0)]; i0 = idx[L0]
plist = d["plist"]
# the free elements must fix L0 (they are elements of P) and preserve distance classes
assert all(p[i0] == i0 for p in plist) and all((dist[p] == dist).all() for p in plist)
def Mp(p, v):  # M(p) v : (M(p) v)[p[L]] = v[L]
    out = np.zeros_like(v); out[p] = v; return out
def lin(c, v): return sum(ck * Mp(p, v) for ck, p in zip(c, plist))
rng = np.random.default_rng(12345)
Y = rng.normal(size=(W, 3)) + 1j * rng.normal(size=(W, 3))   # 3 fresh random scalar inputs per array slot
pairs = [tuple(x) for x in d["pairs"]]; gam = {pr: d["gam"][q] for q, pr in enumerate(pairs)}
tau, b, a = d["tau"], d["b"], d["a"]
worst = 0.0
for col in range(Y.shape[1]):
    y = Y[:, col]
    inputs = [np.zeros(n, complex) for _ in range(W)]
    for l in range(W): inputs[l][i0] = y[l]                  # register l holds delta_{L0} y_l
    w = []
    for i in range(M):
        u = sum(b[i, l] * inputs[l] for l in range(W))         # free mixing of the input registers
        for ip in range(i): u = u + lin(gam[(ip, i)], w[ip])   # free letters + mixing of stored differences
        tu = u[np.argsort(t)]                                  # the move: M(t) u   (M(t) v = v[t^-1])
        w.append(tu - u)
    for l in range(W):
        o = sum(tau[l, k] * inputs[k] for k in range(W))
        for i in range(M):
            o = o + (a[l, i] * w[i] if mode == "hecke" else lin(a[l, i], w[i]))
        tgt = np.zeros(n, complex); tgt[iK] = y[l]
        if mode == "hecke":
            diff = np.array([o[dist == q].sum() - tgt[dist == q].sum() for q in range(m + 1)])
        else:
            diff = o - tgt
        worst = max(worst, float(np.abs(diff).max()))
print("%s: mode=%s m=%d W=%d M=%d; fresh random inputs (3 per slot): max |output - target| = %.2e  %s"
      % (sys.argv[1], mode, m, W, M, worst, "CONFIRMED" if worst < 1e-8 else "FAILED"))
