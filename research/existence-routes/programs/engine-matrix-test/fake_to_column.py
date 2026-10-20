#!/usr/bin/env python3
"""fake_to_column.py: take a saved realised Hecke fake, rebuild its move outputs w_i, and test the one-sided (column)
condition with the FULL free algebra: is delta_{K L0} (x) e_l in V_M = V_0 + sum_i C[P] w_i for every l?
If yes, outputs o_l = (element of V_M) exist, i.e. the same moves and links give a word passing the column component.
usage: fake_to_column.py fake.npz"""
import sys
import numpy as np
from sp import Sp
from modgame import build, closure, orth
d = np.load(sys.argv[1], allow_pickle=True)
m, W, M = int(d["m"]), int(d["W"]), int(d["M"])
n, i0, iK, t, gens = build(m)
S = Sp(m); L0, F = S.frames(); idx = {L: i for i, L in enumerate(F)}
plist = d["plist"]; pairs = [tuple(x) for x in d["pairs"]]; gam = {pr: d["gam"][q] for q, pr in enumerate(pairs)}
b = d["b"]
def Mp(p, v): out = np.zeros_like(v); out[p] = v; return out
def lin(c, v): return sum(ck * Mp(p, v) for ck, p in zip(c, plist))
tinv = np.argsort(t)
w = []
for i in range(M):
    u = np.zeros((n, W)); u[i0] = b[i]
    for ip in range(i): u = u + lin(gam[(ip, i)], w[ip])
    w.append(u[tinv] - u)
V0 = np.zeros((n * W, W))
for l in range(W): V0[l * n + i0, l] = 1
cols = [V0] + [wi.T.reshape(-1)[:, None] for wi in w]   # stack the W blocks: block l = column l of w_i
T = np.zeros((n * W, W))
for l in range(W): T[l * n + iK, l] = 1
for tol in (1e-7, 1e-9, 1e-11):
    Vt = closure(np.hstack(cols), gens, n, W, tol)
    print("   tol %.0e: dim V_M = %d, distance of targets = %.2e" % (tol, Vt.shape[1], np.linalg.norm(T - Vt @ (Vt.T @ T))))
V = closure(np.hstack(cols), gens, n, W)
res = np.linalg.norm(T - V @ (V.T @ T))
# P-averaged (trivial isotypic) part, as a check of the Hecke condition
print("%s: m=%d W=%d M=%d: dim V_M = %d of %d; distance of the W targets to V_M = %.2e -> %s"
      % (sys.argv[1], m, W, M, V.shape[1], n * W, res, "COLUMN CONDITION HOLDS (column fake)" if res < 1e-8 else "column condition fails for these links"))
