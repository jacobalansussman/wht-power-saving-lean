#!/usr/bin/env python3
"""v_vpcheck.py (skeptic v-e3-engineE): check a hecke_vp.py save (z = (b, gamma), plist, pairs) with my own frame
module (v_frames.py, from the actual letter matrices).  The outer couplings (a, tau) are solved here by least squares.
Checks: every p in plist fixes L0 and is an automorphism of the distance-1 graph on frames (so it lies in P by Chow's
theorem, Aut = Sp(2m,2) for q = 2), and the two-sided condition E M(Z) E == E M(K) E (x) I_W evaluated by the word on
every P-invariant input vector (orbit indicators) and on fresh random complex vectors.
usage: v_vpcheck.py save.npz"""
import sys
import numpy as np
from v_frames import main as frames

d = np.load(sys.argv[1], allow_pickle=True)
m, W, M = int(d["m"]), int(d["W"]), int(d["M"])
R = frames(m)
n, i0, pt, pK, dist, F = R["n"], R["i0"], R["pt"], R["pK"], R["dist"], R["F"]
plist = [np.asarray(p, int) for p in d["plist"]]
L0 = F[i0]
# distance-1 graph
Fl = [np.array(sorted(L)) for L in F]
def dd(A, B): return m - (len(A & B).bit_length() - 1)
adj = np.zeros((n, n), bool)
for a in range(n):
    for b in range(a + 1, n):
        if dd(F[a], F[b]) == 1: adj[a, b] = adj[b, a] = True
ok = all(p[i0] == i0 and (adj[np.ix_(p, p)] == adj).all() for p in plist)
print("plist: %d elements; all fix L0 and are automorphisms of the dual polar graph: %s" % (len(plist), ok))
pairs = [tuple(int(v) for v in x) for x in d["pairs"]]; nb = len(plist)
z = d["z"]; b = z[:M * W].reshape(M, W); g = z[M * W:].reshape(len(pairs), nb)
gam = {pr: g[q] for q, pr in enumerate(pairs)}
def Mp(p, v): out = np.zeros_like(v); out[p] = v; return out
def Mt(v): return Mp(pt, v)
def word_w(v):
    """move outputs w_i (n x W arrays, column l' = input slot l') when every input slot holds v."""
    w = []
    for i in range(M):
        u = np.outer(v, b[i])
        for ip in range(i): u = u + sum(c * Mp(p, w[ip]) for c, p in zip(gam[(ip, i)], plist))
        w.append(Mt(u) - u)
    return w
Ind = np.array([(dist == q).astype(float) for q in range(m + 1)]).T
def Eavg(x): return Ind @ ((Ind.T @ x) / Ind.sum(0)[:, None])
# solve (a, tau) from the delta_L0 equations: S_d(o_l) = [d = m] e_l
delta = np.zeros(n); delta[i0] = 1
w0 = word_w(delta)
Sd = Ind.T
Sdel = Sd @ delta
Bm = np.hstack([np.concatenate([Sdel * (1.0 if lp == k else 0.0) for lp in range(W)])[:, None] for k in range(W)]
               + [np.concatenate([Sd @ wi[:, lp] for lp in range(W)])[:, None] for wi in w0])
coef = []
res = 0.0
for l in range(W):
    T = np.concatenate([(np.arange(m + 1) == m) * (1.0 if lp == l else 0.0) for lp in range(W)])
    c, *_ = np.linalg.lstsq(Bm, T, rcond=None); coef.append(c); res = max(res, np.abs(Bm @ c - T).max())
print("outer couplings solved from delta_L0: max residual %.2e; max |coef| %.2e" % (res, max(np.abs(c).max() for c in coef)))
tau = np.array([c[:W] for c in coef]); a = np.array([c[W:] for c in coef])
rng = np.random.default_rng(777)
tests = [Ind[:, q] for q in range(m + 1)] + [Eavg((rng.normal(size=n) + 1j * rng.normal(size=n))[:, None])[:, 0] for _ in range(3)]
worst = 0.0
for v in tests:
    w = word_w(v.astype(complex))
    for l in range(W):
        o = np.outer(v, tau[l]) + sum(a[l, i] * w[i] for i in range(M))      # n x W
        tgt = np.outer(Mp(pK, v), np.eye(W)[l])
        worst = max(worst, np.abs(Eavg(o) - Eavg(tgt)).max())
print("two-sided E M(Z) E v == E M(K) E v (x) e_l on %d P-invariant and %d fresh random vectors: max error %.2e -> %s"
      % (m + 1, 3, worst, "HOLDS" if worst < 1e-9 else "FAILS"))
