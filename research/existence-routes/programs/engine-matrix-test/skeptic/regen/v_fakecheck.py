#!/usr/bin/env python3
"""v_fakecheck.py (skeptic v-e3-engineE): re-check e3-engineE's saved realised Hecke fakes with my own frame module
(built from the actual letter matrices in v_frames.py; nothing of e3-engineE's code is imported; only the saved
numbers in the .npz are read).
 (a) every saved free element p (plist) is the frame permutation of an element of P = image of the actual free letters
     (m <= 3: all of P enumerated by BFS and compared);
 (b) the two-sided condition as an OPERATOR identity on the whole frame module:  E M(Z)_{l l'} E == delta_{ll'} E M(K) E,
     with M(K) the frame permutation of the actual K matrix and E the orbit-averaging projector; evaluated also on
     fresh random complex vectors of the frame module;
 (c) the one-sided condition: delta_{K L0} (x) e_l in V_M = V_0 + sum_i C[P] w_i, with C[P] = span of ALL of P;
     singular values printed so the numerical rank can be judged.
usage: v_fakecheck.py fake.npz"""
import sys
import numpy as np
from v_frames import main as frames

d = np.load(sys.argv[1], allow_pickle=True)
mode, m, W, M = str(d["mode"]), int(d["m"]), int(d["W"]), int(d["M"])
assert mode == "hecke"
R = frames(m)
n, i0, pt, pK, dist, gens = R["n"], R["i0"], R["pt"], R["pK"], R["dist"], R["gens"]
# (a) enumerate P (as frame permutations)
seen = {tuple(range(n))}; todo = [np.arange(n)]
while todo:
    p = todo.pop()
    for g in gens:
        q = g[p]; k = tuple(q)
        if k not in seen: seen.add(k); todo.append(q)
print("|P| as frame permutations: %d" % len(seen))
plist = d["plist"]
inP = [tuple(int(v) for v in p) in seen for p in plist]
print("(a) saved free elements: %d, all in P: %s" % (len(plist), all(inP)))
# (b) operators
def Mperm(p):
    A = np.zeros((n, n)); A[np.asarray(p), np.arange(n)] = 1; return A
Mt, MK = Mperm(pt), Mperm(pK)
Mp = [Mperm(p) for p in plist]
pairs = [tuple(int(v) for v in x) for x in d["pairs"]]; gam = {pr: d["gam"][q] for q, pr in enumerate(pairs)}
tau, b, a = d["tau"], d["b"], d["a"]
print("   coefficient sizes: max|tau| %.2e max|b| %.2e max|a| %.2e max|gam| %.2e"
      % (np.abs(tau).max(), np.abs(b).max(), np.abs(a).max(), max(np.abs(g).max() for g in gam.values())))
I = np.eye(n)
Wop = {}  # Wop[(i, l')] operator n x n: w_i for input slot l'
for lp in range(W):
    for i in range(M):
        U = b[i, lp] * I
        for ip in range(i):
            G = sum(c * P for c, P in zip(gam[(ip, i)], Mp))
            U = U + G @ Wop[(ip, lp)]
        Wop[(i, lp)] = (Mt - I) @ U
Ind = np.array([(dist == q).astype(float) for q in range(m + 1)]).T
E = Ind @ np.diag(1.0 / Ind.sum(0)) @ Ind.T
EKE = E @ MK @ E
worst = 0.0; termsize = 0.0
rng = np.random.default_rng(20261015)
V = rng.normal(size=(n, 5)) + 1j * rng.normal(size=(n, 5))
worstv = 0.0
for l in range(W):
    for lp in range(W):
        Z = tau[l, lp] * I + sum(a[l, i] * Wop[(i, lp)] for i in range(M))
        termsize = max(termsize, max(np.abs(a[l, i] * Wop[(i, lp)]).max() for i in range(M)))
        D = E @ Z @ E - (EKE if l == lp else 0)
        worst = max(worst, np.abs(D).max())
        worstv = max(worstv, np.abs(D @ V).max())
print("(b) two-sided: max |E M(Z) E - E M(K) E (x) I| over the whole frame module = %.2e; on 5 fresh random complex vectors %.2e;"
      " largest single term %.2e  -> %s" % (worst, worstv, termsize, "HOLDS" if worst < 1e-9 else "FAILS"))
# sanity: the same word with K replaced by a wrong target must fail (power of the check)
wrong = np.abs(E @ (tau[0, 0] * I + sum(a[0, i] * Wop[(i, 0)] for i in range(M))) @ E - E @ Mt @ E).max()
print("    control: distance to the wrong target E M(t) E = %.2e (must be large)" % wrong)
# (c) one-sided: vectors w_i(delta_L0) as n x W matrices, flattened (block l' = input slot l')
delta = np.zeros(n); delta[i0] = 1
wvec = [np.concatenate([Wop[(i, lp)] @ delta for lp in range(W)]) for i in range(M)]
Pall = [np.array(p) for p in seen]
cols = []
for l in range(W):
    v = np.zeros(n * W); v[l * n + i0] = 1; cols.append(v)
for wv in wvec:
    blocks = wv.reshape(W, n)
    for p in Pall:
        out = np.zeros_like(blocks); out[:, p] = blocks   # M(p) on each block
        cols.append(out.ravel())
A = np.array(cols).T
U_, s, _ = np.linalg.svd(A, full_matrices=False)
print("(c) one-sided: spanning set %d vectors in dim %d; singular values (largest 3, around the gap, smallest 3):"
      % (A.shape[1], n * W))
print("    top", np.round(s[:3], 3), " ...", ["%.2e" % x for x in s[n * W - 8:n * W]])
for tol in (1e-6, 1e-8, 1e-10):
    r = int((s > tol * s[0]).sum()); Q = U_[:, :r]
    T = np.zeros((n * W, W))
    for l in range(W): T[l * n + int(pK[i0]), l] = 1
    print("    rel tol %.0e: rank %d of %d; distance of the targets delta_{K L0} (x) e_l to V_M = %.2e"
          % (tol, r, n * W, np.linalg.norm(T - Q @ (Q.T @ T))))
