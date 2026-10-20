#!/usr/bin/env python3
"""my_cert_check.py (v-e3-rmbuild): print and re-verify an exclusion certificate for 4 in 15 (one move per direction)
in the no-phase game with the cyclic address group <h>: a set Y of virtual moves (t,k), and for every pair of Y a
quadratic form Q (given by its coefficients q_ij) with Q tight (|D_Q| = W*nodd(Q)) and both members in D_Q.
|Y| > W*n then contradicts Lemma V(c).  Usage: my_cert_check.py cols   (h e_1..h e_4)."""
import sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from my_vtests import max_clique
m = 4; N = 16; W = 4
cols = [int(t) for t in sys.argv[1].split(",")]
def lin(cols):
    img = []
    for x in range(N):
        y = 0
        for k in range(m):
            if x >> k & 1: y ^= cols[k]
        img.append(y)
    return img
g = lin(cols); H = [list(range(N))]
while True:
    nx = [g[v] for v in H[-1]]
    if nx == list(range(N)): break
    H.append(nx)
n = len(H); Hinv = []
for h in H:
    inv = [0] * N
    for x in range(N): inv[h[x]] = x
    Hinv.append(inv)
fixed = [x for x in range(1, N) if g[x] == x]
Z = list(range(1, N)); M = 15
pairs = [(i, j) for i in range(m) for j in range(i, m)]
def Qtab(bits):
    q = [p for b, p in enumerate(pairs) if bits >> b & 1]
    return [sum((z >> i & 1) & (z >> j & 1) for i, j in q) & 1 for z in range(N)]
tight = {}
for bits in range(1 << len(pairs)):
    Q = Qtab(bits)
    D = frozenset((t, k) for t in range(M) for k in range(n) if Q[Hinv[k][Z[t]]])
    nodd = sum(1 for k in range(n) if sum(Q[Hinv[k][1 << j]] for j in range(m)) & 1)
    if D and len(D) == W * nodd:
        tight[bits] = D
adj = {}
for D in tight.values():
    for v in D: adj.setdefault(v, set()).update(D - {v})
Y = max_clique(adj, list(adj))
# re-verify every pair with an explicit witness form
ok = True; witnesses = {}
for a in range(len(Y)):
    for b in range(a + 1, len(Y)):
        w = [bits for bits, D in tight.items() if Y[a] in D and Y[b] in D]
        if not w: ok = False
        else: witnesses[(a, b)] = w[0]
print("h e_j =", cols, "order", n, "fixed nonzero vectors", fixed)
print("tight forms (coefficient bitmask over pairs i<=j):", sorted(tight))
print("clique Y (move t = direction t+1, channel k):", sorted(Y))
print("|Y| = %d, bound W*n = %d, every pair has a tight witness: %s, distinct witness forms used: %d"
      % (len(Y), W * n, ok, len(set(witnesses.values()))))
