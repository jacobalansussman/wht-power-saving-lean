#!/usr/bin/env python3
"""my_vtests.py (v-e3-rmbuild, skeptic): independent re-implementation of Lemma V's deletion and tight-clique tests
for the no-phase game with a cyclic address group H = <sigma> (written from PROOF.md's statement only, not from
vclique.py / vtest5.py).

Virtual move (t,k): physical move t in channel k, direction sigma^-k z_t.  For a binary quadratic form Q:
  D_Q   = {(t,k) : Q(sigma^-k z_t) = 1}
  nodd  = #{k : sum_j Q(sigma^-k e_j) = 1}
  slack = |D_Q| - W*nodd   (deletion: slack >= 0)
  tight = slack == 0 and D_Q nonempty; clique = set of virtual moves pairwise in a common tight D_Q; bound Wn.
Usage: my_vtests.py m n_power poly W dirs   (n_power 0 = slice; dirs 'all' = one move per nonzero direction,
       'e1xM' = M moves all in direction e_1; poly = primitive polynomial as int)
"""
import sys, itertools
from math import comb


def mulx(v, m, poly):
    v <<= 1
    if v >> m & 1:
        v ^= poly
    return v


def singer(m, poly, power):
    tab = list(range(1 << m))
    for _ in range(power):
        tab = [mulx(v, m, poly) for v in tab]
    return tab


def order(tab):
    N = len(tab); cur = list(range(N)); k = 0
    while True:
        cur = [tab[v] for v in cur]; k += 1
        if cur == list(range(N)):
            return k


def forms(m):
    pairs = [(i, j) for i in range(m) for j in range(i, m)]
    N = 1 << m
    for bits in range(1 << len(pairs)):
        q = [p for b, p in enumerate(pairs) if bits >> b & 1]
        yield bits, [sum((z >> i & 1) & (z >> j & 1) for i, j in q) & 1 for z in range(N)]


def max_clique(adj, verts):
    best = []
    def bk(R, P, X):
        nonlocal best
        if not P and not X:
            if len(R) > len(best):
                best = R[:]
            return
        if len(R) + len(P) <= len(best):
            return
        u = max(P | X, key=lambda w: len(adj[w] & P))
        for v in list(P - adj[u]):
            bk(R + [v], P & adj[v], X & adj[v])
            P = P - {v}; X = X | {v}
    bk([], set(verts), set())
    return best


def run(m, power, poly, W, dirs):
    N = 1 << m
    if power == 0:
        H = [list(range(N))]
    else:
        s = singer(m, poly, power); n = order(s)
        H = [list(range(N))]
        for _ in range(n - 1):
            H.append([s[v] for v in H[-1]])
    n = len(H)
    Hinv = []
    for h in H:
        inv = [0] * N
        for x in range(N):
            inv[h[x]] = x
        Hinv.append(inv)
    if dirs == "all":
        Z = list(range(1, N))
    else:
        Z = [1] * int(dirs.split("x")[1])
    M = len(Z)
    nfail = 0; tight = []; minslack = None; fails = []
    for bits, Q in forms(m):
        D = [(t, k) for t in range(M) for k in range(n) if Q[Hinv[k][Z[t]]]]
        nodd = sum(1 for k in range(n) if sum(Q[Hinv[k][1 << j]] for j in range(m)) & 1)
        sl = len(D) - W * nodd
        minslack = sl if minslack is None else min(minslack, sl)
        if sl < 0:
            nfail += 1; fails.append((bits, len(D), W * nodd))
        elif sl == 0 and D:
            tight.append(D)
    adj = {}
    for D in tight:
        S = set(D)
        for v in D:
            adj.setdefault(v, set()).update(S - {v})
    cl = max_clique(adj, list(adj)) if adj else []
    print("m=%d n=%d (power %d, poly %s) W=%d M=%d: deletion failures %d %s; min slack %d; tight forms %d; "
          "max tight clique %d; bound Wn = %d -> %s"
          % (m, n, power, bin(poly), W, M, nfail, fails[:4], minslack, len(tight), len(cl), W * n,
             "EXCLUDED (deletion)" if nfail else ("EXCLUDED (clique)" if len(cl) > W * n else "not excluded")))
    sys.stdout.flush()


if __name__ == "__main__":
    m, power, poly, W = (int(x, 0) for x in sys.argv[1:5])
    run(m, power, poly, W, sys.argv[5])
