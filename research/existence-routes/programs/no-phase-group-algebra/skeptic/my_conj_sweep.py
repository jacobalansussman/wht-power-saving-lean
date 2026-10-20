#!/usr/bin/env python3
"""my_conj_sweep.py (v-e3-rmbuild): the Lemma V tests depend on the pair (H, target I), not on H up to conjugacy
(conjugating H by g changes the target [I] into [g^-1 g^-T]).  Sweep all conjugates g H g^-1 (g in GL(m,2)) of the
Singer powers used by e3-rmbuild and report the spread of verdicts.
Usage: my_conj_sweep.py m power W dirs   (dirs 'all' or 'e1xM'; poly fixed to x^3+x+1 / x^4+x+1)"""
import sys, itertools
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from my_vtests import singer, order, forms, max_clique

m, power, W = (int(x) for x in sys.argv[1:4]); dirs = sys.argv[4]
N = 1 << m
poly = {3: 0b1011, 4: 0b10011}[m]
s = singer(m, poly, power); n = order(s)


def lin(cols):
    img = []
    for x in range(N):
        y = 0
        for k in range(m):
            if x >> k & 1:
                y ^= cols[k]
        img.append(y)
    return img


GL = [lin(c) for c in itertools.product(range(1, N), repeat=m) if len(set(lin(c))) == N]
FORMS = list(forms(m))
seen = {}
for g in GL:
    gi = [0] * N
    for x in range(N):
        gi[g[x]] = x
    gen = [g[s[gi[x]]] for x in range(N)]          # g sigma g^-1
    H = [list(range(N))]
    for _ in range(n - 1):
        H.append([gen[v] for v in H[-1]])
    key = frozenset(tuple(h) for h in H)
    if key in seen:
        continue
    Hinv = []
    for h in H:
        inv = [0] * N
        for x in range(N):
            inv[h[x]] = x
        Hinv.append(inv)
    if dirs == "all":
        # one move per direction; with H not transitive the multiset of orbits is what matters
        Z = list(range(1, N))
    elif dirs == "plain":
        Z = [1 << j for j in range(m)] * W          # the plain schedule: exists in every sub-game
    elif dirs.startswith("list:"):
        Z = [int(t) for t in dirs[5:].split(",")]
    else:
        Z = [1] * int(dirs.split("x")[1])
    M = len(Z)
    nfail = 0; tight = []
    for bits, Q in FORMS:
        D = [(t, k) for t in range(M) for k in range(n) if Q[Hinv[k][Z[t]]]]
        nodd = sum(1 for k in range(n) if sum(Q[Hinv[k][1 << j]] for j in range(m)) & 1)
        sl = len(D) - W * nodd
        if sl < 0:
            nfail += 1
        elif sl == 0 and D:
            tight.append(D)
    adj = {}
    for D in tight:
        S = set(D)
        for v in D:
            adj.setdefault(v, set()).update(S - {v})
    cl = len(max_clique(adj, list(adj))) if adj else 0
    verdict = "EXCL-deletion" if nfail else ("EXCL-clique" if cl > W * n else "not-excluded")
    seen[key] = (nfail, len(tight), cl, verdict)
from collections import Counter
print("m=%d n=%d W=%d dirs=%s: %d distinct conjugate subgroups" % (m, n, W, dirs, len(seen)))
for k, v in sorted(Counter(seen.values()).items(), key=lambda kv: -kv[1]):
    print("  (deletion failures, tight forms, max clique, verdict) =", k, "x", v, " bound Wn =", W * n)
