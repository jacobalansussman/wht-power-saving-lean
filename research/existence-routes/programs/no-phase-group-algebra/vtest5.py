#!/usr/bin/env python3
"""vtest5.py (e3-rmbuild): Lemma V deletion + tight-clique tests, fast version on direction multiplicities
(works at block 5).  Usage: vtest5.py m n power W dirs|uniform:r  (uniform:r = r copies of every nonzero direction)."""
import sys
exec(open("vclique.py").read().split("m, n, power, W =")[0].split('"""', 2)[2])
m, n, power, W = (int(x) for x in sys.argv[1:5])
N = 1 << m
dirs = ([z for z in range(1, N) for _ in range(int(sys.argv[5].split(":")[1]))] if sys.argv[5].startswith("uniform")
        else [int(d) for d in sys.argv[5].split(",")])
sig = tf.singer_table(m, power) if n > 1 else list(range(N))
if n > 1:
    assert tf.order_of(sig) == n
inv = [0] * N
for v in range(N):
    inv[sig[v]] = v
sk = [list(range(N))]
for k in range(1, n):
    sk.append([inv[v] for v in sk[-1]])
c = [0] * N
for z in dirs:
    for k in range(n):
        c[sk[k][z]] += 1
pairs = [(j, k) for j in range(m) for k in range(j, m)]
pz = [[b for b, (j, k) in enumerate(pairs) if (z >> j) & 1 and (z >> k) & 1] for z in range(N)]
zmask = [sum(1 << b for b in pz[z]) for z in range(N)]
basis_k = [[sk[k][1 << j] for j in range(m)] for k in range(n)]
bad = 0; tight = []
for q in range(1 << len(pairs)):
    val = [bin(q & zmask[z]).count("1") & 1 for z in range(N)]
    D = sum(c[z] for z in range(1, N) if val[z])
    nodd = sum(1 for k in range(n) if sum(val[w] for w in basis_k[k]) & 1)
    if D < W * nodd:
        bad += 1
    elif D == W * nodd and nodd:
        tight.append(frozenset(z for z in range(1, N) if val[z]))
live = set().union(*tight) if tight else set()
nb = {z: set() for z in live}
for S in tight:
    for z in S:
        nb[z] |= S
best = [0]
def ext(cl, cand, wt):
    if wt > best[0]:
        best[0] = wt
    if wt + sum(c[z] for z in cand) <= best[0]:
        return
    for i, z in enumerate(cand):
        ext(cl + [z], [u for u in cand[i + 1:] if u in nb[z]], wt + c[z])
ext([], sorted(live), 0)
print("m=%d n=%d W=%d M=%d: deletion failures %d, tight forms %d, live directions %d, max tight clique (virtual) %d, bound Wn = %d"
      % (m, n, W, len(dirs), bad, len(tight), len(live), best[0], W * n))
print("VERDICT:", "EXCLUDED (deletion)" if bad else ("EXCLUDED (tight clique)" if best[0] > W * n else "not excluded"))
