#!/usr/bin/env python3
"""v_forms6.py (v-e3-slice6): my own rebuild, no code from e3-slice6 or ex2-families.
Enumerates all 2^21 binary quadratic forms Q(z)=sum_{j<=k} q_jk z_j z_k on F_2^6 (Q(0)=0) as support bitmasks over
the 63 nonzero points, with trace tr Q = sum_j q_jj. Then:
 (1) checks F = {z1=z2=1} is the support of an EVEN form (z1 z2), |F| = 16;
 (2) for each of the 30 affine hyperplanes H of F (8-point subsets {z in F : l(z)=c}), finds the minimum slack
     |D_G| - W over ODD forms G with H inside supp(G) (U_6, W = 13: |D_G| = |supp G|);
 (3) the near-uniform scan U_6 + Z: for each multiset Z of t extra directions, is there an even weight-16 support F
     avoiding Z such that every hyperplane H of F lies in an odd G with |supp G| + #(Z in G) - W <= 3 ?
"""
import itertools, sys
m = 6
pts = list(range(1, 1 << m))
def bit(z, j): return (z >> j) & 1
mons = []  # (j,k) j<=k
for j in range(m):
    for k in range(j, m):
        mons.append((j, k))
monmask = []
for (j, k) in mons:
    mk = 0
    for z in pts:
        if bit(z, j) and bit(z, k): mk |= 1 << z
    monmask.append(mk)
diag = [1 if j == k else 0 for (j, k) in mons]
# enumerate all 2^21 forms by Gray code
N = len(mons)
supports_odd16 = set(); supports_even16 = set()
allsup = {}
mask = 0; tr = 0
pc = lambda x: bin(x).count('1')
cnt = 0
for g in range(1 << N):
    if g:
        b = (g & -g).bit_length() - 1
        mask ^= monmask[b]; tr ^= diag[b]
    w = pc(mask)
    if w == 16:
        (supports_odd16 if tr else supports_even16).add(mask)
    if mask in allsup and allsup[mask] != tr:
        print("ERROR: same support, different trace", mask); sys.exit(1)
    if w <= 20: allsup[mask] = tr
    cnt += 1
print("forms enumerated", cnt, "; weight-16 supports: odd", len(supports_odd16), "even", len(supports_even16))
ws = sorted(set(pc(s) for s in allsup)); print("weights <= 20 occurring:", ws)
F = 0
for z in pts:
    if bit(z, 0) and bit(z, 1): F |= 1 << z
print("F = {z1=z2=1}: |F| =", pc(F), "; even support?", F in supports_even16, "; odd?", F in supports_odd16)
Fpts = [z for z in pts if (F >> z) & 1]
# affine hyperplanes of F: l a linear functional on F_2^6 not constant on F
hyps = set()
for l in range(1, 1 << m):
    for c in (0, 1):
        H = 0
        for z in Fpts:
            if pc(z & l) % 2 == c: H |= 1 << z
        if pc(H) == 8: hyps.add(H)
print("affine hyperplanes of F:", len(hyps))
W = 13
worst = 0
for H in hyps:
    best = min((pc(G) - W for G in supports_odd16 if H & G == H), default=None)
    if best is None:
        print("hyperplane NOT covered by an odd weight-16 support", bin(H)); worst = 99
    else: worst = max(worst, best)
print("U_6, W=13: max over the 30 hyperplanes of the min slack of an odd form containing it =", worst)
# also: every hyperplane is covered by >= how many odd 16-supports
covs = [sum(1 for G in supports_odd16 if H & G == H) for H in hyps]
print("number of odd 16-supports containing each hyperplane: min", min(covs), "max", max(covs))
