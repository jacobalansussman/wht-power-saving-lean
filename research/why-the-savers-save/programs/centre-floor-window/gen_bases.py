#!/usr/bin/env python3
"""ctr-window / gen_bases.py TOPW N SEED OUT: S-side column generation for lpmwu. Reads the worst subspaces W
(line masks) of the last exact check, and builds N bases S of R^8 (8 P112 lines, exact rank 8 over Q) that meet as
few of those L(W) as possible (randomised greedy plus swap descent). Writes them in pool format (sets only)."""
import sys, os, random
from fractions import Fraction
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../medium-block-networks')))
from designs_cx import design, labels, ip
L = design('P112'); n = len(L)
top = []
for line in open(sys.argv[1]):
    r, d, a, b = line.split(); top.append((float(r), int(d), int(a) | (int(b) << 64)))
N = int(sys.argv[2]); rng = random.Random(int(sys.argv[3])); outp = sys.argv[4]
zmax = max(t[0] for t in top)
D = [m for r, d, m in top if r >= 0.97 * zmax]
lm = [0] * n                                  # for each line: bitmask of the D-members whose L(W) contains it
for i, m in enumerate(D):
    for t in range(n):
        if m >> t & 1: lm[t] |= 1 << i
def indep(basis, v):
    w = [Fraction(x) for x in v]
    for piv, b in basis:
        if w[piv] != 0:
            f = w[piv] / b[piv]; w = [x - f * y for x, y in zip(w, b)]
    nz = next((i for i, x in enumerate(w) if x != 0), None)
    return None if nz is None else (nz, w)
def cost(S):
    m = 0
    for t in S: m |= lm[t]
    return bin(m).count('1')
def build():
    S, basis, m = [], [], 0
    while len(S) < 8:
        cands = list(range(n)); rng.shuffle(cands)
        cands.sort(key=lambda t: bin(m | lm[t]).count('1') + rng.random() * 3)
        for t in cands:
            if t in S: continue
            e = indep(basis, L[t])
            if e is not None:
                S.append(t); basis.append(e); m |= lm[t]; break
    return S
def rank8(S):
    basis = []
    for t in S:
        e = indep(basis, L[t])
        if e is None: return False
        basis.append(e)
    return True
out = set(); costs = []
while len(out) < N:
    S = build(); c = cost(S)
    for _ in range(40):                                    # swap descent
        i = rng.randrange(8); t = rng.randrange(n)
        if t in S: continue
        T = S[:i] + [t] + S[i+1:]
        ct = cost(T)
        if ct <= c and rank8(T): S, c = T, ct
    key = tuple(sorted(S))
    if key not in out: out.add(key); costs.append(c)
with open(outp, 'w') as f:
    for S in out: f.write(' '.join(map(str, S)) + '\n')
print('D size', len(D), 'generated', len(out), 'met D-members: min %d mean %.1f max %d' % (min(costs), sum(costs) / len(costs), max(costs)))
