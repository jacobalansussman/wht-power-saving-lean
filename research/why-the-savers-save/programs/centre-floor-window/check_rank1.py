#!/usr/bin/env python3
"""ctr-window / check_rank1.py: re-check the cheapest term found by a2clique_long.py from root inner products only
(no labels): B_ra = 1 iff r = a or <r,a> = 0; F2 rank of B[R,A]; exact-frame hits recounted by brute force."""
import sys, os, json, itertools
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../medium-block-networks')))
from designs_cx import design, labels, ip
L = design('P112'); n = len(L); labs, k, _ = labels(L)
lab2port = {u: p for p, u in enumerate(labs)}
d = json.load(open('a2clique_long.json'))['A2_search']['1']
def span(vs):
    els = {0}
    for v in vs:
        if v not in els: els |= {x ^ v for x in els}
    return els
R = sorted(lab2port[x] for x in span([int(x, 16) for x in d['W1']]) if x in lab2port)
A = sorted(lab2port[x] for x in span([int(x, 16) for x in d['W2']]) if x in lab2port)
Bm = [[1 if (r == a or ip(L[r], L[a]) == 0) else 0 for a in A] for r in R]
rows = [int(''.join(map(str, row)), 2) for row in Bm]
basis = []
for v in rows:
    for b in basis: v = min(v, v ^ b)
    if v: basis.append(v)
print('R lines', R, 'A lines', A)
print('B[R,A] rows:', [''.join(map(str, row)) for row in Bm])
print('F2 rank', len(basis))
print('R pairwise orthogonal?', all(ip(L[a], L[b]) == 0 for a, b in itertools.combinations(R, 2)))
# frames by brute force
orth = [[ip(L[i], L[j]) == 0 for j in range(n)] for i in range(n)]
fr = []
def ext(c, cand):
    if len(c) == 8: fr.append(set(c)); return
    for i, v in enumerate(cand): ext(c + [v], [w for w in cand[i+1:] if orth[v][w]])
ext([], list(range(n)))
print('frames', len(fr), 'hit on both sides', sum(1 for F in fr if F & set(R) and F & set(A)))
