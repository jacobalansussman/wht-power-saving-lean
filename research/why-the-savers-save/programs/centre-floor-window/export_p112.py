#!/usr/bin/env python3
"""ctr-window / export_p112.py: write P112 labels (verified against root inner products) and the 1260 orthogonal
8-frames (own enumeration) to p112.txt for onesided.c."""
import sys, os
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../medium-block-networks')))
from designs_cx import design, labels, ip
L = design('P112'); n = len(L); labs, k, _ = labels(L)
for i in range(n):
    for j in range(n):
        assert (bin(labs[i] & labs[j]).count('1') & 1) == (1 if (i == j or ip(L[i], L[j]) == 0) else 0)
orth = [[ip(L[i], L[j]) == 0 for j in range(n)] for i in range(n)]
fr = []
def ext(c, cand):
    if len(c) == 8: fr.append(c); return
    for i, v in enumerate(cand): ext(c + [v], [w for w in cand[i+1:] if orth[v][w]])
ext([], list(range(n)))
with open('p112.txt', 'w') as f:
    f.write('%d %d\n' % (n, len(fr)))
    f.write(' '.join(map(str, labs)) + '\n')
    for c in fr: f.write(' '.join(map(str, c)) + '\n')
print(n, len(fr))
