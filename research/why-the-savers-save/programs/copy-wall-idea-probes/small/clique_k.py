"""For a maximum clique C of pairwise NON-orthogonal E8 lines (|C| = 8), how many clique lines is each other line
orthogonal to?  (A rank-8 decoder of the bit pattern puts x_T, y_T in the span of those clique coordinates.)"""
import sys; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from common import *
import collections
cl = []
def rec(cur, cand):
    global cl
    if len(cur) == 8:
        cl.append(tuple(cur)); return
    for i, c in enumerate(cand):
        if len(cur) + len(cand) - i < 8: break
        rec(cur + [c], [d for d in cand[i + 1:] if NONORTH[c][d]])
rec([0], [d for d in range(1, V) if NONORTH[0][d]])
print('8-cliques of non-orth graph through line 0:', len(cl))
c9 = []
for C in cl[:2000]:
    for T in range(V):
        if T not in C and all(NONORTH[T][c] for c in C):
            c9.append(C + (T,))
print('9-cliques found from them:', len(c9))
dist = collections.Counter()
for C in cl[:50]:
    ks = collections.Counter(sum(ORTH[T][c] for c in C) for T in range(V) if T not in C)
    dist[tuple(sorted(ks.items()))] += 1
for k, n in dist.items(): print('k_T distribution (k: #lines)', k, 'cliques', n)
