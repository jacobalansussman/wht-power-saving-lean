#!/usr/bin/env python3
"""e8_minrank_check.py (t-family) -- runs e8_minrank7.solve on (i) every total a and two cliques each (uniform twist:
forced zeros at every non-orthogonal pair), (ii) the NON-UNIFORM relaxation (per-total twists; a non-orthogonal
pair S,T is forced to be known in total a only when no zero-sum class of totals containing a exists, i.e. when S,T
share exactly two coordinates; same-support pairs are released), and (iii) a sanity instance where a completion of
clique size is known to exist (D4 frame decoder total: 6 lines, clique 3, colouring 3).  python3 -I e8_minrank_check.py"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import e8_minrank7 as E
G, L, live = E.G, E.L, E.live
out = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)
supp = {i: frozenset(c for c in range(8) if L[i][c]) for i in range(120)}
for mode in ('uniform', 'nonuniform'):
    for a in range(8):
        A = [S for S in live if G[S][a] != 0]
        if mode == 'uniform':
            forced = lambda S, T: G[S][T] != 0
        else:
            forced = lambda S, T: G[S][T] != 0 and len(supp[S] & supp[T]) == 2
        adj = {S: {T for T in A if T != S and forced(S, T)} for S in A}
        cl = E.all_cliques(A, adj, 7); cl8 = E.all_cliques(A, adj, 8)
        res = []
        if not cl8:
            G_saved = E.G
            # solve() reads E.G[U][S] == 0 for "may be 1" and != 0 for "forced"; emulate the relaxation with a mask
            E.G = [[(1 if (S != T and forced(S, T)) else 0) if (S in A and T in A) else G_saved[S][T]
                    for T in range(120)] for S in range(120)]
            for K in (cl[0], cl[len(cl) // 2]):
                r, n, dt, ds, asg = E.solve(A, K, adj)
                res.append('EXISTS' if r else 'none')
            E.G = G_saved
        say('%-10s total a=%d: forced-known pairs %d, 7-cliques %d, 8-cliques %d, rank-7 completion: %s'
            % (mode, a, sum(len(adj[S]) for S in A) // 2, len(cl), len(cl8), res or 'n/a (clique 8 => rank >= 8)'))
# sanity: D4 lines e_i +- e_j, total a = lines through coordinate 0 (6 lines), rank-3 completion must exist
D4 = [tuple((1 if k == i else (s if k == j else 0)) for k in range(4)) for i, j in itertools.combinations(range(4), 2) for s in (1, -1)]
A = [x for x in D4 if x[0] != 0]
ip = lambda x, y: sum(p * q for p, q in zip(x, y))
idx = {x: i for i, x in enumerate(A)}
Gs = [[ip(x, y) for y in A] for x in A]
G_saved = E.G; E.G = Gs
adj = {i: {j for j in range(6) if j != i and Gs[i][j] != 0} for i in range(6)}
cl = E.all_cliques(range(6), adj, 3)
r, n, dt, ds, asg = E.solve(list(range(6)), cl[0], adj)
E.G = G_saved
say('sanity D4 total (6 lines, clique 3): rank-3 completion', 'EXISTS' if r else 'none', '(expected EXISTS)')
open(os.path.join(HERE, 'out', 'e8_minrank.txt'), 'w').write('\n'.join(out) + '\n')
