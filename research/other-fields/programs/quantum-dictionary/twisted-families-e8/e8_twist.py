#!/usr/bin/env python3
"""e8_twist.py (t-family) -- can ANY uniform twist (any Lagrangian label family) make the hub's complex P112 frame
decoder cheaper?  Standard library only; python3 -I e8_twist.py.

Model (twisted pairing floor, t-sandwich): in a family with fixed twists, port T carries a source vector p_T and a
target vector q_S (pulled back into one symplectic space); N[T,S] = omega(p_T, q_S) in F_2, N[S,S] = 1 (forced).
Target S cancels x_T at no copy cost only if N[T,S] = 0.  A total reading R and serving D costs at least
rank_F2 N[R, D].  The hub's frame decoder delivers B = G (the real Gram of the E8 lines, totals A_a = lines
non-orthogonal to 2e_a, a = 1..8).  B[S,T] != 0 forces N[T,S] = 0, so on A_a the matrix N has unit diagonal and
zeros at every non-orthogonal pair; it may be 1 only at orthogonal pairs.  Hence
     cost(A_a) >= minrank_F2(orthogonality pattern on A_a) >= clique number of the NON-orthogonality graph on A_a
(a set of pairwise non-orthogonal lines gives an identity submatrix), and the hub's labels (N = B0 = 1 + G mod 2,
which is 1 at EVERY orthogonal pair) give 8.  Upper bound for any twist: chromatic number of the non-orthogonality
graph (block-diagonal all-ones N on colour classes)."""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '../lib')))
from copyfloor import e8_lines, rank

L, G = e8_lines()
live = list(range(8, 120))
out = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)

def max_clique(V, adj):
    best = []
    def expand(R, P):
        nonlocal best
        if not P:
            if len(R) > len(best): best = list(R)
            return
        if len(R) + len(P) <= len(best): return
        for v in list(P):
            if len(R) + len(P) <= len(best): return
            expand(R + [v], [w for w in P if w in adj[v]])
            P = [w for w in P if w != v]
    expand([], list(V))
    return best

def dsatur(V, adj):
    col = {}
    order = sorted(V, key=lambda v: -len(adj[v] & set(V)))
    while len(col) < len(V):
        un = [v for v in V if v not in col]
        v = max(un, key=lambda v: (len({col[w] for w in adj[v] if w in col}), len(adj[v] & set(V))))
        used = {col[w] for w in adj[v] if w in col}
        c = 0
        while c in used: c += 1
        col[v] = c
    return 1 + max(col.values())

def exact_chromatic(V, adj, ub):
    """smallest k <= ub with a proper colouring (backtracking, DSATUR order)."""
    V = list(V)
    def colourable(k):
        col = {}
        def rec():
            if len(col) == len(V): return True
            un = [v for v in V if v not in col]
            v = max(un, key=lambda v: (len({col[w] for w in adj[v] if w in col}), len(adj[v])))
            used = {col[w] for w in adj[v] if w in col}
            mx = max(col.values()) + 1 if col else 0
            for c in range(min(k, mx + 1)):
                if c in used: continue
                col[v] = c
                if rec(): return True
                del col[v]
            return False
        return rec()
    k = ub
    while k > 1 and colourable(k - 1): k -= 1
    return k

B0 = [[(1 + G[i][j]) % 2 for j in range(120)] for i in range(120)]
tot_lo = tot_old = tot_hi = 0
for a in range(8):
    A = [S for S in live if G[S][a] != 0]
    adj = {S: {T for T in A if T != S and G[S][T] != 0} for S in A}      # non-orthogonality graph on A_a
    cl = max_clique(A, adj)
    ub = dsatur(A, adj)
    chi = ub if ub == len(cl) else '%d..%d' % (len(cl), ub)
    old = rank([[B0[i][j] for j in A] for i in A], 2)
    qr = rank([[B0[i][j] for j in A] for i in A], 'Q')
    tot_lo += len(cl); tot_old += old; tot_hi += ub
    say('total a=%d: |A_a| = %d; clique (pairwise non-orthogonal) %d %s; chromatic %s; hub labels F2-rank %d (Q-rank %d)'
        % (a, len(A), len(cl), [tuple(L[i]) for i in cl][:3] + ['...'], chi, old, qr))
say('frame decoder, any uniform twist: cst >= %d (cliques), hub labels %d, block-colouring twist %d; v = 112'
    % (tot_lo, tot_old, tot_hi))
open(os.path.join(HERE, 'out', 'e8_twist.txt'), 'w').write('\n'.join(out) + '\n')
