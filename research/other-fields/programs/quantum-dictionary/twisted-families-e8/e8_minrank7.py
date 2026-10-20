#!/usr/bin/env python3
"""e8_minrank7.py (t-family) -- is there an F_2 matrix N on one frame-decoder total A_a (56 E8 lines) with unit
diagonal, N[U,V] = 0 at every non-orthogonal pair, and rank 7 (one less than the hub's labels)?
Clique K of 7 pairwise non-orthogonal lines (e8_twist.py) => N[K,K] = I_7, so the rows p_U and columns q_V of a rank-7
factorisation N = P^T Q are pinned in the basis K: p_U = (N[U,S])_{S in K}, q_V = (N[S,V])_{S in K}, and
N[U,V] = sum_{S in K} N[U,S] N[S,V].  N[U,S] may be 1 only if U, S are orthogonal.  CSP: one variable per U not in K,
value (a_U, b_U) in F_2^7 x F_2^7 supported on the K-lines orthogonal to U; unary a_U.b_U = 1; binary a_U.b_V = 0 and
a_V.b_U = 0 for every non-orthogonal pair U, V.  Exact backtracking with forward checking.  Standard library only.
usage: python3 -I e8_minrank7.py [a]"""
import os, sys, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '../lib')))
from copyfloor import e8_lines
L, G = e8_lines()
live = list(range(8, 120))

def max_clique(V, adj):
    best = []
    def expand(R, P):
        nonlocal best
        if not P:
            if len(R) > len(best): best = list(R)
            return
        for v in list(P):
            if len(R) + len(P) <= len(best): return
            expand(R + [v], [w for w in P if w in adj[v]])
            P = [w for w in P if w != v]
    expand([], list(V))
    return best

def all_cliques(V, adj, k):
    res = []
    def rec(R, P):
        if len(R) == k: res.append(tuple(R)); return
        for i, v in enumerate(P):
            rec(R + [v], [w for w in P[i+1:] if w in adj[v]])
    rec([], sorted(V))
    return res

def solve(A, K, adj, limit=None):
    K = list(K); k = len(K)
    others = [U for U in A if U not in K]
    # allowed bit positions: S in K orthogonal to U
    dom = {}
    for U in others:
        pos = [i for i, S in enumerate(K) if G[U][S] == 0]
        vals = []
        for am in range(1 << len(pos)):
            a = sum(1 << pos[t] for t in range(len(pos)) if am >> t & 1)
            for bm in range(1 << len(pos)):
                b = sum(1 << pos[t] for t in range(len(pos)) if bm >> t & 1)
                if bin(a & b).count('1') & 1:
                    vals.append((a, b))
        dom[U] = vals
    nbr = {U: [V for V in others if V != U and G[U][V] != 0] for U in others}
    assign = {}
    nodes = [0]
    def ok(U, val, V, w):
        return not (bin(val[0] & w[1]).count('1') & 1) and not (bin(w[0] & val[1]).count('1') & 1)
    def rec(doms):
        nodes[0] += 1
        if len(assign) == len(others): return True
        U = min((u for u in others if u not in assign), key=lambda u: len(doms[u]))
        for val in doms[U]:
            nd = dict(doms); bad = False
            for V in nbr[U]:
                if V in assign: continue
                nd[V] = [w for w in doms[V] if ok(U, val, V, w)]
                if not nd[V]: bad = True; break
            if bad: continue
            assign[U] = val
            if rec(nd): return True
            del assign[U]
        return False
    t = time.time()
    r = rec(dom)
    return r, nodes[0], time.time() - t, {U: len(dom[U]) for U in others}, dict(assign)

if __name__ == '__main__':
    a = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    A = [S for S in live if G[S][a] != 0]
    adj = {S: {T for T in A if T != S and G[S][T] != 0} for S in A}
    cl = all_cliques(A, adj, 7)
    print('total a=%d: |A_a|=%d, 7-cliques of pairwise non-orthogonal lines: %d, 8-cliques: %d'
          % (a, len(A), len(cl), len(all_cliques(A, adj, 8))), flush=True)
    K = cl[0]
    r, n, dt, ds, asg = solve(A, K, adj)
    print('clique', [tuple(L[i]) for i in K])
    print('domain sizes', sorted(ds.values()))
    print('rank-7 completion on A_a:', 'EXISTS' if r else 'NONE', '(nodes %d, %.1fs)' % (n, dt))
    if r:
        # verify
        Kl = list(K)
        def vec(U, side):
            if U in Kl: return 1 << Kl.index(U)
            return asg[U][side]
        N = {(U, V): bin(vec(U, 0) & vec(V, 1)).count('1') & 1 for U in A for V in A}
        assert all(N[(U, U)] == 1 for U in A)
        assert all(N[(U, V)] == 0 for U in A for V in A if U != V and G[U][V] != 0)
        print('verified: unit diagonal, zero at all non-orthogonal pairs; ones at %d of %d orthogonal ordered pairs'
              % (sum(N[(U, V)] for U in A for V in A if U != V and G[U][V] == 0),
                 sum(1 for U in A for V in A if U != V and G[U][V] == 0)))
