"""an-synth calc/clique.py -- lower bound on the binary window of a two-field design.
If ports S_1..S_w are pairwise junk (A != 0), their labels must be pairwise F_2-orthogonal with odd weight,
so B restricted to them is the identity and rank_F2(B) >= w. Hence window k >= clique number of the junk graph.
Exact maximum clique by simple branch and bound (graphs here have <= 120 vertices). python3 -I calc/clique.py"""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ade import roots_E8_scaled, lines_from_roots, dot, complex_E8, herm

def max_clique(adj):
    n = len(adj); best = []
    order = sorted(range(n), key=lambda v: -bin(adj[v]).count('1'))
    def expand(R, P):
        nonlocal best
        if not P:
            if len(R) > len(best): best = R[:]
            return
        # bound: |R| + popcount(P)
        if len(R) + bin(P).count('1') <= len(best): return
        # greedy colouring bound
        cols = 0; Q = P; colour_bound = 0
        while Q:
            colour_bound += 1; avail = Q
            while avail:
                v = avail.bit_length() - 1; avail &= ~adj[v] & ~(1 << v); Q &= ~(1 << v)
        if len(R) + colour_bound <= len(best): return
        while P:
            v = P.bit_length() - 1
            expand(R + [v], P & adj[v])
            P &= ~(1 << v)
            if len(R) + bin(P).count('1') <= len(best): return
    expand([], (1 << n) - 1)
    return best

def graph(n, f):
    adj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and f(i, j): adj[i] |= 1 << j
    return adj

L = lines_from_roots(roots_E8_scaled())
adjR = graph(len(L), lambda i, j: abs(dot(L[i], L[j])) == 4)
cR = max_clique(adjR)
print("real E8 (120 lines): clique number of junk graph =", len(cR), "(window >= this; achieved window 9)")
Z = complex_E8()
adjC = graph(len(Z), lambda i, j: herm(Z[i], Z[j]) != (0, 0))
cC = max_clique(adjC)
print("complex E8 (60 lines): clique number of junk graph =", len(cC), "=> window >=", len(cC),
      "=> copy ratio t*k/n >=", round(4 * len(cC) / 60, 3), "(real E8: 0.600)")
