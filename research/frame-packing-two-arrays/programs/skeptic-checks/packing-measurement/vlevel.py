#!/usr/bin/env python3
"""v-fp-measure: (1) level formula in my own exact engine (vcal.Eng): for every pair of subspaces U, S of F_2^m with an
orthonormal basis (m = 2, 3), level(Phi_U Phi_S^-1) = dim U + dim S - 2 dim(U cap S), Phi_U = product of dir u over
an orthonormal basis; (2) Lagrangian check over ALL subspace pairs (degenerate included) at m = 2..4:
d(F_U, F_S) = m - dim(F_U cap F_S) equals the formula, F_U = {(c | a + c) : c in U, a in U^perp}."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcal import Eng

def rank(vs):
    rows = []
    for x in vs:
        for b in rows: x = min(x, x ^ b)
        if x: rows.append(x); rows.sort(reverse=True)
    return len(rows)

def span(vs):
    out = {0}
    for v in vs: out |= {x ^ v for x in out}
    return frozenset(out)

def subspaces(m):
    seen = set()
    for k in range(m + 1):
        for vs in itertools.combinations(range(1, 1 << m), k):
            if rank(vs) == k:
                S = span(vs)
                if S not in seen:
                    seen.add(S); yield S

def dot(a, b): return bin(a & b).count("1") & 1

def onb(S, m):
    """an orthonormal basis of S, or None"""
    d = rank(list(S))
    els = [x for x in S if x]
    for vs in itertools.combinations(els, d):
        if rank(vs) == d and all(dot(a, b) == (1 if a == b else 0) for a in vs for b in vs):
            return vs
    return None

def lag(S, m):
    perp = [a for a in range(1 << m) if all(dot(a, c) == 0 for c in S)]
    return span([(c << m) | (a ^ c) for c in S for a in perp])

bad = n1 = 0
for m in (2, 3):
    E = Eng(m)
    subs = [(S, onb(S, m)) for S in subspaces(m)]
    subs = [(S, b) for S, b in subs if b is not None]
    phi = {}
    for S, b in subs:
        g = E.I
        for u in b: g = E.mul(E.dirc[u], g)
        phi[S] = g
    for (U, bu), (S, bs) in itertools.product(subs, repeat=2):
        g = E.mul(phi[U], E.dag(phi[S]))
        lv = E.level(E.symp(g))
        f = rank(list(U)) + rank(list(S)) - 2 * rank(list(U & S))
        n1 += 1; bad += (lv != f)
print("engine: %d pairs (orthonormal-basis subspaces, m = 2, 3), mismatches %d" % (n1, bad))
bad = n2 = 0
for m in (2, 3, 4):
    subs = list(subspaces(m))
    L = {S: lag(S, m) for S in subs}
    for U, S in itertools.product(subs, repeat=2):
        assert len(L[U]) == 1 << m
        d = m - (len(L[U] & L[S]).bit_length() - 1)
        f = rank(list(U)) + rank(list(S)) - 2 * rank(list(U & S))
        n2 += 1; bad += (d != f)
print("Lagrangian: %d pairs (all subspaces, m = 2..4), mismatches %d" % (n2, bad))
