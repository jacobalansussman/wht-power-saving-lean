#!/usr/bin/env python3
"""v-sandwich-math: independent counts (no import of frames.py).
Vectors of F_2^(2m): int v = x | z<<m ; omega(u,v) = x.z' + x'.z.
Lagrangians enumerated by extension of isotropic subspaces (canonical = frozenset of elements)."""
import sys, itertools
from collections import Counter
def popc(a): return bin(a).count("1")
def om(u, v, m):
    M = (1 << m) - 1
    return (popc((u & M) & (v >> m)) + popc((v & M) & (u >> m))) & 1
def span(vs):
    s = {0}
    for v in vs:
        if v not in s: s |= {x ^ v for x in s}
    return frozenset(s)
def lagrangians(m):
    N = 1 << (2 * m)
    level = {frozenset([0])}
    for k in range(m):
        nxt = set()
        for S in level:
            for v in range(1, N):
                if v in S: continue
                if all(om(v, w, m) == 0 for w in S):
                    nxt.add(frozenset(S | {x ^ v for x in S}))
        level = nxt
    return sorted(level, key=lambda s: sorted(s))
def dim(S): return len(S).bit_length() - 1
def subspaces(m):
    out = set()
    for k in range(m + 1):
        for vs in itertools.combinations(range(1, 1 << m), k):
            out.add(span(vs))
    return out
def perp(U, m): return frozenset(v for v in range(1 << m) if all(popc(u & v) % 2 == 0 for u in U))
def L_css(U, m):   # CSS convention (U|0)+(0|U^perp)
    return frozenset(u | (v << m) for u in U for v in perp(U, m))
def L_ft(U, m):    # frames.py convention {(u|u+v)}
    return frozenset(u | ((u ^ v) << m) for u in U for v in perp(U, m))
for m in (1, 2, 3, 4):
    Ls = lagrangians(m)
    want = 1
    for i in range(1, m + 1): want *= 2 ** i + 1
    L0 = frozenset(z << m for z in range(1 << m)); LF = frozenset(z | (z << m) for z in range(1 << m))
    subs = subspaces(m)
    Sft = {L_ft(U, m) for U in subs}; Scss = {L_css(U, m) for U in subs}
    d = lambda A, B: m - dim(A & B)
    I0F = {Z for Z in Ls if d(L0, Z) + d(Z, LF) == m}
    exc = Counter(d(L0, Z) + d(Z, LF) - m for Z in Ls)
    nb = [Z for Z in Ls if d(L0, Z) == 1]
    nbsub = [Z for Z in nb if Z in Sft]
    nbnon = [Z for Z in nb if Z not in Sft]
    transv = all(d(Z, LF) == m for Z in nbnon)
    print(f"m={m}: Lagrangians {len(Ls)} (formula {want}); subspaces {len(subs)}; frames.py-subspace frames all Lagrangian: {Sft <= set(Ls)}; "
          f"I(L0,LF)==subspace frames(frames.py conv): {I0F == Sft}; I(L0,LF)==CSS frames: {I0F == Scss}")
    print(f"   excess histogram {dict(sorted(exc.items()))}; neighbours of L0: {len(nb)} = 2(2^m-1)={2*(2**m-1)}; subspace {len(nbsub)}, "
          f"non-subspace {len(nbnon)}, non-subspace all transversal to L_F: {transv}")
    odd = [u for u in range(1, 1 << m) if popc(u) & 1]
    orth = sum(1 for a in odd for b in odd if a != b and popc(a & b) % 2 == 0)
    print(f"   odd labels {len(odd)}; ordered distinct pairs {len(odd)*(len(odd)-1)}; orthogonal ordered pairs {orth}; "
          f"claim formula (2^(m-1)-1)2^(m-2) = {(2**(m-1)-1)*2**(m-2) if m>=2 else '-'}")
for m in (5, 6, 7, 8, 9):
    odd = [u for u in range(1, 1 << m) if popc(u) & 1]
    orth = sum(1 for a in odd for b in odd if a != b and popc(a & b) % 2 == 0)
    print(f"m={m}: odd {len(odd)}, all ordered {len(odd)*(len(odd)-1)}, orthogonal {orth}, claim formula {(2**(m-1)-1)*2**(m-2)}, "
          f"ratio {len(odd)*(len(odd)-1)/orth:.4f}; L0-neighbours 2(2^m-1)={2*(2**m-1)}, subspace {2**m-1}, odd-label subspace {len(odd)}")
