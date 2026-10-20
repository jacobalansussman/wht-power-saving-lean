#!/usr/bin/env python3
"""Check two dual-polar facts used in the x-lit report (m = 3, all 135 Lagrangians of F_2^6):
 (1) every Lagrangian K on a geodesic between L and L' (d(L,K)+d(K,L') = d(L,L')) contains L cap L';
 (2) for subspace pairs A<=B, I(A,B) = {(a, a+b): a in A, b in B^perp} is isotropic; a sum of several I(A_i,B_i)
     is isotropic iff sum A_i <= cap B_j; the Lagrangians containing I(A,B) number prod_{i=1..k}(2^i+1), k = dim B - dim A,
     and the subspace frames among them are exactly L_W, A <= W <= B.
"""
import itertools, random, sys
sys.path.insert(0, __file__.rsplit('/',1)[0])
from dictcheck import lagrangians, span, symp, popc, rank
m = 3
w = lambda a, b: symp(a, b, m) == 0 if a != b else True
Ls = sorted(lagrangians(m, w), key=lambda S: sorted(S))
dim = lambda S: len(S).bit_length() - 1
d = lambda A, B: m - dim(A & B)
bad1 = 0; ngeo = 0
for L in Ls:
    for L2 in Ls:
        W = L & L2
        for K in Ls:
            if d(L, K) + d(K, L2) == d(L, L2):
                ngeo += 1
                if not W <= K: bad1 += 1
subs = set()
for k in range(m + 1):
    for vecs in itertools.combinations(range(1, 1 << m), k):
        if rank(list(vecs)) == k: subs.add(span(list(vecs)))
subs = sorted(subs, key=lambda S: (len(S), sorted(S)))
def perp(U): return frozenset(v for v in range(1 << m) if all(popc(v & u) % 2 == 0 for u in U))
def I(A, B):
    Bp = perp(B); return span([a | ((a ^ b) << m) for a in A for b in Bp])
def LU(U): return I(U, U)
def isotropic(S):
    S = list(S); return all(symp(a, b, m) == 0 for a in S for b in S)
pairs = [(A, B) for A in subs for B in subs if A <= B]
bad2 = 0; bad3 = 0
for A, B in pairs:
    IA = I(A, B)
    if not isotropic(IA): bad2 += 1
    cont = [K for K in Ls if IA <= K]
    k = dim(B) - dim(A); pred = 1
    for i in range(1, k + 1): pred *= (2 ** i + 1)
    sub_cont = set(U for U in subs if LU(U) in cont)
    interval = set(U for U in subs if A <= U <= B)
    if len(cont) != pred or sub_cont != interval: bad3 += 1
random.seed(1); bad4 = 0; trials = 0
for _ in range(3000):
    r = random.randint(2, 3)
    ps = [random.choice(pairs) for _ in range(r)]
    S = span([v for (A, B) in ps for v in I(A, B)])
    SA = span([a for (A, B) in ps for a in A])
    capB = frozenset.intersection(*[B for (A, B) in ps])
    trials += 1
    if isotropic(S) != (SA <= capB): bad4 += 1
print(f"(1) geodesic triples {ngeo}, violations of 'K contains L cap L2': {bad1}")
print(f"(2) subspace pairs A<=B: {len(pairs)}; non-isotropic I(A,B): {bad2}; count/interval mismatches: {bad3}")
print(f"(3) random sums of 2-3 I(A_i,B_i): {trials}; mismatches of 'isotropic iff sum A <= cap B': {bad4}")
