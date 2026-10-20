# flatcover.py: independent check (no quadratic-form tables) of the covering behind Corollary 13.
# F = {z : z_1 = z_2 = 1} (an even (m-2)-flat avoiding 0: it holds no basis vector).  Claim: for d <= m-3, every
# d-dimensional affine subflat A of F lies in an (m-2)-flat G = {z : u.z = c1, w.z = c2} with 0 not in G and an ODD
# number of basis vectors in G.  (G's indicator (u.z+c1+1)(w.z+c2+1) is then an odd binary quadratic form with Q(0)=0.)
# Brute force over all subflats A of F of dimension d and all pairs (u, w).
import sys, itertools
m = int(sys.argv[1]); d = int(sys.argv[2])
def dot(a, b): return bin(a & b).count("1") & 1
F = [z for z in range(1 << m) if (z & 3) == 3]
VF = [z for z in range(1 << m) if (z & 3) == 0]           # direction space of F
# all d-dim linear subspaces of VF, via reduced bases
def span(B):
    S = {0}
    for b in B: S |= {s ^ b for s in S}
    return S
subs = set()
for B in itertools.combinations([v for v in VF if v], d):
    S = frozenset(span(B))
    if len(S) == 1 << d: subs.add(S)
pairs = [(u, w) for u in range(1, 1 << m) for w in range(1, 1 << m) if u < w]   # u != w, both nonzero -> independent
bad = 0; nflats = 0
for V in subs:
    Vb = [v for v in V if v]
    for a in F:
        A = frozenset(a ^ v for v in V)
        if min(A) != a: continue                            # each flat once (a = its least element)
        nflats += 1
        found = False
        for (u, w) in pairs:
            if any(dot(u, v) or dot(w, v) for v in Vb): continue
            c1, c2 = dot(u, a), dot(w, a)
            if c1 == 0 and c2 == 0: continue                # 0 in G
            N = sum(1 for j in range(m) if dot(u, 1 << j) == c1 and dot(w, 1 << j) == c2)
            if N % 2 == 1: found = True; break
        if not found: bad += 1
print(f"m={m} d={d}: {nflats} subflats of F of dimension {d}; without an odd (m-2)-flat avoiding 0 around them: {bad}")
