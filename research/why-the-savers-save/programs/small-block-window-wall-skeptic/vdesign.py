#!/usr/bin/env python3
"""v-an2-smallnet / vdesign.py (skeptic's own code, stdlib only; python3 -I).
Rebuilds from scratch, independently of an2-smallnet/calc: the line systems E8 (from the extended Hamming code
construction: D8 roots + half-integer vectors with even minus signs), P112 = E8 minus one orthogonal frame of 8 lines
(4 coordinate pairs e_{2i} +- e_{2i+1}), E7 = lines orthogonal to one root, E6 = lines orthogonal to an A2, D_n, A_n.
For each: n, t = Q-rank of the Gram, label Gram B over F_2 (B_ij = 1 iff lines orthogonal or equal), window
k = rank_F2(B), junk s per port, sigma = rank of a port's partners' labels, omega = max clique of mutually
non-orthogonal labels, f = max number of a port's partners whose labels lie in a hyperplane of S(a),
and the early/late bound n f / (s-1) (or n f/(s-2) when twin labels exist)."""
import itertools, json, sys
from fractions import Fraction

def e8_lines():
    V = set()
    for i, j in itertools.combinations(range(8), 2):
        for a in (2, -2):
            for b in (2, -2):
                v = [0]*8; v[i] = a; v[j] = b; V.add(tuple(v))
    for s in itertools.product((1, -1), repeat=8):
        if s.count(-1) % 2 == 0: V.add(s)
    return canon(V)

def canon(V):
    out = set()
    for v in V:
        w = tuple(-x for x in v)
        out.add(max(v, w))
    return sorted(out)

def dot(a, b): return sum(x*y for x, y in zip(a, b))

def qrank(M):
    M = [[Fraction(x) for x in r] for r in M]; rk = 0
    cols = len(M[0])
    for c in range(cols):
        p = None
        for i in range(rk, len(M)):
            if M[i][c] != 0: p = i; break
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [x - f*y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk

def f2rank(vecs):
    piv = {}
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in piv: v ^= piv[h]
            else: piv[h] = v; break
    return len(piv)

def f2basis(vecs):
    piv = {}
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in piv: v ^= piv[h]
            else: piv[h] = v; break
    return list(piv.values())

def analyse(name, L):
    n = len(L); nn = dot(L[0], L[0])
    G = [[dot(a, b) for b in L] for a in L]
    for i in range(n):
        for j in range(n):
            assert G[i][j] in (0, nn, -nn, nn//2, -nn//2), (name, G[i][j])
    B = [[1 if G[i][j] in (0, nn, -nn) else 0 for j in range(n)] for i in range(n)]
    rows = [int(''.join(map(str, B[i])), 2) for i in range(n)]
    t = qrank(G); k = f2rank(rows)
    partners = [[j for j in range(n) if j != i and G[i][j] != 0] for i in range(n)]
    s = len(partners[0]); assert all(len(p) == s for p in partners)
    sig = [f2rank([rows[j] for j in partners[i]]) for i in range(n)]
    twins = any(rows[i] == rows[j] for i in range(n) for j in range(i+1, n))
    # f: max over nonzero functionals on S(a) of partners in its kernel.  Functionals on span(S) = restrictions of
    # functionals on F_2^n: phi(v) = parity(v & w) for w in F_2^n; enumerate via a basis dual: coordinates.
    fmax = 0
    for i in range(n):
        bas = f2basis([rows[j] for j in partners[i]])
        # express each partner in basis coordinates by Gaussian elimination
        coords = []
        for j in partners[i]:
            # solve: find subset of bas xoring to rows[j]
            v = rows[j]; c = 0
            # build reduced echelon with tracking
            coords.append(v)
        d = len(bas)
        # brute force functionals: pick d independent "test" coordinates positions where bas is full rank
        # simpler: enumerate all 2^d - 1 nonzero combinations of a dual basis by enumerating hyperplanes as spans of
        # (d-1)-subsets is costly; instead enumerate functionals phi_w(v) = popcount(v & w) % 2 over w in a set of
        # positions that separates the span.
        pos = []
        M = list(bas)
        # choose coordinate positions giving an invertible d x d minor
        cur = []
        for p in range(n):
            cand = cur + [p]
            proj = [sum(((b >> (n-1-q)) & 1) << idx for idx, q in enumerate(cand)) for b in bas]
            if f2rank(proj) == len(cand):
                cur = cand
            if len(cur) == d: break
        pos = cur
        best = 0
        for wbits in range(1, 1 << d):
            w = 0
            for idx, q in enumerate(pos):
                if (wbits >> idx) & 1: w |= 1 << (n-1-q)
            cnt = sum(1 for v in coords if bin(v & w).count('1') % 2 == 0)
            best = max(best, cnt)
        fmax = max(fmax, best)
    # omega: max clique in graph B=1 (i != j)
    adj = [set(j for j in range(n) if j != i and B[i][j] == 1) for i in range(n)]
    best = [0]
    def bk(R, P, X):
        if not P and not X:
            best[0] = max(best[0], len(R)); return
        if len(R) + len(P) <= best[0]: return
        u = max(P | X, key=lambda q: len(adj[q] & P))
        for v in list(P - adj[u]):
            bk(R | {v}, P & adj[v], X & adj[v]); P = P - {v}; X = X | {v}
    bk(set(), set(range(n)), set())
    tight = all(sum(Fraction(G[i][l]*G[l][j], nn*nn) for l in range(n)) == Fraction(n, t)*Fraction(G[i][j], nn)
                for i in range(n) for j in range(n))
    den = (s - 2) if twins else (s - 1)
    rec = dict(design=name, n=n, t=t, k=k, s=s, sigma_min=min(sig), sigma_max=max(sig), omega=best[0], f=fmax,
               twins=twins, tight=tight, early_bound=n*fmax/den, hub_charge=t*(k-1)/n,
               floor2=1+2*min(sig), floor3=1+3*min(sig))
    print(rec, flush=True)
    return rec

def main():
    E8 = e8_lines()
    frame = set(canon({tuple((2 if q == 2*i else (2*sg if q == 2*i+1 else 0)) for q in range(8))
                       for i in range(4) for sg in (1, -1)}))
    P112 = [l for l in E8 if l not in frame]
    E7 = [l for l in E8 if dot(l, (2, -2, 0, 0, 0, 0, 0, 0)) == 0]
    E6 = [l for l in E8 if dot(l, (0, 0, 0, 0, 0, 0, 2, -2)) == 0 and dot(l, (0, 0, 0, 0, 0, 2, -2, 0)) == 0]
    def Dn(n):
        V = set()
        for i, j in itertools.combinations(range(n), 2):
            for a in (1, -1):
                for b in (1, -1):
                    v = [0]*n; v[i] = a; v[j] = b; V.add(tuple(v))
        return canon(V)
    def An(n):
        V = set()
        for i, j in itertools.permutations(range(n+1), 2):
            v = [0]*(n+1); v[i] = 1; v[j] = -1; V.add(tuple(v))
        return canon(V)
    designs = [('E8', E8), ('P112', P112), ('E7', E7), ('E6', E6)] + [(f'D{q}', Dn(q)) for q in (4, 5, 6, 7, 8)] + \
              [(f'A{q}', An(q)) for q in (3, 4, 5, 6, 7)]
    out = [analyse(nm, L) for nm, L in designs]
    json.dump(out, open(sys.argv[1], 'w'), indent=1)

if __name__ == '__main__':
    main()
