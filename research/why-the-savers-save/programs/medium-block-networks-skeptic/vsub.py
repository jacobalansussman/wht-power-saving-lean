#!/usr/bin/env python3
"""v-cx-medium / vsub.py (own code, stdlib; python3 -I).  A gap screen: template points cx-medium did not list.
With one pair the block is 1 + 3(k-1): k = 6 -> 16, k = 8 -> 22 (k = 4 -> 10, P = 2 -> 16; P = 3 -> 22).  Root designs
have odd windows (A, D, E: 3, 5, 7, 9), so even windows come only from SUB-designs.  Screen: every sub-design of E8 (and
of P112) whose labels lie in a hyperplane (k = 8) of the label space, and every sub-design of E7 in a hyperplane (k = 6):
its size n, clique number omega, and the 'hub-type' charge omega*(k-1) against the one-pair break-even 2n/3.
(Only a screen: hub-type is neither a floor nor a decoder.)"""
import itertools, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from vclique import design, gram_B, rowint, add_basis, reduce_, max_cliques

def screen(name, L, top=5):
    n = len(L); B = gram_B(L); rows = [rowint(B[i]) for i in range(n)]
    basis = {}
    for r in rows: add_basis(basis, r)
    bv = list(basis.values()); k = len(bv)
    # coordinates of each row in the basis bv (solve by elimination over the dict basis)
    def coords(v):
        c = 0
        for i, b in sorted(enumerate(bv), key=lambda t: -t[1].bit_length()):
            pass
        # brute: express v as xor of a subset of bv using the pivot structure
        piv = {b.bit_length() - 1: (i, b) for i, b in enumerate(bv)}
        while v:
            p = v.bit_length() - 1; i, b = piv[p]; v ^= b; c ^= 1 << i
        return c
    C = [coords(r) for r in rows]
    res = []
    for a in range(1, 1 << k):
        S = [i for i in range(n) if bin(C[i] & a).count('1') % 2 == 0]
        if len(S) < 4: continue
        sub = [L[i] for i in S]
        Bs = gram_B(sub)
        adj = [{j for j in range(len(S)) if j != i and Bs[i][j] == 1} for i in range(len(S))]
        bs = {}
        for i in S: add_basis(bs, rows[i])
        ks = len(bs)
        res.append((len(S), ks, a, adj))
    res.sort(key=lambda t: -t[0])
    out = []
    seen = set()
    for nS, ks, a, adj in res:
        if (nS, ks) in seen: continue
        seen.add((nS, ks))
        om, cl = max_cliques(nS, adj)
        out.append((nS, ks, om, len(cl), om * (ks - 1), 2 * nS / 3))
        if len(out) >= top: break
    for o in out:
        print(f"{name} hyperplane section: n={o[0]} k={o[1]} omega={o[2]} maxcliques={o[3]} hub-type {o[4]} vs break-even {o[5]:.1f}"
              f" -> {'SAVES at hub-type' if o[4] < o[5] else 'no'}")

if __name__ == '__main__':
    import json
    E8 = design('E8')
    screen('E8', E8)
    fr = set()
    for i in range(0, 8, 2):
        for s in (1, -1):
            v = [0] * 8; v[i] = 1; v[i + 1] = s; fr.add(tuple(v))
    from fractions import Fraction as Fr
    P112 = [l for l in E8 if tuple(int(x) if x.denominator == 1 else None for x in l) not in fr]
    print('P112 size', len(P112))
    screen('P112', P112)
    screen('E7', design('E7'))
