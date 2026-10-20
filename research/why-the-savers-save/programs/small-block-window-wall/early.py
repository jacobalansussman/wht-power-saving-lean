#!/usr/bin/env python3
"""an2-smallnet / early.py (own code, stdlib; python3 -I).  Closes the ordering gap of NOTE sec. 2.3 by counting.
A copy of x_b reaching target a is EARLY if it arrives while Y_a's frame does not yet contain S(a), the span of the
labels of a's junk partners; the labels of early arrivals then lie in a proper subspace of S(a), so a target has at
most f(a) early arrivals, f(a) = the largest number of a's partners whose labels lie in one hyperplane of S(a).
A LATE copy is erased (condition A) at a frame containing S(a).  A source with two late copies, at targets with
different ports, takes up S(a1) + S(a2) = the whole window.  The number of sources with at most one late copy is at
most n * max f / (s - 1) (s = junk partners per port, each source sending s copies).  If that is < n, some source of
every amortised call takes up its window, which is all NOTE sec. 2.4 needs (the block is common to all arrays).
Labels are represented by rows of the label Gram B (an injective linear image of the labels)."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import designs as D

def label_gram(L):
    n = len(L); norm = D.ip(L[0], L[0]); unit = norm // 2
    G = [[D.ip(a, b) for b in L] for a in L]
    return G, [[(1 + G[i][j] // unit) % 2 for j in range(n)] for i in range(n)]

def toint(row): return int(''.join(str(x) for x in row), 2)

def reduce_basis(vs):
    """returns (basis list, function giving coordinates of a vector in that basis)"""
    piv = []   # list of (pivotbit, vector, index)
    basis = []
    for v in vs:
        w = v
        for p, b in piv:
            if (w >> p) & 1: w ^= b
        if w:
            p = w.bit_length() - 1
            # keep reduced echelon: eliminate p from existing
            piv = [(q, (c ^ w) if (c >> p) & 1 else c) for q, c in piv]
            piv.append((p, w))
    return piv

def coords(piv, v):
    """coordinates of v (in the span) w.r.t. the reduced echelon basis piv: bit i set if basis i is used."""
    c = 0
    for i, (p, b) in enumerate(piv):
        if (v >> p) & 1:
            v ^= b; c |= (1 << i)
    assert v == 0
    return c

def analyse(name, L):
    G, B = label_gram(L)
    n = len(L)
    rows = [toint(B[i]) for i in range(n)]
    fmax, smin = 0, None
    for a in range(n):
        part = [b for b in range(n) if b != a and G[a][b] != 0]
        piv = reduce_basis([rows[b] for b in part])
        sig = len(piv)
        cs = [coords(piv, rows[b]) for b in part]
        best = 0
        for phi in range(1, 1 << sig):
            cnt = sum(1 for c in cs if bin(c & phi).count('1') % 2 == 0)
            best = max(best, cnt)
        fmax = max(fmax, best)
        smin = len(part) if smin is None else min(smin, len(part))
    bound = n * fmax / (smin - 1)
    ok = bound < n
    print(f"{name:5s} n={n:4d} junk/port s={smin:3d} max early arrivals per target f={fmax:3d}: sources with <= 1 late copy"
          f" <= {bound:7.1f} {'< n: some source takes up the window in every amortised call' if ok else '>= n: the count does not close the gap'}",
          flush=True)
    return dict(design=name, n=n, s=smin, f=fmax, unsafe_bound=bound, closes=ok)

if __name__ == '__main__':
    E8 = D.lines(D.roots_E8())
    frame = set()
    for i in range(0, 8, 2):
        for sj in (1, -1):
            v = [0] * 8; v[i] = 2; v[i + 1] = 2 * sj; frame.add(tuple(v))
    P112 = [l for l in E8 if l not in frame]
    E7 = [l for l in E8 if D.ip(l, (2, 2, 0, 0, 0, 0, 0, 0)) == 0]
    out = [analyse(nm, L) for nm, L in [('E8', E8), ('P112', P112), ('E7', E7), ('D5', D.lines(D.roots_D(5))),
                                         ('D6', D.lines(D.roots_D(6))), ('A4', D.lines(D.roots_A(4))),
                                         ('A5', D.lines(D.roots_A(5)))]]
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'early.json'), 'w'), indent=1)
