#!/usr/bin/env python3
"""exact_hook.py: EXACT (over Q(i), Fractions, no mod p, no floats) values of f for the hook certificate at m=3, f=5:
module V (x) Lambda^4 V, N = span{e_x (x) e_T : x in T}; f(g) = rank of pi g |_N.  Words: d1, d2*d1 (rank-2 block),
d4*d2*d1 (kernel), and d1 * h * d1 for a few free h (second move after a free letter).
Directions doubled: 2 dir(z) = (1+i) I + (1-i) X_z (scaling does not change ranks)."""
import itertools, sys
from fractions import Fraction
n = 8; k = 4
S4 = list(itertools.combinations(range(n), k))
basis = [(x, T) for x in range(n) for T in S4]
def gm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def dir2(z):   # dict col -> list of (row, gaussian)
    return {c: [(c, (1, 1)), (c ^ z, (1, -1))] for c in range(n)}
def mono(pi, q):  # e_x -> i^q(x) e_pi(x)
    ip = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    return {c: [(pi[c], ip[q[c] % 4])] for c in range(n)}
def compose(A, B):  # (A B) as col-dict: first B then A
    out = {}
    for c in range(n):
        acc = {}
        for (r, v) in B[c]:
            for (r2, w) in A[r]:
                t = gm(v, w); o = acc.get(r2, (0, 0)); acc[r2] = (o[0] + t[0], o[1] + t[1])
        out[c] = [(r, v) for r, v in acc.items() if v != (0, 0)]
    return out
def sort_sign(idx):
    s = 1; idx = list(idx)
    for i in range(len(idx)):
        for j in range(len(idx) - 1 - i):
            if idx[j] > idx[j + 1]: idx[j], idx[j + 1] = idx[j + 1], idx[j]; s = -s
    return s, tuple(idx)
def piG_on_N(G):
    cols = []
    for (x, T) in basis:
        if x not in T: continue
        wT = {}
        for choice in itertools.product(*[G[t] for t in T]):
            rows = [r for (r, c) in choice]
            if len(set(rows)) < k: continue
            coef = (1, 0)
            for (r, c) in choice: coef = gm(coef, c)
            s, key = sort_sign(rows); o = wT.get(key, (0, 0)); wT[key] = (o[0] + s * coef[0], o[1] + s * coef[1])
        col = {}
        for (y, cy) in G[x]:
            for T2, cT in wT.items():
                if cT == (0, 0) or y in T2: continue
                t = gm(cy, cT); o = col.get((y, T2), (0, 0)); col[(y, T2)] = (o[0] + t[0], o[1] + t[1])
        cols.append({kk: (Fraction(v[0]), Fraction(v[1])) for kk, v in col.items() if v != (0, 0)})
    return cols
def cmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def cinv(a):
    d = a[0] * a[0] + a[1] * a[1]; return (a[0] / d, -a[1] / d)
def rank_sparse(vecs):
    """exact rank over Q(i) of a list of sparse vectors (dict key -> (Fraction, Fraction))"""
    piv = {}   # key -> reduced vector with leading key
    order = sorted({kk for v in vecs for kk in v})
    pos = {kk: i for i, kk in enumerate(order)}
    r = 0
    for v in vecs:
        v = dict(v)
        while v:
            lead = min(v, key=lambda kk: pos[kk])
            if lead in piv:
                w = piv[lead]; fct = v[lead]
                for kk, val in w.items():
                    t = cmul(fct, val); o = v.get(kk, (Fraction(0), Fraction(0)))
                    nv = (o[0] - t[0], o[1] - t[1])
                    if nv[0] == 0 and nv[1] == 0: v.pop(kk, None)
                    else: v[kk] = nv
            else:
                inv = cinv(v[lead]); piv[lead] = {kk: cmul(val, inv) for kk, val in v.items()}; r += 1; break
    return r
ident = list(range(n))
TOF2 = mono([x ^ (((x & 1) & ((x >> 1) & 1)) << 2) for x in range(n)], [0] * n)
CS01 = mono(ident, [((x & 1) & ((x >> 1) & 1)) for x in range(n)])
S0 = mono(ident, [(x & 1) for x in range(n)])
words = {
    'd1': dir2(1),
    'd2*d1 (rank-2 block)': compose(dir2(2), dir2(1)),
    'd4*d2*d1 (kernel)': compose(dir2(4), compose(dir2(2), dir2(1))),
    'd2*S0*d1': compose(dir2(2), compose(S0, dir2(1))),
    'd2*CS01*d1': compose(dir2(2), compose(CS01, dir2(1))),
    'd2*TOF2*d1': compose(dir2(2), compose(TOF2, dir2(1))),
}
sel = sys.argv[1:] or list(words)
for name in words:
    if not any(name.startswith(s) for s in sel): continue
    print(f"EXACT over Q(i): f({name}) = {rank_sparse(piG_on_N(words[name]))}", flush=True)
