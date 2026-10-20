"""Functions shared with cones.py (copied verbatim)."""
import sys, gzip, json
from fractions import Fraction as F
from collections import Counter, defaultdict
def basis(vs):
    piv = {}
    for a in vs:
        while a:
            bit = a.bit_length()-1
            if bit in piv: a ^= piv[bit]
            else: piv[bit] = a; break
    for k in sorted(piv):
        for j in piv:
            if j > k and piv[j] >> k & 1: piv[j] ^= piv[k]
    return tuple(piv[k] for k in sorted(piv, reverse=True))
def contained(a, b): return basis(list(a)+list(b)) == b
def registers(e): return e[1:3] if e[0] == 'add' else e[1:2]
def apply(e, rows):
    if e[0] == 'add':
        t, s = e[1:3]; q = F(e[3], e[4]); dst = rows[t]
        for k, x in list(rows[s].items()):
            a = dst.get(k, F(0))+q*x
            if a: dst[k] = a
            else: dst.pop(k, None)
    else:
        t = e[1]; q = F(e[2], e[3]); rows[t] = {k: q*x for k, x in rows[t].items()}
def invert(e):
    r = list(e)
    if r[0] == 'add': r[3] = -r[3]
    else: r[2], r[3] = r[3], r[2]
    return r
