#!/usr/bin/env python3
"""fp-measure: check of the position formula used by fm_ledger (pos=True): for every frame of the U8X word and every
port line <p>, the position T_p F_U is a Lagrangian (isotropic, dimension h) and its distance to the raw frame L0 equals
the level dim U + 1 - 2[p in U]; for twist 0 the position F_U has distance dim U.  usage: python3 -I pos_check.py U8X.gcert.json.gz"""
import json, gzip
D = json.load(gzip.open(__import__("sys").argv[1], "rt"))   # the U8X word (see the README)
h, fr, ports = D["h"], D["frames"], D["ports"]
def rank(vs):
    piv = {}
    for x in vs:
        for pb in sorted(piv, reverse=True):
            if x >> pb & 1: x ^= piv[pb]
        if x: piv[x.bit_length() - 1] = x
    return len(piv)
def inspan(basis, u):
    for b in sorted(basis, key=lambda b: -b.bit_length()):
        pass
    return rank(list(basis) + [u]) == rank(basis)
def symp(v, w):
    a1, b1, a2, b2 = v >> h, v & ((1 << h) - 1), w >> h, w & ((1 << h) - 1)
    return (bin(a1 & b2).count("1") + bin(a2 & b1).count("1")) & 1
bad = n = 0
for f, basis in enumerate(fr):
    perp = [a for a in range(1, 1 << h) if all(bin(a & c).count("1") % 2 == 0 for c in basis)]
    pb = {}
    for a in perp:                                   # a basis of U^perp
        for q in sorted(pb, reverse=True):
            if a >> q & 1: a ^= pb[q]
        if a: pb[a.bit_length() - 1] = a
    F = [(c << h) | c for c in basis] + list(pb.values())
    assert rank(F) == h and all(symp(v, w) == 0 for v in F for w in F)
    dU = len(basis)
    for p in [None] + ports:
        if p is None:
            P, lev = F, dU
        else:
            P = []
            for w in F:
                a, b = w >> h, w & ((1 << h) - 1)
                if bin(p & b).count("1") & 1: a ^= p
                P.append((a << h) | b)
            lev = dU + 1 - 2 * (1 if inspan(basis, p) else 0)
        ok = rank(P) == h and all(symp(v, w) == 0 for v in P for w in P)
        dist = rank([w >> h for w in P])          # d(X, L0) = h - dim(X cap L0) = rank of the x parts
        n += 1
        if not ok or dist != lev: bad += 1
print("frames %d x (1 + %d twists): %d positions checked, %d mismatches" % (len(fr), len(ports), n, bad))
