#!/usr/bin/env python3
"""uscan.py (w8-rmexact): the decoupled multi-target rank-cover test (cell partitions) when the virtual direction
counts are UNIFORM (forced when H is transitive on directions: Singer cycle Z_{2^m-1}, GL).  For c_x = t/(2^m-1):
  pass(t)  iff  slack(G) = c|G| - a_G >= 0 for all G  and, for every F,
               sum_{x in F} min(c, ms(x)) + 2 sum_{x<y in F} min(c, ms(x,y)) >= slack(F),
  ms(.) = min of slack(G) over G containing the cell's directions.  Any unit W in M with this H has pass(M/W).
Targets: H = singer:k (h h^T over the cyclic group), GL (every non-alternating invertible symmetric matrix, uniform),
trivial.  Exact fractions.  Scans t over --tlist or a grid; prints pass/fail and the first failing row."""
import sys, argparse, itertools
from fractions import Fraction as Fr
ap_ = argparse.ArgumentParser()
ap_.add_argument('--m', type=int, default=5); ap_.add_argument('--H', default='singer:1')
ap_.add_argument('--tlist', default='')
ap_.add_argument('--float', action='store_true')
args = ap_.parse_args(); m = args.m
poly = {3: 0b1011, 4: 0b10011, 5: 0b100101}[m]
idx = [(i, j) for i in range(m) for j in range(i, m)]; nb = len(idx); nf = 1 << nb
dirs = list(range(1, 1 << m)); nd = len(dirs)
# Q(x) via bit masks: for each x, the mask of form bits b with x_i x_j = 1
xmask = [sum(1 << b for b, (i, j) in enumerate(idx) if (x >> i) & 1 and (x >> j) & 1) for x in range(1 << m)]
def gfmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a >> m: a ^= poly
    return r
def gram(cols):
    S = 0
    for b, (i, j) in enumerate(idx):
        v = 0
        for k in range(m): v ^= ((cols[k] >> i) & 1) & ((cols[k] >> j) & 1)
        if v: S |= 1 << b
    return S
if args.H.startswith('singer'):
    k = int(args.H.split(':')[1]); g = 1
    for _ in range(k): g = gfmul(g, 2)
    T = []; p = 1
    while True:
        T.append(gram([gfmul(p, 1 << j) for j in range(m)])); p = gfmul(p, g)
        if p == 1: break
elif args.H == 'trivial': T = [gram([1 << j for j in range(m)])]
elif args.H == 'GL':
    def symrows(S):
        rows = [0] * m
        for b, (i, j) in enumerate(idx):
            if (S >> b) & 1: rows[i] |= 1 << j; rows[j] |= 1 << i
        return rows
    def rank(rows):
        rows = list(rows); r = 0
        for c in range(m):
            pv = next((i for i in range(r, m) if (rows[i] >> c) & 1), None)
            if pv is None: continue
            rows[r], rows[pv] = rows[pv], rows[r]
            for i in range(m):
                if i != r and (rows[i] >> c) & 1: rows[i] ^= rows[r]
            r += 1
        return r
    diagbits = sum(1 << b for b, (i, j) in enumerate(idx) if i == j)
    T = [S for S in range(nf) if (S & diagbits) and rank(symrows(S)) == m]
n = len(T)
print("m=%d H=%s: %d virtual targets (distinct %d)" % (m, args.H, n, len(set(T))), flush=True)
pc = [bin(v).count('1') & 1 for v in range(1 << 16)]
def par(v): return pc[v & 0xFFFF] ^ pc[(v >> 16) & 0xFFFF]
# a_F numerators and supports
supp = [[x for x in dirs if par(F & xmask[x])] for F in range(nf)]
anum = [sum(par(F & S) for S in T) for F in range(nf)]
# for each cell (x,y) (x<=y): the set of distinct (s, anum) among G containing x, y
from collections import defaultdict
cell = defaultdict(set)
for G in range(1, nf):
    P = supp[G]; key = (len(P), anum[G])
    for i in range(len(P)):
        for j in range(i, len(P)): cell[(P[i], P[j])].add(key)
print("tables built", flush=True)
def test(t):
    c = (float(t) if args.float else Fr(t)) / nd
    sl = (lambda s, an: c * s - an / n) if args.float else (lambda s, an: c * s - Fr(an, n))
    for G in range(1, nf):
        if sl(len(supp[G]), anum[G]) < (-1e-12 if args.float else 0): return ('deletion', G, len(supp[G]), Fr(anum[G], n))
    ms = {k: min(min(sl(s, an) for (s, an) in v), c) for k, v in cell.items()}
    worst = None
    for F in range(1, nf):
        P = supp[F]
        tot = sum(ms[(P[i], P[i])] for i in range(len(P))) + 2 * sum(ms[(P[i], P[j])] for i in range(len(P)) for j in range(i + 1, len(P)))
        gap = tot - sl(len(P), anum[F])
        if gap < (-1e-12 if args.float else 0): return ('cover', F, len(P), Fr(anum[F], n), gap)
        if worst is None or gap < worst[0]: worst = (gap, F, len(P))
    return ('pass', worst)
ts = [Fr(v) for v in args.tlist.split(',')] if args.tlist else [Fr(k, 20) for k in range(40, 101)]
for t in ts:
    print("t = M/W = %s (%.4f): %s" % (t, float(t), test(t)), flush=True)
