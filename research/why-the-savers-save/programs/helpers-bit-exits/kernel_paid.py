#!/usr/bin/env python3
"""d-bit-exits / kernel_paid.py: response-kernel entrances for #589 helpers when donor or pivot runs may be split.

Same model as kernel_z.py, but an inserted 'd += p' at time k may also run at a frame G strictly between a register's
current frame c_r(k) and its next A frame n_r(k): that splits the run c_r -> n_r into c_r -> G -> n_r (one extra paid
child, same rank mass). For every helper p: the least number of splits s (0..3) with which K_p is a sum of available
injection responses (zero-cost vectors span a free subspace; costly vectors searched by meet-in-the-middle mod that span).
Also flags whether p's entrance frame is a coordinate frame (spanned by lines 0..7: the only entrance frames the pricing
chain's residual-isometry family accepts, bank_routes.validate_isometries).
usage: python3 -I kernel_paid.py WORD.json.gz [out.json]"""
import gzip, json, os, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frames_lib import Frames

w = json.load(gzip.open(sys.argv[1], 'rt'))
V, R = w['v'], w['R']; N = 2 * V + R; F = Frames(w)
A, B = w['A'], w['B']; helpers = list(range(2 * V, N))
isz = lambda r: not (V <= r < 2 * V)
touches = defaultdict(list); first_write = {}
for k, (f, t, ss) in enumerate(A):
    for r in [t, *ss]: touches[r].append((k, f))
    first_write.setdefault(t, k)
def cn(r, k, start):
    c = start; n = None
    for kk, f in touches[r]:
        if kk < k: c = f
        else: n = f; break
    return c, n
resp = [0] * N
for t in range(V): resp[V + t] = 1 << t
for f, t, ss in reversed(B):
    for s in ss: resp[s] ^= resp[t]
for t, ks in enumerate(w['scatter']):
    for k in ks: resp[w['ret'][k]] ^= 1 << t
respA = [None] * (len(A) + 1); respA[len(A)] = list(resp)
for k in range(len(A) - 1, -1, -1):
    f, t, ss = A[k]
    for s in ss: resp[s] ^= resp[t]
    respA[k] = list(resp)
K = {h: respA[0][h] for h in helpers}
F1 = {h: touches[h][0][1] for h in helpers}
between = {}
def frames_between(a, b):
    """catalogue frames G with a < G < b strictly (nested)"""
    key = (a, b)
    if key not in between:
        between[key] = [g for g in range(len(F.mask)) if g not in (a, b) and F.le(a, g) and F.le(g, b)]
    return between[key]
def coord(E):
    labs = F.labels(E)
    return F.dim[E] > 0 and all(i < 8 for i in labs) and len(labs) == F.dim[E]

def reduce(basis, x):
    while x:
        b = x.bit_length() - 1
        if b in basis: x ^= basis[b]
        else: break
    # full reduction to a canonical coset representative
    y = x; out = 0
    while y:
        b = y.bit_length() - 1
        if b in basis: y ^= basis[b]
        else: out ^= 1 << b; y ^= 1 << b
    return out
def addb(basis, v):
    x = v
    while x:
        b = x.bit_length() - 1
        if b in basis: x ^= basis[b]
        else: basis[b] = x; return True
    return False

res = {}; stat = Counter()
for p in helpers:
    E = F1[p]; fw = first_write.get(p, len(A))
    free = {}; costly = {}
    for k in range(0, fw + 1):
        cp, np_ = cn(p, k, E)
        Gp = {cp: 0} if np_ is None else {cp: 0, np_: 0, **{g: 1 for g in frames_between(cp, np_)}}
        for d in range(N):
            if d == p or not isz(d): continue
            cd, nd = cn(d, k, w['start'][d])
            if nd is None: Gd = {cd: 0}
            else: Gd = {cd: 0, nd: 0, **{g: 1 for g in frames_between(cd, nd)}}
            for G, a in Gp.items():
                if G in Gd and F.dim[G] > 0:
                    c = a + Gd[G]; v = respA[k][d]
                    if c == 0: addb(free, v)
                    elif c == 1: costly.setdefault(v, (k, d, G))
    t = reduce(free, K[p])
    best = None
    if t == 0: best = 0
    else:
        red = {}
        for v in costly:
            red.setdefault(reduce(free, v), v)
        if t in red: best = 1
        else:
            keys = list(red)
            for a in keys:
                if reduce(free, a ^ t) in red: best = 2; break
    res[p] = dict(splits=best, coordinate=coord(E), entrance_dim=F.dim[E], free_rank=len(free), costly=len(costly))
    stat['splits=%s' % best] += 1
    if coord(E): stat['coordinate entrance, splits=%s' % best] += 1
out = dict(stats=dict(sorted(stat.items())))
print(json.dumps(out, indent=1))
if len(sys.argv) > 2: json.dump(dict(out, per_helper={str(p): v for p, v in res.items()}), open(sys.argv[2], 'w'), indent=1)
