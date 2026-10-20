#!/usr/bin/env python3
"""d-bit-exits / kernel.py: entrances for #589's bit helpers through target-response kernels (the early-reset kernel
of PR435/417, carried to the bit supplier #589).

A pivot helper p starts at its first gate frame E instead of 0. D (frame 0) can no longer read p, so p's target
response column K_p (the targets whose D gate reads p) must be cancelled: donors d with K_p = XOR of K_d, each donor
at E, get 'd += p' at E before A, and 'd += p' at FULL after C. Then the main word sees s'_d = s_d + s_p and the
response of s_p is K_p + sum K_d = 0, so D drops p. A donor whose first frame is E pays nothing extra; a donor whose first
frame strictly contains E splits its first run (0 -> E -> F1).

Within one first-frame class G_E (no split donors), the largest pivot set is |G_E| - rank_F2(K columns of G_E): any
column basis of G_E serves as donors and every other member is a pivot. This script computes that exactly, then builds
the word (pivots start at E, pre-gates, post-gates, D without the pivots), replays it exactly on all 661 columns with
symbolic helper contents, and writes the local paid profile.
usage: python3 -I kernel.py WORD.json.gz OUTDIR [--split]"""
import gzip, json, os, sys, copy
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frames_lib import Frames

w = json.load(gzip.open(sys.argv[1], 'rt')); outdir = sys.argv[2]; os.makedirs(outdir, exist_ok=True)
SPLIT = '--split' in sys.argv
V, R = w['v'], w['R']; N = 2 * V + R
F = Frames(w)
helpers = list(range(2 * V, N))
K = defaultdict(int)                                   # helper -> bitmask of targets whose D gate reads it
for f, t, ss in w['D']:
    assert F.dim[f] == 0 and V <= t < 2 * V
    for s in ss: K[s] ^= 1 << (t - V)
first = {}
cur = list(w['start'])
for f, t, ss in w['A'] + w['B']:
    for r in [t, *ss]:
        if r >= 2 * V and r not in first: first[r] = f
assert all(h in first for h in helpers)

def reduce(basis, x):
    while x:
        p = x.bit_length() - 1
        if p in basis: x ^= basis[p][0]
        else: return x
    return 0

def express(basis, x):
    """basis: pivotbit -> (vector, set of donor helpers); returns the donor set whose XOR is x, or None"""
    used = set()
    while x:
        p = x.bit_length() - 1
        if p not in basis: return None
        v, s = basis[p]; x ^= v; used ^= s
    return used

groups = defaultdict(list)
for h in helpers: groups[first[h]].append(h)
pivots = {}                                          # pivot -> (E, donor set)
for E, G in sorted(groups.items()):
    basis = {}
    for h in G:                                      # greedy in register order: first independent ones are donors
        x = K[h]; s = {h}
        while x:
            p = x.bit_length() - 1
            if p in basis: x ^= basis[p][0]; s ^= basis[p][1]
            else: break
        if x: basis[x.bit_length() - 1] = (x, s)
        else:
            d = s - {h}
            assert d and not (d & set(pivots)), (h, d)
            pivots[h] = (E, sorted(d))
summary = dict(classes=len(groups), class_sizes=dict(sorted(Counter(len(G) for G in groups.values()).items())),
               K_rank_all=None, pivots=len(pivots),
               pivot_entrance_dims=dict(sorted(Counter(F.dim[E] for E, _ in pivots.values()).items())),
               donors=len({d for _, ds in pivots.values() for d in ds}))
# global rank of K
basis = {}
for h in helpers:
    x = K[h]
    while x:
        p = x.bit_length() - 1
        if p in basis: x ^= basis[p]
        else: basis[p] = x; break
summary['K_rank_all'] = len(basis)

# build the word
w2 = copy.deepcopy(w)
P = set(pivots)
w2['D'] = [[f, t, [s for s in ss if s not in P]] for f, t, ss in w['D']]
w2['D'] = [g for g in w2['D'] if g[2]]
pre = []; post = []
for p, (E, ds) in sorted(pivots.items()):
    for d in ds:
        pre.append([E, d, [p]])
    w2['start'][p] = E
FULL = F.FULL
w2['A'] = pre + w['A']
# C: literal reverse of the z-gates of A and B at FULL (pre-gates included, so they are undone last)
z = [g for g in w2['A'] + w2['B'] if g[1] < V or g[1] >= 2 * V]
w2['C'] = [[FULL, t, list(ss)] for f, t, ss in reversed(z)]

# exact replay (same rules as xcheck_p.py Part F, with helper starts allowed at their entrance)
def replay(w):
    cur = list(w['start']); val = [1 << i for i in range(N)]; H = Counter()
    def climb(r, f):
        a = cur[r]
        if a == f: return
        assert F.le(a, f) and F.dim[f] > F.dim[a], ('climb', r, a, f)
        H[F.dim[f] - F.dim[a]] += 1; cur[r] = f
    def gates(sec):
        for f, t, ss in w[sec]:
            for r in [t, *ss]: climb(r, f)
            for s in ss: val[t] ^= val[s]
    for f, t, ss in w['D']:
        assert F.dim[f] == 0 and all(F.dim[cur[r]] == 0 for r in [t, *ss]), 'D operand not at frame 0'
    gates('D'); gates('A')
    tops = [F.dim[cur[q]] for q in w['ret']]; snap = [val[q] for q in w['ret']]
    for t, ks in enumerate(w['scatter']):
        for k in ks: val[V + t] ^= snap[k]
    gates('B'); gates('C')
    for r in range(N): climb(r, w['final'][r])
    for r in range(N):
        want = (1 << r) ^ (1 << (r - V)) if V <= r < 2 * V else 1 << r
        assert val[r] == want, ('endpoint', r)
    return dict(sorted((H + Counter(tops)).items()))

base = replay(w)
assert base == {1: 1323, 2: 1342, 3: 79, 4: 44, 5: 102, 6: 11, 7: 19}
new = replay(w2)
summary['local_profile_base'] = base
summary['local_profile_new'] = new
summary['mass_base'] = sum(r * n for r, n in base.items()); summary['mass_new'] = sum(r * n for r, n in new.items())
summary['ENT_per_call'] = sum(F.dim[w2['start'][h]] for h in helpers)
summary['width_inventory'] = dict(sorted(Counter(8 - F.dim[w2['start'][h]] for h in helpers).items()))
print(json.dumps(summary, indent=1))
json.dump(dict(summary=summary, pivots={str(p): [E, ds] for p, (E, ds) in pivots.items()}),
          open(os.path.join(outdir, 'kernel.json'), 'w'), indent=1)
data = gzip.compress(json.dumps(w2, sort_keys=True, separators=(',', ':')).encode(), compresslevel=9, mtime=0)
open(os.path.join(outdir, 'bitE8mut-p111-ent.json.gz'), 'wb').write(data)
