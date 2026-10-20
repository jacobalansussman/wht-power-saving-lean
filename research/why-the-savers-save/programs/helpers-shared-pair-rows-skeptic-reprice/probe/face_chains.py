#!/usr/bin/env python3
"""v-pair-market-reprice: face-copy sharing once readers sit at their own planes.  A copy read by t1..tk in that order
lifts to P_face + G_t1 + ... ; reader tj needs every earlier reader's frame inside u_tj^perp.  U8's frames already
satisfy this; the pair-market design adds P_own(ti) to each reader's frame, so the NEW condition is
P_own(ti) inside u_tj^perp for all i < j on the copy.  For each of the 168 face planes (8 readers): the minimum
number of copies = minimum cover of the readers by such chains (exact, brute force over orders), capped at q = 4.
Usage: python3 -I face_chains.py WORD DEMAND BLOCKS"""
import sys, json, gzip, itertools
from collections import Counter, defaultdict
from functools import lru_cache
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
d = json.load(open(sys.argv[2])); U8 = json.load(open(sys.argv[3]))['blocks']
v = c['v']; ports = c['ports']; h = c['h']
blk = {s: i for i, x in enumerate(U8) for s in x}
def B_(gens):
    piv = {}
    for a in gens:
        while a:
            bit = a.bit_length() - 1
            if bit in piv: a ^= piv[bit]
            else: piv[bit] = a; break
    return len(piv), piv
def cont(a, b): return B_(list(a) + list(b))[0] == B_(b)[0]
perp = {t: [x for x in range(1, 1 << h) if bin(x & ports[t]).count('1') % 2 == 0] for t in range(v)}
dem = {int(t): {int(s): x for s, x in dd.items()} for t, dd in d['demand'].items()}
readers = defaultdict(list)
for t in range(v):
    per = defaultdict(list)
    for s in dem[t]: per[blk[s]].append(s)
    for b, ss in per.items(): readers[tuple(sorted(ss))].append(t)
own = {t: [ports[s] for s in ss] for ss, ts in readers.items() if len(ts) == 4 for t in ts}
def inperp(vecs, t): return all(bin(x & ports[t]).count('1') % 2 == 0 for x in vecs)
comp = {(a, b): inperp(own[a], b) for a in range(v) for b in range(v) if a != b}
print('ordered pairs with P_own(a) in u_b^perp:', sum(comp.values()), 'of', len(comp))
faces = {ss: ts for ss, ts in readers.items() if len(ts) == 8}
print('face planes:', len(faces))
def chain_ok(seq): return all(comp[(seq[i], seq[j])] for i in range(len(seq)) for j in range(i + 1, len(seq)))
def min_cover(ts, q=4):
    ts = tuple(ts)
    @lru_cache(None)
    def f(mask):
        if mask == 0: return 0
        rest = [ts[i] for i in range(len(ts)) if mask >> i & 1]
        lead = rest[0]; best = 99
        others = rest[1:]
        for k in range(0, min(q, len(rest)) ):
            for comb in itertools.combinations(others, k):
                grp = (lead,) + comb
                if any(chain_ok(p) for p in itertools.permutations(grp)):
                    m2 = mask
                    for x in grp: m2 &= ~(1 << ts.index(x))
                    best = min(best, 1 + f(m2))
        return best
    return f((1 << len(ts)) - 1)
dist = Counter(); maxchain = Counter()
for ss, ts in faces.items():
    dist[min_cover(ts)] += 1
print('min copies per face plane (q<=4, own-plane condition):', dict(sorted(dist.items())))
print('U8 / model face side uses 2 copies per face plane (8 readers, q = 4)')
