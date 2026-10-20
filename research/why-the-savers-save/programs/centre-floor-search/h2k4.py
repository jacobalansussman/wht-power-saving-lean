#!/usr/bin/env python3
"""ctr-search / h2k4.py: closes the last Problem D profiles (m = 11: (1,6^9,7), (1,6^9,8), (1,6^8,7,7), (2,6^9,7);
m = 12: (1,1,6^9,7)).  Each has >= 8 dim-6 directions, and NC_7 with J = {its 1 or 2 smallest directions} plus ANY
three dim-6 directions says: every triple of the dim-6 port sets covers a co-flat up to a fixed set of <= EPS ports
(EPS = 1 or 2).  Such sets hold <= 62 ports, so the co-flat is some Gamma_a (53) or N_i (56).
Relaxation: allow each triple its own missing set (<= EPS ports).  If the hypergraph H_EPS of such triples has no
4 vertices with all 4 triples present, every one of these profiles is impossible.  usage: h2k4.py EPS"""
import sys, os, json, time, itertools, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *
from nd import rref_subspaces, mask
EPS = int(sys.argv[1]) if len(sys.argv) > 1 else 2
t0 = time.time()
S6 = set()
for rows in rref_subspaces(9, 3):
    m = (1 << N) - 1
    for r in rows: m &= mask[r]
    if m: S6.add(m)
S6 = sorted(S6)
print('distinct dim-6 port sets:', len(S6), f'({time.time()-t0:.0f} s)', flush=True)
def words(m): return (m & ((1 << 64) - 1), m >> 64)
W0 = np.array([words(s)[0] for s in S6], dtype=np.uint64); W1 = np.array([words(s)[1] for s in S6], dtype=np.uint64)
POP = np.array([bin(i).count('1') for i in range(256)], dtype=np.uint8)
def popcnt(a):
    return POP[a.view(np.uint8)].reshape(a.shape[0], 8).sum(1)
cof = []
for a in range(N): cof.append((1 << a) | sum(1 << b for b in range(N) if JUNK[a][b]))
for i in range(8): cof.append(sum(1 << b for b in range(N) if LINES[b][i] != 0))
triples = set()
for gi, g in enumerate(cof):
    size = bin(g).count('1')
    g0, g1 = words(g)
    inside = popcnt(W0 & np.uint64(g0)) + popcnt(W1 & np.uint64(g1))
    need1 = size - EPS - 40
    cand = np.nonzero(inside >= need1)[0]
    cw0, cw1 = W0[cand], W1[cand]
    cin = inside[cand]
    order = np.argsort(-cin); cand, cw0, cw1, cin = cand[order], cw0[order], cw1[order], cin[order]
    L = len(cand)
    for x in range(L):
        if cin[x] < (size - EPS) / 3: break                 # the largest of the three holds >= a third
        sx = S6[cand[x]] & g
        for y in range(x + 1, L):
            if cin[x] + cin[y] < size - EPS - 20: break
            u = sx | (S6[cand[y]] & g)
            if bin(u).count('1') < size - EPS - 20: continue
            R = g & ~u
            r0, r1 = words(R)
            miss = popcnt(np.uint64(r0) & ~cw0[y + 1:]) + popcnt(np.uint64(r1) & ~cw1[y + 1:]) if L > y + 1 else np.array([])
            for z in (np.nonzero(miss <= EPS)[0] + y + 1).tolist():
                triples.add(tuple(sorted((int(cand[x]), int(cand[y]), int(cand[z])))))
    if gi % 20 == 0 or gi >= N:
        print(f'  co-flat {gi} (size {size}): candidates {L}; triples so far {len(triples)} ({time.time()-t0:.0f} s)', flush=True)
print('covering triples with defect <=', EPS, ':', len(triples), flush=True)
pairs = collections.defaultdict(set)
for t in triples:
    for p in itertools.combinations(t, 2): pairs[p].add(next(v for v in t if v not in p))
k4 = []
for (a, b), cs in pairs.items():
    for c, d in itertools.combinations(sorted(cs), 2):
        if tuple(sorted((a, c, d))) in triples and tuple(sorted((b, c, d))) in triples:
            k4.append((a, b, c, d))
print('4-sets with all four triples in H_eps:', len(k4), f'({time.time()-t0:.0f} s)')
json.dump(dict(eps=EPS, sets=len(S6), triples=len(triples), K4=len(k4), K4_examples=k4[:5]),
          open(os.path.join(HERE, f'h2k4_eps{EPS}.json'), 'w'))
