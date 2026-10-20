#!/usr/bin/env python3
"""ctr-search / h2k4sym.py: the same hypergraph H_EPS as h2k4.py, built by symmetry: triples for Gamma_0 and N_0 only,
mapped to every Gamma_a / N_i by automorphisms of the design (signed coordinate maps preserving the code; each acts
on ports and maps label-subspace port sets to label-subspace port sets, which is asserted).  Cross-check of h2k4.py."""
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
S6 = sorted(S6); IDX = {s: i for i, s in enumerate(S6)}
def words(m): return (m & ((1 << 64) - 1), m >> 64)
W0 = np.array([words(s)[0] for s in S6], dtype=np.uint64); W1 = np.array([words(s)[1] for s in S6], dtype=np.uint64)
POP = np.array([bin(i).count('1') for i in range(256)], dtype=np.uint8)
def popcnt(a): return POP[a.view(np.uint8)].reshape(a.shape[0], 8).sum(1)
def triples_for(g):
    size = bin(g).count('1'); g0, g1 = words(g); out = set()
    inside = popcnt(W0 & np.uint64(g0)) + popcnt(W1 & np.uint64(g1))
    cand = np.nonzero(inside >= size - EPS - 40)[0]
    cin = inside[cand]; order = np.argsort(-cin); cand = cand[order]; cin = cin[order]
    cw0, cw1 = W0[cand], W1[cand]; L = len(cand)
    for x in range(L):
        if cin[x] * 3 < size - EPS: break
        sx = S6[cand[x]] & g
        for y in range(x + 1, L):
            if cin[x] + cin[y] < size - EPS - 20: break
            u = sx | (S6[cand[y]] & g)
            if bin(u).count('1') < size - EPS - 20: continue
            R = g & ~u; r0, r1 = words(R)
            if L > y + 1:
                miss = popcnt(np.uint64(r0) & ~cw0[y + 1:]) + popcnt(np.uint64(r1) & ~cw1[y + 1:])
                for z in (np.nonzero(miss <= EPS)[0] + y + 1).tolist():
                    out.add(tuple(sorted((int(cand[x]), int(cand[y]), int(cand[z])))))
    return out
G0 = (1 << 0) | sum(1 << b for b in range(N) if JUNK[0][b])
N0 = sum(1 << b for b in range(N) if LINES[b][0] != 0)
TA = triples_for(G0); TN = triples_for(N0)
print(f'EPS {EPS}: triples for Gamma_0: {len(TA)}, for N_0: {len(TN)} ({time.time()-t0:.0f} s)', flush=True)
# group elements: port 0 -> a, and coordinate 0 -> i
perms = code_perms()
def find(pred):
    for p in perms:
        for sg in itertools.product((1, -1), repeat=8):
            r = act(p, sg)
            if r and pred(p, r): return r
    return None
def mapset(pi, s):
    m = 0
    while s:
        b = (s & -s).bit_length() - 1; s &= s - 1; m |= 1 << pi[b]
    return m
allT = set()
for a in range(N):
    r = None
    for p in perms:
        rr = act(p, (1,) * 8)
        if rr and rr[0][0] == a: r = rr; break
    if r is None:
        r = find(lambda p, rr: rr[0][0] == a)
    pi = r[0]
    for t in TA:
        im = []
        for k in t:
            ms = mapset(pi, S6[k]); assert ms in IDX, 'image of a subspace port set is not one'
            im.append(IDX[ms])
        allT.add(tuple(sorted(im)))
for i in range(8):
    p = next(p for p in perms if p[0] == i); pi = act(p, (1,) * 8)[0]
    for t in TN:
        allT.add(tuple(sorted(IDX[mapset(pi, S6[k])] for k in t)))
print('all covering triples:', len(allT), f'({time.time()-t0:.0f} s)', flush=True)
pairs = collections.defaultdict(set)
for t in allT:
    for pq in itertools.combinations(t, 2): pairs[pq].add(next(v for v in t if v not in pq))
k4 = 0
for (a, b), cs in pairs.items():
    for c, d in itertools.combinations(sorted(cs), 2):
        if tuple(sorted((a, c, d))) in allT and tuple(sorted((b, c, d))) in allT: k4 += 1
print('4-sets with all four triples covering:', k4, f'({time.time()-t0:.0f} s)')
json.dump(dict(eps=EPS, triples_G0=len(TA), triples_N0=len(TN), all=len(allT), K4=k4), open(os.path.join(HERE, f'h2k4sym_eps{EPS}.json'), 'w'))
