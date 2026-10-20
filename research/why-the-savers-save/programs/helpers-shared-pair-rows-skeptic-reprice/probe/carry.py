#!/usr/bin/env python3
"""v-pair-market-reprice: the most favourable timing for the shared second-half row.  For a sharing pair (A first,
B second) on an own plane, B reads the shared copy LAST (nothing carried), A reads it right after the face planes it
can read first (those inside u_B^perp; the copy then lifts to G_A + L2, inside u_B^perp).  A then carries P_own(A)
into its remaining face reads.  A carrying reader of a face copy must be that copy's last reader unless P_own(A) is
inside the later readers' u^perp.  Exact search over split and roles per own plane (3 splits x 2 pairs x 2 roles),
minimising sum over face planes of extra copies max(0, minimal copies - 2), by coordinate descent from all starts
of the per-plane minimum; reports carried incidences and the face-copy lower bound (q = 4).
Usage: python3 -I carry.py WORD DEMAND BLOCKS"""
import sys, json, gzip, itertools, random
from collections import Counter, defaultdict
from fractions import Fraction as Fr
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
d = json.load(open(sys.argv[2])); U8 = json.load(open(sys.argv[3]))['blocks']
v = c['v']; ports = c['ports']
blk = {s: i for i, x in enumerate(U8) for s in x}
dem = {int(t): {int(s): Fr(x) for s, x in dd.items()} for t, dd in d['demand'].items()}
readers = defaultdict(list); need = defaultdict(dict)
for t in range(v):
    per = defaultdict(list)
    for s in dem[t]: per[blk[s]].append(s)
    for b, ss in per.items(): readers[tuple(sorted(ss))].append(t); need[t][b] = tuple(sorted(ss))
def inperp(srcs, t): return all(bin(ports[s] & ports[t]).count('1') % 2 == 0 for s in srcs)
ownp = {ss: ts for ss, ts in readers.items() if len(ts) == 4}
faces = {ss: ts for ss, ts in readers.items() if len(ts) == 8}
fidx = {ss: i for i, ss in enumerate(faces)}
own_of = {t: ss for ss, ts in ownp.items() for t in ts}
comp = {(a, b): inperp(own_of[a], b) for a in range(v) for b in range(v) if a != b}
def options(ss):
    ts = ownp[ss]; S = list(ss); out = []
    for a in S[1:]:
        L2 = tuple(x for x in S if x not in (S[0], a))
        # pair targets by their ratio on L2
        key = lambda t: dem[t][L2[1]] / dem[t][L2[0]]
        grp = defaultdict(list)
        for t in ts: grp[key(t)].append(t)
        assert all(len(g) == 2 for g in grp.values()), grp
        pairs = list(grp.values())
        for r1 in (0, 1):
            for r2 in (0, 1):
                As = []
                for (p, r) in zip(pairs, (r1, r2)):
                    A, B = (p[0], p[1]) if r == 0 else (p[1], p[0])
                    As.append((A, B))
                carried = set()
                for A, B in As:
                    for b, fs in need[A].items():
                        if fs in faces and not inperp(fs, B): carried.add((A, fs))
                out.append(frozenset(carried))
    return out
opts = {ss: options(ss) for ss in ownp}
def copies(fs, carriers):
    # lower bound on copies: carriers that are pairwise incompatible in BOTH orders need distinct copies
    # (a carrier can share a copy only as its last reader, or before a reader whose u^perp holds its own plane)
    car = [t for t in faces[fs] if t in carriers]; n = len(car); best = 0
    for m in range(1 << n):
        grp = [car[i] for i in range(n) if m >> i & 1]
        if len(grp) > best and all(not comp[(a, b)] and not comp[(b, a)] for a, b in itertools.combinations(grp, 2)):
            best = len(grp)
    return max(2, best)
def total(choice):
    carriers_by_face = defaultdict(set); inc = 0
    for ss, k in choice.items():
        for (A, fs) in opts[ss][k]: carriers_by_face[fs].add(A); inc += 1
    cp = Counter(copies(fs, carriers_by_face[fs]) for fs in faces)
    return sum(k * x for k, x in cp.items()), inc, cp
random.seed(1)
best = None
for start in range(4):
    choice = {ss: (min(range(len(opts[ss])), key=lambda k: len(opts[ss][k])) if start == 0 else random.randrange(len(opts[ss]))) for ss in ownp}
    cur = total(choice)
    improved = True
    while improved:
        improved = False
        for ss in ownp:
            for k in range(len(opts[ss])):
                if k == choice[ss]: continue
                old = choice[ss]; choice[ss] = k; t2 = total(choice)
                if t2[0] < cur[0]: cur = t2; improved = True
                else: choice[ss] = old
    if best is None or cur[0] < best[0]: best = cur
print('min carried (target, face) incidences per own-plane option:', Counter(min(len(o) for o in opts[ss]) for ss in ownp))
print('best found: total face copies', best[0], '(U8/model: 2 x 168 = 336); carried incidences', best[1], '; copies per face plane', dict(sorted(best[2].items())))
