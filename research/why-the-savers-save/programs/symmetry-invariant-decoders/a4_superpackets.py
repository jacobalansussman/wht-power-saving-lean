#!/usr/bin/env python3
"""g-search a4: (i) cross-block 'super-packets': for two needed signed packets (rows) in different blocks, how many
targets need BOTH with one consistent relative sign (so one 8-source row could serve them); (ii) own-block planes:
are the 4 patterns of each own-block plane the 4 Hadamard rows (so two in-place pairs give all four)?"""
import json, itertools
from fractions import Fraction as Fr
from collections import Counter, defaultdict
d = json.load(open('a0_demand.json')); ports = d['ports']; v = len(ports)
dem = {int(t): {int(s): Fr(x) for s, x in dd.items()} for t, dd in d['demand'].items()}
blocks = json.load(open('a1_blocks.json'))['blocks']; blk = {s: b for b, f in enumerate(blocks) for s in f}
rows_of = defaultdict(list)   # target -> list of (block, canon row, sign)
reader = defaultdict(set)
for t in range(v):
    per = defaultdict(dict)
    for s, x in dem[t].items(): per[blk[s]][s] = x
    for b, m in per.items():
        if len(m) != 4: continue
        ss = sorted(m); x0 = m[ss[0]]
        key = (b, tuple((s, int(m[s]/x0)) for s in ss))
        rows_of[t].append((key, x0)); reader[key].add(t)
face = {k for k, r in reader.items() if len(r) >= 6}
pair_common = Counter()
for t in range(v):
    fr = [(k, x) for k, x in rows_of[t] if k in face]
    for (k1, x1), (k2, x2) in itertools.combinations(sorted(fr), 2):
        pair_common[(k1, k2, x1 / x2)] += 1
hist = Counter(pair_common.values())
print('super-packet (two face rows, one relative sign): readers ->', sorted(hist.items()))
# triples
tri = Counter()
for t in range(v):
    fr = sorted((k, x) for k, x in rows_of[t] if k in face)
    for (k1, x1), (k2, x2), (k3, x3) in itertools.combinations(fr, 3):
        tri[(k1, k2, k3, x1/x2, x1/x3)] += 1
print('three-row super-packets: readers ->', sorted(Counter(tri.values()).items()))
# own-block planes
H = {(1, 1, 1, 1), (1, 1, -1, -1), (1, -1, 1, -1), (1, -1, -1, 1)}
odd = [k for k, r in reader.items() if len(r) < 6]
pats = defaultdict(set)
for k in odd: pats[(k[0], tuple(s for s, _ in k[1]))].add(tuple(e for _, e in k[1]))
def is_hadamard_coset(ps):
    # the 4 patterns, as +-1 vectors, are mutually orthogonal
    ps = list(ps)
    return len(ps) == 4 and all(sum(a*b for a, b in zip(p, q)) == 0 for p, q in itertools.combinations(ps, 2))
print('own-block planes:', len(pats), ' with 4 mutually orthogonal patterns:', sum(1 for p in pats.values() if is_hadamard_coset(p)))
# own-block line splits that give pair rows: for a plane {a,b,c,d}, split {a,b}|{c,d}: patterns need ratios on ab and cd
ok_split = Counter()
for (b, quad), ps in pats.items():
    for i, j in [(0, 1), (0, 2), (0, 3)]:
        rest = [x for x in range(4) if x not in (i, j)]
        r1 = {p[i]*p[j] for p in ps}; r2 = {p[rest[0]]*p[rest[1]] for p in ps}
        ok_split[(len(r1), len(r2))] += 1
print('own-block plane splits: (#ratios on first pair, #ratios on second pair) ->', dict(ok_split))
# face rows: for each line {a,b} inside a block, sign ratios demanded by the face packets containing it
lr = defaultdict(set)
for (b, kk) in face:
    m = dict(kk)
    for s1, s2 in itertools.combinations(sorted(m), 2): lr[(s1, s2)].add(m[s1]*m[s2])
print('lines inside face packets:', len(lr), ' ratio sets:', Counter(len(x) for x in lr.values()))
