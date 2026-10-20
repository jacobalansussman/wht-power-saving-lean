"""d-mcopy-code STEP 1b: group the 112 M copies (components of M-M gates) and list per group: members, entrance frames,
targets, sources/helpers they read, values delivered to targets (x-part, as source combinations), gate counts by phase/frame dim.
Usage: python3 -I mc_groups.py FULL-WORD.json.gz mc_groups.json"""
import sys, json
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, '.')
from wordlib import *
w = load(sys.argv[1]); v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
start, end, chain, paid, runs = frames(w)
vis3 = {q for q in range(n) if any(len(U) == 3 for _, U in chain[q])}
ent = {int(q) for q in w['start_extra']}
M = sorted(q for q in range(2*v, n) if q in ent and q not in vis3); Ms = set(M)
cls = lambda q: 'x' if q < v else ('y' if q < 2*v else ('M' if q in Ms else 'h'))
prog = program(w)
# union-find on M-M gates
par = {q: q for q in M}
def find(a):
    while par[a] != a: par[a] = par[par[a]]; a = par[a]
    return a
for ph, i, e in prog:
    if e[0] == 'add' and e[1] in Ms and e[2] in Ms: par[find(e[1])] = find(e[2])
comps = defaultdict(list)
for q in M: comps[find(q)].append(q)
print('components', len(comps), 'sizes', Counter(len(c) for c in comps.values()),
      'entrance-dim signature', Counter(tuple(sorted(len(start[q]) for q in c)) for c in comps.values()))
# symbolic replay; contributions of M reads to targets
rows = [{q: F(1)} for q in range(n)]
deliv = defaultdict(lambda: defaultdict(dict))  # comp -> target -> {sym: coef} summed over target reads of members
treads = defaultdict(list); feeds = defaultdict(Counter); gph = defaultdict(Counter)
for t, (ph, i, e) in enumerate(prog):
    if e[0] == 'add':
        a, b = e[1], e[2]; k = F(e[3], e[4])
        fd = None if e[-1] is None else len(basis(e[-1]))
        if b in Ms and a not in Ms:
            g = find(b)
            if cls(a) == 'y':
                for s, x in rows[b].items():
                    deliv[g][a][s] = deliv[g][a].get(s, 0) + k*x
                treads[g].append((t, ph, a, b, str(k), fd))
            else: feeds[g]['read-by-' + cls(a) + ':' + ph] += 1
        if a in Ms and b not in Ms: feeds[find(a)]['writes-from-' + cls(b) + ':' + ph] += 1
        for q in (a, b):
            if q in Ms: gph[find(q)][(ph, fd)] += 1
    apply(e, rows)
# check restoration
restored = all(rows[q] == {q: F(1)} for q in M)
print('all M copies restored (identity row at end):', restored)
out = []
for g, mem in sorted(comps.items()):
    D = deliv[g]
    dsum = {}
    for y, d in D.items():
        xs = {s: c for s, c in d.items() if c != 0}
        dsum[y] = dict(xpart={s: str(c) for s, c in xs.items() if s < v}, other={s: str(c) for s, c in xs.items() if s >= v})
    out.append(dict(members=mem, entr_dims=[len(start[q]) for q in mem], entrances=[list(start[q]) for q in mem],
                    targets=sorted(D), target_reads=len(treads[g]), deliveries=dsum, feeds=dict(feeds[g]),
                    gates=Counter({f'{a}@{b}': x for (a, b), x in gph[g].items()})))
json.dump(out, open(sys.argv[2], 'w'), indent=0)
for o in out[:3]:
    print('\nGROUP members', o['members'], 'entr dims', o['entr_dims'])
    print('  entrances', o['entrances'])
    print('  targets', o['targets'], 'target reads', o['target_reads'])
    for y, d in o['deliveries'].items():
        print('   target', y, 'x-part nnz', len(d['xpart']), 'non-x leftover', d['other'])
    print('  feeds', o['feeds'])
    print('  gates', dict(o['gates']))
print('\nsummary: targets per group', Counter(len(o['targets']) for o in out), ' target reads per group', Counter(o['target_reads'] for o in out))
print('non-x leftover in deliveries (should be empty):', sum(1 for o in out for d in o['deliveries'].values() if d['other']))
print('x-part sizes per target delivery', Counter(len(d['xpart']) for o in out for d in o['deliveries'].values()))
