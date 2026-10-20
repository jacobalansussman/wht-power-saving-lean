"""v-mcopy-code-rules probe 3 (group 0): the cut's delivery to each group target, restricted to each super-packet,
per centre read; and the net packet content each target holds right after the cut and at the end of B."""
import sys, json
sys.path.insert(0, '../helpers-m-copies')
from wordlib import *
w = load(sys.argv[1]); v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
G = json.load(open('../helpers-m-copies/mc_groups.json')); allM = set(q for x in G for q in x['members'])
gi = int(sys.argv[2]); g = G[gi]; mem = set(g['members']); T = g['targets']
prog = program(w); rows = [{q: F(1)} for q in range(n)]; pk = {}
snap = {}
for ph, i, e in prog:
    if e[0] == 'add' and ph == 'B' and e[1] in mem and e[2] not in allM:
        xs = sorted(s for s in rows[e[2]] if s < v)
        if len(xs) == 8:
            j = len(set(pk.values()))
            for s in xs: pk[s] = (j, rows[e[2]][s])
    if e[0] == 'add' and ph == 'cut' and e[1] in T:
        d = defaultdict(set)
        for s, c in rows[e[2]].items():
            if s < v: d['all'].add(1)
        snap.setdefault(e[1], []).append((e[2], str(F(e[3], e[4])), {s: c for s, c in rows[e[2]].items() if s < v}))
    if ph == 'B' and 'afterCut' not in snap:
        snap['afterCut'] = {t: dict(rows[t]) for t in T}
    apply(e, rows)
    if ph == 'B': snap['endB'] = {t: dict(rows[t]) for t in T}
def bypk(row):
    out = defaultdict(set)
    for s, c in row.items():
        if s in pk: out[pk[s][0]].add(str(c / pk[s][1]))
    return {j: sorted(x) for j, x in sorted(out.items())}
for t in T:
    for c, q, row in snap[t]:
        print('target', t, 'reads centre', c, 'coef', q, 'centre content by packet (coef relative to packet):', bypk(row), ' own source in centre:', str(row.get(t - v, 0)))
    print('   target', t, 'after cut, by packet:', bypk(snap['afterCut'][t]))
    print('   target', t, 'end of B, by packet:', bypk(snap['endB'][t]))
