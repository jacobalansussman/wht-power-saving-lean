"""v-mcopy-code-rules probe 1: for group 0 targets, sum contributions of packet-source symbols by giver class and phase."""
import sys, json
sys.path.insert(0, '../helpers-m-copies')
from wordlib import *
w = load(sys.argv[1]); v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
G = json.load(open('../helpers-m-copies/mc_groups.json'))
g = G[0]; mem = set(g['members']); T = g['targets']
allM = set(q for x in G for q in x['members'])
prog = program(w)
print('precondition gates', len(w['precondition']), 'cleanup', len(w['cleanup']), 'early', len(w.get('early_cleanup', [])))
print('precondition regs classes', Counter(('x' if q < v else 'y' if q < 2*v else 'h') for e in w['precondition'] for q in (e[1:3] if e[0]=='add' else e[1:2])))
rows = [{q: F(1)} for q in range(n)]
contrib = {t: defaultdict(lambda: defaultdict(F)) for t in T}
for ph, i, e in prog:
    if e[0] == 'add' and e[1] in T:
        t, s = e[1], e[2]; k = F(e[3], e[4])
        key = (ph, 'Mgrp' if s in mem else 'M' if s in allM else 'x' if s < v else 'y' if s < 2*v else 'h')
        for sym, c in rows[s].items():
            if sym < v: contrib[t][key][sym] += k*c
    apply(e, rows)
for t in T:
    print('target', t, 'own source', t - v, 'final row', {k: str(c) for k, c in rows[t].items()})
    for key, d in sorted(contrib[t].items()):
        nz = {s: c for s, c in d.items() if c}
        print('   ', key, 'nnz', len(nz), 'own', str(d.get(t - v, 0)))
