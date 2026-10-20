"""v-mcopy-code-reprice: where do the 64 packet sources of group G enter/leave the group's targets?
Tracks x-content of the group targets through the audit program order, by gate phase and partner class."""
import sys, json
sys.path.insert(0, '.')
from wordlib import *
from fractions import Fraction as F
from collections import Counter, defaultdict
w = load(sys.argv[1]); gi = int(sys.argv[2]); v = w['report']['v']; n = 2*v + w['report']['R']
Gs = json.load(open('mc_groups.json')); allM = set(q for g in Gs for q in g['members'])
A = json.load(open('all_groups.json'))[gi]; T = A['targets']; mem = set(A['members'])
prog = program(w)
rows = [{q: F(1)} for q in range(n)]; pk = set()
for t, (ph, i, e) in enumerate(prog):
    if e[0] == 'add' and ph == 'B' and e[1] in mem and e[2] not in allM:
        xs = [s for s in rows[e[2]] if s < v]
        if len(xs) == 8: pk |= set(xs)
    apply(e, rows)
print('group', gi, 'targets', T, 'packet sources', len(pk), 'own sources of targets', [t - v for t in T], 'in packets?', [t - v in pk for t in T])
cls = lambda q: 'x' if q < v else ('y' if q < 2*v else ('M' if q in mem else ('Mo' if q in allM else 'h')))
rows = [{q: F(1)} for q in range(n)]
flow = defaultdict(Counter)
for t, (ph, i, e) in enumerate(prog):
    if e[0] == 'add' and e[1] in T:
        k = F(e[3], e[4]); src = e[2]
        d = {s: k * x for s, x in rows[src].items() if s in pk}
        if d: flow[e[1]][(ph, cls(src))] += len(d)
    apply(e, rows)
    if ph == 'cut' and (t + 1 == len(prog) or prog[t+1][0] != 'cut'):
        for y in T: print(' after cut: target', y, 'packet-source nnz', sum(1 for s in rows[y] if s in pk))
for y in T:
    print(' target', y, 'entries of packet sources by (phase, partner class):', dict(flow[y]))
    print('   final x-part', {s: str(c) for s, c in rows[y].items() if s < v})
