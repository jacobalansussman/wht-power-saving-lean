#!/usr/bin/env python3
"""d-I1-idle: local re-implementation of the remainder derivation in hub code/audit_allcold.py lines 175-181
(one remainder child per (parent child, port, replica): family f = rep*8V + (2*pair+side)*V + port,
rank = parent rank - hosted cells, pending_offset = hosted) and the equality assert at line 181.
Control: the derivation reproduces the hosted #1042 idle_children exactly. Mutations: (A) merge a rank-48 remainder
with a rank-1 remainder of the same replica and port (different family) into one rank-49 child; (B) the same with the
rank-1 piece taken from a different replica. Both must fail the equality. Also checks the axis structure.
args: PARENT-TEMPLATE HOSTED-TEMPLATE DATA-PHASES OWNERS"""
import json, gzip, sys, copy
from collections import defaultdict
rd = lambda p: json.loads(gzip.decompress(open(p,'rb').read())) if p.endswith('.gz') else json.load(open(p))
parent, hosted, ph, owners = map(rd, sys.argv[1:5])
V, T = 112, 16
pairs, sides, lv = {}, {}, {}
for e in parent['execution_schedule']:
    if e['kind'] == 'idle_cover_sweep':
        for side, i in enumerate(e['child_indices']): pairs[i], sides[i], lv[i] = e['pair'], side, (e['from_level'], e['to_level'])
used = defaultdict(set)
for o in owners: used[o['child'], o['owner_replica'], o['port']] |= set(range(o['offset'], o['offset'] + o['width']))
wanted = []
for old in parent['idle_children']:
    i, port = old['child'], old['port']
    for rep in range(T):
        u = used[i, rep, port]; assert u == set(range(len(u))); r = old['rank'] - len(u)
        if r:
            f = rep*(8*V) + (2*pairs[i]+sides[i])*V + port
            wanted.append(dict(child=i, port=port, replica=rep, replicas=1, family=f, rank=r, hosted=len(u), pending_offset=len(u)))
print('control: derivation == hosted idle_children:', wanted == hosted['idle_children'], len(wanted))
def audit_eq(children): return children == wanted
# find a rank-48 and a rank-1 remainder at the same (replica, port)
by = defaultdict(list)
for k, z in enumerate(hosted['idle_children']): by[z['replica'], z['port']].append(k)
pick = None
for (rep, port), ks in sorted(by.items()):
    a = [k for k in ks if hosted['idle_children'][k]['rank'] == 48]; b = [k for k in ks if hosted['idle_children'][k]['rank'] == 1]
    if a and b: pick = (a[0], b[0]); break
ka, kb = pick; A, B = hosted['idle_children'][ka], hosted['idle_children'][kb]
print('merge pair (same replica+port):', {k: A[k] for k in ('child','replica','port','family','rank')}, '+', {k: B[k] for k in ('child','replica','port','family','rank')})
print('  families differ:', A['family'] != B['family'], ' sweep levels', lv[A['child']], lv[B['child']])
q = copy.deepcopy(hosted['idle_children']); m = dict(A, rank=49); q[ka] = m; del q[kb]
print('(A) merged rank-49 child accepted by line-181 equality?', audit_eq(q))
other = next(k for k, z in enumerate(hosted['idle_children']) if z['rank'] == 1 and z['replica'] != A['replica'])
q = copy.deepcopy(hosted['idle_children']); q[ka] = dict(A, rank=49); del q[other]
print('(B) wrong-replica merge accepted?', audit_eq(q))
# A rank-49 child cannot even be written in the remainder-frame export (audit_allcold lines 241-244):
axes = list(map(int, ph['recipes'][A['port']]['idle_children'][A['child']]['orthogonal_axes']))
print('  axes available to that child in its own family:', len(axes[A['pending_offset']:]), '(a rank-49 tail needs 49)')
# axis structure of the rank-48 children
def bits(v): return [k for k in range(73) if v >> k & 1]
for i in (8, 10, 12, 13):
    pat = set()
    for rec in ph['recipes']:
        ax = list(map(int, rec['idle_children'][i]['orthogonal_axes']))
        pat.add((tuple(tuple(bits(a)) for a in ax[:3]), tuple(sorted(b for a in ax[3:] for b in bits(a)))))
    print(f'child {i}: distinct axis patterns over 112 ports: {len(pat)};', 'example first-3', list(pat)[0][0], 'unit tail', list(pat)[0][1][:3], '...', list(pat)[0][1][-3:], 'len', len(list(pat)[0][1]))
# per family total idle rank by pair
tot = defaultdict(int)
for i in pairs: tot[pairs[i], sides[i]] += {c['child']: c['rank'] for c in parent['idle_children']}.get(i, 0)
print('total idle rank per family (pair, side):', dict(sorted(tot.items())))
