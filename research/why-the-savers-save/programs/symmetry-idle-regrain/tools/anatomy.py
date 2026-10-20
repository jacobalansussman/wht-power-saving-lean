#!/usr/bin/env python3
"""d-I1-idle: anatomy of #1042's idle children and idle remainders (read-only).
args: GLOBAL-TEMPLATE(lift, parent) HOSTED-TEMPLATE(global/hosted) DATA-PHASES PLAN OWNERS"""
import json, gzip, sys
from collections import Counter, defaultdict
rd = lambda p: json.loads(gzip.decompress(open(p,'rb').read())) if p.endswith('.gz') else json.load(open(p))
parent, hosted, ph, plan, owners = map(rd, sys.argv[1:6])
V = 112
prank = {c['child']: c['rank'] for c in parent['idle_children']}
# pair/side/levels per child from the schedule (as audit_allcold.py does)
info = {}
for e in parent['execution_schedule']:
    if e['kind'] == 'idle_cover_sweep':
        for side, i in enumerate(e['child_indices']):
            info[i] = dict(pair=e['pair'], side=side, from_level=e['from_level'], to_level=e['to_level'])
print('child table: child rank pair side from->to  8*(to-from)')
for i in sorted(prank):
    z = info[i]; print(f"  {i:2d} {prank[i]:3d}  p{z['pair']} s{z['side']}  {z['from_level']:2d}->{z['to_level']}  {8*(z['to_level']-z['from_level'])}")
# orthogonal axes structure
def bits(v): return [k for k in range(73) if v >> k & 1]
unit = Counter(); blocks = Counter(); overlap = 0
fam_axes = defaultdict(list)
for rec in ph['recipes']:
    port = rec['port']
    for i, ch in enumerate(rec['idle_children']):
        if i not in prank: continue
        ax = list(map(int, ch['orthogonal_axes']))
        assert len(ax) == ch['rank'] == prank[i]
        unit[i, all(bin(a).count('1') == 1 for a in ax)] += 1
        fam_axes[port, info[i]['pair'], info[i]['side']].append((i, ax))
print('unit-vector axes? (child, all unit):', sorted(unit.items())[:40])
# rank of union within a family (pair, side, port): are children of one family disjoint/independent?
def rank2(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)
fam_rank = Counter()
for (port, pair, side), L in fam_axes.items():
    allax = [a for i, ax in L for a in ax]
    fam_rank[pair, side, tuple(i for i, _ in L), len(allax), rank2(allax)] += 1
print('family (pair, side, children, sum ranks, rank of union): ports')
for k, n in sorted(fam_rank.items()): print(' ', k, n)
# 48 = six blocks of 8? show child 8 / 12 of port 0 as bit positions
rec0 = ph['recipes'][0]
for i in (8, 10, 12, 13):
    ax = list(map(int, rec0['idle_children'][i]['orthogonal_axes']))
    print(f'port {rec0["port"]} child {i} axes (first 12 as bit lists):', [bits(a) for a in ax[:12]])
# remainder pieces by child
cap = {}; fill = defaultdict(int)
for o in owners: fill[o['owner_replica'], o['child'], o['port']] += o['width']
rem = Counter(); fam_rem = defaultdict(dict)
for k, (rank, times) in plan['windows'].items():
    r, c = map(int, k.split(','))
    for p in range(V):
        left = rank - fill.get((r, c, p), 0)
        if left:
            rem[c, rank, left] += 1
            fam_rem[r, info[c]['pair'], info[c]['side'], p][c] = left
print('remainders (child, child rank, remainder): count')
for k, n in sorted(rem.items()): print(' ', k, n)
# do rank-1/3 remainders share a family (replica, pair, side, port) with a rank-48 remainder?
share = Counter(); samerep_port = Counter()
fullidle48 = {(r, c, p) for k, (rank, t) in plan['windows'].items() for r, c in [tuple(map(int, k.split(',')))] for p in range(V) if rank == 48 and fill.get((r, c, p), 0) == 0}
for fam, d in fam_rem.items():
    small = [c for c, x in d.items() if x in (1, 3)]
    big = [c for c, x in d.items() if x == 48]
    if small and big: share['same family has 1/3 and 48'] += 1
for (r, c, p) in fullidle48:
    for c2 in prank:
        x = fam_rem.get((r, info[c2]['pair'], info[c2]['side'], p), {}).get(c2)
    # same replica and port (any family)
for r in range(16):
    for p in range(V):
        ks = {}
        for c in prank:
            x = fam_rem.get((r, info[c]['pair'], info[c]['side'], p), {}).get(c)
            if x: ks[c] = x
        has48 = any(x == 48 for x in ks.values()); has13 = [c for c, x in ks.items() if x in (1, 3)]
        if has48 and has13: samerep_port[tuple(sorted(c for c, x in ks.items() if x == 48)), tuple(sorted(has13))] += 1
print('same family (replica,pair,side,port) holding a 1/3 remainder AND a 48 remainder:', dict(share) or 0)
print('same replica+port, different family: (48-children, 1/3-children): count  [top 12]')
for k, n in samerep_port.most_common(12): print(' ', k, n)
print('total (replica,port) with both:', sum(samerep_port.values()))
# hosted remainder children as the audit writes them
hc = hosted['idle_children']
print('hosted idle remainder children:', len(hc), 'keys', sorted(hc[0]))
print('hosted remainder ranks:', dict(sorted(Counter(z['rank'] for z in hc).items())))
print('family formula check f = rep*8V + (2*pair+side)*V + port:', all(z['family'] == z['replica']*8*V + (2*info[z['child']]['pair']+info[z['child']]['side'])*V + z['port'] for z in hc))
