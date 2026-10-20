"""d-mcopy-code: exact spec of one M group as a small case.
Replaces each packet load 'M += holder' by a fresh packet symbol pi_j and replays; returns
 C (targets x packets) as delivered at the end, holders, packet frames, target chains, pair frames.
Usage: python3 -I group_spec.py FULL-WORD.json.gz GROUP_INDEX"""
import sys, json
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, '.')
from wordlib import *
w = load(sys.argv[1]); v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
gi = int(sys.argv[2]); G = json.load(open('mc_groups.json'))[gi]
mem = set(G['members']); T = sorted(G['targets'])
start, end, chain, paid, runs = frames(w)
prog = program(w)
rows = [{q: F(1)} for q in range(n)]
loads = []
for t, (ph, i, e) in enumerate(prog):
    if e[0] == 'add' and ph == 'B' and e[1] in mem and e[2] not in mem:
        # packet load
        pid = ('pi', len(loads)); loads.append(dict(t=t, i=i, M=e[1], holder=e[2], k=str(F(e[3], e[4])), frame=list(basis(e[-1])),
                                                    holder_val={s: str(x) for s, x in rows[e[2]].items()}))
        rows[e[1]][pid] = rows[e[1]].get(pid, 0) + F(e[3], e[4])
        continue
    apply(e, rows)
C = [[str(rows[y].get(('pi', j), 0)) for j in range(len(loads))] for y in T]
ports = w['ports']
# x-content of each packet & block membership
blocks = None
print('group', gi, 'members', sorted(mem), 'targets', T, 'ports', [ports[y - v] for y in T])
print('C (targets x packets, as delivered at the end, after target butterflies):')
for y, r in zip(T, C): print('  ', y, r)
for j, L in enumerate(loads):
    xs = sorted(int(s) for s in L['holder_val'] if int(s) < v)
    print(' packet', j, 'M', L['M'], 'holder', L['holder'], 'cls', 'x' if L['holder'] < v else 'h', 'frame', L['frame'], 'sources', xs, 'labels', [ports[s] for s in xs])
    print('     holder chain', [(str(tag), len(U)) for tag, U in chain[L['holder']]][:12])
for y in T:
    print(' target', y, 'end', end[y], 'chain', [(str(tag), list(U)) for tag, U in chain[y]])
json.dump(dict(group=gi, members=sorted(mem), targets=T, C=C, loads=loads,
               target_chain={y: [[str(tag), list(U)] for tag, U in chain[y]] for y in T},
               holder_chain={L['holder']: [[str(tag), list(U)] for tag, U in chain[L['holder']]] for L in loads},
               member_chain={q: [[str(tag), list(U)] for tag, U in chain[q]] for q in mem},
               end={y: list(end[y]) for y in T}), open(f'group{gi}_spec.json', 'w'), indent=0)
