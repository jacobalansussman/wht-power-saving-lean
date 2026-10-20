"""d-mcopy-code STEP 2a: giver inventory for one group. For every non-M register g and every frame G of g's (fixed) chain
with G inside U_t for some group target t, record: dim, time window [arrive, leave), and the max number of the group's
needed source-coefficients (over the 8 packets) that g's value carries inside that window.  Needed content of the group:
the 64 sources of its 8 packets.  Usage: python3 -I givers.py FULL-WORD.json.gz GROUP"""
import sys, json
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, '.')
from wordlib import *
w = load(sys.argv[1]); v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
gi = int(sys.argv[2])
A = json.load(open('all_groups.json'))[gi]
Gs = json.load(open('mc_groups.json')); allM = set(q for g in Gs for q in g['members'])
start, end, chain, paid, runs = frames(w)
T = A['targets']; Ut = {t: tuple(end[t]) for t in T}
prog = program(w)
# time of each chain step in program index
tagtime = {}
for t, (ph, i, e) in enumerate(prog): tagtime[(ph, i)] = t
def ttime(tag):
    if tag is None: return -1
    ph, i = tag
    if ph == 'final': return 10**9
    return tagtime[(ph, i)]
pk = [p for p in A['packets'] if p['nsrc'] == 8]
need = set()
rows = [{q: F(1)} for q in range(n)]
# packet source sets
grp_src = []
for p in pk: pass
# recompute packet sources from holder value at load time
loads_t = {}
for t, (ph, i, e) in enumerate(prog):
    if e[0] == 'add' and ph == 'B' and e[1] in set(A['members']) and e[2] not in allM:
        loads_t[e[1]] = (t, sorted(s for s in rows[e[2]] if s < v))
    apply(e, rows)
for q, (t, xs) in loads_t.items():
    if len(xs) == 8: need |= set(xs)
print('group', gi, 'targets', T, 'needed sources', len(need))
# windows
win = []
for g in range(n):
    if g in allM: continue
    ch = chain[g]
    for k, (tag, U) in enumerate(ch):
        U = tuple(U)
        if not U: continue
        ins = [t for t in T if contained(U, Ut[t])]
        if not ins: continue
        a = ttime(tag); b = ttime(ch[k+1][0]) if k + 1 < len(ch) else 10**9
        win.append((g, k, len(U), a, b, ins))
# replay to get values in windows
rows = [{q: F(1)} for q in range(n)]
best = defaultdict(int); nonempty = defaultdict(set)
open_w = defaultdict(list)
for (g, k, d, a, b, ins) in win: open_w[g].append((a, b, k, d, ins))
for t, (ph, i, e) in enumerate(prog):
    apply(e, rows)
    for q in (e[1:3] if e[0] == 'add' else e[1:2]):
        for (a, b, k, d, ins) in open_w.get(q, []):
            if a <= t < b:
                c = len([s for s in rows[q] if s in need])
                if c > best[(q, k)]: best[(q, k)] = c
cls = lambda q: 'x' if q < v else ('y' if q < 2*v else 'h')
tab = Counter()
rec = []
for (g, k, d, a, b, ins) in win:
    c = best[(g, k)]
    if c:
        tab[(cls(g), d, c)] += 1
        rec.append(dict(g=g, cls=cls(g), dim=d, window=[a, b], inside_U_of=ins, max_needed_coeffs=c))
print('givers with needed content (class, frame dim, max #needed source coeffs held): count')
for k2, x in sorted(tab.items()): print('  ', k2, x)
json.dump(rec, open(f'givers_g{gi}.json', 'w'), indent=0)
