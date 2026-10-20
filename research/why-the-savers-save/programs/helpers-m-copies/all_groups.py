"""d-mcopy-code STEP 1c: for all 14 M groups: packet loads, F_j frames, visitors of each F_j, delivery matrix C pattern,
reader pairs, target-pair frames W, whether F_j lies in W, entrance frames, timing windows.
Usage: python3 -I all_groups.py FULL-WORD.json.gz all_groups.json"""
import sys, json
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, '.')
from wordlib import *
w = load(sys.argv[1]); v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
start, end, chain, paid, runs = frames(w)
Gs = json.load(open('mc_groups.json'))
allM = set(q for g in Gs for q in g['members'])
cls = lambda q: 'M' if q in allM else ('x' if q < v else ('y' if q < 2*v else 'h'))
byframe = defaultdict(list)
for r in range(n):
    for tag, U in chain[r]: byframe[tuple(U)].append(r)
def span(b):
    s = {0}
    for x in b: s |= {y ^ x for y in s}
    return s
def inter(a, b):
    A = span(a); return basis([x for x in span(b) if x in A])
prog = program(w)
memb = {q: gi for gi, g in enumerate(Gs) for q in g['members']}
rows = [{q: F(1)} for q in range(n)]
loads = defaultdict(list)
for t, (ph, i, e) in enumerate(prog):
    if e[0] == 'add' and ph == 'B' and e[1] in allM and e[2] not in allM:
        gi = memb[e[1]]; j = len(loads[gi]); pid = ('pi', gi, j)
        loads[gi].append(dict(t=t, M=e[1], holder=e[2], frame=basis(e[-1]), k=F(e[3], e[4]),
                              xs=sorted(s for s in rows[e[2]] if s < v), other=sorted(s for s in rows[e[2]] if s >= v)))
        rows[e[1]][pid] = rows[e[1]].get(pid, 0) + F(e[3], e[4]); continue
    apply(e, rows)
out = []; summary = Counter()
xs_used = Counter(); hsrc = Counter()
for gi, g in enumerate(Gs):
    T = sorted(g['targets'])
    C = [[rows[y].get(('pi', gi, j), F(0)) for j in range(len(loads[gi]))] for y in T]
    # readers of each packet
    rd = [tuple(T[a] for a in range(4) if C[a][j] != 0) for j in range(len(loads[gi]))]
    Wp = {}
    for a in range(4):
        for b in range(a + 1, 4):
            Wp[(T[a], T[b])] = inter(end[T[a]], end[T[b]])
    # target pairs = pairs whose intersection is in both chains
    tp = [p for p, W in Wp.items() if any(tuple(U) == W for tag, U in chain[p[0]]) and any(tuple(U) == W for tag, U in chain[p[1]])]
    rec = dict(group=gi, targets=T, target_pairs=tp, members=g['members'], entr_dims=g['entr_dims'],
               entrances=sorted(set(tuple(start[q]) for q in g['members'])))
    pk = []
    for j, L in enumerate(loads[gi]):
        Fj = L['frame']; vis = byframe[Fj]
        inW = [p for p in tp if contained(Fj, Wp[p])]
        inU = [t for t in T if contained(Fj, end[t])]
        pk.append(dict(M=L['M'], holder=L['holder'], holder_cls=cls(L['holder']), frame=list(Fj), nsrc=len(L['xs']), other=L['other'],
                       visitors=[(r, cls(r)) for r in vis], readers=rd[j], coefs=[str(C[a][j]) for a in range(4)], inW=inW, inU=inU))
        summary['visitors=' + ''.join(sorted(cls(r) for r in vis))] += 1
        summary['readers=%d' % len(rd[j])] += 1
        summary['readers_in_one_pair=%s' % any(set(rd[j]) == set(p) for p in tp)] += 1
        summary['F_in_some_W=%s' % bool(inW)] += 1
        summary['F_in_U_exactly_readers=%s' % (set(inU) == set(rd[j]))] += 1
        for s in L['xs']: xs_used[s] += 1
        for r in vis:
            if r < v: hsrc[r] += 1
    rec['packets'] = pk
    rec['C_pattern'] = [[('+' if c > 0 else '-' if c < 0 else '0') for c in row] for row in C]
    out.append(rec)
    summary['target_pairs=%d' % len(tp)] += 1
    summary['entr=' + str(sorted(Counter(g['entr_dims']).items()))] += 1
print(json.dumps(dict(summary), indent=0))
print('sources visiting some F_j (as holder/partner):', len(hsrc), ' max visits per source', max(hsrc.values()))
print('distinct F frames', len(set(tuple(p['frame']) for o in out for p in o['packets'])))
json.dump(out, open(sys.argv[2], 'w'), indent=0, default=str)
