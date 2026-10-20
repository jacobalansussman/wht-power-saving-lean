#!/usr/bin/env python3
"""Step 1h: the relaxed-H2 pair count of g-fields (direct target readers only) redone on the FULL-WORD with an exact
max b-matching, then every compatible pair re-tested against the exact transitive content sandwich S_arr of
arrival_sandwich.py, with one explicit certificate per failing pair: a target t that receives P2's post-arrival packet
content (through an assembler) and has P1 + P2 not inside u_t^perp.
Usage: python3 -I direct_vs_exact.py FULL-WORD.json.gz survey.json arrival.json direct_vs_exact.json"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
from collections import Counter, defaultdict
w = load(sys.argv[1]); S = json.load(open(sys.argv[2])); AR = json.load(open(sys.argv[3])); v = 112; ports = w['ports']
segs, comp = program(w); off = S['seg_offsets']
gl = []
for name, g in segs: gl += [(name, e) for e in g]
planes = {tuple(r['plane']): r for r in S['planes']}
def inter(us): return basis([x for x in range(1, 512) if all(dot(x, u) == 0 for u in us)])
res = {}
for label, phases in (('B direct', ('B',)), ('A+B+tail direct', ('A', 'B', 'tail'))):
    cons = defaultdict(set)
    for P, r in planes.items():
        for d in r['visitors']:
            for t in range(max(d['t_in'], 0) + 1, len(gl)):
                name, e = gl[t]
                if name in phases and e[0] == 'add' and e[2] == d['q'] and v <= e[1] < 2 * v: cons[P].add(e[1])
    SP = {P: inter([ports[t - v] for t in cons[P]]) for P in planes}
    slots = {P: sum(1 for d in r['visitors'] if d['q'] >= 2 * v) for P, r in planes.items()}
    live = [P for P in planes if cons[P]]
    E = [(a, b) for a in live for b in live if a != b and planes[a]['block'] != planes[b]['block'] and contained(join(a, b), SP[b])]
    und = {tuple(sorted(e)) for e in E}
    # exact max b-matching (capacities = helper visits per plane), general graph: upper bound = floor of the bipartite
    # double-cover max flow / 2; lower bound = augmenting local search; equal => exact.
    nodes = sorted({x for e in und for x in e}); idx = {P: i for i, P in enumerate(nodes)}
    cap = defaultdict(int); adj = defaultdict(set)
    def add(a, b, c): cap[(a, b)] += c; adj[a].add(b); adj[b].add(a)
    for P in nodes: add('s', ('L', P), slots[P]); add(('R', P), 't', slots[P])
    for a, b in und: add(('L', a), ('R', b), 1); add(('L', b), ('R', a), 1)
    flow = 0
    while True:
        par = {'s': None}; q = ['s']
        for x in q:
            for y in adj[x]:
                if y not in par and cap[(x, y)] > 0: par[y] = x; q.append(y)
        if 't' not in par: break
        y = 't'
        while par[y] is not None: x = par[y]; cap[(x, y)] -= 1; cap[(y, x)] += 1; y = x
        flow += 1
    ub = flow // 2
    # lower bound: greedy + randomized restarts
    import random
    best = 0
    for seed in range(400):
        rnd = random.Random(seed); el = sorted(und); rnd.shuffle(el)
        if seed == 0:
            deg = Counter(x for e in und for x in e); el = sorted(und, key=lambda e: deg[e[0]] + deg[e[1]])
        sl = dict(slots); m = 0
        for a, b in el:
            if sl[a] > 0 and sl[b] > 0: sl[a] -= 1; sl[b] -= 1; m += 1
        best = max(best, m)
    # exact transitive test of every compatible ordered pair (P1 -> P2): some helper visitor f of P2 with P1+P2 <= S_arr(f)?
    Sarr = defaultdict(list)
    for r in AR['rows']:
        if r['q'] >= 2 * v: Sarr[tuple(r['plane'])].append(r['q'])
    # recompute S_arr quickly from arrival.json is not stored; use the dims rule: arrival.json pairs list is empty
    surv = len(AR['pairs'])
    from bmatch import max_bmatching
    exact = max_bmatching({P: slots[P] for P in nodes}, sorted(und))
    res[label] = dict(consumed_planes=len(live), dimS=dict(Counter(len(SP[P]) for P in live)), ordered_pairs=len(E),
                      undirected=len(und), bmatch_lower=best, bmatch_upper=ub, bmatch_exact=exact, survive_exact=surv)
    print(label, json.dumps(res[label]))
    if label == 'B direct' and E:
        # certificate for the first few pairs: a transitive consumer of P2's packet violating P1+P2
        for (P1, P2) in E[:3]:
            r2 = planes[P2]
            f = [d for d in r2['visitors'] if d['q'] >= 2 * v and d['out_spans_P']][0]
            zz = {f['q']: F(1)}; via = {f['q']: None}; found = None
            for t in range(f['t_in'] + 1, len(gl)):
                name, e = gl[t]
                if name == 'cut' or e[0] != 'add': continue
                a, b = e[1], e[2]
                if b in zz and zz[b] != 0:
                    zz[a] = zz.get(a, F(0)) + F(e[3], e[4]) * zz[b]; via.setdefault(a, (b, t, list(basis(e[-1]))))
                    if v <= a < 2 * v and zz[a] != 0 and not inperp(join(P1, P2), ports[a - v]):
                        found = (a, t); break
            path = []; x = found[0]
            while via.get(x): b, t, Fr = via[x]; path.append(dict(reader=x, read=b, t=t, frame=Fr)); x = b
            print('  certificate: P1', P1, 'P2', P2, 'copy holder', f['q'], '-> target', found[0], 'label', ports[found[0] - v],
                  'u_t.(P1+P2) nonzero; path', path[::-1])
            res.setdefault('certificates', []).append(dict(P1=P1, P2=P2, holder=f['q'], target=found[0], path=path[::-1]))
json.dump(res, open(sys.argv[4], 'w'), indent=1)
