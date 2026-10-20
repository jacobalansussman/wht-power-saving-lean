"""Global re-listing: every legal order is a linear extension of the dependency DAG (pairs of gates on a common
register that do not commute, plus the cut for gates that cannot cross it).  Kahn list-scheduling with a priority
(lower current frame first, etc.), then exact global frame MILP.  A cheap 'extraction' over the whole trace class."""
import sys, os, heapq, time
sys.path.insert(0, '.')
import fw, search, ilp, reorder2
from ilp_wide import all_between
W = search.Word(fw.load(os.environ.get('WORD', 'K8-FULL-WORD.json.gz'))); O = reorder2.OW2(W)
N = O.N; CUT = O.CUT
order0 = list(O.order)
# edges
succ = {g: set() for g in order0}; indeg = {g: 0 for g in order0}
pos0 = {g: i for i, g in enumerate(order0)}
for r, gs in O.chain.items():
    for i in range(len(gs)):
        for j in range(i + 1, len(gs)):
            a, b = gs[i], gs[j]
            if not O.commute(a, b) and b not in succ[a]: succ[a].add(b)
for g in range(N):
    if O.sec0[g] == 'T': continue
    if not O.cross_ok[g]:
        if pos0[g] < pos0[CUT]: succ[g].add(CUT)
        else: succ[CUT].add(g)
# tail gates after the cut and after every A/B gate that shares a register (already via chains) -- also keep tail after CUT
for g in range(N):
    if O.sec0[g] == 'T': succ[CUT].add(g)
for a in succ:
    for b in succ[a]: indeg[b] += 1
def schedule(prio):
    ind = dict(indeg); heap = [(prio(g), g) for g in order0 if ind[g] == 0]; heapq.heapify(heap); out = []
    while heap:
        _, g = heapq.heappop(heap); out.append(g)
        for b in succ[g]:
            ind[b] -= 1
            if ind[b] == 0: heapq.heappush(heap, (prio(b), b))
    assert len(out) == len(order0), 'cycle'
    return out
def evaluate(order, label):
    O.order = order; O.refresh()
    if not O.feasible:
        print(label, 'infeasible (some P > Q)', sum(1 for g in O.order if not search.le(O.P[g], O.Q[g])), flush=True); return None
    # build a word with these gates/order and current frames clipped, then exact MILP
    for g in O.order:
        if g != CUT and not (search.le(O.P[g], O.f[g]) and search.le(O.f[g], O.Q[g])): O.f[g] = O.P[g]
    w2 = O.to_word() if all(search.le(O.P[g], O.f[g]) and search.le(O.f[g], O.Q[g]) for g in O.order) else None
    W2 = search.Word(w2)
    P, Q = ilp.bounds_PQ(W2)
    extra = {k: all_between(P[k], Q[k]) for k in range(len(W2.seq)) if P[k] != Q[k] and len(Q[k]) - len(P[k]) <= 3}
    fr, info = ilp.solve(W2, extra_pool=extra, verbose=False)
    if fr is None: print(label, 'MILP failed', info); return None
    W2.f = fr; c = W2.total_cost()
    moved = sum(1 for i, g in enumerate(order) if g != order0[i])
    print(label, 'gates displaced', moved, 'optimal cost', round(c, 4), 'delta', round(c - 20241.35972374547, 4), dict(sorted(W2.profile().items())), flush=True)
    return W2, c
f0 = {g: O.f[g] for g in range(N)}
fd = lambda g: len(f0[g]) if g != CUT else 99
res = {}
res['original'] = evaluate(schedule(lambda g: pos0[g]), 'original-priority')
O.f = dict(f0); O.f[CUT] = ()
res['dim'] = evaluate(schedule(lambda g: (fd(g), pos0[g])), 'lowest-frame-first')
O.f = dict(f0); O.f[CUT] = ()
res['dimrev'] = evaluate(schedule(lambda g: (-pos0[g],)), 'reverse-index (ALAP-ish)')
O.f = dict(f0); O.f[CUT] = ()
res['asap'] = evaluate(schedule(lambda g: (0, pos0[g])), 'asap')

import json, random
d = json.load(open('lbdep.json')); Pd = {int(k): len(v) for k, v in d['Pd'].items()}; Qd = {int(k): len(v) for k, v in d['Qd'].items()}
for label, pr in (('Pd-first', lambda g: (Pd.get(g, 99), pos0[g])), ('Qd-first', lambda g: (Qd.get(g, 99), pos0[g])),
                  ('Pd+Qd', lambda g: (Pd.get(g, 99) + Qd.get(g, 99), pos0[g])), ('frame-then-Qd', lambda g: (fd(g), Qd.get(g, 99), pos0[g])),
                  ('frame-then-rev', lambda g: (fd(g), -pos0[g]))):
    O.f = dict(f0); O.f[CUT] = ()
    evaluate(schedule(pr), label)
for s in range(4):
    rnd = random.Random(s); key = {g: rnd.random() for g in order0}
    O.f = dict(f0); O.f[CUT] = ()
    evaluate(schedule(lambda g: (fd(g), key[g])), f'frame-then-random{s}')
