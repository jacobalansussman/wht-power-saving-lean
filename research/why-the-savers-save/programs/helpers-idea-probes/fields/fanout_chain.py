#!/usr/bin/env python3
"""g-fields idea I1 (Steiner fan-out / Dilworth): can ONE copy register serve all B-phase target consumers of a packet plane
in a single nested chain? Necessary frame condition (time order ignored, i.e. assuming a legal re-timing exists):
order the consumers t_1..t_k so that every earlier consumer's actual gate frame f_l lies inside u_{t_j}^perp of every
later consumer (the copy carries f_1+...+f_{j-1} as fill into t_j's frame, which must stay inside u_{t_j}^perp).
Exact over all orders (k <= 8). Also reports the same for the dirt reads (target reads of the copies before the packet).
Usage: python3 -I fanout_chain.py   (needs timeline.json, consumers.json)"""
import json, itertools
from collections import Counter
T = json.load(open('timeline.json')); R = json.load(open('consumers.json'))
ev = {int(q): E for q, E in T['ev'].items()}; cls = {int(q): k for q, k in T['cls'].items()}
dims = T['dims']; Bf = T['frames']; ports = T['ports']; v = 112
dot = lambda a, b: bin(a & b).count('1') & 1
def perp_ok(f, t):   # frame f inside u_t^perp
    return all(not dot(x, ports[t - v]) for x in Bf[f])
out = Counter(); ex = []
for r in R:
    cons = [(t, f, q, p) for t, f, q, p in [(c[0], None, None, c[3]) for c in r['consumers']]]
    # recover actual gate frames from the timeline (consumers.json stored dims only)
    holders = [q for q, k, tv in r['visitors']]
    E = []
    for q, k, tv in r['visitors']:
        for (t, ph, role, p, f) in ev[q]:
            if t > tv and role == 'r' and cls[p] == 'y' and 11003 <= t < 12921: E.append((t, f, q, p))
    if not E: out['no B consumers'] += 1; continue
    tg = sorted({p for _, _, _, p in E})
    fr = {}
    for t, f, q, p in E: fr.setdefault(p, []).append(f)
    k = len(tg); best = None
    if k > 8: out['k>8'] += 1; continue
    ok = False
    for order in itertools.permutations(tg):
        if all(perp_ok(f, order[j]) for j in range(1, k) for l in range(j) for f in fr[order[l]]):
            ok = True; break
    out[(k, len(holders and {q for _, _, q, _ in E}), 'one chain feasible' if ok else 'needs >1')] += 1
    if ok and len(ex) < 3: ex.append(dict(plane=r['plane'], targets=tg, order=order, frames={p: [dims[f] for f in fr[p]] for p in tg}))
for k, x in sorted(out.items(), key=str): print(x, k)
print(json.dumps(ex))
