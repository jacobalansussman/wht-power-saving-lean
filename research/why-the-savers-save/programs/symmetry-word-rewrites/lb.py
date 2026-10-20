"""Per-register decoupled lower bound for the fixed order: each register alone picks a monotone chain with f(k) in
[P(k),Q(k)]; DP over positions with frames from the full interval (all gaps <= 3 here).  Sum <= any coupled cost."""
import sys; sys.path.insert(0,'.')
import fw, search, ilp
from search import le, join, meet, FC
from collections import Counter
W = search.Word(fw.load('K8-FULL-WORD.json.gz'))
P, Q = ilp.bounds_PQ(W)
from ilp_wide import all_between
tot = 0.0; act = 0.0; slack = Counter(); per = {}
for r in range(W.n):
    ps = W.pos.get(r, [])
    # DP: states = frame value at current position
    states = {W.start[r]: 0.0}
    for k in ps:
        cands = all_between(P[k], Q[k]) if P[k] != Q[k] else {P[k]}
        new = {}
        for U in cands:
            best = None
            for V, c in states.items():
                if le(V, U):
                    cc = c + (FC[len(U)-len(V)] if U != V else 0.0)
                    if best is None or cc < best: best = cc
            if best is not None: new[U] = best
        states = new
    E = W.end[r]
    lbr = min(c + (FC[len(E)-len(V)] if V != E else 0.0) for V, c in states.items() if le(V, E))
    g = lambda k: W.f[k]
    a = W.chain_cost(r, g)
    tot += lbr; act += a
    role = 'src' if r < W.v else ('tgt' if r < 2*W.v else 'hlp')
    slack[role] += a - lbr
    per[r] = a - lbr
print('decoupled LB', round(tot,3), 'actual', round(act,3), 'gap', round(act-tot,3))
print({k: round(v,3) for k,v in slack.items()})
print('registers with slack', sum(1 for x in per.values() if x > 1e-9))

import json
top = sorted(per.items(), key=lambda x: -x[1])[:12]
print('top slack', [(r, round(s,3)) for r,s in top])
from collections import Counter as C2
print('slack values', C2(round(x,3) for x in per.values() if x > 1e-9))
json.dump({str(r): s for r, s in per.items() if s > 1e-9}, open('lb_slack.json','w'))
