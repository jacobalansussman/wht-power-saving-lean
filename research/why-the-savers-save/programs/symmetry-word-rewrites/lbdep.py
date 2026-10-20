"""Optimistic reorder potential: intervals [Pd,Qd] from the dependency DAG only (frames propagate only along
non-commuting gate pairs, which every legal reorder keeps in order), then the per-register decoupled chain DP in the
current order of each register's gates.  Not a proof of anything; a gauge of how much reordering could possibly buy."""
import sys; sys.path.insert(0,'.')
import fw, search, reorder2
from search import le, join, meet, FC
from ilp_wide import all_between
from collections import Counter
W = search.Word(fw.load('K8-FULL-WORD.json.gz')); O = reorder2.OW2(W)
order = O.order; N = O.N
# dependency predecessors: for gate g, for each register r in g, the latest earlier gate on r that does NOT commute with g
# (transitively this captures all non-commuting predecessors that share r; non-commuting pairs always share a register)
Pd = {}; Qd = {}
last_on = {}  # register -> list of gates so far
for g in order:
    if g == O.CUT: continue
    U = ()
    for r in O.regs[g]:
        U = join(U, W.start[r])
        for h in reversed(last_on.get(r, [])):
            if not O.commute(g, h): U = join(U, Pd[h]); break
            # commuting h: skip, but a non-commuting one further back may exist
    # also cut ordering: B-gates that cannot cross the cut depend on A side? (ignored: optimistic)
    if g in O.frozen: U = join(U, O.frozen[g])
    Pd[g] = U
    for r in O.regs[g]: last_on.setdefault(r, []).append(g)
next_on = {}
for g in reversed(order):
    if g == O.CUT: continue
    U = W.full
    for r in O.regs[g]:
        U = meet(U, W.end[r])
        if O.sec[g] != 'T': U = meet(U, W.ub[r])
        if O.sec[g] == 'A' and r in W.cframe: U = meet(U, W.cframe[r])
        for h in reversed(next_on.get(r, [])):
            if not O.commute(g, h): U = meet(U, Qd[h]); break
    if g in O.frozen: U = meet(U, O.frozen[g])
    Qd[g] = U
    for r in O.regs[g]: next_on.setdefault(r, []).append(g)
bad = sum(1 for g in Pd if not le(Pd[g], Qd[g])); print('Pd<=Qd violations', bad)
tot = 0.0; act = 0.0; gaps = Counter()
for r in range(W.n):
    states = {W.start[r]: 0.0}
    for k in O.chain.get(r, []):
        gap = len(Qd[k]) - len(Pd[k])
        cands = all_between(Pd[k], Qd[k]) if gap <= 3 else {Pd[k], Qd[k], O.f[k]}
        gaps[gap] += 1
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
    tot += min(c + (FC[len(E)-len(V)] if V != E else 0.0) for V, c in states.items() if le(V, E))
    act += O.reg_cost(r, O.f)
print('interval gaps', sorted(gaps.items()))
print('dependency-only decoupled estimate', round(tot, 3), 'actual', round(act, 3), 'potential', round(act - tot, 3))

# per-register potential and phase breakdown
per = {}
for r in range(W.n):
    states = {W.start[r]: (0.0, None)}
    for k in O.chain.get(r, []):
        gap = len(Qd[k]) - len(Pd[k])
        cands = all_between(Pd[k], Qd[k]) if gap <= 3 else {Pd[k], Qd[k], O.f[k]}
        new = {}
        for U in cands:
            best = None
            for V, (c, _) in states.items():
                if le(V, U):
                    cc = c + (FC[len(U)-len(V)] if U != V else 0.0)
                    if best is None or cc < best: best = cc
            if best is not None: new[U] = (best, None)
        states = new
    E = W.end[r]
    lbr = min(c + (FC[len(E)-len(V)] if V != E else 0.0) for V, (c, _) in states.items() if le(V, E))
    per[r] = O.reg_cost(r, O.f) - lbr
role = lambda r: 'src' if r < W.v else ('tgt' if r < 2*W.v else 'hlp')
from collections import defaultdict
agg = defaultdict(float)
for r, p in per.items(): agg[role(r)] += p
print({k: round(x, 2) for k, x in agg.items()})
top = sorted(per.items(), key=lambda x: -x[1])[:15]
print([(r, role(r), round(p, 2)) for r, p in top])
import json; json.dump({'Pd': {str(g): list(U) for g, U in Pd.items()}, 'Qd': {str(g): list(U) for g, U in Qd.items()}, 'per': per}, open('lbdep.json', 'w'))
