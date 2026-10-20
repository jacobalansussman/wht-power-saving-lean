"""Rewrite 'sink to cleanup': a main-section gate whose two registers both end full and which commutes with every
later main/tail gate is moved to the start of the (free, unclimbed) cleanup.  Frames re-optimised by the exact
fixed-order MILP on the new word; the result is written and audited."""
import sys, copy, json, subprocess
sys.path.insert(0, '.')
import fw, search, ilp, reorder2
from ilp_wide import all_between
def sink(w, gids_A):
    """gids_A: indices into w['A'] to sink.  Returns a new word (frames not yet re-optimised)."""
    w2 = copy.deepcopy(w)
    full = list(fw.basis([1 << i for i in range(w['report']['h'])]))
    moved = [w['A'][i][:-1] + [full] for i in sorted(gids_A)]
    w2['A'] = [e for i, e in enumerate(w['A']) if i not in set(gids_A)]
    m = len(moved)
    w2['cleanup'] = moved + w['cleanup']
    w2['early_cleanup'] = [[i + m, U] for i, U in w['early_cleanup']]
    return w2
def reframe(w2):
    W = search.Word(w2)
    P, Q = ilp.bounds_PQ(W)
    extra = {k: all_between(P[k], Q[k]) for k in range(len(W.seq)) if P[k] != Q[k] and len(Q[k]) - len(P[k]) <= 3}
    fr, info = ilp.solve(W, extra_pool=extra, verbose=False)
    if fr is None: return None, info
    W.f = fr
    return W.to_word(), info
if __name__ == '__main__':
    w = fw.load('K8-FULL-WORD.json.gz')
    W = search.Word(w); O = reorder2.OW2(W)
    cand = []
    for p, g in enumerate(O.order):
        if g == O.CUT or O.sec0[g] == 'T': continue
        if all(W.end[r] == W.full for r in O.regs[g]) and all(O.commute(g, h) for h in O.order[p + 1:]):
            cand.append(W.seq[g][1])
    print('sinkable A indices', cand)
    w2 = sink(w, cand)
    w3, info = reframe(w2)
    print(info)
    if w3:
        print('profile', w3['report']['local_profile'])
        fw.save(w3, 'W-sink10.json.gz')
