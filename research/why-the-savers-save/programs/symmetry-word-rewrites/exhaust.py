"""Exhaustive small sub-words: for a window of consecutive positions, enumerate every order of its gates that keeps
non-commuting pairs in order, up to swaps of register-disjoint gates (one representative per trace class, the
lexicographic normal form), and re-frame each with the local MILP.  Reports the best delta per window."""
import sys, time, json
sys.path.insert(0, '.')
import fw, search, reorder2
import os; W = search.Word(fw.load(os.environ.get('WORD', 'K8-FULL-WORD.json.gz'))); O = reorder2.OW2(W)
L = int(sys.argv[1]); maxtr = int(sys.argv[2]); starts = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else None
def traces(gs):
    n = len(gs); rs = [set(O.regs[g]) for g in gs]
    dep = [[j for j in range(i) if not O.commute(gs[i], gs[j])] for i in range(n)]
    out = []
    def rec(placed, seq, last):
        if len(out) > maxtr: return
        if len(seq) == n: out.append([gs[i] for i in seq]); return
        for i in range(n):
            if i in placed or any(j not in placed for j in dep[i]): continue
            if last is not None and not (rs[i] & rs[last]) and i < last: continue
            placed.add(i); seq.append(i); rec(placed, seq, i); seq.pop(); placed.discard(i)
    rec(set(), [], None)
    return out
res = []
t0 = time.time()
cands = starts if starts else list(range(0, len(O.order) - L, L))
for a in cands:
    win = O.order[a:a + L]
    if O.CUT in win: continue
    secs = {O.sec[g] for g in win}
    if len(secs) > 1: continue
    trs = traces(win)
    if len(trs) > maxtr or len(trs) <= 1:
        res.append((a, len(trs), None)); print(a, 'traces', len(trs), '(skipped)', flush=True); continue
    best = 0.0; nev = 0
    base_order = list(O.order); base_f = dict(O.f); oldP, oldQ = O.P, O.Q
    for tr in trs:
        if tr == win: continue
        O.order = base_order[:a] + tr + base_order[a + L:]
        O.refresh()
        if not O.feasible: continue
        changed = {h for h in O.order if O.P[h] != oldP.get(h) or O.Q[h] != oldQ.get(h)} | set(win)
        r = O.local_milp(changed)
        if r is None: continue
        newf, obj, regsF = r
        newfr = dict(base_f); newfr.update(newf)
        newc = sum(O.reg_cost(q, newfr) for q in regsF)
        O.order = base_order; O.refresh_chains_only()
        oldc = sum(O.reg_cost(q, base_f) for q in regsF)
        nev += 1; best = min(best, newc - oldc)
    O.order = base_order; O.f = base_f; O.refresh()
    res.append((a, len(trs), round(best, 4), nev))
    print(a, O.sec[win[0]], 'traces', len(trs), 'evaluated', nev, 'best delta', round(best, 4), f'{time.time()-t0:.0f}s', flush=True)
json.dump(res, open(f'exhaust_L{L}.json', 'w'))
print('windows', len(res), 'exhaustively searched', sum(1 for x in res if x[2] is not None), 'improving', sum(1 for x in res if x[2] is not None and x[2] < -1e-9))
