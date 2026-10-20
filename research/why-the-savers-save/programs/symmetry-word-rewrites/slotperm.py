"""Exhaustive accumulation orders: for register r, take a maximal run of its consecutive chain gates that pairwise
commute (e.g. all adds into a target between two butterflies), and try every permutation of those gates over their
slots; a permutation is legal if every gate commutes with every non-run gate between its old and new slot.  Each legal
permutation is re-framed with the local MILP.  Exhaustive per run (runs of length <= KMAX)."""
import sys, time, itertools, json
sys.path.insert(0, '.')
import fw, search, reorder2
import os; W = search.Word(fw.load(os.environ.get('WORD', 'K8-FULL-WORD.json.gz'))); O = reorder2.OW2(W)
regs = [int(x) for x in sys.argv[1].split(',')]; KMAX = int(sys.argv[2])
out = []
for r in regs:
    ch = O.chain[r]
    runs = []; cur = [ch[0]]
    for g in ch[1:]:
        if O.sec[g] == O.sec[cur[0]] and all(O.commute(g, h) for h in cur): cur.append(g)
        else:
            runs.append(cur); cur = [g]
    runs.append(cur)
    for run in runs:
        k = len(run)
        if k < 2 or k > KMAX: continue
        slots = sorted(O.posof[g] for g in run)
        base_order = list(O.order); base_f = dict(O.f); oldP, oldQ = O.P, O.Q
        runset = set(run)
        best = 0.0; legal = 0; t0 = time.time()
        for perm in itertools.permutations(run):
            if list(perm) == [base_order[s] for s in slots]: continue
            ok = True
            for g, s in zip(perm, slots):
                a, b = sorted((O.posof[g], s))
                if any(not O.commute(g, h) for h in base_order[a:b + 1] if h not in runset): ok = False; break
            if not ok: continue
            neworder = list(base_order)
            for g, s in zip(perm, slots): neworder[s] = g
            O.order = neworder; O.refresh()
            if not O.feasible: continue
            changed = {h for h in O.order if O.P[h] != oldP.get(h) or O.Q[h] != oldQ.get(h)} | runset
            res = O.local_milp(changed)
            if res is None: continue
            newf, obj, regsF = res
            newfr = dict(base_f); newfr.update(newf)
            newc = sum(O.reg_cost(q, newfr) for q in regsF)
            O.order = base_order; O.refresh_chains_only()
            oldc = sum(O.reg_cost(q, base_f) for q in regsF)
            legal += 1; best = min(best, newc - oldc)
        O.order = base_order; O.f = base_f; O.refresh()
        out.append(dict(reg=r, run_len=k, legal_perms=legal, best_delta=round(best, 4)))
        print(r, 'run', k, 'legal perms', legal, 'best delta', round(best, 4), f'{time.time()-t0:.0f}s', flush=True)
json.dump(out, open('slotperm_%s.json' % sys.argv[1].replace(',', '_')[:40], 'w'))
