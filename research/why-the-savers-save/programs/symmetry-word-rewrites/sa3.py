"""t-egraph SA v3: window moves (incl. across the cut) with a secondary objective that steers the plateau walk:
contamination S = sum_g (dim P(g) - dim Pd(g)) + (dim Qd(g) - dim Q(g)), i.e. how far each gate's legal interval is
pinched by frame coupling through COMMUTING gate pairs (Pd/Qd from the dependency DAG, lbdep.json).  Accept: primary
cost down; or primary equal and S down (or up with prob exp(-dS/T))."""
import sys, time, random, math, pickle, json
sys.path.insert(0, '.')
import fw, search, reorder2
seed = int(sys.argv[1]); minutes = float(sys.argv[2]); TS = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
rng = random.Random(seed)
W = search.Word(fw.load('K8-FULL-WORD.json.gz'))
O = reorder2.OW2(W)
d = json.load(open('lbdep.json')); Pd = {int(k): len(v) for k, v in d['Pd'].items()}; Qd = {int(k): len(v) for k, v in d['Qd'].items()}
def contam():
    return sum((len(O.P[g]) - Pd[g]) + (Qd[g] - len(O.Q[g])) for g in range(O.N))
cost = O.total_cost(); best = cost; c0 = cost; S = contam(); S0 = S
print('seed', seed, 'cost0', cost, 'contamination0', S, flush=True)
t0 = time.time(); it = 0; acc = 0
while time.time() - t0 < minutes * 60:
    it += 1
    mv = O.random_window_move(rng)
    if mv is None: continue
    snap = (list(O.order), dict(O.f))
    dd = O.try_move(mv[0], mv[1], accept=lambda x: x < 1e-9)
    if dd is None: continue
    S2 = contam()
    if dd < -1e-9 or S2 <= S or (TS > 0 and rng.random() < math.exp(-(S2 - S) / TS)):
        acc += 1; cost += dd; S = S2
    else:
        O.order, O.f = snap; O.refresh(); continue
    if cost < best - 1e-9:
        best = cost
        print(f'{time.time()-t0:7.1f}s it {it} acc {acc} NEW BEST {best:.4f} (gain {c0-best:.4f}) S {S} profile {dict(sorted(O.profile().items()))}', flush=True)
        pickle.dump(dict(order=O.order, f=O.f, cost=best), open(f'sa3_best_{seed}.pkl', 'wb'))
    if it % 500 == 0:
        print(f'{time.time()-t0:7.1f}s it {it} acc {acc} cost {cost:.4f} best {best:.4f} S {S} (S0 {S0})', flush=True)
print('done', it, acc, 'best', best, 'gain', c0 - best, 'S', S, flush=True)
pickle.dump(dict(order=O.order, f=O.f, cost=cost), open(f'sa3_final_{seed}.pkl', 'wb'))
