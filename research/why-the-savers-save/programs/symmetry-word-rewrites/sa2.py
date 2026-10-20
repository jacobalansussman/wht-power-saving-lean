"""t-egraph SA v2: long-range window moves (incl. across the cut), plateau walk + optional uphill; frames by local MILP."""
import sys, time, random, math, pickle
sys.path.insert(0, '.')
import fw, search, reorder2
seed = int(sys.argv[1]); minutes = float(sys.argv[2]); T0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
rng = random.Random(seed)
W = search.Word(fw.load(sys.argv[4] if len(sys.argv) > 4 else 'K8-FULL-WORD.json.gz'))
O = reorder2.OW2(W)
cost = O.total_cost(); best = cost; c0 = cost
print('seed', seed, 'cost0', cost, 'cross-capable gates', sum(O.cross_ok), flush=True)
t0 = time.time(); it = 0; acc = 0; up = 0; crossed = 0
while time.time() - t0 < minutes * 60:
    it += 1
    mv = O.random_window_move(rng)
    if mv is None: continue
    T = T0 * max(0.0, 1 - (time.time() - t0) / (minutes * 60))
    secb = O.sec[mv[0]]
    d = O.try_move(mv[0], mv[1], accept=lambda d: d < 1e-9 or (T > 0 and rng.random() < math.exp(-d / max(T, 1e-9))))
    if d is None: continue
    acc += 1; cost += d
    if O.sec[mv[0]] != secb: crossed += 1
    if d > 1e-9: up += 1
    if cost < best - 1e-9:
        best = cost
        print(f'{time.time()-t0:7.1f}s it {it} acc {acc} NEW BEST {best:.4f} (gain {c0-best:.4f}) profile {dict(sorted(O.profile().items()))}', flush=True)
        pickle.dump(dict(order=O.order, f=O.f, cost=best), open(f'sa2_best_{seed}.pkl', 'wb'))
    if it % 500 == 0:
        print(f'{time.time()-t0:7.1f}s it {it} acc {acc} up {up} crossed {crossed} cost {cost:.4f} best {best:.4f}', flush=True)
print('done', it, acc, 'crossed', crossed, 'best', best, 'gain', c0 - best, flush=True)
pickle.dump(dict(order=O.order, f=O.f, cost=cost), open(f'sa2_final_{seed}.pkl', 'wb'))
