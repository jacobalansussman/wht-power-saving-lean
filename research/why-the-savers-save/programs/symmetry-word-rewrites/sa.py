"""t-egraph SA over commuting reorders (plateau walk + rare uphill), frames re-optimised by local MILP per move.
Saves the best state (order of gids + frames) whenever the cost improves."""
import sys, time, random, json, math, pickle
sys.path.insert(0, '.')
import fw, search, reorder
seed = int(sys.argv[1]); minutes = float(sys.argv[2]); T0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
random.seed(seed)
W = search.Word(fw.load('K8-FULL-WORD.json.gz'))
O = reorder.OW(W)
cost = O.total_cost(); best = cost; c0 = cost
print('seed', seed, 'cost0', cost, flush=True)
t0 = time.time(); it = 0; acc = 0; up = 0
regs = [r for r in O.chain if len(O.chain[r]) > 1]
def sample_move():
    r = random.choice(regs); gs = O.chain[r]
    i = random.randrange(len(gs) - 1); g1, g2 = gs[i], gs[i + 1]
    if O.sec[g1] != O.sec[g2] or not O.commute(g1, g2): return None
    a, b = O.posof[g1], O.posof[g2]
    between = O.order[a + 1:b]
    opts = []
    if all(O.commute(g2, h) for h in between): opts.append((g2, a))
    if all(O.commute(g1, h) for h in between): opts.append((g1, b))
    return random.choice(opts) if opts else None
while time.time() - t0 < minutes * 60:
    it += 1
    mv = sample_move()
    if mv is None: continue
    T = T0 * max(0.0, 1 - (time.time() - t0) / (minutes * 60))
    d = O.try_move(mv[0], mv[1], accept=lambda d: d < 1e-9 or (T > 0 and random.random() < math.exp(-d / max(T, 1e-9))))
    if d is None: continue
    acc += 1; cost += d
    if d > 1e-9: up += 1
    if cost < best - 1e-9:
        best = cost
        print(f'{time.time()-t0:7.1f}s it {it} acc {acc} up {up} NEW BEST {best:.4f} (gain {c0-best:.4f}) profile {dict(sorted(O.profile().items()))}', flush=True)
        pickle.dump(dict(order=O.order, f=O.f, cost=best), open(f'sa_best_{seed}.pkl', 'wb'))
    if it % 2000 == 0:
        print(f'{time.time()-t0:7.1f}s it {it} acc {acc} up {up} cost {cost:.4f} best {best:.4f}', flush=True)
print('done', it, acc, 'best', best, 'gain', c0 - best, flush=True)
