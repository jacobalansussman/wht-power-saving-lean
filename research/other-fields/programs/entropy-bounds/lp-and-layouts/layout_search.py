# layout_search.py -- look for W-pair layouts on Gamma_m whose directed cut LP is as SMALL as possible (where cuts fail
# to exclude a coherent saver). Random sampling, then local search: replace one pair by a random opposite pair whose
# endpoints are within distance 1 of the other endpoints (or a random pair), accept if the cut LP does not increase.
import sys, json, random, time
import dpg
from layout_cut2 import cut_lp

def main():
    m, W, samples, steps, seed = (int(x) for x in sys.argv[1:6])
    verts, idx, adj, dist = dpg.build(m)
    n = len(verts)
    arcs = [(u, v) for u in range(n) for v in adj[u]]
    opp = [[t for t in range(n) if dist[s][t] == m] for s in range(n)]
    rng = random.Random(seed); t0 = time.time()
    def rand_pair(D):
        used = [x for p in D for x in p]
        if used and rng.random() < 0.7:
            s = rng.choice(adj[rng.choice(used)]) if rng.random() < 0.7 else rng.choice(used)
            cands = [t for t in opp[s] if any(dist[t][x] <= 1 for x in used)] or opp[s]
            return (s, rng.choice(cands))
        s = rng.randrange(n); return (s, rng.choice(opp[s]))
    def evalD(D):
        v, npat, r = cut_lp(n, arcs, D, dist=dist, stop_above=W * m - 1)
        return v if v is not None else 1e9
    best = None; log = []
    for i in range(samples):
        D = [(0, opp[0][0])]
        while len(D) < W:
            p = rand_pair(D)
            if p not in D: D.append(p)
        v = evalD(D)
        log.append(('sample', v)); print(f"   sample {i}: {v:.6f} ({time.time()-t0:.0f}s)"); sys.stdout.flush()
        if best is None or v < best[0]: best = (v, D)
    print(f"m={m} W={W}: best of {samples} samples {best[0]:.6f} (routing {W*m}) ({time.time()-t0:.0f}s)"); sys.stdout.flush()
    cur = best
    for st in range(steps):
        D = list(cur[1]); j = rng.randrange(1, W) if W > 1 else 0
        p = rand_pair(D[:j] + D[j + 1:])
        if p in D: continue
        D[j] = p
        v = evalD(D)
        if v <= cur[0] + 1e-9:
            cur = (v, D)
            if v < best[0] - 1e-9:
                best = (v, D)
                print(f"   step {st}: cut LP {v:.6f} ({time.time()-t0:.0f}s)"); sys.stdout.flush()
    print(f"m={m} W={W} (early stop once > {W*m-1}, i.e. no saving certified): minimum found {best[0]:.6f} vs routing {W*m}; layout {best[1]}  ({time.time()-t0:.0f}s)")
    json.dump(dict(m=m, W=W, best=best[0], layout=best[1]), open(f"lsearch_m{m}_W{W}_s{seed}.json", "w"))

if __name__ == '__main__':
    main()
