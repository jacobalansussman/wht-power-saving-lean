# layout_cut2.py -- the directed cut LP for a FIXED layout D (same LP as layout_cut.py), solved by constraint
# generation: LP over the arc uses c only; separation by max-flow / min-cut for each pattern P of endpoints
# (integer max-flow on floor(1e6 c), so a reported "no violated cut" is rigorous up to the 1e-6 scaling:
# every cut has c(in Q) >= floor-cut >= need - tol).
import sys, json, random, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
from scipy.sparse.csgraph import maximum_flow, breadth_first_order
import dpg

SCALE = 10 ** 6

def cut_lp(n, arcs, D, max_rounds=400, verbose=False, dist=None, stop_above=None):
    """returns (value, patterns, rounds). With stop_above=T it may stop early once the partial LP (a valid LOWER
    bound on the cut LP) exceeds T; the value returned is then that certified lower bound."""
    X = sorted({v for p in D for v in p})
    pats = []
    for mask in range(1, (1 << len(X)) - 1):
        P = frozenset(X[i] for i in range(len(X)) if mask >> i & 1)
        need = sum(1 for (s, t) in D if t in P and s not in P)
        if need > 0: pats.append((P, need))
    na = len(arcs)
    tail = np.array([u for (u, v) in arcs]); head = np.array([v for (u, v) in arcs])
    rows = []; rhs = []
    # seed cuts: for each pair, the balls around its sink (each arc enters at most one ball of a family)
    if dist is not None:
        Dm = np.asarray(dist)
        for (s0, t0) in D:
            for k in range(1, int(Dm[s0, t0]) + 1):
                inQ = Dm[t0] <= k - 1
                row = np.zeros(na); row[(~inQ[tail]) & inQ[head]] = 1.0
                rows.append(row); rhs.append(sum(1 for (s, t) in D if inQ[t] and not inQ[s]))
    c = np.zeros(na)
    S_, T_ = n, n + 1   # super source, super sink
    for rnd in range(max_rounds):
        if rows:
            A = sp.csr_matrix(np.array(rows))
            res = linprog(np.ones(na), A_ub=-A, b_ub=-np.array(rhs, dtype=float), bounds=(0, None), method='highs')
            c = res.x; val = res.fun
        else:
            val = 0.0
        if stop_above is not None and val > stop_above + 1e-7:
            return val, len(pats), rnd
        added = 0
        cap = np.floor(c * SCALE).astype(np.int64)
        for (P, need) in pats:
            src = [v for v in X if v not in P]; snk = list(P)
            big = SCALE * 10 ** 3
            r_ = np.concatenate([tail, np.full(len(src), S_), np.array(snk)])
            c_ = np.concatenate([head, np.array(src), np.full(len(snk), T_)])
            w_ = np.concatenate([cap, np.full(len(src), big), np.full(len(snk), big)])
            G = sp.csr_matrix((w_.astype(np.int32) if w_.max() < 2 ** 31 else w_, (r_, c_)), shape=(n + 2, n + 2))
            G.sum_duplicates()
            fl = maximum_flow(G, S_, T_)
            if fl.flow_value < need * SCALE - 1:   # violated: extract the min cut (source side) from the residual
                F = fl.flow.tocsr() if hasattr(fl, 'flow') else fl.residual.tocsr()
                Res = (G - F).tocsr()
                Res.data[Res.data < 0] = 0; Res.eliminate_zeros()
                # residual also has reverse arcs with capacity = flow
                Rev = F.T.tocsr(); Rev.data[Rev.data < 0] = 0; Rev.eliminate_zeros()
                Rg = (Res + Rev).tocsr(); Rg.data = (Rg.data > 0).astype(np.int8); Rg.eliminate_zeros()
                reach = set(breadth_first_order(Rg, S_, directed=True, return_predecessors=False).tolist())
                inQ = np.array([v not in reach for v in range(n)])      # Q = sink side
                row = np.zeros(na); row[(~inQ[tail]) & inQ[head]] = 1.0
                need_Q = sum(1 for (s, t) in D if inQ[t] and not inQ[s])
                if row @ c < need_Q - 1e-7:
                    rows.append(row); rhs.append(need_Q); added += 1
        if verbose: print(f"   round {rnd}: LP {val:.6f}, cuts {len(rows)}, added {added}")
        if added == 0:
            return val, len(pats), rnd
    return None, len(pats), max_rounds

def main():
    m = int(sys.argv[1]); W = int(sys.argv[2]); trials = int(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    verts, idx, adj, dist = dpg.build(m)
    n = len(verts)
    arcs = [(u, v) for u in range(n) for v in adj[u]]
    opp = [[t for t in range(n) if dist[s][t] == m] for s in range(n)]
    rng = random.Random(seed)
    results = []; t0 = time.time()
    def layout(kind):
        D = [(0, opp[0][0])]
        guard = 0
        while len(D) < W and guard < 1000:
            guard += 1
            used = [x for p in D for x in p]
            if kind == 'random':
                s = rng.randrange(n); t = rng.choice(opp[s])
            elif kind == 'near':
                s = rng.choice(adj[rng.choice(used)])
                cands = [t for t in opp[s] if any(dist[t][x] <= 1 for x in used)] or opp[s]
                t = rng.choice(cands)
            else:  # 'disjoint-near': new endpoints, each adjacent to an old endpoint
                us = set(used)
                cs = [v for x in used for v in adj[x] if v not in us]
                s = rng.choice(cs)
                cands = [t for t in opp[s] if t not in us and any(dist[t][x] == 1 for x in used)]
                if not cands: continue
                t = rng.choice(cands)
            if (s, t) not in D: D.append((s, t))
        return D
    kinds = ['random', 'near', 'disjoint-near']
    worst = None
    for tr in range(trials):
        kind = kinds[tr % 3]
        D = layout(kind)
        val, npat, rounds = cut_lp(n, arcs, D)
        results.append(dict(kind=kind, D=D, cut=val, rounds=rounds))
        if val is not None and (worst is None or val < worst[0] - 1e-6):
            worst = (val, D, kind)
            print(f"  m={m} W={W} trial {tr} ({kind}): cut LP {val:.6f} (routing {W*m}), patterns {npat}, rounds {rounds} ({time.time()-t0:.0f}s) NEW WORST")
            sys.stdout.flush()
    vals = [r['cut'] for r in results if r['cut'] is not None]
    print(f"m={m} W={W}: {len(vals)} layouts solved; cut LP min {min(vals):.6f}, max {max(vals):.6f}, routing {W*m}; "
          f"below routing: {sum(1 for v in vals if v < W*m - 1e-5)} ({time.time()-t0:.0f}s)")
    json.dump(dict(m=m, W=W, results=results), open(f"layout_cut2_m{m}_W{W}_s{seed}.json", "w"))

if __name__ == '__main__':
    main()
