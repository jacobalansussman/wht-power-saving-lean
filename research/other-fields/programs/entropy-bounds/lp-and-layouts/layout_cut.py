# layout_cut.py -- for a FIXED layout D of W opposite pairs on Gamma_m: the directed cut LP
#   min sum_arcs c   s.t. for every node set Q: c(arcs into Q) >= #{a : T_a in Q, S_a not in Q}
# (the cut-set rows of the entropy LP; a code with arc uses c satisfies them). Solved exactly as an LP by writing,
# for each pattern P of endpoints inside Q, a flow from (endpoints outside) to P of the required value under c
# (max-flow = min-cut). If the value is W*m, no code for D beats routing, so no coherent unit with that layout saves.
import sys, json, random, time, itertools, collections
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
import dpg

def cut_lp(n, arcs, D, verbose=False):
    X = sorted({v for p in D for v in p})
    xi = {v: i for i, v in enumerate(X)}
    na = len(arcs)
    pats = []
    for mask in range(1, (1 << len(X)) - 1):
        P = {X[i] for i in range(len(X)) if mask >> i & 1}
        need = sum(1 for (s, t) in D if t in P and s not in P)
        if need > 0: pats.append((P, need))
    # variables: c (na), then per pattern: f (na) + source arcs (|X\P|) + sink arcs (|P|)
    nvar = na; blocks = []
    for (P, need) in pats:
        srcs = [v for v in X if v not in P]; snks = [v for v in X if v in P]
        blocks.append((nvar, srcs, snks, need)); nvar += na + len(srcs) + len(snks)
    rows_I = []; rows_J = []; rows_V = []; b_ub = []; r = 0
    eq_I = []; eq_J = []; eq_V = []; b_eq = []; re = 0
    tail = np.array([u for (u, v) in arcs]); head = np.array([v for (u, v) in arcs])
    for (off, srcs, snks, need) in blocks:
        # capacity f_a - c_a <= 0
        rows_I += list(range(r, r + na)) * 2; rows_J += list(range(off, off + na)) + list(range(na))
        rows_V += [1.0] * na + [-1.0] * na; b_ub += [0.0] * na; r += na
        # conservation at each node: inflow - outflow + (source injection) - (sink extraction) = 0
        so = off + na; sk = so + len(srcs)
        for a in range(na):
            eq_I += [re + head[a], re + tail[a]]; eq_J += [off + a, off + a]; eq_V += [1.0, -1.0]
        for j, v in enumerate(srcs): eq_I.append(re + v); eq_J.append(so + j); eq_V.append(1.0)
        for j, v in enumerate(snks): eq_I.append(re + v); eq_J.append(sk + j); eq_V.append(-1.0)
        b_eq += [0.0] * n; re += n
        # value: sum of source injections >= need  ->  -sum <= -need
        for j in range(len(srcs)): rows_I.append(r); rows_J.append(so + j); rows_V.append(-1.0)
        b_ub.append(-float(need)); r += 1
    A_ub = sp.csr_matrix((rows_V, (rows_I, rows_J)), shape=(r, nvar))
    A_eq = sp.csr_matrix((eq_V, (eq_I, eq_J)), shape=(re, nvar))
    cost = np.zeros(nvar); cost[:na] = 1.0
    res = linprog(cost, A_ub=A_ub, b_ub=np.array(b_ub), A_eq=A_eq, b_eq=np.array(b_eq), bounds=(0, None), method='highs')
    return res.fun if res.status == 0 else None, len(pats)

def main():
    m = int(sys.argv[1]); W = int(sys.argv[2]); trials = int(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    verts, idx, adj, dist = dpg.build(m)
    n = len(verts)
    arcs = [(u, v) for u in range(n) for v in adj[u]]
    opp = [[t for t in range(n) if dist[s][t] == m] for s in range(n)]
    nbr = adj
    rng = random.Random(seed)
    results = []
    t0 = time.time()
    def layout(kind):
        D = [(0, opp[0][0])]
        while len(D) < W:
            if kind == 'random':
                s = rng.randrange(n); t = rng.choice(opp[s])
            elif kind == 'near':      # new pair with endpoints near existing endpoints (dense interleaving)
                a = rng.choice(D); s = rng.choice(nbr[rng.choice(a)])
                cands = [t for t in opp[s] if any(dist[t][x] <= 1 for p in D for x in p)] or opp[s]
                t = rng.choice(cands)
            elif kind == 'disjoint-near':  # endpoints distinct from all previous, at distance 1 from some
                used = {x for p in D for x in p}
                s = rng.choice([v for p in D for x in p for v in nbr[x] if v not in used] or list(range(n)))
                cands = [t for t in opp[s] if t not in used and any(dist[t][x] <= 1 for x in used)] or [t for t in opp[s] if t not in used]
                t = rng.choice(cands)
            if (s, t) not in D: D.append((s, t))
        return D
    kinds = ['random', 'near', 'disjoint-near']
    worst = None
    for tr in range(trials):
        kind = kinds[tr % 3]
        D = layout(kind)
        val, npat = cut_lp(n, arcs, D)
        results.append(dict(kind=kind, D=D, cut=val))
        if worst is None or val < worst[0] - 1e-9:
            worst = (val, D, kind)
            print(f"  m={m} W={W} trial {tr} ({kind}): cut LP {val:.6f} (routing {W*m}), patterns {npat}  ({time.time()-t0:.0f}s) NEW WORST")
            sys.stdout.flush()
    vals = [r['cut'] for r in results]
    print(f"m={m} W={W}: {trials} layouts; cut LP min {min(vals):.6f}, max {max(vals):.6f}, routing {W*m}; "
          f"layouts below routing: {sum(1 for v in vals if v < W*m - 1e-7)}")
    json.dump(dict(m=m, W=W, results=results), open(f"layout_cut_m{m}_W{W}_s{seed}.json", "w"))

if __name__ == '__main__':
    main()
