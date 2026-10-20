# cells_lp.py -- the entropy LP on an actual contraction of Gamma_m (Lemma 2): cells as nodes, crossing arcs grouped
# per ordered pair of adjacent cells, sessions grouped per ordered pair of cells (rate = number of opposite pairs).
# Reports routing on the quotient, the quotient's cut LP, the entropy LP (closure-reduced full Shannon, plus
# Ingleton/DFZ cutting planes), and the per-session bound each gives after averaging.
import sys, json, collections, itertools
import numpy as np
from scipy.optimize import linprog
import dpg, nclp

def quotient(m, part):
    verts, idx, adj, dist = dpg.build(m)
    n = len(verts)
    cells = sorted(set(part)); ci = {c: i for i, c in enumerate(cells)}; k = len(cells)
    p = [ci[part[v]] for v in range(n)]
    E = collections.Counter(); 
    for u in range(n):
        for v in adj[u]:
            if p[u] != p[v]: E[(p[u], p[v])] += 1      # arcs u->v
    dem = collections.Counter((p[s], p[t]) for s in range(n) for t in range(n) if dist[s][t] == m and p[s] != p[t])
    nA = sum(len(a) for a in adj); nU = sum(1 for s in range(n) for t in range(n) if dist[s][t] == m)
    return k, E, dem, nA, nU

def run(m, part, label, linear_rank=True):
    k, E, dem, nA, nU = quotient(m, part)
    cross = sum(E.values())
    arcs = sorted(E)                         # one ground element per ordered adjacent cell pair
    sess = sorted(dem)
    # quotient distances
    import collections as C
    nb = {i: set() for i in range(k)}
    for (a, b) in arcs: nb[a].add(b)
    def qd(s, t):
        seen = {s}; fr = [s]; d = 0
        while fr:
            if t in fr: return d
            d += 1; nf = []
            for x in fr:
                for y in nb[x]:
                    if y not in seen: seen.add(y); nf.append(y)
            fr = nf
        return None
    R = sum(c * qd(a, b) for (a, b), c in dem.items())
    scale = nA / (cross * nU)
    # cut LP on the quotient (arc groups have no capacity limit; cost per use 1)
    rows = []; rhs = []
    for q in range(1, 1 << k):
        Q = {t for t in range(k) if q >> t & 1}
        need = sum(c for (a, b), c in dem.items() if a not in Q and b in Q)
        if need: rows.append([-1.0 if (a not in Q and b in Q) else 0.0 for (a, b) in arcs]); rhs.append(-need)
    cut = linprog(np.ones(len(arcs)), A_ub=np.array(rows), b_ub=np.array(rhs), bounds=(0, None), method='highs').fun
    print(f"{label}: {k} cells, {len(arcs)} arc groups, {len(sess)} session groups -> ground {len(arcs)+len(sess)}")
    print(f"   routing on quotient {R} -> {R*scale:.5f} per session;  cut LP {cut} -> {cut*scale:.5f}")
    sys.stdout.flush()
    out = dict(label=label, cells=k, ground=len(arcs) + len(sess), routing=R, cut=cut, scale=scale)
    if len(arcs) + len(sess) <= 16:
        net = nclp.Net(list(range(k)), arcs, [(a, [b], dem[(a, b)]) for (a, b) in sess])
        lp = nclp.build_closed_lp(net)
        Rtot = sum(dem.values())
        primal, LB, res = nclp.solve_pdlp(lp, Rtot)
        out['shannon'] = primal; out['shannon_rigorous_lb'] = LB
        print(f"   entropy LP (full Shannon) {primal:.6f} (rigorous >= {LB:.6f}) -> {LB*scale:.5f} per session"); sys.stdout.flush()
        if linear_rank:
            h = res.x; canon = lp['canon']
            letters = sorted({1 << i for i in range(net.n)} | {(1 << i) | (1 << j) for i in range(net.n) for j in range(i + 1, net.n)})
            ni, vi = nclp.find_cuts(h, canon, net.n, letters, 'ingleton', tries=200000)
            nd, vd = nclp.find_cuts(h, canon, net.n, letters, 'dfz', tries=200000)
            print(f"   Ingleton violated {len(ni)} (max {vi:.1e}), DFZ violated {len(nd)} (max {vd:.1e})")
            if ni or nd:
                primal, LB, res = nclp.solve_pdlp(lp, Rtot, extra_ub=ni + nd); out['linear_rank'] = LB
                print(f"   with them: {primal:.6f} (rigorous >= {LB:.6f})")
    return out

if __name__ == '__main__':
    results = []
    # block 2: best 3-cell contraction from quot.py's search (a line, and a split of the rest), rebuilt here by search
    for m, fn in ((2, None), (3, 'quot2_m3.json')):
        if fn:
            d = json.load(open(fn))
            for kk in ('5', '6'):
                part = d[kk]['part']
                results.append(run(m, part, f"m={m} annealed k={kk} (rho {d[kk]['rho']:.4f})"))
        else:
            verts, idx, adj, dist = dpg.build(2)
            v = 0
            line = [v] + [x for x in adj[v] if any(y in adj[x] for y in adj[v] if y != x)][:2]
            part = [0 if u in line else (1 if dist[v][u] == 2 else 2) for u in range(len(verts))]
            results.append(run(2, part, "m=2 three cells (a line / far / rest)"))
    json.dump(results, open('cells_lp.json', 'w'), indent=1)
