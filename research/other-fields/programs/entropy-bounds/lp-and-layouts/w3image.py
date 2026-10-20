# w3image.py -- for three opposite pairs that are NOT compatible: cut value of (a) Gamma_m itself, (b) the image graph
# of psi(v) = (d(u1,v), d(u2,v), d(u3,v)) (u_i one end of pair i; vertices = images, edges = images of edges), (c) the
# full king graph {0..m}^3 (l-infinity adjacency). Pull-backs: (c) <= (b) <= (a).
import sys, json, random, itertools, time, numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
from scipy.sparse.csgraph import maximum_flow, breadth_first_order
import w3cut, w3compat
SCALE = 10 ** 6
def cut_value(n, edges, pairs, seeds=(), max_rounds=300):
    eu = np.array([u for u, v in edges]); ev = np.array([v for u, v in edges]); E = len(edges)
    X = sorted({v for p in pairs for v in p})
    pats = []
    for mask in range(0, 1 << (len(X) - 1)):
        A = {X[0]} | {X[i + 1] for i in range(len(X) - 1) if mask >> i & 1}
        if len(A) == len(X): continue
        s_ = sum(1 for (s, t) in pairs if (s in A) != (t in A))
        if s_ > 0: pats.append((sorted(A), [x for x in X if x not in A], s_))
    rows, rhs = [], []
    def sep(inQ): return sum(1 for (s, t) in pairs if inQ[s] != inQ[t])
    for q in seeds: rows.append((q[eu] != q[ev]).astype(float)); rhs.append(sep(q))
    for x in X:
        q = np.zeros(n, bool); q[x] = True; rows.append((q[eu] != q[ev]).astype(float)); rhs.append(sep(q))
    S_, T_ = n, n + 1
    for rnd in range(max_rounds):
        res = linprog(np.ones(E), A_ub=-sp.csr_matrix(np.array(rows)), b_ub=-np.array(rhs, float), bounds=(0, None), method='highs')
        c = res.x; cap = np.floor(c * SCALE).astype(np.int64); added = 0
        for (A, B, s_) in pats:
            big = SCALE * 1000
            r_ = np.concatenate([eu, ev, np.full(len(A), S_), np.array(B)]); c_ = np.concatenate([ev, eu, np.array(A), np.full(len(B), T_)])
            w_ = np.concatenate([cap, cap, np.full(len(A), big), np.full(len(B), big)])
            Cap = sp.csr_matrix((w_.astype(np.int32), (r_, c_)), shape=(n + 2, n + 2)); Cap.sum_duplicates()
            fl = maximum_flow(Cap, S_, T_)
            if fl.flow_value < s_ * SCALE - 1:
                R = (Cap - fl.flow.tocsr()).tocsr(); R.data = (R.data > 0).astype(np.int8); R.eliminate_zeros()
                reach = breadth_first_order(R, S_, directed=True, return_predecessors=False)
                q = np.zeros(n, bool); q[[v for v in reach.tolist() if v < n]] = True
                row = (q[eu] != q[ev]).astype(float)
                if row @ c < sep(q) - 1e-6: rows.append(row); rhs.append(sep(q)); added += 1
        if added == 0: return res.fun, -np.array(res.ineqlin.marginals), rows, rhs
    # not converged: the last dual is still feasible for the full dual (a packing of the cuts generated so far),
    # so it is still a valid lower bound once checked exactly
    return res.fun, -np.array(res.ineqlin.marginals), rows, rhs
def image_graph(g, us):
    m = g.m
    prof = np.stack([g.D[u] for u in us], axis=1)   # n x k
    keys = {}; lab = np.zeros(g.n, int)
    for v in range(g.n):
        k = tuple(prof[v].tolist()); lab[v] = keys.setdefault(k, len(keys))
    E = set()
    for (a, b) in g.edges:
        x, y = lab[a], lab[b]
        if x != y: E.add((min(x, y), max(x, y)))
    return len(keys), sorted(E), lab, keys
def king(m, k=3):
    pts = list(itertools.product(range(m + 1), repeat=k)); idx = {p: i for i, p in enumerate(pts)}
    E = []
    for p in pts:
        for d in itertools.product((-1, 0, 1), repeat=k):
            q = tuple(a + b for a, b in zip(p, d))
            if q in idx and idx[q] > idx[p]: E.append((idx[p], idx[q]))
    return len(pts), E, idx
if __name__ == '__main__':
    m = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3]); full = len(sys.argv) > 4 and sys.argv[4] == 'full'
    t0 = time.time(); g = w3cut.G(m); rng = random.Random(seed); D = g.D
    kn, kE, kidx = king(m)
    done = 0; out = []
    while done < trials:
        P = w3cut.rand_3k2(g, rng)
        if w3compat.any_compat(D, m, P): continue
        done += 1
        us3 = [p[0] for p in P]
        nH, EH, lab, keys = image_graph(g, us3)
        vH = cut_value(nH, EH, [(lab[s], lab[t]) for s, t in P])[0]
        us6 = [x for p in P for x in p]
        nH6, EH6, lab6, keys6 = image_graph(g, us6)
        vH6 = cut_value(nH6, EH6, [(lab6[s], lab6[t]) for s, t in P])[0]
        vK = cut_value(kn, kE, [(kidx[tuple(int(D[u, s]) for u in us3)], kidx[tuple(int(D[u, t]) for u in us3)]) for s, t in P])[0]
        rec = dict(P=P, image3=(nH, vH), image6=(nH6, vH6), king=vK)
        if full:
            rec['gamma'] = cut_value(g.n, g.edges, P)[0]
        out.append(rec)
        print(f"{P}: image3 {nH} pts cut {vH:.4f}; image6 {nH6} pts cut {vH6:.4f}; king cut {vK:.4f}" + (f"; Gamma cut {rec['gamma']:.4f}" if full else '') + f" ({time.time()-t0:.0f}s)", flush=True)
    json.dump(out, open(f"w3image_m{m}_s{seed}.json", "w"))
