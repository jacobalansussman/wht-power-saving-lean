# layout_cut_pdlp.py -- the directed cut LP of a fixed layout D in its compact flow form (one flow per endpoint
# pattern, all under the same arc uses c), solved by HiGHS PDLP, with a RIGOROUS lower bound from the approximate
# dual: every variable is boxed in [0, W*m] (an optimal c never needs more, since routing costs W*m), so
# LB = b.y + sum_j min(0, reduced cost_j) * W*m is valid for any multipliers y (<=-rows clipped to y <= 0).
import sys, json, random, time
import numpy as np
import scipy.sparse as sp
import highspy
import dpg

def build(n, arcs, D, Wm):
    X = sorted({v for p in D for v in p})
    pats = []
    for mask in range(1, (1 << len(X)) - 1):
        P = {X[i] for i in range(len(X)) if mask >> i & 1}
        need = sum(1 for (s, t) in D if t in P and s not in P)
        if need > 0: pats.append((P, need))
    na = len(arcs); tail = np.array([u for (u, v) in arcs]); head = np.array([v for (u, v) in arcs])
    nvar = na; offs = []
    for (P, need) in pats:
        srcs = [v for v in X if v not in P]; snks = [v for v in X if v in P]
        offs.append((nvar, srcs, snks, need)); nvar += na + len(srcs) + len(snks)
    I = []; J = []; V = []; lo = []; hi = []; r = 0
    for (off, srcs, snks, need) in offs:
        idx = np.arange(na)
        I.append(r + idx); J.append(off + idx); V.append(np.ones(na))
        I.append(r + idx); J.append(idx); V.append(-np.ones(na))
        lo.append(np.full(na, -np.inf)); hi.append(np.zeros(na)); r += na
        so = off + na; sk = so + len(srcs)
        I.append(r + head); J.append(off + idx); V.append(np.ones(na))
        I.append(r + tail); J.append(off + idx); V.append(-np.ones(na))
        I.append(r + np.array(srcs)); J.append(so + np.arange(len(srcs))); V.append(np.ones(len(srcs)))
        I.append(r + np.array(snks)); J.append(sk + np.arange(len(snks))); V.append(-np.ones(len(snks)))
        lo.append(np.zeros(n)); hi.append(np.zeros(n)); r += n
        I.append(np.full(len(srcs), r)); J.append(so + np.arange(len(srcs))); V.append(np.ones(len(srcs)))
        lo.append(np.array([float(need)])); hi.append(np.array([np.inf])); r += 1
    A = sp.csc_matrix((np.concatenate(V), (np.concatenate(I), np.concatenate(J))), shape=(r, nvar))
    c = np.zeros(nvar); c[:na] = 1.0
    return A, c, np.concatenate(lo), np.concatenate(hi), nvar, len(pats)

def solve(A, c, lo, hi, nvar, U, time_limit=1800):
    h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('solver', 'pdlp')
    h.setOptionValue('time_limit', float(time_limit)); h.setOptionValue('threads', 1)
    model = highspy.HighsLp(); model.num_col_ = nvar; model.num_row_ = A.shape[0]
    model.col_cost_ = c; model.col_lower_ = np.zeros(nvar); model.col_upper_ = np.full(nvar, float(U))
    INF = highspy.kHighsInf
    model.row_lower_ = np.where(np.isinf(lo), -INF, lo); model.row_upper_ = np.where(np.isinf(hi), INF, hi)
    model.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    model.a_matrix_.start_ = A.indptr; model.a_matrix_.index_ = A.indices; model.a_matrix_.value_ = A.data
    h.passModel(model); h.run()
    sol = h.getSolution(); y = np.array(sol.row_dual); x = np.array(sol.col_value)
    # rigorous bound: rows lo <= Ax <= hi. Use y_i >= 0 only on rows with finite lo, y_i <= 0 only with finite hi.
    yl = np.where(np.isfinite(lo), np.maximum(y, 0), 0.0); yu = np.where(np.isfinite(hi), np.minimum(y, 0), 0.0)
    yy = yl + yu
    red = c - A.T @ yy
    LB = float(np.where(np.isfinite(lo), lo, 0) @ yl + np.where(np.isfinite(hi), hi, 0) @ yu + np.minimum(red, 0).sum() * U)
    return float(c @ x), LB, h.modelStatusToString(h.getModelStatus())

if __name__ == '__main__':
    m = int(sys.argv[1]); W = int(sys.argv[2]); trials = int(sys.argv[3]); seed = int(sys.argv[4])
    verts, idx, adj, dist = dpg.build(m); n = len(verts)
    arcs = [(u, v) for u in range(n) for v in adj[u]]
    opp = [[t for t in range(n) if dist[s][t] == m] for s in range(n)]
    rng = random.Random(seed); t0 = time.time(); out = []
    layouts = json.loads(sys.argv[5]) if len(sys.argv) > 5 else []
    for tr in range(trials):
        if tr < len(layouts): D = [tuple(p) for p in layouts[tr]]
        else:
            D = [(0, opp[0][0])]
            while len(D) < W:
                s = rng.randrange(n); p = (s, rng.choice(opp[s]))
                if p not in D: D.append(p)
        A, c, lo, hi, nvar, npat = build(n, arcs, D, W * m)
        primal, LB, st = solve(A, c, lo, hi, nvar, W * m)
        out.append(dict(D=D, primal=primal, lb=LB))
        print(f"m={m} W={W} layout {D}: PDLP {st}, primal {primal:.5f}, rigorous LB {LB:.5f} (routing {W*m}); vars {nvar}, patterns {npat}  ({time.time()-t0:.0f}s)")
        sys.stdout.flush()
    json.dump(out, open(f"lpd_m{m}_W{W}_s{seed}.json", "w"))
