# vfull4.py -- v-nc-gain: cut LP of one three-pair layout on ALL of Gamma_4 (my graph, saved by vprep4.py in dpg's
# indexing), my own flow-form model (vcut.cut_flow model_only), HiGHS PDLP, and my own rigorous Lagrangian bound:
#   x in [0, U] (U = W), any y_ub <= 0, any y_eq:  c.x >= y_ub.b_ub + y_eq.b_eq + sum_j min(0, red_j) U.
# usage: vfull4.py '<json layout>' directed|undirected [time_limit]
import sys, json, time
import numpy as np, scipy.sparse as sp
import highspy
import vcut

def solve(c, Aub, bub, Aeq, beq, U, tl):
    A = sp.vstack([Aub, Aeq]).tocsc()
    lo = np.concatenate([np.full(Aub.shape[0], -np.inf), beq]); hi = np.concatenate([bub, beq])
    h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('solver', 'pdlp')
    h.setOptionValue('time_limit', float(tl)); h.setOptionValue('threads', 1)
    lp = highspy.HighsLp(); nv = len(c); lp.num_col_ = nv; lp.num_row_ = A.shape[0]
    lp.col_cost_ = c; lp.col_lower_ = np.zeros(nv); lp.col_upper_ = np.full(nv, float(U))
    INF = highspy.kHighsInf
    lp.row_lower_ = np.where(np.isinf(lo), -INF, lo); lp.row_upper_ = np.where(np.isinf(hi), INF, hi)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
    h.passModel(lp); h.run()
    sol = h.getSolution(); y = np.array(sol.row_dual); x = np.array(sol.col_value)
    nu = Aub.shape[0]
    yu = np.minimum(y[:nu], 0.0); ye = y[nu:]
    red = c - Aub.T @ yu - Aeq.T @ ye
    LB = float(yu @ bub + ye @ beq + np.minimum(red, 0).sum() * U)
    # same bound with the opposite sign convention for the duals (in case HiGHS reports -y); keep the larger, both valid
    yu2 = np.minimum(-y[:nu], 0.0); ye2 = -y[nu:]
    red2 = c - Aub.T @ yu2 - Aeq.T @ ye2
    LB2 = float(yu2 @ bub + ye2 @ beq + np.minimum(red2, 0).sum() * U)
    return float(c @ x), max(LB, LB2), h.modelStatusToString(h.getModelStatus())

if __name__ == '__main__':
    P = [tuple(p) for p in json.loads(sys.argv[1])]; directed = sys.argv[2] == 'directed'
    tl = float(sys.argv[3]) if len(sys.argv) > 3 else 5400
    E = np.load('gamma4_edges_dpgidx.npy'); D = np.load('gamma4_D_dpgidx.npy')
    edges = [tuple(map(int, e)) for e in E]; n = D.shape[0]
    assert all(D[s, t] == 4 for s, t in P)
    t0 = time.time()
    c, Aub, bub, Aeq, beq, U, npat = vcut.cut_flow(n, edges, P, directed=directed, model_only=True)
    print(f"{P} {'directed' if directed else 'undirected'}: {len(c)} vars, {Aub.shape[0] + Aeq.shape[0]} rows, {npat} patterns", flush=True)
    pr, LB, st = solve(c, Aub, bub, Aeq, beq, U, tl)
    print(f"{P} {'directed' if directed else 'undirected'}: PDLP {st}, primal {pr:.6f}, rigorous LB {LB:.6f}  ({time.time()-t0:.0f}s)", flush=True)
