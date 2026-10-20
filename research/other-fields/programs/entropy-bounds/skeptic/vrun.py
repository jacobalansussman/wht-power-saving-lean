# vrun.py -- v-nc-gain Colab runner: my cut LPs (vcut model, my bound) on my saved Gamma_m (dpg indexing).
# usage: vrun.py <task file .json>  ; task = list of dicts {m, P, mode: directed|undirected|image, orient (image),
#        solver: pdlp|ipm, claim, tl}
import sys, json, time, itertools
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
import vcut

def bound(c, Aub, bub, Aeq, beq, U, yub, yeq):
    yu = np.minimum(yub, 0.0)
    red = c - Aub.T @ yu - Aeq.T @ yeq
    return float(yu @ bub + yeq @ beq + np.minimum(red, 0).sum() * U)

def solve_ipm(model):
    c, Aub, bub, Aeq, beq, U, npat = model
    res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=(0, U), method='highs-ipm')
    LB = bound(c, Aub, bub, Aeq, beq, U, np.array(res.ineqlin.marginals), np.array(res.eqlin.marginals))
    return float(res.fun), LB, res.message

def solve_pdlp(model, tl):
    import highspy
    c, Aub, bub, Aeq, beq, U, npat = model
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
    sol = h.getSolution(); y = np.array(sol.row_dual); x = np.array(sol.col_value); nu = Aub.shape[0]
    LB = max(bound(c, Aub, bub, Aeq, beq, U, y[:nu], y[nu:]), bound(c, Aub, bub, Aeq, beq, U, -y[:nu], -y[nu:]))
    return float(c @ x), LB, h.modelStatusToString(h.getModelStatus())

def image(n, edges, D, us):
    prof = np.stack([D[u] for u in us], axis=1)
    keys = {}; lab = np.zeros(n, int)
    for v in range(n): lab[v] = keys.setdefault(tuple(prof[v].tolist()), len(keys))
    E = sorted({(min(lab[a], lab[b]), max(lab[a], lab[b])) for a, b in edges if lab[a] != lab[b]})
    return len(keys), E, lab

if __name__ == '__main__':
    tasks = json.load(open(sys.argv[1])); graphs = {}; t0 = time.time()
    for T in tasks:
        m = T['m']
        if m not in graphs:
            E = np.load(f'gamma{m}_edges_dpgidx.npy'); D = np.load(f'gamma{m}_D_dpgidx.npy').astype(np.int16)
            graphs[m] = ([tuple(map(int, e)) for e in E], D)
        edges, D = graphs[m]; n = D.shape[0]
        P = [tuple(p) for p in T['P']]
        assert all(D[s, t] == m for s, t in P), 'pair not opposite'
        mode = T['mode']
        if mode == 'image':
            vals = []
            for o in ([T['orient']] if T.get('orient') is not None else itertools.product((0, 1), repeat=len(P))):
                us = [P[i][o[i]] for i in range(len(P))]
                nH, EH, lab = image(n, edges, D, us)
                model = vcut.cut_flow(nH, EH, [(int(lab[s]), int(lab[t])) for s, t in P], model_only=True)
                pr, LB, st = solve_ipm(model) if T.get('solver', 'ipm') == 'ipm' else solve_pdlp(model, T.get('tl', 600))
                vals.append((list(o), nH, pr, LB))
            best = max(vals, key=lambda v: v[3])
            print(f"m={m} {P} image (orient {best[0]}, {best[1]} pts): primal {best[2]:.6f} LB {best[3]:.6f}; claim {T.get('claim')}  ({time.time()-t0:.0f}s)", flush=True)
        else:
            model = vcut.cut_flow(n, edges, P, directed=(mode == 'directed'), model_only=True)
            pr, LB, st = solve_ipm(model) if T.get('solver', 'ipm') == 'ipm' else solve_pdlp(model, T.get('tl', 5400))
            print(f"m={m} {P} {mode}: {st}; primal {pr:.6f} LB {LB:.6f}; vars {len(model[0])}; claim {T.get('claim')}  ({time.time()-t0:.0f}s)", flush=True)
