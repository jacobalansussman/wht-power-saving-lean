# os_pdlp.py -- the OS calibration LP (same as os_lp.py), solved with highspy directly: PDLP (first-order) for a fast
# approximate value, then dual simplex warm from nothing with presolve, both with time limits.
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # modules of ../
import sys, time, numpy as np, scipy.sparse as sp
import highspy, nclp
L, R = ['a', 'b'], ['x', 'y', 'z']
nodes = L + R
arcs = []
for u in L:
    for v in R: arcs += [(u, v), (v, u)]
sessions = [('x', ['y'], 1), ('y', ['z'], 1), ('z', ['x'], 1), ('a', ['b'], 1)]
net = nclp.Net(nodes, arcs, sessions)
rot = {'a': 'a', 'b': 'b', 'x': 'y', 'y': 'z', 'z': 'x'}
perm = [None] * net.n
for i, (s, t, r) in enumerate(sessions):
    perm[i] = next(j for j, (s2, t2, r2) in enumerate(sessions) if s2 == rot[s] and t2 == [rot[u] for u in t])
for a, (u, v) in enumerate(arcs): perm[net.ns + a] = net.ns + arcs.index((rot[u], rot[v]))
G = nclp.close_group([tuple(perm)], net.n)
lp = nclp.build_lp(net, perms=G)
A = sp.vstack([lp['A_ub'], lp['A_eq']]).tocsc()
nu = lp['A_ub'].shape[0]; ne = lp['A_eq'].shape[0]
lo = np.concatenate([np.full(nu, -highspy.kHighsInf), lp['b_eq']]); hi = np.concatenate([lp['b_ub'], lp['b_eq']])
def run(solver, tl):
    h = highspy.Highs()
    h.setOptionValue('output_flag', True); h.setOptionValue('solver', solver); h.setOptionValue('time_limit', float(tl))
    h.setOptionValue('threads', 1)
    model = highspy.HighsLp()
    model.num_col_ = lp['nv']; model.num_row_ = A.shape[0]
    model.col_cost_ = lp['c']; model.col_lower_ = np.zeros(lp['nv']); model.col_upper_ = np.full(lp['nv'], highspy.kHighsInf)
    model.row_lower_ = lo; model.row_upper_ = hi
    model.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    model.a_matrix_.start_ = A.indptr; model.a_matrix_.index_ = A.indices; model.a_matrix_.value_ = A.data
    h.passModel(model)
    t0 = time.time(); h.run()
    info = h.getInfo()
    print(f"== {solver}: status {h.modelStatusToString(h.getModelStatus())}, objective {info.objective_function_value:.9f}, {time.time()-t0:.0f}s"); sys.stdout.flush()
    return h
run('pdlp', 1800)
run('simplex', 5400)
