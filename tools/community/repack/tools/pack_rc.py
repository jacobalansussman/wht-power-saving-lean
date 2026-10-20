"""x-pack: provably optimal bank contents for Shea's C588 residual list, by LP duality + reduced-cost pruning.

For any dual prices pi feasible for the LP relaxation (c_p - sum_w p_w pi_w = rc_p >= 0 for every pattern p),
EVERY bank table has cost  P = LP* + sum_p x_p rc_p.  So, given an incumbent of cost P_inc, an optimal table uses
only patterns with rc_p <= P_inc - LP*; the MIP restricted to those patterns (few) is solved to gap 0, and its
optimum is the optimum of the full problem.  Costs are g_a(pad)/a with g_a(p) = p((73/p)^a - 1).
Usage: pack_rc.py <C588Sched.lean> <out.schedule.json.gz> <a> [incumbent_cost]"""
import sys, os, json, gzip
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pack_opt import all_patterns, g, M, CALLS, table
from pack_exact import validate, read_res
import highspy, numpy as np


def lp(ws, need, pats, cost):
    h = highspy.Highs(); h.setOptionValue('output_flag', False)
    for j in range(len(pats)):
        h.addVar(0, highspy.kHighsInf); h.changeColCost(j, cost[j])
    for i in range(len(ws)):
        idx = [j for j, p in enumerate(pats) if p[i]]
        h.addRow(need[i], need[i], len(idx), np.array(idx, dtype=np.int32), np.array([pats[j][i] for j in idx], dtype=np.double))
    h.run()
    return h.getInfo().objective_function_value, list(h.getSolution().row_dual)


def mip(ws, need, pats, cost, tl=900):
    h = highspy.Highs(); h.setOptionValue('output_flag', False)
    h.setOptionValue('mip_rel_gap', 0.0); h.setOptionValue('mip_abs_gap', 0.0); h.setOptionValue('time_limit', float(tl))
    for j in range(len(pats)):
        h.addVar(0, highspy.kHighsInf); h.changeColCost(j, cost[j]); h.changeColIntegrality(j, highspy.HighsVarType.kInteger)
    for i in range(len(ws)):
        idx = [j for j, p in enumerate(pats) if p[i]]
        h.addRow(need[i], need[i], len(idx), np.array(idx, dtype=np.int32), np.array([pats[j][i] for j in idx], dtype=np.double))
    h.run()
    info = h.getInfo()
    return str(h.getModelStatus()), info.objective_function_value, info.mip_dual_bound, [int(round(v)) for v in h.getSolution().col_value]


if __name__ == '__main__':
    sched_lean, out, a = sys.argv[1], sys.argv[2], float(sys.argv[3])
    inc = float(sys.argv[4]) if len(sys.argv) > 4 else None
    res = read_res(sched_lean)
    cnt = Counter(res); ws = sorted(cnt, reverse=True); need = [CALLS * cnt[w] for w in ws]
    pats = all_patterns(ws, M)
    loads = [sum(p[i] * ws[i] for i in range(len(ws))) for p in pats]
    cost = [g(M - l, a) / a for l in loads]
    LP, pi = lp(ws, need, pats, cost)
    rc = [cost[j] - sum(pats[j][i] * pi[i] for i in range(len(ws))) for j in range(len(pats))]
    print('LP* =', LP, 'duals', dict(zip(ws, pi)), 'min rc', min(rc))
    if inc is None:
        inc = LP + 40
    gap = inc - LP + 1e-6
    keep = [j for j in range(len(pats)) if rc[j] <= gap]
    print('incumbent', inc, 'allowed rc <=', gap, 'patterns kept', len(keep), 'of', len(pats))
    st, obj, bnd, x = mip(ws, need, [pats[j] for j in keep], [cost[j] for j in keep])
    print('restricted MIP', st, 'objective', obj, 'bound', bnd)
    contents = []
    for k, j in enumerate(keep):
        for _ in range(x[k]):
            contents.append(Counter({ws[i]: pats[j][i] for i in range(len(ws)) if pats[j][i]}))
    contents.sort(key=lambda c: -sum(w * n for w, n in c.items()))
    pads = [M - sum(w * n for w, n in c.items()) for c in contents]
    print('banks', len(contents), 'pads', sorted(Counter(p for p in pads if p).items()), 'P/a', sum(g(p, a) for p in pads) / a)
    print('patterns', sorted(((dict(k), v) for k, v in Counter(tuple(sorted(c.items())) for c in contents).items()), key=lambda t: -t[1]))
    rows, pads = table(res, contents)
    validate(rows, pads, res)
    print('VALID (Shea checks: 72 rows, injective rows, load+pad=73, duplicate-owner and padding controls rejected)')
    open(out, 'wb').write(gzip.compress(json.dumps(dict(rows=rows, pads=pads, widths=res), separators=(',', ':')).encode(), mtime=0))
    json.dump(dict(a=a, LP=LP, duals=dict(zip(map(str, ws), pi)), incumbent=inc, kept=len(keep), status=st, objective=obj, bound=bnd,
                   banks=len(pads), pads=dict(sorted(Counter(pads).items()))), open(out.replace('.schedule.json.gz', '') + '.rc.json', 'w'), indent=1)
