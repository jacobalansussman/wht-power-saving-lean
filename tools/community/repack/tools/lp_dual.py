"""LP relaxation of the bank-contents programme (pack_opt.py) with its dual prices: the dual gives a lower bound
for EVERY bank table of Shea's C588 residual list: P_a >= sum_w 72 n_w pi_w."""
import sys, os, json
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pack_opt import all_patterns, g, M, CALLS
import highspy, numpy as np
cnt = {9: 388, 8: 140, 7: 32, 5: 8, 4: 48}
a = float(sys.argv[1]) if len(sys.argv) > 1 else 0.00124
ws = sorted(cnt, reverse=True)
need = [CALLS * cnt[w] for w in ws]
pats = all_patterns(ws, M)
h = highspy.Highs(); h.setOptionValue('output_flag', False)
for j, p in enumerate(pats):
    load = sum(p[i] * ws[i] for i in range(len(ws)))
    h.addVar(0, highspy.kHighsInf); h.changeColCost(j, g(M - load, a) / a)
for i in range(len(ws)):
    idx = [j for j, p in enumerate(pats) if p[i]]
    h.addRow(need[i], need[i], len(idx), np.array(idx, dtype=np.int32), np.array([pats[j][i] for j in idx], dtype=np.double))
h.run()
sol = h.getSolution()
print('LP objective (P_a/a)', h.getInfo().objective_function_value)
print('duals pi_w', dict(zip(ws, sol.row_dual)))
for j, v in enumerate(sol.col_value):
    if v > 1e-9:
        p = pats[j]; load = sum(p[i] * ws[i] for i in range(len(ws)))
        print('  %10.3f x %s load %d pad %d' % (v, dict((ws[i], p[i]) for i in range(len(ws)) if p[i]), load, M - load))
