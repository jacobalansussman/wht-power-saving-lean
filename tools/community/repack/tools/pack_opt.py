"""x-pack: OPTIMAL bank packing for Ryan Shea's C588 certificate (Kothe's PR588 atm3 + 112 source copies,
#401/#377 dedicated banks, Lean checker GX.schedK).

Model.  A bank table = 72 rows (calls); each row maps the 616 residual helpers injectively to banks; bank b's load
(sum of the widths it receives over all rows) must be 73 - pad_b with 0 <= pad_b < 73 (GX.schedK).  The invoice
is 72 x local + 112 x bankPrice377 + one block of rank pad_b per bank with pad_b > 0; the deficit stays 2560 and
the moment equation  sum_r r n_r ((73/r)^a - 1) = 2560  does not involve W otherwise, so the best table is the
one minimising  P_a = sum_b g_a(pad_b),  g_a(p) = p ((73/p)^a - 1)  (at the final a; checked by a fixed point).

Reduction.  Bank CONTENTS (a multiset of widths per bank) determine everything: given contents that use exactly
72 n_w slots of width w, deal slots to the helpers of width w (72 each) and colour the bipartite multigraph
banks x helpers with 72 colours (Koenig: max degree 72), colour = row.  So the optimum is the integer programme
   min sum_p g_a(73 - load p) x_p   s.t.  sum_p p_w x_p = 72 n_w (w in {4,5,7,8,9}),  x_p in Z>=0,
over all width patterns p with load <= 73.  Solved exactly by HiGHS (MIP gap 0).
Usage: pack_opt.py <C588Sched.lean> <out.schedule.json.gz> [a]
"""
import sys, os, json, gzip, math
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pack_exact import edge_colour, validate, read_res

M, CALLS = 73, 72


def all_patterns(ws, cap):
    out = []
    def rec(i, rem, cur):
        if i == len(ws):
            if rem < cap:   # load >= 1
                out.append(tuple(cur))
            return
        for k in range(rem // ws[i] + 1):
            rec(i + 1, rem - k * ws[i], cur + [k])
    rec(0, cap, [])
    return out


def g(p, a):
    return 0.0 if p == 0 else p * ((M / p) ** a - 1)


def solve(cnt, a, time_limit=600):
    import highspy, numpy as np
    ws = sorted(cnt, reverse=True)
    need = [CALLS * cnt[w] for w in ws]
    pats = all_patterns(ws, M)
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    h.setOptionValue('mip_rel_gap', 0.0)
    h.setOptionValue('mip_abs_gap', 0.0)
    h.setOptionValue('time_limit', float(time_limit))
    inf = highspy.kHighsInf
    for j, p in enumerate(pats):
        load = sum(p[i] * ws[i] for i in range(len(ws)))
        h.addVar(0, inf)
        # cost scaled by 1/a (first-order units ~ p ln(73/p)); tiny tie-break towards fewer banks
        h.changeColCost(j, g(M - load, a) / a + 1e-7)
        h.changeColIntegrality(j, highspy.HighsVarType.kInteger)
    for i in range(len(ws)):
        idx = [j for j, p in enumerate(pats) if p[i]]
        h.addRow(need[i], need[i], len(idx), np.array(idx, dtype=np.int32), np.array([pats[j][i] for j in idx], dtype=np.double))
    h.run()
    st = h.getModelStatus()
    info = h.getInfo()
    x = [int(round(v)) for v in h.getSolution().col_value]
    contents = []
    for j, p in enumerate(pats):
        for _ in range(x[j]):
            contents.append(Counter({ws[i]: p[i] for i in range(len(ws)) if p[i]}))
    contents.sort(key=lambda c: -sum(w * n for w, n in c.items()))     # full banks first
    return contents, dict(status=str(st), objective=info.objective_function_value, bound=info.mip_dual_bound,
                          gap=info.mip_gap, npatterns=len(pats))


def table(res, contents):
    hw = defaultdict(list)
    for q, w in enumerate(res):
        hw[w].append(q)
    ptr = {w: 0 for w in hw}
    banks_items = []
    for c in contents:
        items = []
        for w, n in sorted(c.items()):
            for _ in range(n):
                items.append(hw[w][ptr[w] % len(hw[w])]); ptr[w] += 1
        banks_items.append(items)
    assert all(ptr[w] == CALLS * len(hw[w]) for w in hw)
    edges, col = edge_colour(banks_items, list(range(len(res))))
    rows = [[None] * len(res) for _ in range(CALLS)]
    for (b, q), c in zip(edges, col):
        assert rows[c][q] is None
        rows[c][q] = b
    assert all(x is not None for row in rows for x in row)
    pads = [M - sum(w * n for w, n in c.items()) for c in contents]
    return rows, pads


if __name__ == '__main__':
    sched_lean, out = sys.argv[1], sys.argv[2]
    a = float(sys.argv[3]) if len(sys.argv) > 3 else 0.00124
    res = read_res(sched_lean)
    cnt = Counter(res)
    contents, info = solve(cnt, a)
    print('MIP', info)
    pads = [M - sum(w * n for w, n in c.items()) for c in contents]
    print('banks', len(contents), 'pads', sorted(Counter(p for p in pads if p).items()), 'pad total', sum(pads),
          'P_a/a = %.4f' % (sum(g(p, a) for p in pads) / a))
    print('patterns', sorted(Counter(tuple(sorted(c.items())) for c in contents).items(), key=lambda t: -t[1]))
    rows, pads = table(res, contents)
    validate(rows, pads, res)
    print('VALID (Shea checks: 72 rows, injective rows, load+pad=73, duplicate-owner and padding controls rejected)')
    open(out, 'wb').write(gzip.compress(json.dumps(dict(rows=rows, pads=pads, widths=res), separators=(',', ':')).encode(), mtime=0))
    json.dump(dict(info, a=a, banks=len(pads), pads=dict(sorted(Counter(pads).items())),
                   patterns=[[dict(k), v] for k, v in Counter(tuple(sorted(c.items())) for c in contents).items()]),
              open(out.replace('.schedule.json.gz', '') + '.mip.json', 'w'), indent=1)
