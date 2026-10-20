"""x-pack: a pad-minimal explicit bank table for Ryan Shea's C588 certificate (Kothe's PR588 atm3 + 112 source
copies, #401/#377 dedicated banks).  Replaces Shea's cyclic 70-bin schedule (5,040 banks) by the 4,999-bank optimum.

Ceiling: the 72 rows carry 72 * 5068 = 364896 residual width, a bank holds at most 73, so at least
ceil(364896/73) = 4999 banks with total pad 4999*73 - 364896 = 31.  Pads enter the invoice as blocks p with price
p ((73/p)^a - 1) >= 0, concave in p, so ONE pad of rank 31 (4998 full banks) is the cheapest possible pad multiset,
and every extra bank adds 73 to the pad total.  Nothing below this exists for this residual list.

Construction:
 1. bank CONTENTS by an integer programme (HiGHS): x_p >= 0 banks of each width pattern p over {4,5,7,8,9} with
    load exactly 73, one bank of load 42, using exactly 72 * n_w slots of width w;
 2. width-w slots are dealt round-robin to the n_w helpers of width w (72 slots each; consecutive slots of one bank
    go to distinct helpers when possible);
 3. the call (row) of every slot = a proper 72-edge-colouring of the bipartite multigraph banks x helpers
    (helper degree 72, bank degree <= 18 <= 72; Koenig's theorem), found by alternating-path recolouring.
Every row is then injective (proper colouring at the banks), every helper appears in every row exactly once
(proper at the helpers with degree 72), and every bank's load is its pattern load.
Usage: pack_ilp.py <C588Sched.lean (Shea; for schedRes order)> <out.schedule.json.gz>
"""
import sys, os, re, json, gzip
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tight_pack import patterns
from pack_exact import edge_colour, validate, read_res

M, CALLS = 73, 72


def solve_contents(cnt, last):
    import highspy, numpy as np
    ws = sorted(cnt, reverse=True)
    need = [CALLS * cnt[w] for w in ws]
    pf = patterns(ws, M)
    pl = patterns(ws, last) if last else []
    cols = [(p, 'full') for p in pf] + [(p, 'last') for p in pl]
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    inf = highspy.kHighsInf
    for j, (p, kind) in enumerate(cols):
        # objective: prefer few items per bank (smaller bank degree) -- any feasible point is fine
        h.addVar(0, inf if kind == 'full' else 1)
        h.changeColCost(j, sum(p))
        h.changeColIntegrality(j, highspy.HighsVarType.kInteger)
    for i, w in enumerate(ws):
        idx = [j for j, (p, _) in enumerate(cols) if p[i]]
        h.addRow(need[i], need[i], len(idx), np.array(idx, dtype=np.int32), np.array([cols[j][0][i] for j in idx], dtype=np.double))
    if last:
        idx = [j for j, (_, k) in enumerate(cols) if k == 'last']
        h.addRow(1, 1, len(idx), np.array(idx, dtype=np.int32), np.ones(len(idx)))
    h.run()
    st = h.getModelStatus()
    if st != highspy.HighsModelStatus.kOptimal:
        return None, str(st)
    x = [int(round(v)) for v in h.getSolution().col_value]
    banks = []
    for j, (p, kind) in enumerate(cols):
        for _ in range(x[j]):
            banks.append((kind, Counter({ws[i]: p[i] for i in range(len(ws)) if p[i]})))
    banks.sort(key=lambda b: b[0] == 'last')        # the padded bank last
    return [b[1] for b in banks], 'optimal'


def build(res, nbanks=None):
    cnt = Counter(res)
    tot = CALLS * sum(res)
    if nbanks is None:
        nbanks = -(-tot // M)
    last = M - (nbanks * M - tot)          # load of the single padded bank (pad = nbanks*M - tot)
    assert 0 < last <= M
    contents, st = solve_contents(cnt, last if last < M else 0)
    if contents is None:
        raise SystemExit('integer programme: ' + st)
    assert len(contents) == nbanks
    # deal slots round-robin to helpers of each width
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
    return rows, pads, contents


if __name__ == '__main__':
    sched_lean, out = sys.argv[1], sys.argv[2]
    nb = int(sys.argv[3]) if len(sys.argv) > 3 else None
    res = read_res(sched_lean)
    print('widths', sorted(Counter(res).items()), 'sum', sum(res), 'ceil banks', -(-CALLS * sum(res) // M))
    rows, pads, contents = build(res, nb)
    validate(rows, pads, res)
    print('VALID (Shea checks: 72 rows, injective rows, load+pad=73, duplicate-owner and padding controls rejected)')
    print('banks', len(pads), 'pads', sorted(Counter(pads).items()), 'max items per bank', max(sum(c.values()) for c in contents))
    print('patterns', sorted(Counter(tuple(sorted(c.items())) for c in contents).items(), key=lambda t: -t[1]))
    open(out, 'wb').write(gzip.compress(json.dumps(dict(rows=rows, pads=pads, widths=res), separators=(',', ':')).encode(), mtime=0))
