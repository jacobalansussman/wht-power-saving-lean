# bulk.py -- refusal pipeline on every exported a2-nag relaxation: exact NulLA certificates of low degree, compared with
# a2-nag's Singular verdicts.  Usage: bulk.py out.jsonl file1.json file2.json ...
import sys, json, time, os
from certnum import *
out = open(sys.argv[1], 'a')
for fn in sys.argv[2:]:
    t0 = time.time(); F = PolySys.from_a2nag_json(fn); info = F.info or {}; na = info.get('nalpha_free', 0)
    sing = None; so = fn.replace('.json', '.sout')
    if os.path.exists(so):
        txt = open(so).read().split()
        try:
            i = txt.index('size:'); gbsize = int(txt[i + 1]); j = txt.index('dim:'); dim = int(txt[j + 1])
            sing = 'infeasible' if dim == -1 else 'feasible(dim %d)' % dim
        except Exception: sing = '?'
    rec = {'file': fn, 'n': F.n, 'N': F.N, 'nalpha': na, 'singular': sing, 'nulla': None}
    for d, sub in [(0, None), (1, list(range(na))), (2, list(range(na))), (1, None)]:
        t1 = time.time(); found, msg = nulla(F, d, var_subset=sub)
        rec.setdefault('tries', []).append([d, 'alpha' if sub is not None else 'all', found, msg, round(time.time() - t1, 1)])
        if found:
            rec['nulla'] = 'degree %d (%s multipliers)' % (d, 'alpha-only' if sub is not None else 'all-variable'); break
        if found is None and 'too large' in msg: continue
    rec['time'] = round(time.time() - t0, 1)
    out.write(json.dumps(rec) + '\n'); out.flush()
