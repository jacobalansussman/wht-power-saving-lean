# dim A (span of the free group) in VV and VL2, engine vs widened (CS added), mod p.
import sys, os, json, numpy as np; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common, jac
from common import rep
for game in ('engine', 'widened'):
    p = 2000029; free, Wr = common.letters(2, p); I, h = common.setup(p)
    if game == 'widened': free = free + [np.diag([I if x == 3 else 1 for x in range(4)]).astype(np.int64)]
    Bel = jac.group_elems(free, p)
    for r in ('VV', 'VL2'):
        print(json.dumps({'game': game, 'rep': r, '|B|': len(Bel), 'dimA': jac.span_basis([rep(r, g, p) for g in Bel], p).shape[0]}), flush=True)
