"""Exact frame re-optimisation (fixed order, every subspace in [P,Q]) of any word; writes the re-framed word."""
import sys
sys.path.insert(0, '.')
import fw, search, ilp
from ilp_wide import all_between
w = fw.load(sys.argv[1])
W = search.Word(w)
c0 = W.total_cost(); p0 = dict(sorted(W.profile().items()))
P, Q = ilp.bounds_PQ(W)
gaps = [len(Q[k]) - len(P[k]) for k in range(len(W.seq))]
extra = {k: all_between(P[k], Q[k]) for k in range(len(W.seq)) if 0 < gaps[k] <= 4}
print('max gap', max(gaps), 'positions with gap > 4:', sum(1 for g in gaps if g > 4))
fr, info = ilp.solve(W, extra_pool=extra, verbose=False)
print('word', sys.argv[1], 'cost0', round(c0, 4), p0)
print(info)
if fr:
    W.f = fr; c1 = W.total_cost(); print('cost1', round(c1, 4), 'gain', round(c0 - c1, 4), dict(sorted(W.profile().items())))
    if len(sys.argv) > 2 and c1 < c0 - 1e-9: fw.save(W.to_word(), sys.argv[2]); print('saved', sys.argv[2])
