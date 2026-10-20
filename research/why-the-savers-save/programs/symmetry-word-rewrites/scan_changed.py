"""How many legal single reorders change any frame interval (P or Q) at all, vs purely cosmetic ones."""
import sys, random
sys.path.insert(0, '.')
import fw, search, reorder
from collections import Counter
W = search.Word(fw.load('K8-FULL-WORD.json.gz')); O = reorder.OW(W)
mv = O.candidate_moves(); random.seed(5); random.shuffle(mv)
c = Counter()
for g, p in mv[:int(sys.argv[1])]:
    oldorder = list(O.order); oldP, oldQ = O.P, O.Q
    O.order.remove(g); O.order.insert(p, g); O.refresh()
    if not O.feasible: c['infeasible (P>Q)'] += 1
    else:
        ch = sum(1 for h in O.order if O.P[h] != oldP[h] or O.Q[h] != oldQ[h])
        fr = sum(1 for h in O.order if not (search.le(O.P[h], O.f[h]) and search.le(O.f[h], O.Q[h])))
        c[('changes P/Q' if ch else 'cosmetic', 'current frames still legal' if fr == 0 else 'frames must move')] += 1
    O.order = oldorder; O.refresh()
print(c)
