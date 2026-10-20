"""export_progs.py (v-pk2-complete): run pk2-complete's OWN program generator (t_units_cases.py, same seeds as
t_units.py) and dump the frame-game programs + their readouts as plain JSON. Nothing else of pk2's code is used
downstream: verification, lift and evaluation are in vcheck.py (my code)."""
import sys, json, random, os
PK = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'frame-game-completion')
sys.path.insert(0, PK)
import warnings; warnings.filterwarnings("ignore")
from lift import *
from t_units_cases import cases
m = int(sys.argv[1])
E = Eng(m); rng = random.Random(7 * m + 1)
out = []
for name, P in cases(E, rng):
    st = P.solve(); assert st == 'yes', (name, st)
    ops = []
    for op in P.ops:
        if op[0] == 'in': ops.append(['in', op[1]])
        elif op[0] == 'mix': ops.append(['mix', [[str(x), list(p), r] for x, p, r in op[1]]])
        else: ops.append(['move', op[1], op[2]])
    read = [[[r, list(p), str(x)] for r, p, x in P.read[l]] for l in range(P.W)]
    out.append(dict(name=name, W=P.W, ops=ops, read=read))
json.dump(out, open(f'progs_m{m}.json', 'w'))
print(m, [(o['name'], sum(1 for op in o['ops'] if op[0] == 'move'), sum(len(r) for r in o['read'])) for o in out])
