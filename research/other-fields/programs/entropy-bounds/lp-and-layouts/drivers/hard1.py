# hard1.py -- one hill-climbed block-4 layout with cross deficiency e = 3 whose three-end images all fall short
# (best 10.75 <= 11). Escalate: images to 4 and 5 ends, all 6 ends, then Gamma_4 itself (undirected cut LP).
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # modules of ../
import sys, time, itertools, json
import w3cut, w3image, edef
from w3screen import exact_value
t0 = time.time(); g = w3cut.G(4)
P = [(331, 751), (1319, 596), (837, 1071)]
print("e =", edef.edef(g.D, 4, P), "terminal distances", [[int(g.D[a, b]) for b in [x for p in P for x in p]] for a in [x for p in P for x in p]], flush=True)
ends = [x for p in P for x in p]; out = {}
for k in (6,):
    best = 0; bestset = None
    for S in itertools.combinations(range(6), k):
        if not all(any(i in S for i in (2 * j, 2 * j + 1)) for j in range(3)): continue
        us = [ends[i] for i in S]
        nH, EH, lab, keys = w3image.image_graph(g, us)
        v = exact_value(w3image.cut_value(nH, EH, [(lab[s], lab[t]) for s, t in P]))
        if v > best: best, bestset = v, (S, nH)
        print(f"  ends {S}: {nH} points, exact value {float(v):.5f} ({time.time()-t0:.0f}s)", flush=True)
        if v > 11: break
    out[k] = (float(best), bestset); json.dump(out, open("hard1.json", "w"), default=str)
    print(f"k={k}: best {float(best):.5f}", flush=True)
    if best > 11: break
if False:  # Gamma_4 itself is left to job_hardpd (directed PDLP)
    v = exact_value(w3image.cut_value(g.n, g.edges, P)); out['gamma4'] = float(v); json.dump(out, open("hard1.json", "w"), default=str)
    print(f"Gamma_4: exact value {float(v):.5f} ({time.time()-t0:.0f}s)", flush=True)
