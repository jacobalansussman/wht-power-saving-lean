# full4.py -- block 4, three disjoint opposite pairs that are NOT compatible: fractional cut value on the image graphs
# and on Gamma_4 itself (undirected cut LP by constraint generation, exact rational dual). Is it exactly 12?
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # modules of ../
import sys, time, random, itertools, json
import w3cut, w3compat, w3image
from w3screen import exact_value
t0 = time.time(); g = w3cut.G(4); rng = random.Random(int(sys.argv[2])); out = []
k = 0
while k < int(sys.argv[1]):
    P = w3cut.rand_3k2(g, rng)
    if w3compat.any_compat(g.D, 4, P): continue
    k += 1; rec = dict(P=P)
    vf = exact_value(w3image.cut_value(g.n, g.edges, P)); rec['gamma4'] = float(vf)
    print(f"{P}: Gamma_4 exact {float(vf):.6f} ({time.time()-t0:.0f}s)", flush=True)
    out.append(rec); json.dump(out, open(f"full4_s{sys.argv[2]}.json", "w"))
