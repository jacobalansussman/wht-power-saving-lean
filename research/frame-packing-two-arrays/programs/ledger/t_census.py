#!/usr/bin/env python3
"""t_census.py (sl-combinatorics): census of ALL exact units among words with exactly `nm` paid moves and up to
`k` other letters (free letters, two-register mixes, 2x2 gates), block m, W = 2, R registers.  For each unit: the
ledger in the trivial gauge, the unit-form identity, and the coherence class.  usage: python3 t_census.py m R nm k"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ledger import *

def main():
    m = int(sys.argv[1]); R = int(sys.argv[2]); nm = int(sys.argv[3]); k = int(sys.argv[4])
    W = 2
    E = Engine(m)
    moves = [("move", r, z) for r in range(R) for z in range(1, 1 << m)]
    frees = ["X%d" % ((1 << m) - 1), "P%d_0" % ((1 << m) - 1), "P%d_1" % ((1 << m) - 1)]
    other = [("free", r, n) for r in range(R) for n in frees]
    coefs = [(1, 1), (1, -1), (1, GQ(0, 1)), (1, GQ(0, -1)), (-2, 1), (Fr(-1, 2), Fr(1, 2)), (Fr(1, 2), Fr(1, 2)),
             (0, 1), (GQ(0, 1), 0)]
    for r in range(R):
        for s in range(R):
            if s == r: continue
            for (a, b) in coefs:
                terms = ([(r, a, None)] if a != 0 else []) + ([(s, b, None)] if b != 0 else [])
                other.append(("mix", r, terms))
    h = Fr(1, 2)
    gates = [[[1, 1], [1, -1]], [[h, h], [h, -h]], [[1, GQ(0, 1)], [1, GQ(0, -1)]], [[h, h], [GQ(0, -h), GQ(0, h)]],
             [[0, 1], [1, 0]]]
    for i in range(R):
        for j in range(i + 1, R):
            for A in gates: other.append(("gate", [i, j], A))
    t0 = time.time(); n = 0; bad = 0; units = []
    for L in range(nm, nm + k + 1):
        for pos in itertools.combinations(range(L), nm):
            for mv in itertools.product(moves, repeat=nm):
                for ot in itertools.product(other, repeat=L - nm):
                    word = []; mi = iter(mv); oi = iter(ot)
                    for p in range(L): word.append(next(mi) if p in pos else next(oi))
                    res = run(E, W, R, word)
                    n += 1
                    if not res["raw_ok"]: bad += 1
                    if res["unit"]:
                        if not res["unit_ok"]: bad += 1
                        gc, why = find_coherent_gauge(E, W, R, word)
                        units.append(dict(word=str(word), M=res["M"], P=res["P"], B=res["B"], D=res["D"], Z=res["Z"],
                                          coherence=why.split(" (")[0]))
        print("L=%d done: %d runs, %d bad, %d units, %.0f s" % (L, n, bad, len(units), time.time() - t0), flush=True)
    cls = {}
    for u in units: cls[(u["M"], u["coherence"], u["P"], u["B"])] = cls.get((u["M"], u["coherence"], u["P"], u["B"]), 0) + 1
    print("alphabet: %d moves, %d other letters; runs %d; identity violations %d; units %d" % (len(moves), len(other), n, bad, len(units)))
    for kk in sorted(cls): print("  M=%d %-16s P=%d B=%d : %d units" % (kk[0], kk[1], kk[2], kk[3], cls[kk]))
    inco = [u for u in units if u["coherence"] != "coherent"]
    for u in inco[:6]: print("  incoherent example:", u)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "t_census_m%d_R%d_nm%d_k%d.json" % (m, R, nm, k)), "w") as f:
        json.dump(dict(runs=n, bad=bad, nunits=len(units), classes={str(a): b for a, b in cls.items()}, incoherent=inco[:50]), f, indent=1)

if __name__ == "__main__":
    main()
