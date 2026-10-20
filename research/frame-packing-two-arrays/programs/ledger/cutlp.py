#!/usr/bin/env python3
"""cutlp.py (sl-combinatorics): the cut-condition LP on the block-2 frame graph Gamma_2 (15 Lagrangians).
For a configuration of W opposite (source, sink) pairs: minimise sum_e c_e subject to c(delta(Q)) >= #pairs split by Q
for every vertex set Q.  PROOF.md section 4: every coherent unit has M >= this LP value.  W = 2: Hu's theorem says the
value is 2m = 4 for every configuration (checked here on all 120 configurations with a fixed first pair).  W >= 3:
the value can drop below 2 per pair (flow-cut gap), so cuts alone cannot exclude coherent savers.
usage: python3 cutlp.py [samples]"""
import sys, random, itertools
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from cutconst import setup

def main():
    samples = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    m = 2
    labs, L0, D = setup(m)
    V = len(labs)
    edges = [(i, j) for i in range(V) for j in range(i + 1, V) if D[i, j] == 1]
    opp = [(i, j) for i in range(V) for j in range(V) if D[i, j] == m]
    masks = list(range(1, 1 << (V - 1)))
    Qs = np.array([[(mk >> v) & 1 for v in range(V)] for mk in masks], dtype=bool)
    # constraint matrix rows: -(cut edges)  <= -(demand)
    Acut = np.array([[Qs[k, i] != Qs[k, j] for (i, j) in edges] for k in range(len(masks))], dtype=float)
    def lp(pairs):
        dem = np.array([sum(Qs[k, s] != Qs[k, t] for (s, t) in pairs) for k in range(len(masks))], dtype=float)
        r = linprog(np.ones(len(edges)), A_ub=-Acut, b_ub=-dem, bounds=[(0, None)] * len(edges), method="highs")
        assert r.status == 0
        return r.fun
    p0 = opp[0]
    vals2 = sorted(set(round(lp([p0, p]), 9) for p in opp))
    print("W=2: LP values over all %d configurations with a fixed first pair: %s (2m = 4)" % (len(opp), vals2))
    rng = random.Random(3)
    for W in (3, 4, 6, 8, 12):
        best = None
        for _ in range(samples):
            pairs = [rng.choice(opp) for _ in range(W)]
            v = lp(pairs)
            if best is None or v < best[0]: best = (v, pairs)
        print("W=%d: smallest cut-LP value over %d random configurations: %.4f = %.4f per pair (plain 2; engine floor 2)"
              % (W, samples, best[0], best[0] / W), flush=True)
    # the symmetric configuration: all 120 ordered opposite pairs at once
    v = lp(opp)
    print("all %d ordered opposite pairs: LP %.4f = %.4f per pair (c_2 = 1.5)" % (len(opp), v, v / len(opp)))

if __name__ == "__main__":
    main()
