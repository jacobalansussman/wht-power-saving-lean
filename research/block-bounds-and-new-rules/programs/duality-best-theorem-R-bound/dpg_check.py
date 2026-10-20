# dpg_check.py -- at the optimum of the whole-graph block-3 LP (cuts of dpg.py), evaluate EVERY DFZ instance whose
# letters are seven neighbours of L0, one for each point of the Fano plane L0^* (2 choices per point, 2^7 choices),
# under all 168 ordered bases (A,B,C) of the plane, with T = {} and T = {L0}; and Ingleton instances on random
# 4-sets of frames with T in {}, {L0}.  Reports the most violated values (negative = a cut the search missed).
import sys, os, json, itertools, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from dpg import DPG, omega
from lpm import DFZ, ING, inst_row, solve
from cutsearch import Evaluator
I = DPG(5)
data = json.load(open(sys.argv[1]))
cuts = [(nm, {int(k): v for k, v in r.items()}) for nm, r in data["cuts"]]
rows = I.polymatroid_rows(); eqs = I.uniform_eqs()
out = solve(I, rows, eqs, cuts, uniform=True)
print("re-solved optimum:", out["opt"]); x = out["x"]
L0 = I.frame_of([8, 16, 32]); S0 = I.SUBS[L0]
nbr = [i for i in range(I.NF) if I.D[L0, i] == 1]
# quotient point of a neighbour: the 1-dim (L+L0)/L0, identified by the x-part of any vector of L outside L0
def qpoint(i):
    for u in I.SUBS[i]:
        if u not in S0: return u & 7
byp = {}
for i in nbr: byp.setdefault(qpoint(i), []).append(i)
print("neighbours per point:", {p: len(v) for p, v in byp.items()})
E = Evaluator(I); worst = []
for A in range(1, 8):
    for B in range(1, 8):
        for C in range(1, 8):
            if len({0, A, B, C, A ^ B, A ^ C, B ^ C, A ^ B ^ C}) < 8: continue
            pts = [A, B, C, A ^ B, A ^ C, B ^ C, A ^ B ^ C]          # A,B,C,W=A+B,X=A+C,Y=B+C,Z=A+B+C
            for choice in itertools.product(range(2), repeat=7):
                letters = [byp[p][c] for p, c in zip(pts, choice)]
                for T in ((), (L0,)):
                    v = E.value(DFZ, letters, T, x)
                    worst.append((v, tuple(letters), T))
worst.sort()
print("DFZ over L0-neighbours: instances", len(worst), " min value", worst[0][0], " max", worst[-1][0])
rng = random.Random(5); ws = []
for t in range(200000):
    L = rng.sample(range(I.NF), 4); T = rng.choice([(), (L0,)])
    v = E.value(ING, L, T, x)
    if v is not None: ws.append(v)
print("Ingleton random 200000: min value", min(ws))
