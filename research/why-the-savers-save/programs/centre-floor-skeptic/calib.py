#!/usr/bin/env python3
"""v-centre / calib.py: own exactness and one-sided cost of the two calibration decoders.
 hub64: directions e_i (8 terms), gather x_b(i) = r_b[i];  coord72: directions e_2a +- e_2a+1 (scaled), gather
 x_b(k) = (l_k . r_b)/2.  Exact check: sum_k x(k) x(k)^T == Gram (integers / Fractions), on all 112x112 entries,
 hence on the diagonal and every orthogonal pair.  Cost: sum_k dim span{labels of ports with x_b(k) != 0}."""
import json, os
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'data.json')))
L = [tuple(v) for v in D['lines']]; lab = D['lab']
def rank2(vs):
    b = {}
    for x in vs:
        for p in sorted(b, reverse=True):
            if x >> p & 1: x ^= b[p]
        if x: b[x.bit_length() - 1] = x
    return len(b)
ip = lambda a, b: sum(x * y for x, y in zip(a, b))
def run(name, dirs, scale):
    X = [[Fr(ip(l, r)) * scale for l in dirs] for r in L]   # gather factor
    Y = [[Fr(ip(l, r)) for l in dirs] for r in L]           # scatter factor (asymmetric, same supports)
    bad = 0
    for t in range(112):
        for s in range(112):
            v = sum(Y[t][k] * X[s][k] for k in range(len(dirs)))
            if v != ip(L[t], L[s]): bad += 1
    costs = [rank2([lab[b] for b in range(112) if X[b][k] != 0]) for k in range(len(dirs))]
    sizes = [sum(1 for b in range(112) if X[b][k] != 0) for k in range(len(dirs))]
    print(name, 'entries off Gram:', bad, 'terms', len(dirs), 'gather sizes', sizes, 'one-sided costs', costs, 'total', sum(costs))
E = [tuple(1 if j == i else 0 for j in range(8)) for i in range(8)]
run('hub64  ', E, Fr(1))
bf = []
for a in range(4):
    for sg in (1, -1):
        v = [0] * 8; v[2 * a] = 1; v[2 * a + 1] = sg; bf.append(tuple(v))
run('coord72', bf, Fr(1, 2))   # sum_k l_k l_k^T = 2I
