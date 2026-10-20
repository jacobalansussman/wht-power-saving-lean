# dualcert.py m K4 -- exact LP dual certificate for the cover rows of mycover.py (same row builder).
# Solves the dual LP with HiGHS, rationalises y, rescales so that A^T y <= 1 holds EXACTLY in fractions, and prints the
# certified bound b^T y. Since the excesses are integers, b_m >= ceil(certified bound).
import sys, math
from fractions import Fraction as Fr
sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2], "0"]
src = open("mycover.py").read().split("n = len(edges)")[0]
exec(src)
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
n = len(edges)
A = lil_matrix((len(rows), n))
for i, r in enumerate(rows):
    for c, v in r.items(): A[i, c] = v
A = A.tocsr()
res = linprog(-np.array(lb, float), A_ub=A.T, b_ub=np.ones(n), bounds=(0, None), method="highs")
print("dual LP value", -res.fun)
y = [Fr(v).limit_denominator(10000) if v > 1e-12 else Fr(0) for v in res.x]
col = [Fr(0)] * n
for i, r in enumerate(rows):
    if y[i]:
        for c, v in r.items(): col[c] += v * y[i]
M = max(col); assert M > 0
cert = sum(Fr(l) * yi for l, yi in zip(lb, y)) / max(M, Fr(1))
print("max column sum %s (%.9f); certified exact LP lower bound = %s = %.6f -> b >= %d" %
      (M, float(M), cert, float(cert), math.ceil(cert)))
