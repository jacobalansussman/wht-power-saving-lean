from symsolve import *
from symlib import gl_mul
def order(P, m=4):
    Q = P[:]; k = 1
    while Q != [1, 2, 4, 8]:
        Q = gl_mul(P, Q, m); k += 1
        if k > 100: return None
    return k
# find an order-5 element: search small
import itertools
cand = None
for P in itertools.product(range(1, 16), repeat=4):
    span = {0}; ok = True
    for a in P:
        if a in span: ok = False; break
        span |= {x ^ a for x in span}
    if ok and order(list(P)) == 5:
        cand = list(P); break
print("order-5 P:", cand)
import symsearch
src = open('symsearch.py').read().replace("P = [8, 3, 6, 12]; gens = [(P, cyc(1))]", f"P = {cand}; gens = [(P, cyc(1))]")
open('symsearch.py', 'w').write(src)
