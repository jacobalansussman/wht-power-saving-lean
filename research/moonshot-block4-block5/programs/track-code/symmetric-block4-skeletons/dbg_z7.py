from symsolve import *
import sys
# order of candidate P's for Z6
from symlib import gl_mul
m = 4
def order(P):
    Q = P[:]; k = 1
    while Q != [1, 2, 4, 8]:
        Q = gl_mul(P, Q, m); k += 1
        if k > 100: return None
    return k
print("order of Z6rot P:", order([2, 4, 1, 15]), " Z7 P:", order([2, 4, 3, 8]))
R = Run(4, 7, 'Z7rot', 1)
print("|H| =", len(R.Hc))
rng = random.Random(3)
fd = fixed_dirs(R, rng, 3000)
print("fixed dirs Z7rot:", len(fd))
