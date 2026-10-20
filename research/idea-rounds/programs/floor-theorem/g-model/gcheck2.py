# g-model check 2: weight of an ARBITRARY column-uniform paid operator P (f = 2) against the kernel's.
import random, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from gcheck_lib import *
random.seed(5)
for q in (2, 3, 4):
    n = q * q; Pd = [1 if a // q == a % q else 0 for a in range(n)]
    best = 0
    for _ in range(6):
        P = [[complex(random.randint(-3, 3), random.randint(-3, 3)) for _ in range(q)] for _ in range(q)]
        best = max(best, rank_exact(comm(kron(P, P), Pd)))
    print(" q=%d: largest weight of a random paid operator P(x)P: %d exact; ceiling 2q = %d; kernel (q a power of 2) weighs 2(q-1) = %d" % (q, best, 2 * q, 2 * (q - 1)))
