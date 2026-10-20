"""check to_matrix lifts: symplectic of the exact matrix equals symp_from(P, Bm)"""
import sys, os, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
from sp_check import symp_of
import sp
rng = random.Random(7)
for t in range(40):
    S, P, Bm = sp.letter_full(rng)
    M = sp.to_matrix(P, Bm, rng.randrange(16), rng.randrange(16))
    assert symp_of(M) == S, t
print("lifts ok (40)")
