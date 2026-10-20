"""own construction of the two symmetry groups at block 4: Z7 (order-7 element of GL(4,2) fixing x = e4, lifted into
H_max) and Z2^3 (transvections y -> y + f(y) x, f(x) = 0, lifted).  Checks: lifts in H_max, homomorphism up to scalars."""
import sys, itertools, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from grp import Grp
G = Grp(4); m, n = 4, 16
x = 8  # e4
def glmat_mul(A, B): return [G.apply_gl(A, b) for b in B]
def gl_order(A):
    I = [1, 2, 4, 8]; P = A; k = 1
    while P != I: P = glmat_mul(A, P); k += 1
    return k
def scal(M):  # is M a scalar multiple of I (unit modulus)?
    d = M[0, 0]
    return abs(abs(d) - 1) < 1e-9 and np.allclose(M, d * np.eye(n))
def z7():
    # Singer cycle on <e1,e2,e3>: companion of t^3 + t + 1, e4 fixed
    A = [2, 4, 1 ^ 2, 8]
    assert gl_order(A) == 7
    c = G.hmax_lift(A)
    els = [np.eye(n, dtype=complex)]
    for h in range(1, 7): els.append(G.mul(els[-1], c))
    return els, G.mul(els[-1], c)
def z2cube():
    # f in dual with f(e4) = 0: f = (f1,f2,f3,0); transvection T_f(y) = y + f(y) e4
    lifts = {}
    for j in range(3):
        A = [1, 2, 4, 8]; A[j] ^= 8
        lifts[j] = G.hmax_lift(A)
    els = []
    for bits in itertools.product([0, 1], repeat=3):
        M = np.eye(n, dtype=complex)
        for j in range(3):
            if bits[j]: M = G.mul(M, lifts[j])
        els.append(M)
    return els, lifts
if __name__ == '__main__':
    e7, c7 = z7()
    print('Z7: c^7 scalar?', scal(c7), 'c in H_max (K^-1 c^-1 K monomial):', G.is_mono(G.mul(G.Ki, G.inv(e7[1]), G.K)))
    e8, L = z2cube()
    comm = [[scal(G.mul(L[i], L[j], G.inv(L[i]), G.inv(L[j]))) for j in range(3)] for i in range(3)]
    sq = [scal(G.mul(L[i], L[i])) for i in range(3)]
    print('Z2^3: generators commute up to scalars', comm, 'squares scalar', sq)
