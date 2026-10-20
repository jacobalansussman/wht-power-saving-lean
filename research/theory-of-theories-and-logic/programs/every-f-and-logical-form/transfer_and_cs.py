"""transfer_and_cs.py -- v-ml-logic (skeptic of meta3/ml-logic), 14 Oct 2026.
Part 1. Counterexample to ml-logic PROOF.md 4.7(b) / NOTE.md section 4 (5): "An approximate family whose residual
  tends to 0 transfers only if the residual is exactly 0 in the ultrapower"; "Approximate families with a nonzero
  infinitesimal error do not transfer".
  Shape (block 1, one register): y -> D y (one paid move), then r := mu * r.  Unit condition: mu D^(x)f = D^(x)f for
  all f, i.e. mu = 1.  Family mu_k = 1 + 1/k: residual nonzero for every k and -> 0, so in the ultrapower the
  residual is a NONZERO infinitesimal; yet the family is bounded, its standard part mu = 1 is an exact unit of the
  same shape.  A less trivial bounded family: a 2-parameter shape r := mu*D y + nu*(X D X) y with mu + nu = 1 a
  line of exact units (X D X = D since D commutes with X); the family (mu, nu) = (1/2 + 1/k, 1/2) has nonzero
  residual -> 0 and converges to the exact unit (1/2, 1/2).  General fact: the zero set of the unit equations is
  closed, so any BOUNDED family with residual -> 0 has an exact limit point (no ultrapower needed); only unbounded
  (border-type) families can fail to give a unit, and stability theorems are a second transfer route.
Part 2. The widened group <engine letters, CS> at block 2 is infinite: exhibit a word of infinite order
  (an eigenvalue that is not a root of unity: its minimal polynomial over Q is not cyclotomic).
Usage: python3 -I transfer_and_cs.py   (tiny)
"""
import sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, itertools
a, b = (1 + 1j) / 2, (1 - 1j) / 2
D = np.array([[a, b], [b, a]]); X = np.array([[0, 1], [1, 0]], complex)
def kp(M, f):
    out = np.ones((1, 1), complex)
    for _ in range(f): out = np.kron(out, M)
    return out
print("Part 1: family mu_k = 1 + 1/k on the shape 'D then scale by mu'")
for k in (1, 10, 100, 1000, 10**6):
    mu = 1 + 1 / k
    res = [np.linalg.norm(mu * kp(D, f) - kp(D, f)) for f in (1, 2, 3)]
    print(f"  k = {k:>7}: residual at f = 1,2,3: " + ", ".join(f"{r:.3e}" for r in res) + "  (nonzero, -> 0)")
print("  limit mu = 1: residual", max(np.linalg.norm(1 * kp(D, f) - kp(D, f)) for f in (1, 2, 3)), "(exact unit)")
print("  X D X == D:", np.allclose(X @ D @ X, D))
for k in (10, 1000, 10**6):
    mu, nu = 0.5 + 1 / k, 0.5
    res = [np.linalg.norm(mu * kp(D, f) + nu * kp(X @ D @ X, f) - kp(D, f)) for f in (1, 2, 3)]
    print(f"  2-parameter family k = {k:>7}: residual " + ", ".join(f"{r:.3e}" for r in res))

print("Part 2: an element of infinite order in <engine letters, CS> at block 2")
I2 = np.eye(2)
L = {"D0": np.kron(I2, D), "D1": np.kron(D, I2), "CS": np.diag([1, 1, 1, 1j]),
     "S0": np.kron(I2, np.diag([1, 1j])), "CX": np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]], complex)}
names = list(L)
found = None
for n in range(1, 5):
    for w in itertools.product(names, repeat=n):
        g = np.eye(4, dtype=complex)
        for s in w: g = L[s] @ g
        ev = np.linalg.eigvals(g)
        # finite order <=> all eigenvalues roots of unity; for a matrix over Z[i,1/2] of finite order the order
        # divides a bounded number; test g^N = I for N up to 10^4 is replaced by an exact char-poly test below
        cp = np.poly(g)            # characteristic polynomial (complex coefficients in Z[i,1/2])
        if np.max(np.abs(np.abs(ev) - 1)) > 1e-9: continue
        args = np.angle(ev) / (2 * np.pi)
        # rational with small denominator?
        ok = all(min(abs(x * q - round(x * q)) for q in range(1, 241)) < 1e-9 for x in args)
        if not ok:
            found = (w, ev, cp); break
    if found: break
w, ev, cp = found
print("  word:", w, " eigenvalue angles / 2pi:", np.round(np.angle(ev) / (2 * np.pi), 6))
print("  char poly coefficients:", np.round(cp, 6))
# exact: a 4x4 unitary over Q(i) of finite order has eigenvalues that are roots of unity of degree <= 8 over Q,
# i.e. of order n with phi(n) <= 8, all dividing 5040; so finite order <=> g^5040 = I.  Exact Gaussian-dyadic powers.
import engine_group_exact as E
def to_exact(M):
    N = [[(int(round(M[r, c].real * 4)), int(round(M[r, c].imag * 4))) for c in range(4)] for r in range(4)]
    assert np.allclose(np.array([[complex(*e) for e in row] for row in N]) / 4, M)
    return E.reduce_((N, 2))
Ge = to_exact(np.eye(4))
for s_ in w: Ge = E.matmul(to_exact(L[s_]), Ge)
def mpow(A, n):
    R = (E.ident(4), 0)
    while n:
        if n & 1: R = E.matmul(R, A)
        A = E.matmul(A, A); n >>= 1
    return R
P = mpow(Ge, 5040)
isI = E.key(P) == E.key((E.ident(4), 0))
print("  exact: g^5040 == I ?", isI, "; 2-adic denominator exponent of g^5040:", P[1])
print("  => infinite order:", not isI)
