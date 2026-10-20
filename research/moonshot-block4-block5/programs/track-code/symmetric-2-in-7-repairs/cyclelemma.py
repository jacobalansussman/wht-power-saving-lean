"""cyclelemma.py -- mw2-absorb4: exact checks (sympy) of the two group-algebra facts behind the near miss.

(1) CYCLE LEMMA.  Let w generate a cyclic group of order n acting on n cells (the n frames of the filter's terms), and
    let the filter be phi = sum_j a_j w^j (a_0 = 1).  If the junk that the last move leaves is a cyclic shift of the
    filter, J = w phi, and the only junk cancellers are scalar multiples of phi on the same cells, then J = c phi forces
    a_j = c^-j and c^n = 1, so phi = sum_j (c^-1 w)^j = n * (idempotent), and the K-frame block phi(u) (u the image of w
    at the K frame, of order n) is a ZERO DIVISOR: (1 - c^-1 u) phi(u) = 0.  This holds for EVERY order n, odd included.
    n = 2 is the near miss (junk forces a = +-1, K block 1 +- u).

(2) SIGMA SEPARATION (v-ms-fano).  Two DIFFERENT pieces G1 = x1 + x2, G2 = x1 - s x2 with s of order n: Z1 G1 + Z2 G2
    = K x1 needs Z1 = Z2 s (junk) and Z2 (1 + s) = K.  1 + s is invertible in C[Z_n] iff n is odd; then
    (1 + s)^-1 = (1 - s + s^2 - ... + s^(n-1)) / 2.  The junk equation has a FREE unknown (Z2), so it does not constrain
    the coefficients of the pieces: this is what the near miss lacks.

Both checked in the regular representation of Z_n (exact rationals / roots of unity)."""
import sympy as sp
import numpy as np

def shift(n):
    M = sp.zeros(n, n)
    for j in range(n):
        M[(j + 1) % n, j] = 1
    return M

print("(1) cycle lemma: filter phi = sum a_j w^j, junk w*phi must equal c*phi")
for n in range(2, 7):
    a = sp.symbols(f"a1:{n}")
    c = sp.symbols("c")
    vec = sp.Matrix([1, *a])
    eqs = list(shift(n) * vec - c * vec)
    sols = sp.solve(eqs, [*a, c], dict=True)
    sing = []; form = []
    for s in sols:
        cv = complex(sp.N(s[c], 40))
        coeffs = [complex(sp.N(vec[j].subs(s), 40)) for j in range(n)]
        form.append(max(abs(coeffs[j] - cv ** (-j)) for j in range(n)) < 1e-12 and abs(cv ** n - 1) < 1e-12)
        Kb = sum(coeffs[j] * np.linalg.matrix_power(np.roll(np.eye(n), 1, axis=0), j) for j in range(n))
        sing.append(np.linalg.matrix_rank(Kb, tol=1e-9) < n)
    print(f"  n = {n}: {len(sols)} solutions, all of the form a_j = c^-j with c^n = 1: {all(form)}; K block singular (rank 1) for all: {all(sing)}")

print("(2) sigma separation: is 1 + s invertible in C[Z_n], and the inverse")
for n in range(2, 8):
    S = shift(n); I = sp.eye(n)
    d = (I + S).det()
    if d == 0:
        print(f"  n = {n}: 1 + s singular (det 0)")
    else:
        inv = sum(((-1) ** j * S ** j for j in range(n)), sp.zeros(n, n)) / 2
        print(f"  n = {n}: det(1 + s) = {d}; (1 - s + ... + s^(n-1))/2 is the inverse: {(I + S) * inv == I}")
