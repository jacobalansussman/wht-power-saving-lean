# Direct test of the star lemma on structured matrices (e3-slice6): for E = bent circulant on AG(4,2) (special position,
# hyperplane ranks 4), the Singer permutation on AG(3,2), and random integer matrices, check for flats Psi and every
# nonsingular E[R,C] with R, C inside Psi:  (i) rank E[Phi_j,Phi_j] = k + rank T[Phi_j\R, Phi_j\C] (Schur rank formula),
# (ii) 2(|Psi| - k) <= sum_j (rank E[Phi_j,Phi_j] - k) - (g - 2) a.
import itertools
import numpy as np
from fractions import Fraction
rk = lambda X: np.linalg.matrix_rank(X, tol=1e-8) if X.size else 0
def flats(n, d):
    N = 1 << n; out = set()
    for basis in itertools.combinations(range(1, N), d):
        span = {0}
        for b in basis: span |= {x ^ b for x in span}
        if len(span) != 1 << d: continue
        for t in range(N): out.add(frozenset(x ^ t for x in span))
    return [sorted(f) for f in out]
def test(E, n, d, kmax, label, limit=4000):
    N = 1 << n; X = list(range(N)); viol1 = viol2 = checked = 0; tight = 0
    up = flats(n, d + 1)
    for Psi in flats(n, d):
        Phis = [P for P in up if set(Psi) <= set(P)]
        for k in range(0, kmax + 1):
            for R in itertools.combinations(Psi, k):
                for C in itertools.combinations(Psi, k):
                    if checked >= limit: break
                    if k and abs(np.linalg.det(E[np.ix_(R, C)])) < 1e-9: continue
                    rows = [x for x in X if x not in R]; cols = [x for x in X if x not in C]
                    if k:
                        T = E[np.ix_(rows, cols)] - E[np.ix_(rows, C)] @ np.linalg.solve(E[np.ix_(R, C)], E[np.ix_(R, cols)])
                    else:
                        T = E.copy()
                    ri = {x: i for i, x in enumerate(rows)}; ci = {x: i for i, x in enumerate(cols)}
                    sig = []
                    for P in Phis:
                        s_full = rk(E[np.ix_(P, P)])
                        s_T = rk(T[np.ix_([ri[x] for x in P if x not in R], [ci[x] for x in P if x not in C])])
                        if s_full != k + s_T: viol1 += 1
                        sig.append(s_full)
                    a = rk(T[np.ix_([ri[x] for x in Psi if x not in R], [ci[x] for x in Psi if x not in C])])
                    g = len(Phis)
                    gap = sum(sj - k for sj in sig) - (g - 2) * a - 2 * (len(Psi) - k)
                    if gap < 0: viol2 += 1
                    if gap == 0: tight += 1
                    checked += 1
    print(f"{label}: n={n} d={d}: checked {checked} pivots; Schur-rank mismatches {viol1}; star-lemma violations {viol2}; tight cases {tight}")
phi = lambda x: ((x & 1) * (x >> 1 & 1) + (x >> 2 & 1) * (x >> 3 & 1)) % 2
bent = np.array([[(-1.0) ** phi(x ^ y) for y in range(16)] for x in range(16)])
test(bent, 4, 2, 3, "bent circulant (hyperplane ranks 4)")
test(bent, 4, 1, 2, "bent circulant")
sig = (0, 2, 4, 6, 3, 1, 7, 5)
P = np.array([[1.0 if sig[x] == y else 0.0 for y in range(8)] for x in range(8)])
test(P, 3, 1, 2, "Singer permutation (hyperplane ranks 2)")
test(P, 3, 0, 1, "Singer permutation")
rng = np.random.default_rng(3)
for t in range(3):
    E = rng.integers(-1, 2, (16, 16)) * (rng.random((16, 16)) < 0.25) + np.eye(16) * 0
    if abs(np.linalg.det(E)) < 1e-6: E = E + np.eye(16)
    test(E.astype(float), 4, 2, 2, f"random sparse #{t}", limit=1500)
