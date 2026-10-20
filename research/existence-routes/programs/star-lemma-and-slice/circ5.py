# Circulant candidates for the local problem L(5, s) (block 7): E = H diag((-1)^psi) H / 32 (circulant, +-1 spectrum).
# Rank on a hyperplane of direction h = #{pairs {a, a+h}: psi(a) = psi(a+h)} = (32 + C(h))/4, C = autocorrelation.
# Need max_{h != 0} C(h) <= 4s - 32.  s = 7 is impossible for circulants (sum of C(h) over h != 0 is S^2 - 32 >= -32).
# Search (annealing) for psi with max C(h) <= 0 (s = 8).
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
n = 5; N = 32
X = np.arange(N)
def autoc(t):
    s = 1 - 2 * t
    return np.array([int(np.dot(s, s[X ^ h])) for h in range(1, N)])
best = None
for run in range(200):
    t = rng.integers(0, 2, N)
    cur = autoc(t).max(); T = 4.0
    for it in range(4000):
        j = rng.integers(N); t[j] ^= 1
        new = autoc(t).max()
        if new <= cur or rng.random() < np.exp((cur - new) / T): cur = new
        else: t[j] ^= 1
        T *= 0.999
        if cur <= 0: break
    if best is None or cur < best[0]: best = (cur, t.copy())
    if cur <= 0: break
c, t = best
print("best max autocorrelation", c, "-> hyperplane rank", (32 + c) // 4, "truth table", "".join(map(str, t)))
if c <= 0:
    Hd = np.array([[(-1.0) ** bin(a & x).count("1") for x in range(N)] for a in range(N)])
    E = Hd @ np.diag((-1.0) ** t) @ Hd / N          # circulant with spectrum (-1)^psi: orthogonal, nonsingular
    E = np.round(E * N).astype(int)                 # integer multiple, exact ranks below
    from itertools import product
    worst = 0
    for h in range(1, N):
        for b in (0, 1):
            H = [x for x in range(N) if bin(h & x).count("1") % 2 == b]
            worst = max(worst, np.linalg.matrix_rank(E[np.ix_(H, H)]))
    print("check: rank E =", np.linalg.matrix_rank(E), "max hyperplane rank =", worst)
    # exact ranks over Q
    from fractions import Fraction
    def qrank(M):
        A = [[Fraction(int(v)) for v in row] for row in M]; r = 0
        for c in range(len(A[0])):
            p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
            if p is None: continue
            A[r], A[p] = A[p], A[r]
            for i in range(len(A)):
                if i != r and A[i][c] != 0:
                    f = A[i][c] / A[r][c]; A[i] = [a - f * b for a, b in zip(A[i], A[r])]
            r += 1
        return r
    ex = [qrank(E[np.ix_([x for x in range(N) if bin(h & x).count("1") % 2 == b], [x for x in range(N) if bin(h & x).count("1") % 2 == b])])
          for h in range(1, N) for b in (0, 1)]
    print("EXACT: rank E =", qrank(E), "; hyperplane ranks over Q: min", min(ex), "max", max(ex), "(62 hyperplanes)")
