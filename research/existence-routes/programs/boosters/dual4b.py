# e3-boost dual4b.py: smallest fixed sets of rigorous certificates (plus the plain balance) whose cone contains
# phi_theta for EVERY theta in [0, 1] at registers m = 3, 4.  With exactly m generators the coefficients are
# lambda(theta) = G^-1 phi_theta, exponential sums in theta; printed exactly (Fractions) for the exact proof.
import sys, itertools
from fractions import Fraction as F
import numpy as np
import certlib as C

def prof(name, b, m):
    if name == "frame": return [F(1)] * m
    if name == "plain": return [F(d) for d in range(1, m + 1)]
    P = next(P for (n, M, P, s) in C.PROFILES if n == name)
    return [F(P[b * d - 1]) for d in range(1, m + 1)]

short = {3: [("frame", 1), ("plain", 1), ("A3", 1), ("A6", 2), ("A9", 3), ("B3", 1), ("D3", 1), ("C3", 1), ("hook3", 1),
             ("liftE162", 1), ("E3", 1)],
         4: [("frame", 1), ("plain", 1), ("A4", 1), ("A8", 2), ("A12", 3), ("B4", 1), ("C4k2", 1), ("C4k3", 1),
             ("hook4", 1), ("liftg4c", 1), ("proj3to4", 1), ("Dist4S1", 1), ("Dist4S2", 1), ("Dist4S01", 1),
             ("Dist4S12", 1)]}
grid = [i / 400 for i in range(0, 401)]
for m in (3, 4):
    found = []
    for S in itertools.combinations(short[m], m):
        if ("plain", 1) not in S or ("frame", 1) not in S: continue
        Gm = np.array([[float(x) for x in prof(n, b, m)] for n, b in S])
        if abs(np.linalg.det(Gm)) < 1e-12: continue
        inv = np.linalg.inv(Gm.T)
        ok = all((inv @ (np.arange(1, m + 1, dtype=float) ** t) >= -1e-12).all() for t in grid)
        if ok: found.append(S)
    print(f"m={m}: {len(found)} generator sets of size {m} contain the whole curve phi_theta, theta in [0,1] (grid of 401):")
    for S in found[:10]: print("   ", S)
    if found:
        S = found[0]
        G = [prof(n, b, m) for n, b in S]
        # exact inverse of G^T
        Mx = [[G[j][i] for j in range(m)] for i in range(m)]   # columns = generators
        # Gauss-Jordan
        A = [row[:] + [F(int(i == k)) for k in range(m)] for i, row in enumerate(Mx)]
        for c in range(m):
            p = next(r for r in range(c, m) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
            A[c] = [x / A[c][c] for x in A[c]]
            for r in range(m):
                if r != c and A[r][c] != 0: A[r] = [x - A[r][c] * y for x, y in zip(A[r], A[c])]
        inv = [row[m:] for row in A]
        print(f"   exact: lambda_j(theta) = sum_d c_jd d^theta for the set {S}:")
        for j, (n, b) in enumerate(S):
            coeffs = inv[j]
            signs = [1 if x > 0 else (-1 if x < 0 else 0) for x in coeffs]
            ch = sum(1 for a, b2 in zip([s for s in signs if s], [s for s in signs if s][1:]) if a != b2)
            print(f"     {n}x{b}: coefficients on d^theta, d = 1..{m}: {[str(x) for x in coeffs]}; sign changes {ch}")
