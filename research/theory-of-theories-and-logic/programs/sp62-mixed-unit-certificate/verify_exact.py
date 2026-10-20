# verify_exact.py -- independent exact check of exact_cert.py's ranks: plain Gaussian elimination over Q with Python
# Fractions (no fraction-free tricks) on the saved integer matrices [N, gN]_R, N_R and the basis B_R (whose full rank makes
# the row restriction injective on X).   Usage: verify_exact.py exact_mats_lam3.json
import json, sys
from fractions import Fraction
def rank_Q(rows):
    A = [[Fraction(x) for x in r] for r in rows]; n = len(A); m = len(A[0]) if n else 0; r = 0
    for c in range(m):
        p = next((i for i in range(r, n) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]; pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(n):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        r += 1
        if r == n: break
    return r
S = json.load(open(sys.argv[1]))
for name, D in S.items():
    rB = rank_Q(D['B_R']); rN = rank_Q([list(x) for x in zip(*D['N'])])
    out = [f"{name}: rank B_R = {rB} (must equal dim X), dim N = {rN}"]
    psi = []
    for g in ('h', 'h2', 'k'):
        M = D[g]; r = rank_Q([list(x) for x in zip(*M)]); psi.append(r - rN)
    out.append(f"profile (psi1, psi2, psi3) = {tuple(psi)}; psi3 + 4psi2 - 10psi1 = {psi[2] + 4*psi[1] - 10*psi[0]}; 7psi3 + 38psi2 - 88psi1 = {7*psi[2] + 38*psi[1] - 88*psi[0]}")
    print("; ".join(out), flush=True)
