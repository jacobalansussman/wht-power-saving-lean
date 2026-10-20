"""exact_plain.py (rm-singer): EXACT check (rational arithmetic) that the canonical Singer system accepts the plain
schedule and refuses broken variants.  Plain schedule for one array at block m, all moves in direction 1 with free
rotations: Atil_t = P^{r_t}, Btil_t = P^{-r_t}, Gamma_{s,t} = 2 P^{r_t - r_s} (s > t), where alpha^{r_t} = w_t runs
over an orthonormal basis of the chosen coordinates (the self-dual basis in 'sd').  Every fibre must give
Atil[:,D] (I + Gamma[D,D])^-1 Btil[D,:] = P_J exactly.
Controls: (a) dropping one move; (b) a non-orthonormal basis (w_1 replaced by w_1 + w_2): both must fail somewhere."""
import sys
from fractions import Fraction as Fr
from sfield import Singer

def solve(Lm, Bm):
    """exact L^-1 B by Gauss-Jordan (L square)."""
    n = len(Lm); k = len(Bm[0]) if Bm else 0
    A = [row[:] + brow[:] for row, brow in zip(Lm, Bm)]
    for c in range(n):
        p = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[p] = A[p], A[c]
        piv = A[c][c]; A[c] = [x / piv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]

def check(m, coords, rexps, label):
    S = Singer(m, coords); F = S.F; n = S.n; M = len(rexps)
    shift = lambda r: [[Fr(int((k - kk - r) % n == 0)) for kk in range(n)] for k in range(n)]  # P^r: |kk> -> |kk+r>
    # big matrices: Atil (n x Mn), Btil (Mn x n), Gamma (Mn x Mn)
    At = [[Fr(0)] * (M * n) for _ in range(n)]
    Bt = [[Fr(0)] * n for _ in range(M * n)]
    Gt = [[Fr(0)] * (M * n) for _ in range(M * n)]
    for t, r in enumerate(rexps):
        Pr, Pm = shift(r), shift(-r)
        for k in range(n):
            for kk in range(n):
                At[k][t * n + kk] = Pr[k][kk]; Bt[t * n + k][kk] = Pm[k][kk]
        for s in range(t + 1, M):
            Q_ = shift(r - rexps[s])
            for k in range(n):
                for kk in range(n):
                    Gt[s * n + k][t * n + kk] = 2 * Q_[k][kk]
    seen = set(); bad = 0; nf = 0
    for Q in S.forms():
        if not any(Q): continue
        v = tuple(S.Qval(Q, x) for x in range(1, F.q))
        if v in seen: continue
        seen |= {tuple(S.Qval(Q, F.mul(F.a_exp(-e), x)) for x in range(1, F.q)) for e in range(n)}
        nf += 1
        D = [t * n + k for t in range(M) for k in range(n) if S.Qval(Q, F.a_exp(-k))]
        J = {k for k in range(n) if S.pair(Q, S.Tk(k))}
        T = [[Fr(int(i == j and i in J)) for j in range(n)] for i in range(n)]
        if not D:
            X = [[Fr(0)] * n for _ in range(n)]
        else:
            L = [[Fr(int(i == j)) + Gt[D[i]][D[j]] for j in range(len(D))] for i in range(len(D))]
            Y = solve(L, [Bt[i] for i in D])
            X = [[sum(At[i][D[a]] * Y[a][j] for a in range(len(D))) for j in range(n)] for i in range(n)]
        if X != T: bad += 1
    print(f"m={m} {coords} {label}: moves {M}, fibres {nf}, fibres violated (exact) {bad}")
    return bad

for m in [int(a) for a in sys.argv[1:]] or [3, 4]:
    for coords in ('sd', 'poly'):
        S = Singer(m, coords); F = S.F
        r = [F.log[w] for w in S.w]
        assert check(m, coords, r, "plain (orthonormal basis)") == 0
        assert check(m, coords, r[1:], "plain minus one move") > 0
        w1 = S.w[0] ^ S.w[1]
        assert check(m, coords, [F.log[w1]] + r[1:], "non-orthonormal basis") > 0
print("ALL AS EXPECTED")
