# Dual polar graphs with classical parameters (d, b, 0, beta): build the intersection array, check the sphere
# sizes, check exactly that theta_j = beta*[d-j]_b - [j]_b are the eigenvalues of the intersection matrix,
# and print the Theorem S'-type cut cap 2/(1 - theta1/k).
from fractions import Fraction as F
def qint(i, b): return sum(b**j for j in range(i))
def det(M):
    M = [row[:] for row in M]; n = len(M); s = F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return F(0)
        if p != c: M[c], M[p] = M[p], M[c]; s = -s
        s *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for j in range(c, n): M[r][j] -= f * M[c][j]
    return s
for name, d, b, beta in [("DSp(4,2), engine block 2", 2, 2, 2), ("DSp(6,2), engine block 3", 3, 2, 2),
                         ("DSp(8,2), engine block 4", 4, 2, 2),
                         ("H(3,4), chi game block 2 (27 frames)", 2, 4, 2), ("H(5,4), chi game block 3 (891 frames)", 3, 4, 2)]:
    bi = [(qint(d, b) - qint(i, b)) * beta for i in range(d)]
    ci = [qint(i, b) for i in range(1, d + 1)]
    k = bi[0]
    L = [[F(0)] * (d + 1) for _ in range(d + 1)]
    for i in range(d + 1):
        if i < d: L[i][i + 1] = F(bi[i])
        if i > 0: L[i][i - 1] = F(ci[i - 1])
        L[i][i] = F(k - (bi[i] if i < d else 0) - (ci[i - 1] if i > 0 else 0))
    ks = [1]
    for i in range(d): ks.append(ks[-1] * bi[i] // ci[i])
    th = [beta * qint(d - j, b) - qint(j, b) for j in range(d + 1)]
    ok = all(det([[L[r][c] - (th[j] if r == c else 0) for c in range(d + 1)] for r in range(d + 1)]) == 0 for j in range(d + 1))
    cap = F(2) / (1 - F(th[1], k))
    print(f"{name}: vertices={sum(ks)} spheres={ks} eigenvalues={th} exact={ok} theta1/k={F(th[1], k)} cut cap={cap} = {float(cap):.4f}")
