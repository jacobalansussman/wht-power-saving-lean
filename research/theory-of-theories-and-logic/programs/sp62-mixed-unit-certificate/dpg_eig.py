# dpg_eig.py -- eigenvalues of the frame graph DSp(2m,2) for m = 2, 3 (check of theta_j = 2(2^(m-j) - 1) - (2^j - 1)).
import itertools, numpy as np
for m in (2, 3):
    NV = 1 << (2 * m); MM = (1 << m) - 1
    def Bf(u, v): return (bin((u & MM) & (v >> m)).count('1') + bin((v & MM) & (u >> m)).count('1')) & 1
    def span(vs):
        S = {0}
        for v in vs: S |= {x ^ v for x in S}
        return frozenset(S)
    F = sorted({span(vs) for vs in itertools.combinations(range(1, NV), m) if len(span(vs)) == 1 << m and all(Bf(a, b) == 0 for a in span(vs) for b in span(vs))}, key=sorted)
    Adj = np.array([[1.0 if len(a & b) == 1 << (m - 1) else 0.0 for b in F] for a in F])
    ev = np.round(np.linalg.eigvalsh(Adj), 6)
    vals, cnt = np.unique(ev, return_counts=True)
    print(f"m={m}: {len(F)} frames, eigenvalues {list(zip(vals.tolist(), cnt.tolist()))}; formula {[2*(2**(m-j)-1)-(2**j-1) for j in range(m+1)]}")
