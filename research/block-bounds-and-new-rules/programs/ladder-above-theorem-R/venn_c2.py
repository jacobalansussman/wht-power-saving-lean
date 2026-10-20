# pf-ladder: test of the guess "c2 <= 1.8 c1" on two-factor Venn certificates at block 3 (widened game kept):
# Lambda^a V (x) Lambda^b V, N_J = span{ e_S (x) e_T : |S cap T| in J }.  c_r = rank of the (not N) x N corner of B_r.
# Float ranks (relative tol 1e-8); a guess test, not a proof.
import sys  # needs numpy
import itertools, numpy as np
n = 8
def dirm(z):
    M = np.zeros((n, n), complex)
    for x in range(n): M[x, x] += (1 + 1j) / 2; M[x ^ z, x] += (1 - 1j) / 2
    return M
def wedge(g, k):
    subs = list(itertools.combinations(range(n), k))
    W = np.array([[np.linalg.det(g[np.ix_(S, T)]) if k else 1 for T in subs] for S in subs], complex)
    return W, subs
mats = {1: dirm(1), 2: dirm(1) @ dirm(2), 3: dirm(1) @ dirm(2) @ dirm(4)}
best = []
for a in range(1, 4):
    for b in range(a, 7 - a):
        Wa = {r: wedge(g, a) for r, g in mats.items()}; Wb = {r: wedge(g, b) for r, g in mats.items()}
        sa, sb = Wa[1][1], Wb[1][1]
        inter = np.array([len(set(S) & set(T)) for S in sa for T in sb])
        Gs = {r: np.kron(Wa[r][0], Wb[r][0]) for r in mats}
        vals = sorted(set(inter))
        for size in range(1, len(vals)):
            for J in itertools.combinations(vals, size):
                inN = np.isin(inter, J); Nc = np.nonzero(inN)[0]; Oc = np.nonzero(~inN)[0]
                c = []
                for r in (1, 2, 3):
                    Cr = Gs[r][np.ix_(Oc, Nc)]
                    s = np.linalg.svd(Cr, compute_uv=False)
                    c.append(int((s > 1e-8 * max(1.0, s.max())).sum()))
                if c[0] > 0:
                    best.append((c[1] / c[0], c[2] / c[0], a, b, J, c))
best.sort(reverse=True)
print("largest c2/c1 among two-factor Venn certificates Lambda^a (x) Lambda^b, a <= b, a + b <= 6 (block 3):")
for row in best[:8]:
    print(f"  c2/c1 = {row[0]:.4f}, c3/c1 = {row[1]:.4f}, (a, b) = ({row[2]}, {row[3]}), J = {row[4]}, (c1, c2, c3) = {row[5]}")
print("largest c3/c1:")
for row in sorted(best, key=lambda r: -r[1])[:4]:
    print(f"  c3/c1 = {row[1]:.4f}, c2/c1 = {row[0]:.4f}, (a, b) = ({row[2]}, {row[3]}), J = {row[4]}, (c1, c2, c3) = {row[5]}")
print(f"certificates scanned: {len(best)}; any with c2/c1 > 1.8 + 1e-9: {any(r[0] > 1.8 + 1e-9 for r in best)}")
