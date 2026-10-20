"""s2_index_code.py -- the crossing-slack centre read as a linear index code.

For a simply-laced root system: ports = root lines (N), decoder = Gram matrix / 2 (rank d), the label of a port is a
binary vector u_a with u_a.u_b = 1 + <a,b> (mod 2) for a != b and u_a.u_a = 1 (Alon's complementarity: F_2-orthogonal
<=> Euclidean NON-orthogonal, so the junk falls on cheap pairs). The label Gram U = J + G (mod 2); its F_2-rank is the
window w. Copy ratio = d * w / N.

Checks: (1) reproduce R35's copy ratios (E8 0.600, E7 0.778, E6 1.167, D_n and A_n >= 1.0);
        (2) alpha = largest set of pairwise orthogonal lines = independence number of the side-information graph
            (cheap = non-orthogonal pairs); Haemers: minrank >= alpha, so the Gram decoder (rank d) is optimal
            exactly when alpha = d ("perfect index code");
        (3) for E8 and E7 the side-information graph is strongly regular; Lovasz theta from its eigenvalues.
"""
import itertools
import numpy as np

def roots_A(n):  # A_n in R^{n+1}
    R = []
    for i, j in itertools.permutations(range(n + 1), 2):
        v = [0] * (n + 1); v[i] = 1; v[j] = -1; R.append(v)
    return np.array(R, dtype=float)

def roots_D(n):
    R = []
    for i, j in itertools.combinations(range(n), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * n; v[i] = si; v[j] = sj; R.append(v)
    return np.array(R, dtype=float)

def roots_E8():
    R = list(roots_D(8))
    for s in itertools.product((0.5, -0.5), repeat=8):
        if sum(1 for x in s if x < 0) % 2 == 0: R.append(list(s))
    return np.array(R, dtype=float)

def sub_orth(R, fixed):
    keep = [r for r in R if all(abs(r @ f) < 1e-9 for f in fixed)]
    return np.array(keep)

def lines_of(R):
    L = []
    seen = set()
    for r in R:
        key = tuple(np.round(r, 6)); nkey = tuple(np.round(-r, 6))
        if key in seen or nkey in seen: continue
        seen.add(key); L.append(r)
    return np.array(L)

def f2_rank(M):
    M = (np.array(M) % 2).astype(np.uint8).copy(); r = 0; rows, cols = M.shape
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i, c]), None)
        if p is None: continue
        M[[r, p]] = M[[p, r]]
        for i in range(rows):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
    return r

def max_clique(adj):
    n = len(adj); best = [0]
    nbrs = [set(np.nonzero(adj[i])[0]) for i in range(n)]
    def expand(R, P):
        if not P:
            best[0] = max(best[0], len(R)); return
        if len(R) + len(P) <= best[0]: return
        for v in sorted(P):
            expand(R | {v}, P & nbrs[v]); P = P - {v}
            if len(R) + len(P) <= best[0]: return
    expand(set(), set(range(n)))
    return best[0]

E8 = roots_E8()
E7 = sub_orth(E8, [E8[-1]])
E6 = sub_orth(E8, [E8[-1], E7[0]])
# make sure E6 is E6 (72 roots) -- pick the second fixed root inside an A2 with the first
if len(E6) != 72:
    for r in E7:
        if abs(r @ E8[-1]) < 1e-9:
            cand = sub_orth(E8, [E8[-1], r])
            # E6 = orthogonal complement of an A2, not of an A1 x A1
    a = E8[-1]
    b = next(r for r in E8 if abs(r @ a + 1) < 1e-9)
    E6 = sub_orth(E8, [a, b])
systems = [("E8", E8), ("E7", E7), ("E6", E6)] + [("D%d" % n, roots_D(n)) for n in (4, 5, 6, 7, 8)] + \
          [("A%d" % n, roots_A(n)) for n in (3, 4, 5, 6, 7, 8)]
r35 = {"E8": 0.600, "E7": 0.778, "E6": 1.167}
print("%-4s %5s %3s %6s %5s %5s %8s %8s %s" % ("sys", "roots", "d", "lines", "alpha", "w", "ratio", "N/(a w)", "note"))
rows = []
for name, R in systems:
    L = lines_of(R); N = len(L)
    L2 = np.rint(2 * L).astype(np.int64)
    G = (L2 @ L2.T) // 4   # exact integer Gram (roots doubled to integers)
    assert ((L2 @ L2.T) % 4 == 0).all()
    d = np.linalg.matrix_rank(G.astype(float))
    U = (np.ones((N, N), dtype=int) + G) % 2
    w = f2_rank(U)
    orth = (G == 0).astype(int); np.fill_diagonal(orth, 0)
    alpha = max_clique(orth)
    ratio = d * w / N
    junk = int((np.abs(G) == 1).sum(axis=1)[0])
    note = ("perfect (alpha = d)" if alpha == d else "alpha < d") + ", junk/port %d" % junk
    if name in r35: note += ", R35 says %.3f" % r35[name]
    print("%-4s %5d %3d %6d %5d %5d %8.3f %8.3f %s" % (name, len(R), d, N, alpha, w, ratio, N / (alpha * w), note))
    rows.append((name, N, d, alpha, w, ratio))
    if name in ("E8", "E7"):
        A = (np.abs(G) == 1).astype(int)   # side-information (cheap) graph: non-orthogonal pairs
        ev = np.round(np.linalg.eigvalsh(A), 6)
        vals, counts = np.unique(ev, return_counts=True)
        k = A.sum(axis=1)[0]; smin = vals.min()
        theta = -N * smin / (k - smin)
        print("     side-information graph: %d-regular, eigenvalues %s ; Lovasz theta = %.4f (alpha = %d)" %
              (k, dict(zip(vals.tolist(), counts.tolist())), theta, alpha))
for name, N, d, alpha, w, ratio in rows:
    if name in r35: assert abs(ratio - r35[name]) < 0.001, (name, ratio)
print("R35 copy ratios reproduced by d*w/N for E8, E7, E6: OK")
