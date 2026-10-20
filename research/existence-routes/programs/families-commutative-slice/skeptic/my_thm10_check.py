# Theorem 10 and the rank-cover inequalities on GENUINE sharing units built from my own U_3 solution:
#   block 3: U_3 (2 in 7);  block 4: U_3 + e4 on each array (2 in 9);  block 5: U_3 + e4 + e5 (2 in 11).
import numpy as np, itertools, sys
import geo
d = np.load("u3_solution.npz")
A1, B1, G1, dirs1, W = d["A"], d["B"], d["Gam"], list(d["dirs"]), int(d["W"])
def compose(A1, B1, G1, dirs1, newdir):
    M1 = len(dirs1)
    A = np.hstack([A1, np.eye(W)])
    B = np.vstack([B1, np.eye(W)])
    G = np.zeros((M1 + W, M1 + W), complex)
    G[:M1, :M1] = G1
    G[M1:, :M1] = 2 * A1
    return A, B, G, dirs1 + [newdir] * W
units = {3: (A1, B1, G1, dirs1)}
units[4] = compose(*units[3], 8)
units[5] = compose(*units[4], 16)
def svals(Mx):
    return np.linalg.svd(Mx, compute_uv=False) if Mx.size else np.array([])
tol = 1e-8
for m in (3, 4, 5):
    A, B, G, dirs = units[m]
    M = len(dirs)
    K = B @ A
    E = np.eye(M) + G - K
    forms = geo.all_form_supports(m)
    worst_star = 0; worst_even = np.inf; worst_odd_extra = 0; rank_bad = 0; ntight = 0; zero_bad = 0
    tight = []
    for S, tr in forms:
        D = [i for i, z in enumerate(dirs) if z in S]
        if D:
            X = A[:, D] @ np.linalg.solve(np.eye(len(D)) + G[np.ix_(D, D)], B[D, :])
        else:
            X = np.zeros((W, W))
        worst_star = max(worst_star, np.abs(X - tr * np.eye(W)).max())
        if not D:
            continue
        sv = svals(E[np.ix_(D, D)])
        if tr == 0:
            worst_even = min(worst_even, sv.min())
            if sv.min() < tol: rank_bad += 1
        else:
            slack = len(D) - W
            r = int((sv > tol).sum())
            if r != slack: rank_bad += 1
            worst_odd_extra = max(worst_odd_extra, sv[slack] if slack < len(sv) else 0)
            if slack == 0:
                ntight += 1; tight.append(D)
    for D in tight:
        zero_bad = max(zero_bad, np.abs(E[np.ix_(D, D)]).max())
    print("block %d: %d in %d, forms %d: max |X_Q - [odd]I| = %.2e; Thm10(b) min sing. value on even supports %.3f;"
          " Thm10(c) rank mismatches %d, largest singular value beyond slack %.2e; tight supports %d, max |E| on them %.2e"
          % (m, W, M, len(forms), worst_star, worst_even, rank_bad, worst_odd_extra, ntight, zero_bad))
    # rank-cover check with ACTUAL ranks: for every even (m-2)-flat, rank E_F = |D_F| <= sum of piece ranks,
    # each piece rank <= min(size, slack of every odd (m-2)-flat around it)
    if m >= 4:
        ev, od = geo.even_odd_flats(m, m - 2)
        viol = 0; nchk = 0
        def rk(R, C):
            Rm = [i for i, z in enumerate(dirs) if z in R]; Cm = [i for i, z in enumerate(dirs) if z in C]
            if not Rm or not Cm: return 0
            return int((svals(E[np.ix_(Rm, Cm)]) > tol).sum())
        for F in ev:
            DF = [i for i, z in enumerate(dirs) if z in F]
            tot = sum(rk({x}, {y}) for x in F for y in F)
            if len(DF) > tot: viol += 1
            nchk += 1
        for x in range(1, 1 << m):
            for y in range(1, 1 << m):
                r = rk({x}, {y})
                for G_ in od:
                    if x in G_ and y in G_:
                        nchk += 1
                        slack = sum(1 for z in dirs if z in G_) - W
                        if r > slack: viol += 1
        print("   rank-cover pieces/rows checked %d, violations %d; c-vector sum = %.3f (LP floor %s)"
              % (nchk, viol, M / W, "4" if m == 4 else "155/36=4.306"))
