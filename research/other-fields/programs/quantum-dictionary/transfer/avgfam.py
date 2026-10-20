"""x-transfer, 13 Oct 2026 (python3 -I).  Two quantum-natural families for the Averaging Lemma (w-bounds.md sec. 3),
m = 3: delta(U,V) = dim Y_U - dim(Y_U cap Y_V); ratio a_3/a_1 (best known lift: 105/41 = 2.561; 8/3 needed).
 (i)  Y_U = span{Pauli P : P in U} (the stabiliser group's span; = span of the frame's basis projectors):
      Paulis are linearly independent, so dim(Y_U cap Y_V) = |U cap V| exactly.
 (ii) Y_U = span{e_V : d(U,V) <= 1} in Q^135 (closed neighbourhood indicators)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import frames_geom as fg
m = 3
Ls = list(fg.all_lagrangians(m)); n = len(Ls)
D = [[fg.dist(m, A, B) for B in Ls] for A in Ls]
for name, delta in [("(i) stabiliser-group span", lambda a, b: len(Ls[a]) - len(Ls[a] & Ls[b])),
                    ("(ii) closed neighbourhood", None)]:
    if delta is None:
        N = [set(j for j in range(n) if D[i][j] <= 1) for i in range(n)]
        delta = lambda a, b: len(N[a]) - len(N[a] & N[b])
    avg = {}
    for a in range(n):
        for b in range(n):
            d = D[a][b]
            if d:
                s, c = avg.get(d, (0, 0)); avg[d] = (s + delta(a, b), c + 1)
    a = {d: s / c for d, (s, c) in avg.items()}
    print(f"{name}: a_1 = {a[1]:.4f}, a_2 = {a[2]:.4f}, a_3 = {a[3]:.4f}; a_3/a_1 = {a[3]/a[1]:.4f}; 3 arrays need >= {3*a[3]/a[1]:.3f} moves")
