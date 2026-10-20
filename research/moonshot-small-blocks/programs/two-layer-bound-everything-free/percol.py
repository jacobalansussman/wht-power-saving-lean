"""percol.py -- mw2-xwlayers.  Per-column (column-uniform) layer test at block 2.
For column-uniform words g = g1^(x)f the layer is lambda(g) = log2(max row support of the 4x4 per-column matrix).
The free relation P_r = a + b Zs_r, sandwiched as x P_r y, gives the every-f relation
    x P_r y = a (x y) + b (x Zs_r y)
among three group elements; it is layer-INhomogeneous iff lambda(x y) != lambda(x Zs_r y).
We enumerate x = D_w c D_z (c a column-uniform free letter) and y = D_v, or y free, and report every case.
"""
import numpy as np, itertools
a, b = (1 + 1j) / 2, (1 - 1j) / 2
I4 = np.eye(4, dtype=complex)
def X(z):
    M = np.zeros((4, 4), complex)
    for x in range(4): M[x ^ z, x] = 1
    return M
def diag(fn): return np.diag([fn(x & 1, x >> 1) for x in range(4)]).astype(complex)
def dirz(z): return a * I4 + b * X(z)
Z = {0: diag(lambda p, q: (-1) ** p), 1: diag(lambda p, q: (-1) ** q)}
CS = diag(lambda p, q: 1j ** (p * q)); CZ = diag(lambda p, q: (-1) ** (p * q))
S0 = diag(lambda p, q: 1j ** p); S1 = diag(lambda p, q: 1j ** q)
C01 = np.zeros((4, 4), complex)          # x1 ^= x0
for x in range(4):
    p, q = x & 1, x >> 1; C01[p | ((q ^ p) << 1), x] = 1
FREE = {"1": I4, "CS": CS, "CSi": CS.conj(), "CZ": CZ, "S0": S0, "S1": S1, "C01": C01, "CS.C01": C01 @ CS}
def supp(M): return int(max((np.abs(M) > 1e-9).sum(axis=1)))
def lam(M): return np.log2(supp(M))
rows = []
for w, z, v in itertools.product((1, 2, 3), repeat=3):
    for cn, c in FREE.items():
        x = dirz(w) @ c @ dirz(z)          # x = D_w c D_z  (matrix product: rightmost acts first)
        for r in (0, 1):
            for ynm, y in (("D%d" % v, dirz(v)),):
                l1 = lam(x @ y); l2 = lam(x @ Z[r] @ y)
                if l1 != l2:
                    rows.append((w, cn, z, r, ynm, l1, l2))
print("x = D_w c D_z, y = D_v (3 paid letters), relation x P_r y = a x y + b x Zs_r y")
print("cases with lambda(x y) != lambda(x Zs_r y): %d" % len(rows))
seen = set()
for (w, cn, z, r, ynm, l1, l2) in rows:
    print("  D%d . %-6s . D%d . [P%d] . %s   lambda(x y) = %.0f   lambda(x Zs%d y) = %.0f" % (w, cn, z, r, ynm, l1, r, l2))
# controls: Clifford c never gives inhomogeneity
cl = [r for r in rows if r[1] in ("1", "CZ", "S0", "S1", "C01")]
print("Clifford middle letters among them: %d (expected 0: Clifford x maps Zs to a Pauli)" % len(cl))
# with only 2 paid letters: x = D_w c, y = D_z  and x = D_w, y = c D_z
cnt2 = 0
for w, z in itertools.product((1, 2, 3), repeat=2):
    for cn, c in FREE.items():
        for r in (0, 1):
            for x, y in ((dirz(w) @ c, dirz(z)), (dirz(w), c @ dirz(z))):
                if lam(x @ y) != lam(x @ Z[r] @ y): cnt2 += 1
print("2 paid letters (x = D_w c, y = D_z or x = D_w, y = c D_z): inhomogeneous cases = %d (expected 0)" % cnt2)

print()
print("3 paid letters, other placements: x = D_w c D_z, y = c2 D_v; and x = D_w, y = c D_z c2 D_v")
cnt3 = 0; ex3 = []
for w, z, v in itertools.product((1, 2, 3), repeat=3):
    for (cn, c), (c2n, c2) in itertools.product(FREE.items(), repeat=2):
        for r in (0, 1):
            for tag, x, y in (("A", dirz(w) @ c @ dirz(z), c2 @ dirz(v)), ("B", dirz(w), c @ dirz(z) @ c2 @ dirz(v)),
                              ("C", dirz(w) @ c @ dirz(z) @ c2, dirz(v))):
                l1, l2 = lam(x @ y), lam(x @ Z[r] @ y)
                if l1 != l2:
                    cnt3 += 1; ex3.append((tag, w, cn, z, c2n, v, r, l1, l2))
print("inhomogeneous cases with 3 paid letters: %d" % cnt3)
for e in ex3[:10]: print("  ", e)
print()
print("4 paid letters: x = D_w c D_z, y = D_v c2 D_u  (includes h P h^-1 with h = D_w CS D_z)")
cnt4 = 0; ex4 = []
for w, z, v, u in itertools.product((1, 2, 3), repeat=4):
    for (cn, c), (c2n, c2) in itertools.product(FREE.items(), repeat=2):
        for r in (0, 1):
            x = dirz(w) @ c @ dirz(z); y = dirz(v) @ c2 @ dirz(u)
            l1, l2 = lam(x @ y), lam(x @ Z[r] @ y)
            if l1 != l2:
                cnt4 += 1
                if len(ex4) < 12: ex4.append((w, cn, z, r, v, c2n, u, l1, l2))
print("inhomogeneous cases with 4 paid letters: %d" % cnt4)
for (w, cn, z, r, v, c2n, u, l1, l2) in ex4:
    print("  D%d.%s.D%d.[P%d].D%d.%s.D%d   lambda(xy) = %.0f, lambda(x Zs y) = %.0f" % (w, cn, z, r, v, c2n, u, l1, l2))
h = dirz(2) @ CS @ dirz(1)
for r in (0, 1):
    n = h @ Z[r] @ np.linalg.inv(h)
    print("h = D2.CS.D1: lambda(h Zs%d h^-1) = %.0f   (h P h^-1 = a + b h Zs h^-1: layers {%.0f, 0})" % (r, lam(n), lam(n)))
