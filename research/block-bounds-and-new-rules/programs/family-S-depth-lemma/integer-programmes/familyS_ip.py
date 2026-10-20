# b4-depth4: integer program for the depth constraints of a 2-array scheme whose two K-paths share one mixing
# event e (single-term generators; junk cancelled by nodes of the pure trees T1, T2).  All quantities are move
# counts / depths.  Constraints (derived in the report, section 3):
#  a = depth(h_e) <= dx (node x of T1), a' = depth(h'_e) <= dy (node y of T2); M_e: common prefix c >= 1, branches r1, r2
#  K-paths:      r1 + c + a >= m,   r2 + c + a' >= m
#  junk J1 (x2-part of final 1) cancelled at T2 node z:  dist_T2(z, y) + a >= m,  depth(z) >= dJ1,
#      dJ1 >= m - (r1 + r2),  dJ1 >= m - (a + a')
#  junk J2 symmetric:  dist_T1(z', x) + a' >= m,  depth(z') >= dJ2,  dJ2 >= m - (r1 + r2), dJ2 >= m - (a + a')
#  sizes: |T1| >= dx + dz' - L (Steiner tree of root, x, z'; L = depth of lca), dist = dx + dz' - 2L; same for T2.
import sys, itertools
m = int(sys.argv[1]) if len(sys.argv) > 1 else 4
best = None; sols = []
R = range(0, 2 * m + 1)
for a, a2, c, r1, r2 in itertools.product(range(m + 1), range(m + 1), range(1, m + 1), range(m + 1), range(m + 1)):
    if r1 + c + a < m or r2 + c + a2 < m: continue
    dJ = max(0, m - (r1 + r2), m - (a + a2))
    # T1: x at depth dx >= a, z' at depth dz >= dJ, lca depth L; dist >= m - a2
    def tree(aa, other):
        bestT = None
        for dx in range(aa, 2 * m + 1):
            for dz in range(dJ, 2 * m + 1):
                for L in range(0, min(dx, dz) + 1):
                    if dx + dz - 2 * L < m - other: continue
                    s = dx + dz - L
                    if bestT is None or s < bestT[0]: bestT = (s, dx, dz, L)
        return bestT
    t1 = tree(a, a2); t2 = tree(a2, a)
    n = t1[0] + t2[0] + c + r1 + r2
    if best is None or n < best: best = n; sols = []
    if n == best: sols.append((a, a2, c, r1, r2, dJ, t1, t2))
print(f"m={m}: least n allowed by these constraints = {best} (plain {2*m}); {len(sols)} parameter sets, e.g.")
for s in sols[:12]: print("  a,a',c,r1,r2,dJ =", s[:6], " T1 (size,dx,dz',L) =", s[6], " T2 =", s[7])
