# b4-depth4: relaxed integer program for family S when the junk of the shared lineage e* may be cancelled by
# OTHER mixed lineages (the residual case).  Every constraint is a necessary condition (derived in the report,
# section 2); tree sizes are lower-bounded generously (relaxation), so the minimum is a valid LOWER bound on n.
#   a, a2   : depth of e*'s x1 / x2 parents (nodes x in T1 at depth dx >= a, y in T2 at depth dy >= a2)
#   c, r1, r2 : common prefix and the two branches of e*'s lineage
#   K-paths : a + c + r1 >= m ; a2 + c + r2 >= m
#   junk depth dJ >= max(0, m - (r1+r2), m - (a+a2))   (for both J1 and J2)
#   J1 cover, pure : node z of T2: depth dz >= dJ, dist(z,y) + a >= m   (exact Steiner sizes)
#   J1 cover, lineage f : q >= 1 moves in M_f, node y_f of T2 at depth dyf:  q + dyf >= dJ,  q + dist(y_f,y) + a >= m
#                         T2 size >= Steiner(root, y, y_f) with that distance (the pure canceller, if any, is
#                         not added: relaxation)
#   J2 symmetric (cover in T1, or lineage g, possibly g = f: then cost max(q, q'))
import sys, itertools
m = int(sys.argv[1]) if len(sys.argv) > 1 else 4
H = 2 * m
def steiner(d1, d2, mind):   # least size of a rooted tree holding nodes at depths d1, d2 with dist >= mind
    best = None
    for L in range(0, min(d1, d2) + 1):
        if d1 + d2 - 2 * L >= mind:
            s = d1 + d2 - L
            best = s if best is None or s < best else best
    return best
best = None; arg = None
for a, a2, c, r1, r2 in itertools.product(range(m + 1), range(m + 1), range(1, m + 1), range(m + 1), range(m + 1)):
    if a + c + r1 < m or a2 + c + r2 < m: continue
    dJ = max(0, m - (r1 + r2), m - (a + a2))
    base = c + r1 + r2
    # options for J1 (side T2, uses a for the K relation) and J2 (side T1, uses a2)
    def options(d_par, a_other):
        """yield (tree_size, lineage_q) for covering the junk on the tree whose parent node has depth d_par"""
        out = []
        for dpar in range(d_par, H + 1):
            # pure cover
            for dz in range(dJ, H + 1):
                for dist in range(0, dpar + dz + 1):
                    if dist + a_other < m: continue
                    s = steiner(dpar, dz, dist)
                    if s is not None: out.append((s, 0))
            # lineage cover
            for q in range(1, H + 1):
                for dyf in range(0, H + 1):
                    if q + dyf < dJ: continue
                    for dist in range(0, dpar + dyf + 1):
                        if q + dist + a_other < m: continue
                        st = steiner(dpar, dyf, dist)
                        if st is not None: out.append((st, q))
        # keep Pareto-minimal
        out = sorted(set(out))
        pareto = []
        for s, q in out:
            if not any(s2 <= s and q2 <= q for s2, q2 in pareto): pareto.append((s, q))
        return pareto
    o1 = options(a2, a)     # J1: tree T2 (parent y at depth >= a2), K-relation uses a
    o2 = options(a, a2)     # J2: tree T1 (parent x at depth >= a), K-relation uses a2
    for (s2, q1), (s1, q2) in itertools.product(o1, o2):
        lin = max(q1, q2) if (q1 and q2) else q1 + q2      # same lineage may cover both (relaxed)
        n = s1 + s2 + base + lin
        if best is None or n < best: best = n; arg = (a, a2, c, r1, r2, dJ, (s1, q2), (s2, q1))
print(f"m={m}: relaxed lower bound for family S with both kernel paths through one lineage = {best} (plain {2*m})")
print("  argmin a,a',c,r1,r2,dJ,(T1 size, lineage q for J2),(T2 size, lineage q for J1):", arg)
