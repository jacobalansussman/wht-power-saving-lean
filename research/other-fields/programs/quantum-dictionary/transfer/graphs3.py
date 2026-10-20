"""x-transfer, 13 Oct 2026 (python3 -I).  The 8 graph states on 3 qubits (frames (x | Gamma x), Gamma = adjacency):
our classification (problem symmetry GL(3,2) = stabiliser of L_0 and L_F, acting by congruence on Gamma + I,
invariant rank(Gamma + I) = d(frame, L_F)) versus local-complementation (LC) classes."""
import itertools
def rank_f2(rows):
    rows = list(rows); r = 0
    while rows:
        p = rows.pop()
        if p == 0: continue
        r += 1; lb = p & -p
        rows = [q ^ p if q & lb else q for q in rows]
    return r
m = 3
edges = [(0, 1), (0, 2), (1, 2)]
def lc(G, v):   # local complementation at v: complement the edges among neighbours of v
    nb = [u for u in range(m) if frozenset((u, v)) in G]
    G = set(G)
    for a, b in itertools.combinations(nb, 2):
        e = frozenset((a, b)); G ^= {e}
    return frozenset(G)
graphs = [frozenset(frozenset(e) for e, b in zip(edges, bits) if b) for bits in itertools.product((0, 1), repeat=3)]
def rk(G):
    rows = []
    for i in range(m):
        r = 1 << i
        for j in range(m):
            if frozenset((i, j)) in G: r |= 1 << j
        rows.append(r)
    return rank_f2(rows)
name = lambda G: {0: "empty", 1: "one edge", 2: "path", 3: "triangle"}[len(G)]
# simpler: group by LC orbit membership
seen = []; classes = []
for G in graphs:
    if any(G in c for c in classes): continue
    S = {G}; st = [G]
    while st:
        H = st.pop()
        for v in range(m):
            K = lc(H, v)
            if K not in S: S.add(K); st.append(K)
        for p in itertools.permutations(range(m)):
            K = frozenset(frozenset(p[x] for x in e) for e in H)
            if K not in S: S.add(K); st.append(K)
    classes.append(S)
print("LC classes (with relabelling):", [sorted({name(G) for G in c}) for c in classes])
by_rank = {}
for G in graphs: by_rank.setdefault(rk(G), set()).add(name(G))
print("our classes, rank(Gamma + I) = d(frame, L_F):", {r: sorted(v) for r, v in sorted(by_rank.items())})
