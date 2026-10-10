#!/usr/bin/python3
"""subspace_metrics.py -- 10 October 2026.  Standard library only, exact arithmetic modulo a prime.

Checks of the subspace labels of statement A in ../for-network-coding-readers.md (written by an AI agent; not reviewed by a person):

 B1  K_{2,3}: the five subspaces {0, three lines, k^2} of k^2 have subspace distance
        d_S(U,V) = dim(U+V) - dim(U cap V)
     equal to the shortest-path distance of K_{2,3}.
 B2  Gamma_{3,3} (the 16-vertex graph of arXiv:2608.06070, section VI.B: K_{3,3} with every edge a_i b_j
     subdivided by x_ij, plus a vertex o adjacent to all nine x_ij): the 16 subspaces
        a_i -> A_i, b_j -> B_j   (three + three lines of PG(3,k) forming a 3 x 3 grid: A_i cap B_j a point,
                                  the A_i pairwise skew, the B_j pairwise skew)
        x_ij -> A_i + B_j,   o -> k^4
     have d_S equal to the shortest-path distance of Gamma_{3,3}.  Also: on the six a, b alone d_S = 2 * d_{K_{3,3}}.
 B3  the identity  delta_{Y^perp}(U,V) = delta_Y(V,U), where delta_Y(U,V) = dim Y_U - dim(Y_U cap Y_V); so the
     directed inequality for the family Y and for the family of annihilators add up to the symmetric d_S.
 B4  G_sub(3), the 16 "subspace frames" of block size 3 = Hasse graph of the subspaces of GF(2)^3 = the Heawood graph
     (7 points, 7 lines of the Fano plane) plus two apex vertices; number of ordered pairs at distance 3.

Run:  python3 -I -B subspace_metrics.py
"""
import itertools, random, sys


def rank_mod(rows, p):
    rows = [list(r) for r in rows if any(x % p for x in r)]
    rk = 0
    ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = None
        for i in range(rk, len(rows)):
            if rows[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], p - 2, p)
        rows[rk] = [(x * inv) % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % p:
                f = rows[i][c]
                rows[i] = [(x - f * y) % p for x, y in zip(rows[i], rows[rk])]
        rk += 1
        if rk == len(rows):
            break
    return rk


def dS(U, V, p):
    """U, V lists of spanning vectors (possibly empty).  dim(U+V) - dim(U cap V) = 2 dim(U+V) - dim U - dim V."""
    n = len((U + V)[0]) if (U + V) else 1
    ru = rank_mod(U, p) if U else 0
    rv = rank_mod(V, p) if V else 0
    rs = rank_mod(U + V, p) if (U + V) else 0
    return 2 * rs - ru - rv


def delta(U, V, p):
    """dim U - dim(U cap V) = dim(U+V) - dim V."""
    rv = rank_mod(V, p) if V else 0
    rs = rank_mod(U + V, p) if (U + V) else 0
    return rs - rv


def bfs_all(n, edges):
    adj = {i: set() for i in range(n)}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    D = [[None] * n for _ in range(n)]
    for s in range(n):
        D[s][s] = 0
        frontier = [s]
        while frontier:
            nxt = []
            for u in frontier:
                for w in adj[u]:
                    if D[s][w] is None:
                        D[s][w] = D[s][u] + 1
                        nxt.append(w)
            frontier = nxt
    return D


def compare(name, spaces, edges, p):
    n = len(spaces)
    D = bfs_all(n, edges)
    bad = [(i, j) for i in range(n) for j in range(n) if dS(spaces[i], spaces[j], p) != D[i][j]]
    hist = {}
    for i in range(n):
        for j in range(i + 1, n):
            hist[D[i][j]] = hist.get(D[i][j], 0) + 1
    print("  %-34s over GF(%d): %d vertices, %d edges, distance histogram %s ; pairs where d_S != graph distance: %d"
          % (name, p, n, len(edges), dict(sorted(hist.items())), len(bad)))
    return not bad


ok = True
print("B1  K_{2,3} as the subspaces of k^2")
for p in (2, 3, 5, 7):
    spaces = [[], [[1, 0], [0, 1]], [[1, 0]], [[0, 1]], [[1, 1]]]      # a = 0, b = k^2, p1, p2, p3
    edges = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4)]
    ok &= compare("K_{2,3}", spaces, edges, p)

print("B2  Gamma_{3,3} through a 3 x 3 grid of lines in projective 3-space")
for p in (2, 3, 5, 7, 1000003):
    # k^4 = k^2 (x) k^2 with basis e_s (x) e_t -> coordinate 2s + t.  u_1 = e0, u_2 = e1, u_3 = e0 + e1.
    u = [(1, 0), (0, 1), (1, 1)]

    def ten(x, y):
        return [x[0] * y[0], x[0] * y[1], x[1] * y[0], x[1] * y[1]]
    A = [[ten(u[i], (1, 0)), ten(u[i], (0, 1))] for i in range(3)]     # A_i = u_i (x) k^2
    B = [[ten((1, 0), u[j]), ten((0, 1), u[j])] for j in range(3)]     # B_j = k^2 (x) u_j ; A_i cap B_j = u_i (x) u_j
    X = [[A[i] + B[j] for j in range(3)] for i in range(3)]
    O = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    # vertex order: a1 a2 a3 b1 b2 b3 x11 .. x33 o
    spaces = A + B + [X[i][j] for i in range(3) for j in range(3)] + [O]
    edges = []
    for i in range(3):
        for j in range(3):
            x = 6 + 3 * i + j
            edges += [(i, x), (3 + j, x), (x, 15)]
    ok &= compare("Gamma_{3,3}", spaces, edges, p)
    # K_{3,3} on the six a, b: d_S = 2 d
    e33 = [(i, 3 + j) for i in range(3) for j in range(3)]
    D = bfs_all(6, e33)
    bad = [(i, j) for i in range(6) for j in range(6) if dS(spaces[i], spaces[j], p) != 2 * D[i][j]]
    print("  %-34s over GF(%d): pairs where d_S != 2 * graph distance: %d" % ("K_{3,3} on the a's and b's", p, len(bad)))
    ok &= not bad

print("B3  delta for the annihilators is delta with the arguments exchanged (random subspaces)")
rnd = random.Random(20261010)


def annihilator(U, n, p):
    """a basis of {y : y.u = 0 for all u in U}, by brute force over GF(p)^n (n, p small)."""
    out = []
    for y in itertools.product(range(p), repeat=n):
        if all(sum(a * b for a, b in zip(y, uvec)) % p == 0 for uvec in U):
            out.append(list(y))
    return out


for p, n in ((2, 5), (3, 4), (5, 3)):
    bad = 0
    trials = 300
    for _ in range(trials):
        U = [[rnd.randrange(p) for _ in range(n)] for _ in range(rnd.randrange(0, n + 1))]
        V = [[rnd.randrange(p) for _ in range(n)] for _ in range(rnd.randrange(0, n + 1))]
        Up, Vp = annihilator(U, n, p), annihilator(V, n, p)
        if delta(Up, Vp, p) != delta(V, U, p):
            bad += 1
        if delta(U, V, p) + delta(V, U, p) != dS(U, V, p):
            bad += 1
    print("  GF(%d)^%d: %d random pairs, violations: %d" % (p, n, trials, bad))
    ok &= bad == 0

print("B4  the 16 subspace frames at block size 3")
vecs = [v for v in itertools.product(range(2), repeat=3) if any(v)]
points = [[list(v)] for v in vecs]
planes = []
for w in vecs:                       # the plane w^perp
    planes.append([list(v) for v in vecs if sum(a * b for a, b in zip(v, w)) % 2 == 0][:2])
spaces = [[]] + points + planes + [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]]
n = len(spaces)
edges = [(i, j) for i in range(n) for j in range(i + 1, n) if dS(spaces[i], spaces[j], 2) == 1]
D = bfs_all(n, edges)
bad = [(i, j) for i in range(n) for j in range(n) if D[i][j] != dS(spaces[i], spaces[j], 2)]
deg = sorted(len([e for e in edges if i in e]) for i in range(n))
flags = [(i, j) for (i, j) in edges if 1 <= i <= 7 and 8 <= j <= 14]
opp = [(i, j) for i in range(n) for j in range(n) if D[i][j] == 3]
print("  vertices %d, edges %d (of which point-plane incidences, i.e. edges of the Heawood graph: %d), degrees %s"
      % (n, len(edges), len(flags), deg))
print("  graph distance = d_S on all pairs: %s ; ordered pairs at distance 3: %d ; of them zero<->full: %d, point<->plane not through it: %d"
      % (not bad, len(opp), len([1 for i, j in opp if {i, j} == {0, 15}]), len([1 for i, j in opp if {i, j} != {0, 15}])))
ok &= (not bad) and len(opp) == 58 and len(flags) == 21

print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
sys.exit(0 if ok else 1)
