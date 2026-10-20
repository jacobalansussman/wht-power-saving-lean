"""tightpairs.py: the tight-pair test (PROOF.md Theorem 1) for the uniform-direction family U_m.

For a commutative-slice program, K = B^T A (M x M, rank <= W). For every odd-trace form Q whose support holds
exactly W moves, K[D_Q, D_Q] = I + H_{D_Q} is unit upper triangular in program order. So every pair of moves that
lies in a common such support has K = 0 in its lower position, and a set Y of moves all of whose pairs are covered
gives a unit-triangular principal block K[Y, Y]: rank K >= |Y|, hence |Y| <= W.

For U_m^(r) (r copies of each nonzero point, W = r 2^(m-2)) the tight odd forms are the indicators of affine
(m-2)-flats not through 0 with an odd number of basis vectors. This program checks, for each m, that every pair of
nonzero points (and every single point) lies in a common odd flat; then Y = all r(2^m - 1) moves and the program
cannot exist (r(2^m - 1) > r 2^(m-2)).

Flats: P(L, a) = {x : l(x) = a(l) for l in L}, L a 2-dim subspace of linear forms, a a nonzero functional on L.
"""
import itertools, sys


def dot(a, b):
    return bin(a & b).count("1") & 1


def run(m):
    n = 1 << m
    pts = list(range(1, n))
    basis = [1 << j for j in range(m)]
    # all 2-dim subspaces of the dual (forms are also bitmasks): {0, u, v, u^v}
    subs = set()
    for u in range(1, n):
        for v in range(u + 1, n):
            subs.add(tuple(sorted((u, v, u ^ v))))
    flats = []
    for (u, v, w) in subs:
        for (au, av) in [(1, 0), (0, 1), (1, 1)]:
            P = [x for x in pts if dot(u, x) == au and dot(v, x) == av]
            par = sum(1 for e in basis if e in P) & 1
            flats.append((P, par))
    assert all(len(P) == 1 << (m - 2) for P, _ in flats)
    nodd = sum(p for _, p in flats)
    cov = [[False] * n for _ in range(n)]
    for P, par in flats:
        if par:
            for a in P:
                for b in P:
                    cov[a][b] = True
    uncovered = [(a, b) for a in pts for b in pts if a <= b and not cov[a][b]]
    print(f"m={m}: min-weight flats {len(flats)} (odd {nodd}); uncovered pairs (incl. a=b) {len(uncovered)}"
          + (f" e.g. {uncovered[:5]}" if uncovered else "") +
          f"; => U_m^(r) {'EXCLUDED for every r' if not uncovered else 'not decided by this test'}"
          f" (clique {n-1 if not uncovered else '?'} > W = r*{1 << (m-2)} per copy)", flush=True)
    return uncovered


if __name__ == "__main__":
    for m in map(int, sys.argv[1:]):
        run(m)
