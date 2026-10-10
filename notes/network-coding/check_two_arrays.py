#!/usr/bin/env python3
"""check_two_arrays.py: exact, stand-alone check of the table behind
"at block size 3, two arrays cannot be served with 5 moves" (two-arrays.md, next to this file).
usage:   python3 -B check_two_arrays.py TABLE.json        (standard library only; integers and fractions only)

TABLE.json:  "frames": a list of frames, each a list of labels.  A label is an integer: bits 0..m-1 are x, bits
                       m..2m-1 are z; it stands for the Pauli matrix X^x Z^z.   "m": block size (default 3).
             "rows":   for each frame, in the same order, a few integer rows; they span its subspace Y_U of Q^n.
             "abar":   the declared averages (optional, strings).
Checked, from nothing (distance of two frames = m - dimension of their intersection):
  1. the frames of the file are exactly the Lagrangian subspaces of F_2^(2m), each once
     (sets of 2^m labels, closed under addition, on which the form  x.z' + x'.z  is zero; 135 of them for m = 3);
  2. the maps of labels that keep the form have exactly ONE orbit on ordered pairs of frames at each distance
     (shown with the maps  t_a: v -> v + form(v, a) * a  for a = (z|0) and a = (0|z)); and, for the concrete machine
     whose paid letter acts on labels as t_(z|0): t_(0|z) keeps the frame L0 = {(0|z)}, t_(z|0) moves L0 to
     distance 1, and t_(e_1|0) ... t_(e_m|0) moves L0 to distance m;
  3. delta(U, V) = dim Y_U - dim(Y_U & Y_V) for every ordered pair of frames, by exact integer elimination;
     abar_d = the average of delta over the ordered pairs at distance d;  abar_d <= d * abar_1;
  4. the inequality  (2m - 1) * abar_1 < 2 * abar_m :  2m - 1 moves supply less score than two arrays need.
Exit code 0 only if 1 to 4 hold and the declared averages, if present, are reproduced."""
import hashlib
import json
import sys
from fractions import Fraction
from itertools import combinations


def rank(rows):
    """exact rank over the rationals of an integer matrix (fraction-free elimination; every division is exact)"""
    A = [[int(x) for x in r] for r in rows]
    rk, prev = 0, 1
    for c in range(len(A[0])):
        p = next((r for r in range(rk, len(A)) if A[r][c]), None)
        if p is None:
            continue
        A[rk], A[p] = A[p], A[rk]
        for r in range(rk + 1, len(A)):
            new = [divmod(A[r][k] * A[rk][c] - A[r][c] * A[rk][k], prev) for k in range(len(A[r]))]
            assert all(rem == 0 for _, rem in new)
            A[r] = [q for q, _ in new]
        prev = A[rk][c]
        rk += 1
    return rk


def span(gens):
    S = {0}
    for g in gens:
        S |= {s ^ g for s in S}
    return frozenset(S)


def fail(why):
    print("FAIL:", why)
    sys.exit(1)


def main(path):
    raw = open(path, "rb").read()
    d = json.loads(raw.decode())
    m = int(d.get("m", 3))
    mask = (1 << m) - 1
    par = lambda a: bin(a).count("1") & 1
    form = lambda v, w: par(v & mask & (w >> m)) ^ par(w & mask & (v >> m))
    print("table: %s   sha256 %s   block size m = %d" % (path, hashlib.sha256(raw).hexdigest()[:16], m))

    # 1. the frames, from nothing, against the frames of the file
    lag = set()
    for gens in combinations(range(1, 1 << (2 * m)), m):
        if all(form(a, b) == 0 for a, b in combinations(gens, 2)) and len(span(gens)) == 1 << m:
            lag.add(span(gens))
    frames = [frozenset(f) for f in d["frames"]]
    n = len(frames)
    if len(set(frames)) != n or set(frames) != lag:
        fail("the frames of the file are not exactly the %d Lagrangian subspaces, each once" % len(lag))
    rows = d["rows"]
    if len(rows) != n or any(not r or len({len(x) for x in r}) != 1 or not r[0] for r in rows):
        fail("the file must give every frame a non-empty list of rows of one length")
    if any(type(x) is not int for r in rows for row in r for x in row):
        fail("a row entry is not an integer")
    idx = {f: i for i, f in enumerate(frames)}
    lg = {1 << k: k for k in range(m + 1)}
    dist = [[m - lg[len(frames[i] & frames[j])] for j in range(n)] for i in range(n)]
    npairs = [sum(row.count(k) for row in dist) for k in range(m + 1)]
    print("1. frames: %d, exactly the Lagrangian subspaces; ordered pairs at distance 0..%d: %s" % (n, m, npairs))

    # 2. orbits on ordered pairs, and the three facts about the concrete machine
    t = lambda a: (lambda v: v ^ a if form(v, a) else v)
    act = lambda fn: [idx[frozenset(fn(v) for v in f)] for f in frames]        # a label map as a permutation of frames
    L0 = idx[frozenset(z << m for z in range(1 << m))]
    free = [act(t(z << m)) for z in range(1, 1 << m)]
    dirs = [act(t(z)) for z in range(1, 1 << m)]
    kern = L0
    for k in range(m):
        kern = dirs[(1 << k) - 1][kern]
    if not (all(p[L0] == L0 for p in free) and all(dist[L0][p[L0]] == 1 for p in dirs) and dist[L0][kern] == m):
        fail("a fact about the maps t_(0|z), t_(z|0) or the kernel map and the frame L0 does not hold")
    seen = [[-1] * n for _ in range(n)]
    reps = []
    for a in range(n):
        for b in range(n):
            if seen[a][b] < 0:
                seen[a][b] = len(reps)
                reps.append(dist[a][b])
                stack = [(a, b)]
                while stack:
                    x, y = stack.pop()
                    for p in free + dirs:
                        if seen[p[x]][p[y]] < 0:
                            seen[p[x]][p[y]] = seen[a][b]
                            stack.append((p[x], p[y]))
    if sorted(reps) != list(range(m + 1)):
        fail("the label maps do not have exactly one orbit of ordered pairs per distance: " + str(sorted(reps)))
    print("2. one orbit of ordered pairs of frames per distance (under %d maps t_a); the three facts about L0 hold" % len(free + dirs))

    # 3. delta and its averages
    dim = [rank(r) for r in rows]
    tot = [0] * (m + 1)
    hist = [dict() for _ in range(m + 1)]
    for i in range(n):
        for j in range(i + 1, n):
            r = rank(rows[i] + rows[j])
            for dl in (r - dim[j], r - dim[i]):                 # delta(U_i, U_j) and delta(U_j, U_i)
                tot[dist[i][j]] += dl
                hist[dist[i][j]][dl] = hist[dist[i][j]].get(dl, 0) + 1
    abar = [Fraction(tot[k], npairs[k]) for k in range(m + 1)]
    print("3. dimensions of the Y_U: %s ; row length: %s ; largest entry in absolute value: %d"
          % (sorted(set(dim)), sorted({len(r[0]) for r in rows}), max(abs(x) for r in rows for row in r for x in row)))
    for k in range(1, m + 1):
        print("   distance %d: delta -> number of ordered pairs %s ; average abar_%d = %s = %.6f"
              % (k, dict(sorted(hist[k].items())), k, abar[k], float(abar[k])))
    if "abar" in d and [str(x) for x in d["abar"]] != [str(a) for a in abar[1:]]:
        fail("the declared averages %s are not reproduced" % d["abar"])
    a1, am = abar[1], abar[m]
    if a1 <= 0 or any(abar[k] > k * a1 for k in range(2, m + 1)):
        fail("abar_d <= d * abar_1 fails, which no table can do (triangle inequality): an error somewhere")

    # 4. the inequality
    print("4. THE TWO NUMBERS:  score per move  abar_1 = %s   ;   score per array served  abar_%d = %s" % (a1, m, am))
    print("   every correct run for W arrays uses at least  W * %s = %.6f * W  moves" % (am / a1, float(am / a1)))
    print("   two arrays: %d moves supply at most %d * %s = %s = %.4f ; two arrays need 2 * %s = %s"
          % (2 * m - 1, 2 * m - 1, a1, (2 * m - 1) * a1, float((2 * m - 1) * a1), am, 2 * am))
    if not (2 * m - 1) * a1 < 2 * am:
        fail("(2m - 1) * abar_1 < 2 * abar_m does NOT hold: this table does not exclude %d moves for two arrays" % (2 * m - 1))
    print("PASS: two arrays cannot be served with %d moves at block size %d; they need %d" % (2 * m - 1, m, 2 * m))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
