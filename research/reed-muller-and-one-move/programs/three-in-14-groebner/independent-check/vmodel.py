"""vmodel.py (v-r3-3in14): two checks of MY equation generator (vpoly.equations), exact rational arithmetic.
(1) Against the program itself: registers = W arrays + 3 scratch; a move s at form Q is the reflection
    I - 2 u_s w_s^T when Q acts on it.  The product over D in time order, restricted to the arrays' block, must be
    I - 2 X_D with a_s = (u_s on the arrays), b_s = (w_s on the arrays), G[s][t] = 2 w_s . u_t (s later than t).
    I evaluate my polynomials (+ diag p) at that point and compare with the direct product, for random u, w, orders.
(2) The plain program (a_s = b_s = e_l(s), G[s][t] = 2 [same array]) makes every one of my equations vanish, for
    each of the 44 plain 3-in-15 multisets (deleted move restored), with U = the certificate's 5 moves PLUS the
    restored move, in random orders.  Output: out_vmodel.txt."""
import json, random, itertools
from fractions import Fraction
from vpoly import VRing, equations
from vcore import Inst, dir_vec, tgt_vec, EXP
import numpy as np

random.seed(16102026)
out = []

def evalp(p, pt):
    s = Fraction(0)
    for m, c in p.items():
        t = c
        for v in m: t *= pt[v]
        s += t
    return s

# (1) generic check of the convention, W = 3, |D| up to 5, scratch 3
W = 3; n = W + 3; bad = 0; tot = 0
for trial in range(40):
    k = random.randint(1, 5); U = list(range(k)); order = U[:]; random.shuffle(order)
    p = tuple(random.randint(0, 1) for _ in range(W))
    R = VRing(W, U)
    polys_raw = equations(R, order, [(U, (0, 0, 0))])         # p = 0: polys are exactly the entries of X_D
    u = {s: [Fraction(random.randint(-3, 3), random.randint(1, 3)) for _ in range(n)] for s in U}
    w = {s: [Fraction(random.randint(-3, 3), random.randint(1, 3)) for _ in range(n)] for s in U}
    pt = {}
    for l in range(W):
        for s in U: pt[R.a[l, s]] = u[s][l]; pt[R.b[l, s]] = w[s][l]
    for s in U:
        for t in U:
            if s != t: pt[R.g[s, t]] = 2 * sum(w[s][i] * u[t][i] for i in range(n))
    # direct product, later moves on the left
    P = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for s in order:
        Rf = [[Fraction(int(i == j)) - 2 * u[s][i] * w[s][j] for j in range(n)] for i in range(n)]
        P = [[sum(Rf[i][t] * P[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
    # my polys, in the order l, l2 (entries that are identically zero are dropped by equations(): rebuild all 9)
    from vpoly import matmul  # noqa
    vals = []
    # recompute all 9 entries (equations() drops zero polys; with p = 0 none is zero here since |D| >= 1)
    assert len(polys_raw) == W * W
    for idx, f in enumerate(polys_raw):
        l, l2 = divmod(idx, W)
        X = evalp(f, pt)
        direct = (Fraction(int(l == l2)) - P[l][l2]) / 2
        tot += 1; bad += (X != direct)
out.append(f"(1) my X_D vs the direct product of reflections with 3 scratch registers: {tot} entries, {bad} mismatches")

# (2) plain programs
sets = json.load(open("vsets.json"))
nz = 0; npolys = 0
for s in sets:
    ks = s["ks"]; m14 = s["moves"]; d = s["delete"]
    # rebuild the plain 15: the frames used, with the array of each move
    import vscan_frames as VF
    fr = [VF.frames[k][f] for k, f in zip(ks, s["frames"])]
    moves = [x for F in fr for x in F]; arr = [l for l, F in enumerate(fr) for x in F]
    I = Inst([tgt_vec(EXP[(2 * k) % 31]) for k in ks], [dir_vec(x) for x in moves])
    # certificate moves of the 14-set mapped into the 15-set: same directions, same multiplicities;
    # take U = all copies of the directions in s['U_dirs'] plus every copy of the deleted direction
    Udirs = set(s["U_dirs"]) | {d}
    U = [i for i, x in enumerate(moves) if x in Udirs]
    Um = sum(1 << i for i in U)
    sups = {}
    for q in range(1, 1 << 15):
        dm = int(I.Dmask[q])
        if dm and dm & ~Um == 0:
            sups[dm] = tuple(int(I.Pmask[q]) >> l & 1 for l in range(3))
    sups = [([i for i in U if dm >> i & 1], p) for dm, p in sups.items()]
    R = VRing(3, U)
    for _ in range(3):
        order = U[:]; random.shuffle(order); pos = {x: k for k, x in enumerate(order)}
        pt = {}
        for l in range(3):
            for i in U:
                pt[R.a[l, i]] = Fraction(int(arr[i] == l)); pt[R.b[l, i]] = Fraction(int(arr[i] == l))
        for i in U:
            for j in U:
                if i != j: pt[R.g[i, j]] = Fraction(2 * int(arr[i] == arr[j]))
        for f in equations(R, order, sups):
            npolys += 1; nz += (evalp(f, pt) != 0)
out.append(f"(2) plain programs substituted into my equations (44 plain multisets, U = certificate directions plus "
           f"the restored one, 3 random orders each): {npolys} polynomials, {nz} nonzero")
open("out_vmodel.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
