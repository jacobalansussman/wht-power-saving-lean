# mylift.py LIFTFILE... -- skeptic's exact re-check of interval lifts over Q(omega) (v-meta2-hermitian-compute).
# Own interval (geo.py), own matching of each lift frame to N_s by its reduction mod 2, own covering edges, and
# Q(omega)-ranks by REALIFICATION: rank_Q of the 2n x 12 integer matrix (each a + b w -> its 2x2 multiplication
# matrix in the basis 1, w; w^2 = -1 - w) equals 2 * rank_{Q(w)}; Q-rank by exact fraction-free (Bareiss) elimination.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geo import *
from fractions import Fraction as Fr

def rank_int(M):            # exact rank over Q with Fractions
    M = [[Fr(x) for x in row] for row in M]; n = len(M); m = len(M[0]) if M else 0
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(r + 1, n):
            if M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
        if r == n: break
    return r

def realify(rows):            # rows: list of lists of (a, b)
    out = []
    for row in rows:
        r1 = []; r2 = []
        for a, b in row:
            # z*1 = a + b w ; z*w = -b + (a - b) w.  Row vector v (over Q(w)) times scalar: represent each
            # Q(w)-row v as two Q-rows: coordinates of v and of w*v.
            r1 += [a, b]
            r2 += [-b, a - b]
        out.append(r1); out.append(r2)
    return out
def rankQw(rows):
    rk2 = rank_int(realify(rows)); assert rk2 % 2 == 0, rk2
    return rk2 // 2

def red(a, b): return (a & 1) | ((b & 1) << 1)
def span6(gens):
    S = {(0,) * 6}
    for g in gens:
        S = {vadd(v, vscale(c, g)) for v in S for c in range(4)}
    return frozenset(S)
HP = [herm_perp(i) for i in range(44)]
BP = [bil_perp(i) for i in range(44)]
def Nset(i, perp):
    return frozenset(tuple(x) + tuple(z) for x in SUBS[i] for z in SUBS[perp[i]])
NH = {Nset(i, HP): i for i in range(44)}
NB = {Nset(i, BP): i for i in range(44)}
# isotropy check of N_s (Hermitian) for h(u,v) = sum conj(u_x) v_z + conj(u_z) v_x
def h(u, v):
    t = 0
    for k in range(3): t ^= fmul(fconj(u[k]), v[3 + k]) ^ fmul(fconj(u[3 + k]), v[k])
    return t
for N, i in NH.items():
    assert all(h(u, v) == 0 for u in N for v in N), i
print("own interval: 44 frames, covering edges", len(COVER), "; Hermitian N_s totally isotropic: yes")
absolute = [i for i in range(44) if DIM[i] == 1 and SUBS[i] <= SUBS[HP[i]]]
print("absolute points (p in p^perp, Hermitian):", len(absolute))

for fn in sys.argv[1:]:
    L = open(fn).read().split("\n")
    nf, ne, best = map(int, L[0].split())
    B = []
    for s in range(nf):
        v = list(map(int, L[1 + s].split())); assert len(v) == 36
        B.append([[(v[2 * (6 * r + c)], v[2 * (6 * r + c) + 1]) for c in range(6)] for r in range(3)])
    who = []
    for s in range(nf):
        R = [tuple(red(a, b) for a, b in row) for row in B[s]]
        S = span6(R); assert len(S) == 64, f"{fn} frame {s}: reduction not rank 3"
        if S in NH: who.append(('H', NH[S]))
        elif S in NB: who.append(('B', NB[S]))
        else: raise SystemExit(f"{fn} frame {s}: reduction is no N_s")
        assert rankQw(B[s]) == 3
    kinds = {k for k, _ in who}; ids = [i for _, i in who]
    assert sorted(ids) == list(range(44)), "frames do not cover the interval exactly once"
    pos = {i: s for s, i in enumerate(ids)}
    E = 0; hist = {}; bad = []
    for i, j in COVER:
        e = rankQw(B[pos[i]] + B[pos[j]]) - 4; assert e >= 0
        E += e; hist[e] = hist.get(e, 0) + 1
        if e: bad.append((i, j, e))
    # all pairs: delta >= d, opposite pairs delta = 3
    mind = 0; opp_ok = True
    for i in range(44):
        for j in range(44):
            if i == j: continue
            dl = rankQw(B[pos[i]] + B[pos[j]]) - 3
            assert dl >= D[i][j]
            if D[i][j] == 3 and dl != 3: opp_ok = False
    tang = sum(1 for i, j, e in bad if DIM[i] == 1 and DIM[j] == 2 and i in absolute and j == HP[i])
    levels = {}
    for i, j, e in bad: levels[(DIM[i], DIM[j])] = levels.get((DIM[i], DIM[j]), 0) + 1
    ratio = Fr(3) / (1 + Fr(E, 147))
    print(f"{fn}: perp kind {kinds}; file says E={best}; exact E = {E}; hist {dict(sorted(hist.items()))}; "
          f"bad-edge levels {levels}; tangent flags (p absolute, l = p^perp) among bad: {tang}/{len(bad)}; "
          f"opposite pairs all delta=3: {opp_ok}; 3/(1+E/147) = {ratio} = {float(ratio):.6f}")
