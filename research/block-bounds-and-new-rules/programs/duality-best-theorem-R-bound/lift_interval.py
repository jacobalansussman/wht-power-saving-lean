# lift_interval.py -- lower bounds for the interval relaxation Q_m by LIFTS (pf-duality, 13 Oct 2026).
#
# Interval at block m: frames N_s = s (+) s^0 inside F_2^(2m) = L (+) M, one for each subspace s of L = F_2^m
# (s^0 = annihilator of s in M = dual of L).  A lift gives each frame an integer m x 2m matrix B_s whose rows
# reduce mod 2 to a basis of N_s; Y_s = its row space over Q.  Then delta = m on every complementary pair
# (mod-2 rank already 2m), and on an edge delta = 1 + excess, excess = rank_Q[B_s;B_t] - (m+1) >= 0.
# ratio = m / (1 + E/#edges), E = total excess over unordered edges (delta symmetric: equal dimensions).
# Simulated annealing on the entries (each entry moves by +-2, so the mod-2 reduction is fixed).
# Final answers are re-checked with exact rational ranks.
import sys, random, math
from fractions import Fraction
sys.path.insert(0, __import__('os').path.dirname(__file__))
from interval import subspaces, dimof, rank_frac, rank_mod

def basis_of(S, nbits):
    """a basis (list of ints) of the F_2-space S (frozenset of ints)"""
    basis = []; span = {0}
    for v in sorted(S):
        if v not in span:
            basis.append(v); span |= {x ^ v for x in span}
    return basis

def frames(m):
    subs = subspaces(m)
    out = []
    full = (1 << m) - 1
    for s in subs:
        # annihilator in M: xi with popcount(xi & x) even for all x in s
        ann = frozenset(xi for xi in range(1 << m) if all(bin(xi & x).count('1') % 2 == 0 for x in s))
        vecs = [x for x in basis_of(s, m)] + [xi << m for xi in basis_of(ann, m)]
        assert len(vecs) == m
        out.append((s, vecs))
    return out

def edges(subs):
    E = []
    for i, s in enumerate(subs):
        for j, t in enumerate(subs):
            if i < j:
                ds, dt = dimof(s), dimof(t)
                inter = s & t
                d = ds + dt - 2 * dimof(inter)
                if d == 1: E.append((i, j))
    return E

def rank_float(rows):
    return rank_mod(rows)
def _unused(rows):
    import numpy as np
    return int(np.linalg.matrix_rank(np.array(rows, dtype=float), tol=1e-7))

def run(m, seed, steps, T0=1.0, init='zero01'):
    rnd = random.Random(seed)
    fr = frames(m); subs = [s for s, _ in fr]
    E = edges(subs)
    nb = [[] for _ in subs]
    for k, (i, j) in enumerate(E):
        nb[i].append(k); nb[j].append(k)
    # integer lifts: start from the 0/1 vectors
    B = []
    for s, vecs in fr:
        B.append([[(v >> b) & 1 for b in range(2 * m)] for v in vecs])
    def exc(k):
        i, j = E[k]
        return rank_float(B[i] + B[j]) - (m + 1)
    ex = [exc(k) for k in range(len(E))]
    cur = sum(ex); best = cur; bestB = [ [row[:] for row in Bi] for Bi in B]
    for step in range(steps):
        T = T0 * (1 - step / steps) + 1e-3
        i = rnd.randrange(len(B)); r = rnd.randrange(m); c = rnd.randrange(2 * m)
        old = B[i][r][c]
        new = old + rnd.choice([-2, 2])
        if abs(new) > 3: continue
        B[i][r][c] = new
        newex = {k: exc(k) for k in nb[i]}
        delta = sum(newex[k] - ex[k] for k in nb[i])
        if delta <= 0 or rnd.random() < math.exp(-delta / T):
            for k, v in newex.items(): ex[k] = v
            cur += delta
            if cur < best:
                best = cur; bestB = [[row[:] for row in Bi] for Bi in B]
        else:
            B[i][r][c] = old
    # exact re-check
    Ex = 0
    for (i, j) in E:
        assert rank_frac(bestB[i]) == m and rank_frac(bestB[j]) == m
        Ex += rank_frac(bestB[i] + bestB[j]) - (m + 1)
    # complementary pairs: delta must be m
    nE = len(E)
    return best, Ex, nE, bestB

if __name__ == "__main__":
    m = int(sys.argv[1]); seeds = int(sys.argv[2]); steps = int(sys.argv[3])
    res = []
    for sd in range(seeds):
        b, Ex, nE, BB = run(m, sd, steps)
        ratio = Fraction(m) / (1 + Fraction(Ex, nE))
        print(f"m={m} seed={sd}: E(float)={b} E(exact)={Ex} edges={nE} ratio={ratio} = {float(ratio):.5f}")
        sys.stdout.flush()
        res.append((Ex, sd, BB))
    res.sort(key=lambda t: t[0])
    Ex, sd, BB = res[0]
    with open(f"lift_interval_m{m}_best.txt", "w") as f:
        f.write(f"# m={m} seed={sd} E={Ex}\n")
        for Bi in BB: f.write(repr(Bi) + "\n")
