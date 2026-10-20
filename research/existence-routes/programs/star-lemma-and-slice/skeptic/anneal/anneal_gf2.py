#!/usr/bin/env python3
"""anneal_gf2.py (v-e3-slice6): falsification attempt for Theorem C. Simulated annealing over GF(2) (the star lemma
claims to hold over every field) for a NONSINGULAR 2^n x 2^n matrix with rank <= s on every affine hyperplane of
AG(n,2). Calibration: L(3,2) and L(4,4) are feasible over GF(2) (Singer permutations); L(3,1), L(4,3) are claimed
infeasible. Cost = sum_H max(0, rank_H - s) + 4 * (2^n - rank E). Rows as bitmasks. Usage: n s seed steps"""
import random, sys, math
n, s, seed, steps = map(int, sys.argv[1:5]); N = 1 << n; random.seed(seed)
pc = lambda x: bin(x).count('1')
hyps = []
for l in range(1, N):
    for c in (0, 1): hyps.append([x for x in range(N) if pc(x & l) % 2 == c])
def rank_rows(rows):
    r = 0; rows = list(rows)
    while rows:
        p = rows.pop()
        if p == 0: continue
        r += 1; lb = p & -p
        rows = [q ^ p if q & lb else q for q in rows]
    return r
def sub_rows(E, H):
    out = []
    for x in H:
        v = 0
        for j, y in enumerate(H):
            if (E[x] >> y) & 1: v |= 1 << j
        out.append(v)
    return out
def cost(E):
    c = 4 * (N - rank_rows(E))
    for H in hyps: c += max(0, rank_rows(sub_rows(E, H)) - s)
    return c
best_overall = None
for restart in range(1000):
    perm = list(range(N)); random.shuffle(perm)
    E = [1 << perm[x] for x in range(N)] if restart % 2 == 0 else [random.getrandbits(N) & random.getrandbits(N) & random.getrandbits(N) for _ in range(N)]
    cur = cost(E); T = 2.0; best = cur
    for it in range(steps):
        x = random.randrange(N); y = random.randrange(N)
        E[x] ^= 1 << y
        new = cost(E)
        if new <= cur or random.random() < math.exp((cur - new) / T): cur = new
        else: E[x] ^= 1 << y
        best = min(best, cur); T = max(0.05, T * 0.9995)
        if cur == 0: break
    best_overall = cur if best_overall is None else min(best_overall, cur)
    print("n=%d s=%d restart %d final cost %d best %d" % (n, s, restart, cur, best), flush=True)
    if cur == 0:
        print("FOUND feasible matrix rows:", [format(r, '0%db' % N) for r in E], flush=True); break
    if restart >= int(sys.argv[5]) - 1: break
print("done; best final cost", best_overall)
