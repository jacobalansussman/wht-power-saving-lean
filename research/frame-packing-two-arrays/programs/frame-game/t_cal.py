"""t_cal.py: calibration of fpa.py (frame game): symplectic checks, K, levels, plain units pass, truncations fail."""
import sys, random
from fpa import *
for m in (2, 3, 4, 5):
    S = Sp(m); rng = random.Random(m)
    assert S.is_symp(S.K) and S.level(S.K) == m
    for _ in range(50):
        p = S.rand_P(rng); assert S.is_symp(p) and S.in_P(p) and S.level(p) == 0
        z = rng.randrange(1, 1 << m); assert S.is_symp(S.T(z)) and S.level(S.T(z)) == 1
        A = S.mul(S.T(z), p); assert S.mul(A, S.inv(A)) == S.I
        # level changes by at most 1 under a move
        g = S.I
        for k in range(6):
            g2 = S.mul(S.T(rng.randrange(1, 1 << m)), S.mul(S.rand_P(rng), g))
            assert abs(S.level(g2) - S.level(g)) <= 1; g = g2
    # plain W = 2: states = inputs and the chain outputs
    for W in (1, 2):
        states = [inp(S, l) for l in range(W)]
        for l in range(W):
            c = inp(S, l)
            for j in range(m):
                c = lmul(S, S.T(1 << j), c); states.append(c)
        res = [end_solve(S, states, {(l, S.K): Fr(1)})[0] for l in range(W)]
        # truncated: drop the last move of array 0
        st2 = [s for i, s in enumerate(states) if i != W + m - 1]
        res2 = end_solve(S, st2, {(0, S.K): Fr(1)})
        print(f"m={m} W={W}: plain end solves {res}; array 0 without its last move: {res2[0]} ({res2[1].get('reason')})", flush=True)
    # a scrambled plain: random free letters between moves, a random direction order
    c = inp(S, 0); states = [c]; g_free = S.I
    perm = list(range(m)); rng.shuffle(perm)
    for j in perm:
        c = lmul(S, S.rand_P(rng, levi=False), c)     # phases only: then the next T must be re-aimed
        states.append(c)
        c = lmul(S, S.T(1 << j), c); states.append(c)
    print(f"m={m}: phases inserted between the K directions (generally NOT K): end solve", end_solve(S, states, {(0, S.K): Fr(1)})[0])
