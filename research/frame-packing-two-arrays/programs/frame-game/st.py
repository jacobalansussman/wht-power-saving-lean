"""st.py: single-trunk two-array programs in the frame game (2 in 2m - 1), with designed packing.
Lineages: y1 plain e1..em; the trunk move (shared) is y1's move a1+1; y2 enters the trunk with a random free chart
p2 (relative chart); branch B completes y2's lineage greedily to K.  Injections: extra walkers (stored states with
chosen charts) added to the trunk or to a branch register before a given move.  Then exact end solves."""
import sys, random, time
from fpa import *

def complete_to_K(S, g, rng):
    """random geodesic completion: paid letters T_u (u not in L0, as p T_e1 p^-1) with free letters, so that w g in P K.
    Returns list of (u) left factors: w = T_uk ... T_u1 (each T_u is one paid move times free letters)."""
    m = S.m; word = []
    d = S.level(S.mul(S.K, S.inv(g)))
    allu = [u for u in range(1, 1 << (2 * m)) if (u & ((1 << m) - 1))]   # x-part nonzero: not in L0
    while d > 0:
        rng.shuffle(allu)
        for u in allu:
            Tu = S.transv(u)
            g2 = S.mul(Tu, g)
            if S.level(S.mul(S.K, S.inv(g2))) == d - 1:
                word.append(u); g = g2; d -= 1; break
        else:
            return None
    return word

def transv(S, u):
    m = S.m; n2 = 2 * m
    # T_u(x) = x + w(x,u) u ; w(x,u) = x_a.u_b + x_b.u_a ; row i of matrix: e_i + u_i * (J u)^T
    ua, ub = u & ((1 << m) - 1), u >> m
    Ju = ub | (ua << m)            # w(x,u) = parity(x & Ju)
    rows = []
    for i in range(n2):
        r = 1 << i
        if (u >> i) & 1: r ^= Ju
        rows.append(r)
    return tuple(rows)
Sp.transv = transv

def build(S, m, a1, a2, rng, inj=(), p2=None):
    """returns states list, the two delivering state indices, and walker bookkeeping."""
    y1, y2 = inp(S, 0), inp(S, 1)
    states = [y1, y2]
    # y1 prefix
    c1 = y1
    for k in range(a1):
        c1 = lmul(S, S.T(1 << k), c1); states.append(c1)
    pre1 = list(range(2, len(states)))
    # y2 prefix: plain directions in a random order
    order2 = list(range(m)); rng.shuffle(order2)
    c2 = y2
    for k in range(a2):
        c2 = lmul(S, S.T(1 << order2[k]), c2); states.append(c2)
    # trunk: y1 state + p2 * y2 state ; p2 random until y2's step ascends and stays geodesic to K
    zt = 1 << a1
    for tries in range(2000):
        p = p2 if p2 is not None else S.rand_P(rng)
        g2 = S.mul(S.T(zt), S.mul(p, next(iter(c2)) [1]))
        if S.level(g2) == a2 + 1 and S.level(S.mul(S.K, S.inv(g2))) == m - a2 - 1:
            break
    else:
        return None
    trunk_in = cadd((1, c1), (1, lmul(S, p, c2)))
    for (where, src, chart, co) in inj:
        if where == "trunk": trunk_in = cadd((1, trunk_in), (co, lmul(S, chart, states[src])))
    t1 = lmul(S, S.T(zt), trunk_in); states.append(t1); it = len(states) - 1
    # branch A: y1 continues plain
    cA = t1
    for k in range(a1 + 1, m):
        for (where, src, chart, co) in inj:
            if where == ("A", k): cA = cadd((1, cA), (co, lmul(S, chart, states[src])))
        cA = lmul(S, S.T(1 << k), cA); states.append(cA)
    iA = len(states) - 1
    # branch B: y2 completion
    wB = complete_to_K(S, g2, rng)
    if wB is None or len(wB) != m - 1 - a2: return None
    cB = t1
    for k, u in enumerate(wB):
        for (where, src, chart, co) in inj:
            if where == ("B", k): cB = cadd((1, cB), (co, lmul(S, chart, states[src])))
        cB = lmul(S, S.transv(u), cB); states.append(cB)
    iB = len(states) - 1
    return states, iA, iB, p

if __name__ == "__main__":
    m = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    S = Sp(m); rng = random.Random(seed)
    tally = {}
    t0 = time.time()
    for tr in range(trials):
        a1 = rng.randrange(m); a2 = rng.randrange(m)
        b = build(S, m, a1, a2, rng)
        if b is None: tally["build-fail"] = tally.get("build-fail", 0) + 1; continue
        states, iA, iB, p = b
        r1 = end_solve(S, states, {(0, S.K): Fr(1)})
        r2 = end_solve(S, states, {(1, S.K): Fr(1)})
        key = (a1, a2, r1[0], r2[0])
        tally[key] = tally.get(key, 0) + 1
        if r1[0] == "yes" and r2[0] == "yes":
            print("BOTH YES", a1, a2, p, flush=True)
    print(f"m={m}: {trials} trials in {time.time()-t0:.1f}s", flush=True)
    for k in sorted(tally, key=str): print("  ", k, tally[k])
