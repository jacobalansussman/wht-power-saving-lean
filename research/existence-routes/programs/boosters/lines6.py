# e3-boost lines6.py: incidence certificates on ISOTROPIC LINES at block m (brute force): V = Q^(isotropic lines),
# Y_U = span{ lines W : dim(W cap U) in T }, T a subset of {0, 1, 2}.  The stabiliser of U acts on lines by these
# dimensions, so f(rank-d job) = #{W : a_U(W) in T} - #{W : a_U(W) in T, a_V(W) in T} for dist(U, V) = d [derived
# here, as for C_k; T = {2} is C_2].  Prints the profiles and tests the block-6 survivor S6.
import sys
m = int(sys.argv[1]) if len(sys.argv) > 1 else 6
M = (1 << m) - 1
def omega(a, b):
    ax, az, bx, bz = a >> m, a & M, b >> m, b & M
    return bin((ax & bz) ^ (az & bx)).count("1") & 1
def inL(p, d):  # L_d: x-bits only in positions < d, z-bits only in positions >= d
    x, z = p >> m, p & M
    return (x & ~((1 << d) - 1)) == 0 and (z & ((1 << d) - 1)) == 0
n = 1 << (2 * m)
mem = [[inL(p, d) for p in range(n)] for d in range(m + 1)]
hist = [dict() for _ in range(m + 1)]   # hist[d][(a0, ad)] = count
for p in range(1, n):
    for q in range(p + 1, n):
        r = p ^ q
        if r < q or omega(p, q): continue
        a = [mem[d][p] + mem[d][q] + mem[d][r] for d in range(m + 1)]
        a = [0 if c == 0 else (1 if c == 1 else 2) for c in a]    # 3 points in L <=> line inside L (dim 2)
        for d in range(1, m + 1):
            key = (a[0], a[d]); hist[d][key] = hist[d].get(key, 0) + 1
import itertools
profiles = {}
for r in range(1, 3):
    for T in itertools.combinations(range(3), r):
        prof = []
        for d in range(1, m + 1):
            tot = sum(c for (a0, ad), c in hist[d].items() if a0 in T)
            both = sum(c for (a0, ad), c in hist[d].items() if a0 in T and ad in T)
            prof.append(tot - both)
        profiles[T] = prof
        print(f"m={m} T={T}: profile {prof}", flush=True)
S6 = {2: -1, 5: -4, 3: 2, 4: 1, 6: 2}
if m == 6:
    for T, P in profiles.items():
        print(f"S6 margin under T={T}: {sum(c * P[d - 1] for d, c in S6.items())}")
