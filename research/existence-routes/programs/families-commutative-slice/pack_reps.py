# falsification of Corollaries 13/15 by plain-type units: W disjoint representations of I_m (multisets whose classes
# (z_j z_k)_{j<=k} sum to the class of the identity) packed inside U_m^(r) would be an EXISTING unit with that multiset
# (unused moves idle), so the corollaries would be false.  Search: m = 4, all representations of size <= 6; can W of
# them be packed into r copies of each direction?  (exact ILP-free search: DFS with counts)
import sys, itertools
m, r, W = map(int, sys.argv[1:4])
P = [(j, k) for j in range(m) for k in range(j, m)]
def cls(z): return tuple(((z >> j) & 1) * ((z >> k) & 1) for (j, k) in P)
I = tuple(1 if j == k else 0 for (j, k) in P)
n = (1 << m) - 1
reps = []
for size in range(m, m + 3):
    for combo in itertools.combinations_with_replacement(range(1, n + 1), size):
        s = [0] * len(P)
        for z in combo:
            c = cls(z); s = [(a + b) % 2 for a, b in zip(s, c)]
        if tuple(s) == I: reps.append(combo)
print(f"m={m}: representations of I of size {m}..{m+2}: {len(reps)}; size {m}: {sum(1 for x in reps if len(x)==m)}")
cap = [r] * (n + 1)
best = [0]
def dfs(k, start, used):
    if k > best[0]: best[0] = k
    if k == W: return True
    if used + (W - k) * m > r * n: return False
    for i in range(start, len(reps)):
        R = reps[i]
        ok = all(cap[z] >= R.count(z) for z in set(R))
        if not ok: continue
        for z in R: cap[z] -= 1
        if dfs(k + 1, i, used + len(R)): return True
        for z in R: cap[z] += 1
    return False
found = dfs(0, 0, 0)
print(f"pack {W} disjoint representations into U_{m}^({r}) ({r*n} moves): {'FOUND' if found else 'impossible'}; best packing {best[0]}")
