import sys, itertools
from vfield import Model
m = int(sys.argv[1]) if len(sys.argv) > 1 else 5
M = Model(m); N = (1 << m) - 1
print("self-dual basis", M.B)
# check symmetric and self-dual
for k in range(N):
    T = M.target(k); assert all(T[i][j] == T[j][i] for i in range(m) for j in range(m))
DM = {z: M.dmask(z) for z in range(1, N + 1)}
bucket = {}
for S in itertools.combinations(range(1, N + 1), m):
    v = 0
    for z in S: v ^= DM[z]
    bucket.setdefault(v, []).append(frozenset(S))
cnt = {k: len(bucket.get(M.smask(M.target(k)), [])) for k in range(N)}
print("frames per T_k:", cnt)
# also: any multiset decomposition with fewer than m distinct directions? rank argument says no; check sets of size < m
small = 0
for r in range(1, m):
    for S in itertools.combinations(range(1, N + 1), r):
        v = 0
        for z in S: v ^= DM[z]
        if any(v == M.smask(M.target(k)) for k in range(N)): small += 1
print("decompositions with fewer than m directions:", small)
# distinct target matrices
print("distinct T_k:", len({M.smask(M.target(k)) for k in range(N)}))
# T_0 = I ?
print("T_0 == I:", M.target(0) == [[int(i == j) for j in range(m)] for i in range(m)])
