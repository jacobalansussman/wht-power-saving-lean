# monomial E = permutation matrix of sigma: rank E[H,H] = #{i in H : sigma(i) in H}. Search sigma with max over H <= s.
import itertools, sys
import numpy as np
from loc import hyperplanes
n = int(sys.argv[1]); N = 1 << n
Hs = [set(H.tolist()) for H in hyperplanes(n)]
best = None
cnt = {}
for p in itertools.permutations(range(N)):
    mx = max(sum(1 for i in H if p[i] in H) for H in Hs)
    cnt[mx] = cnt.get(mx, 0) + 1
    if best is None or mx < best[0]:
        best = (mx, p)
print("n", n, "distribution of max_H count:", dict(sorted(cnt.items())), "best", best)
