# kron_check.py (meta/ex-laser): Kronecker (block-sum) powers of a mixed rule never help (PROOF.md Prop. 4).
# N-th Kronecker power of a rule with target ranks k_w (W arrays) and M rank-1 moves: W^N arrays, target rank of
# a tuple = sum of its ranks, cost N W^(N-1) M (run the rule on each factor for every slice).  tau_N solves
# sum_tuples (sum k)^tau = N W^(N-1) M.  Prediction: tau_N increases to 1.
from math import comb, log
from itertools import product
def tau_of(vals, cost):
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if sum(n * k ** mid for k, n in vals.items()) > cost: hi = mid
        else: lo = mid
    return lo
for name, ranks, M in (("1+4 in 10 (block 3)", [3, 2, 2, 2, 2], 10), ("7+38 in 88 (block 3)", [3] * 7 + [2] * 38, 88),
                       ("5 in 14 (block 3)", [3] * 5, 14)):
    W = len(ranks); cnt = {}
    for k in ranks: cnt[k] = cnt.get(k, 0) + 1
    row = []
    for N in (1, 2, 3, 5, 8, 12):
        dist = {0: 1}
        for _ in range(N):
            nd = {}
            for s, c in dist.items():
                for k, ck in cnt.items(): nd[s + k] = nd.get(s + k, 0) + c * ck
            dist = nd
        t = tau_of(dist, N * W ** (N - 1) * M)
        row.append(f"N={N}: tau {t:.4f}")
    print(name + ": " + "; ".join(row))
