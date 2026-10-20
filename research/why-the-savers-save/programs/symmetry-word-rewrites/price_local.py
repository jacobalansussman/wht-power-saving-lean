"""Exact hub price of a changed local profile on the #1042 control (144 calls x local + idle; mass change moved to
or from the idle rank-48 remainder, a remainder paid as one child) -- same decomposition as k8-pairs/price_k8.py."""
import sys, json
sys.path.insert(0, '.')
from pricer import price_b, kappa_from, C_PR589, dec
from fractions import Fraction as F
from collections import Counter
P1042 = {1: 349577, 2: 60193, 3: 62608, 4: 20448, 5: 11808, 7: 3, 8: 1152, 16: 47, 48: 4850}; M, W, D = 73, 14336, 5120
LOC = {1: 2412, 2: 418, 3: 432, 4: 142, 5: 82, 8: 8}; CALLS = 144
B0 = F(418747253984418082117, 250000000000000000000000); K0 = F(1672188119205421, 10**18)
def price_local(loc):
    loc = {int(k): int(v) for k, v in loc.items()}
    G = Counter(P1042)
    for r in set(LOC) | set(loc): G[r] += CALLS * (loc.get(r, 0) - LOC.get(r, 0))
    dm = CALLS * (sum(r * x for r, x in LOC.items()) - sum(r * x for r, x in loc.items()))
    if dm:
        q, rem = divmod(dm, 48); G[48] += q
        if rem: G[rem] += 1
    G = {r: x for r, x in G.items() if x}
    assert min(G.values()) > 0 and sum(r * x for r, x in G.items()) == M * W - D
    b = price_b(G, M, W); k = kappa_from(b, C_PR589)[0]
    return dict(profile=dict(sorted(G.items())), b=dec(b, 25), kappa=dec(k, 19), db=dec(b - B0, 6), dkappa=dec(k - K0, 6))
if __name__ == '__main__':
    print(json.dumps(price_local(json.loads(sys.argv[1])), indent=1))
