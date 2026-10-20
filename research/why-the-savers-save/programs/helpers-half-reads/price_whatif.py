"""d-halfread STEP 3 (WHAT-IF only: no U8 second copy can be retired by half reads, see obstruction.json).
Prices on #942, with m-anatomy's pricer and u-helpers-hub/price_removals.py's two hosting brackets:
 (1) one actual U8 second copy retired in all 333 calls: the entering face copy h581 of plane (265,137,6)
     (chain dims 1,3,6,7,8: runs {2,3,1,1}, width 7);
 (2) the same plus the small case's spend if its two consumers could each add a rank-2 stop at L_a: each splits one
     run (modelled as a 1->6 run of 5 becoming 1->2 + 2->6, i.e. {5} -> {1,4});
 (3) the per-plane ceiling of the device in the small case (one register at P saved per plane, 168 planes) = the
     typical-plane-helper what-if for N = 168 is NOT claimed; only (1)-(2) are priced here.
Usage: python3 -I price_whatif.py"""
import sys, os, json
sys.path.insert(0, __import__('os').environ.get('PRICER_DIR', '.'))  # folder holding pricer.py (it ships with the lever package's tools)
from pricer import price_b, kappa_from, C_PR589, dec
from fractions import Fraction as F
from collections import Counter
P942 = {1: 848490, 2: 112554, 3: 148518, 4: 53946, 5: 23976, 8: 2664, 48: 11086}; M, W, CALLS = 73, 33152, 333
B942 = F(1664727258740678855772, 10**24)
def apply(delta_local, holes_ranks, mode):
    """delta_local: per-call change of the local profile (negative = removed runs). The net freed width per call goes
    back to the idle stock: optimistic = re-pack (fill bins smallest-capacity first), hole = one child per call at each
    rank in holes_ranks."""
    G = Counter(P942)
    for r, n in delta_local.items(): G[r] += CALLS*n
    freed = -CALLS*sum(r*n for r, n in delta_local.items())
    if mode == 'optimistic':
        G[48] -= 11086; G[1] -= 6
        hosted = 37*9*3778 - freed; res = Counter()
        for cap, nb in ((8, 672*37), (16, 448*37), (32, 448*37), (48, 448*37)):
            full = min(nb, hosted//cap); hosted -= full*cap; left = nb - full
            if left and hosted: res[cap-hosted] += 1; left -= 1; hosted = 0
            res[cap] += left
        assert hosted == 0
        G.update(res)
    else:
        assert sum(holes_ranks) * CALLS == freed
        for r in holes_ranks: G[r] += CALLS
    G = {r: x for r, x in G.items() if x}
    assert all(x > 0 for x in G.values()) and sum(r*x for r, x in G.items()) == M*W - 11840
    b = price_b(G, M, W); k = kappa_from(b, C_PR589)[0]
    return b, k
b0 = price_b(P942, M, W); assert b0 == B942
out = {'control_b': dec(b0, 25)}
cases = {
  'one entering second copy (h581, width 7) retired': (Counter({2: -1, 3: -1, 1: -2}), [7]),
  'same + 2 consumers each add a rank-2 stop ({5}->{1,4})': (Counter({2: -1, 3: -1, 1: -2}) + Counter() , [7]),
}
d2 = Counter({2: -1, 3: -1, 1: -2}); d2[5] -= 2; d2[1] += 2; d2[4] += 2
cases['same + 2 consumers each add a rank-2 stop ({5}->{1,4})'] = (d2, [7])
for name, (dl, holes) in cases.items():
    r = {}
    for mode in ('hole', 'optimistic'):
        b, k = apply(dict(dl), holes, mode)
        r[mode] = dict(b=dec(b, 13), b_gain_pct=round(float((b/b0-1)*100), 4), kappa=dec(k, 16))
    out[name] = r; print(name, json.dumps(r))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'price_whatif.json'), 'w'), indent=1)
