# mixed_room.py (meta/literature, 14 Oct 2026).  Transplant T1 (Schoenhage tau-theorem): how much exponent room does
# Theorem R leave for MIXED units at block 3 (3 rows, widened game) that deliver the kernel (rank 3) to W1 arrays and a
# rank-2 target (one spare row) to W2 arrays with M rank-1 moves?
#   Lower bound (Theorem R, any certificate c):  M >= W1*r3(c) + W2*r2(c);  best over c (a direct sum is a mediant for
#   any FIXED linear objective, so single certificates suffice, pf-duality Lemma 1.2).
#   Exponent of a mixed unit (tau-analog, tested in tau_schedule_lp.py):  W1*3^tau + W2*2^tau = M.
# Room = the smallest tau not excluded = min over ratio r = W2/W1 of the tau solving 3^tau + r*2^tau = max_c (r3 + r*r2).
import json, sys
from fractions import Fraction as F
from math import log, ceil
pts = {}
for line in open(sys.argv[1] if len(sys.argv) > 1 else 'out_venn_sweep_m3.jsonl'):
    d = json.loads(line); key = (d['a'], d['b'], tuple(d['J']))
    pts[key] = (F(d['f3'], d['cap']), F(d['f2'], d['cap']))
print(f"certificates read: {len(pts)}")
best3 = max(pts.items(), key=lambda kv: kv[1][0]); best2 = max(pts.items(), key=lambda kv: kv[1][1])
print("best on the kernel  :", best3[0], [float(x) for x in best3[1]])
print("best on rank-2 (1 spare row):", best2[0], [float(x) for x in best2[1]])
# upper hull (Pareto front) of (r3, r2)
front = sorted({v for v in pts.values()}, key=lambda v: (-v[0], -v[1]))
pareto = []; top2 = -1
for v in front:
    if v[1] > top2: pareto.append(v); top2 = v[1]
print("Pareto front (r3, r2):", [(round(float(a), 4), round(float(b), 4)) for a, b in pareto])
def bound(r): return max(float(a) + r * float(b) for a, b in pareto)
def tau_of(r):
    lo, hi = 0.0, 1.0; B = bound(r)
    for _ in range(60):
        mid = (lo + hi) / 2
        if 3 ** mid + r * 2 ** mid >= B: hi = mid
        else: lo = mid
    return hi
grid = [i / 100 for i in range(0, 3001)]
taus = [(tau_of(r), r) for r in grid]
tmin, rmin = min(taus)
t3 = log(float(best3[1][0])) / log(3); t2 = log(float(best2[1][1])) / log(2)
print(f"pure kernel room: tau >= {t3:.4f} (saving <= {1-t3:.4f});  pure rank-2 room: tau >= {t2:.4f} (saving <= {1-t2:.4f})")
print(f"MIXED room: tau >= {tmin:.4f} (saving <= {1-tmin:.4f}) at r = W2/W1 = {rmin}")
# smallest integer mixed units not excluded by these certificates, ranked by saving
rows = []
for W1 in range(1, 13):
    for W2 in range(0, 40):
        B = max(W1 * a + W2 * b for a, b in pareto); Mmin = ceil(B - 1e-12) if B != int(B) else int(B)
        for M in range(Mmin, 3 * W1 + 2 * W2):
            lo, hi = 0.0, 1.0
            for _ in range(60):
                mid = (lo + hi) / 2
                if W1 * 3 ** mid + W2 * 2 ** mid >= M: hi = mid
                else: lo = mid
            rows.append((1 - hi, W1, W2, M, W1 + W2))
rows.sort(reverse=True)
print("best not-excluded mixed units (saving, W1 kernels, W2 rank-2, M moves):")
seen = 0
for s, W1, W2, M, W in rows:
    if W1 == 0: continue
    print(f"   saving {s:.4f}: {W1} kernel + {W2} rank-2 in {M} moves (plain {3*W1+2*W2})"); seen += 1
    if seen == 8: break
small = [x for x in rows if x[4] <= 14 and x[1] >= 1]
small.sort(reverse=True)
print("best with at most 14 arrays in all:", [(round(s, 4), W1, W2, M) for s, W1, W2, M, W in small[:6]])
