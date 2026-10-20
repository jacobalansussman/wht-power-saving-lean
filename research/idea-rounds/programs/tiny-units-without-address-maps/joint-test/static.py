#!/usr/bin/python3
"""joint/cert/static.py W K -- agent "cert".  STATIC cuts: the same set U at every time (U_0 = ... = U_K = U, and the
best B is then U itself).  The inequality a unit must satisfy is then linear in the numbers n_w of moves per direction:
    sum over w of n_w * c_w(U)  >=  W * D(U),
c_w(U) = labels dir(w) carries into U from outside, D(U) = labels the kernel carries into U from outside.
For every multiset of K directions: is it rejected by some static cut?  Plain sums, no flow."""
import sys, itertools, math
W = int(sys.argv[1]); K = int(sys.argv[2])
v = [list(map(int, l.split())) for l in open("model_cert.txt")]
cap = {w: v[15 * (w - 1):15 * w] for w in range(1, 8)}; dem = v[105:120]
def nm(b): return "(0,0)" if b == 0 else ("(%d,0)" % b if b < 8 else "(%d,1)" % (b - 7))
cuts = []
for comp in (list(range(8)), list(range(8, 15))):
    for r in range(1, len(comp)):
        for U in itertools.combinations(comp, r):
            Us = set(U); out = [s for s in comp if s not in Us]
            D = sum(dem[s][d] for s in out for d in U)
            c = tuple(sum(cap[w][s][d] for s in out for d in U) for w in range(1, 8))
            cuts.append((U, D, c))
print("static cuts: %d sets U" % len(cuts))
ineq = {}
for U, D, c in cuts: ineq.setdefault((c, D), []).append(U)
print("distinct inequalities (coefficients c_1..c_7, D): %d" % len(ineq))
tot = rej = 0; totp = rejp = 0; used = {}; surv = []
for ms in itertools.combinations_with_replacement(range(1, 8), K):
    n = [ms.count(w) for w in range(1, 8)]
    orderings = math.factorial(K)
    for x in n: orderings //= math.factorial(x)
    tot += 1; totp += orderings
    best = None
    for (c, D), Us in ineq.items():
        lhs = sum(n[i] * c[i] for i in range(7))
        if lhs < W * D:
            if best is None or W * D - lhs > best[0]: best = (W * D - lhs, c, D)
    if best: rej += 1; rejp += orderings; used[(best[1], best[2])] = used.get((best[1], best[2]), 0) + 1
    else: surv.append((ms, orderings))
print("W=%d K=%d: multisets of directions %d (ordered patterns %d); rejected by a static cut: %d multisets (%d ordered patterns); NOT rejected by any static cut: %d multisets (%d ordered patterns)" % (W, K, tot, totp, rej, rejp, tot - rej, totp - rejp))
print("static inequalities used (largest margin per multiset), most used first:")
for (c, D), cnt in sorted(used.items(), key=lambda kv: -kv[1])[:12]:
    print("  %4d multisets: sum n_w*c_w >= %d*W with c = %s ; example U = %s (%d sets U give it)" % (cnt, D, list(c), " ".join(nm(b) for b in ineq[(c, D)][0]), len(ineq[(c, D)])))
print("first multisets not rejected by any static cut:", [ms for ms, o in surv[:12]])
with open("static_survivors_W%d_K%d.txt" % (W, K), "w") as fh:
    for ms, o in surv: fh.write(" ".join(map(str, ms)) + "  %d\n" % o)
