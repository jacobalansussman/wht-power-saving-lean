"""x-pack: an exact LOWER BOUND on the pad price of every bank table of Shea's C588 residual list, from one
invariant, and the search for the cheapest pad multiset meeting it.

Invariant (any table): with n_b items and pad s_b in bank b,  sum_b [73 (9 - n_b) - 9 s_b] = 9*T - 73*I = 46368,
where T = 72*5068 (total width) and I = 72*616 (total items); per bank n_b >= ceil((73 - s_b)/9) (widths <= 9).
Pad-0 banks contribute 73 (9 - n) <= 0.  So the pad banks must produce w-sum 46368 + 73 K (K >= 0).
Cheapest per unit of w by far: (9 x 8), pad 1, w = 64.  This script enumerates N1 pad-1 banks plus up to
three other pad banks (any s, any admissible n) and prints the cheapest combinations (price sum g_a(s))."""
import sys, math
M = 73; T = 72 * 5068; I = 72 * 616
TARGET = 9 * T - 73 * I
a = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0012383
def g(s): return 0.0 if s == 0 else s * ((M / s) ** a - 1)
print('invariant target', TARGET)
opts = []   # (w, s, n) for pad banks other than (s=1, n=8)
for s in range(1, M):
    nmin = -(-(M - s) // 9)
    for n in range(max(nmin, 1), nmin + 12):
        if (s, n) == (1, 8): continue
        opts.append((M * (9 - n) - 9 * s, s, n))
best = []
for N1 in range(690, 730):
    rem0 = TARGET - 64 * N1
    base = N1 * g(1)
    # zero, one, two, three extra pad banks; remaining must be 73K, K >= 0 (absorbed by pad-0 banks with n > 9)
    def ok(r): return r <= 0 and (-r) % M == 0
    if ok(rem0): best.append((base, N1, ()))
    for i, (w1, s1, n1) in enumerate(opts):
        c1 = base + g(s1)
        if ok(rem0 - w1): best.append((c1, N1, ((s1, n1),)))
        for j in range(i, len(opts)):
            w2, s2, n2 = opts[j]
            c2 = c1 + g(s2)
            if ok(rem0 - w1 - w2): best.append((c2, N1, ((s1, n1), (s2, n2))))
best.sort()
LP = TARGET / 64 * g(1)
print('LP bound (P/a) %.4f' % (LP / a))
for c, N1, ex in best[:12]:
    print('P/a %.4f  (+%.4f over LP)  pad-1 banks %d, other pad banks (s,n) %s' % (c / a, (c - LP) / a, N1, ex))
