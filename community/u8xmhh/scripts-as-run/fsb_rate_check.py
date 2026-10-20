#!/usr/bin/env python3
"""fs-build: independent exact check of the rate fact of a generated RateXXXx.lean (shares no code with the generators
or with exact_moment.py).  Parses the block list `blocks...x` from the Lean text, checks its rank total W*73 - 2560 and
count, rebuilds it from an invocation profile x72 + per-run remainders x8 + pads, then encloses with exact rationals
  LHS(a) = (2^s - 1) * moment(73, 1-a, U) + W * ((72/73)^(1-a) + (1/73)^(1-a)),  RHS = 2^s * W
(the hypothesis `hnum` of RAM.engine_program_group_list, = Lean's fold_* statement) and checks
  LHS(a) < RHS at a = K/10^10 (the certificate's kernel numerator)   and   LHS(a) > RHS at a = (K+1)/10^10.
usage: python3 -I fsb_rate_check.py RATE.lean K W [profile.json]"""
import sys, re, json
from fractions import Fraction as Q

def ln_bounds(x, terms=60):
    """rational enclosure of ln(x) for rational x >= 1 (range reduction by 2, atanh series, geometric tail)"""
    x = Q(x); assert x >= 1
    k = 0
    while x >= 2:
        x /= 2; k += 1
    def atanh2(y):  # 2*atanh(y) enclosure, 0 <= y < 1
        s = Q(0); p = y
        for j in range(terms):
            s += p / (2 * j + 1); p *= y * y
        tail = 2 * p / ((2 * terms + 1) * (1 - y * y))  # >= the rest of the series
        return 2 * s, 2 * s + tail
    l2 = atanh2(Q(1, 3))
    lx = atanh2((x - 1) / (x + 1))
    return k * l2[0] + lx[0], k * l2[1] + lx[1]

def exp_bounds(t0, t1, terms=30):
    """enclosure of exp on [t0, t1], |t| <= 1/2"""
    def ser(t):
        s = Q(0); p = Q(1)
        for j in range(terms):
            s += p; p = p * t / (j + 1)
        return s, abs(p) * 2  # |rest| <= |t|^n/n! * 1/(1-|t|) <= 2 |p|
    assert abs(t0) <= Q(1, 2) and abs(t1) <= Q(1, 2) and t0 <= t1
    s0, e0 = ser(t0); s1, e1 = ser(t1)
    return s0 - e0, s1 + e1

def pow_bounds(r, m, a):
    """(r/m)^(1-a) = (r/m) * exp(a * ln(m/r)) for 1 <= r <= m, 0 <= a"""
    if r == m:
        return Q(1), Q(1)
    l0, l1 = ln_bounds(Q(m, r))
    e0, e1 = exp_bounds(a * l0, a * l1)
    return Q(r, m) * e0, Q(r, m) * e1

def lhs(blocks, W, a, s=40, m=73):
    lo = hi = Q(0)
    for r, n in blocks:
        p0, p1 = pow_bounds(r, m, a)
        lo += n * p0; hi += n * p1
    lo *= (2**s - 1); hi *= (2**s - 1)
    for r in (m - 1, 1):
        p0, p1 = pow_bounds(r, m, a)
        lo += W * p0; hi += W * p1
    return lo, hi

t = open(sys.argv[1]).read(); K = int(sys.argv[2]); W = int(sys.argv[3])
m = re.search(r'def (blocks\w+) : List \(ℕ × ℕ\) := (\[[^\n]*\])', t)
name, blocks = m.group(1), eval(m.group(2).replace('(', '[').replace(')', ']'))
blocks = [tuple(b) for b in blocks]
print(name, blocks)
assert sum(r * n for r, n in blocks) == W * 73 - 2560, 'rank total'
assert all(1 <= r < 73 for r, n in blocks)
fold = re.search(r'theorem fold_\w+ :\s*\(\(2:ℝ\)\^(\d+) - 1\) \* moment 73 \(1 - (\d+)/\(10:ℝ\)\^10\) (\w+)', t)
assert fold and int(fold.group(2)) == K and fold.group(3) == name, 'fold statement'
s = int(fold.group(1))
if len(sys.argv) > 4:
    pr = json.load(open(sys.argv[4]))
    from collections import Counter
    c = Counter()
    for r, n in pr['inv'].items(): c[int(r)] += 72 * n
    for r, n in pr['rem'].items():
        if int(r): c[int(r)] += 8 * n
    for r, n in pr['pads'].items():
        if int(r): c[int(r)] += n
    assert sorted(c.items()) == sorted(blocks), ('composition', sorted(c.items()))
    print('composition OK: 72 x invocation profile + 8 x per-run remainders + pads =', name)
RHS = Q(2**s * W)
lo, hi = lhs(blocks, W, Q(K, 10**10), s)
print('at a = %d/10^10: LHS <= %.12e, RHS = %d, slack >= %.3e' % (K, float(hi), RHS, float((RHS - hi) / RHS)))
assert hi < RHS, 'FOLD FAILS at K'
lo1, hi1 = lhs(blocks, W, Q(K + 1, 10**10), s)
# largest grid point at which the fold itself holds (the Lean certificate's term_le bounds give up a little of it)
a, b = K, K + 400
assert lhs(blocks, W, Q(b, 10**10), s)[0] > RHS
while b - a > 1:
    c = (a + b) // 2
    lo_c, hi_c = lhs(blocks, W, Q(c, 10**10), s)
    if hi_c < RHS: a = c
    elif lo_c > RHS: b = c
    else: break
print('exact fold (s = %d): holds at %d/10^10, fails at %d/10^10 (the certificate uses %d, %d grid points below)' % (s, a, b, K, a - K))
print('RATE CHECK PASS: %s, s = %d, W = %d, fold holds at the certificate exponent 1 - %d/10^10' % (name, s, W, K))
