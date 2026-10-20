#!/usr/bin/env python3
"""Exact rational certificate of the mixed-target rate inequality (Theorem 1G of PROOF.md), and its refusal one grid
step higher.  Standard library only (fractions, decimal).

The inequality, for a word whose paid blocks of one unit have the profile {r: n_r} (children), with W arrays per
unit of which H take the cut job kc and the other W - H the full job m:

    sum_r n_r (r/m)^z  <=  (W - H) + H (kc/m)^z                                              (R)

At the DFT figure N/10^10 the Lean rate fact is stated at the slightly smaller exponent
z1 = 1 - (100 N + 2)/10^12 (two hundredths of a grid step stronger); this program checks (R) at z1 with exact
rational enclosures (upper bounds for every child term, a lower bound for the cut-target term), and checks that
(R) is FALSE at the grid exponent 1 - (N + 1)/10^10 (lower bounds for the children, an upper bound for the cut
target).  Both enclosures are the ones the Lean module uses:
  * exp(x) <= 1 + x + x^2/2 + x^3/6 + (5/96) x^4 for 0 <= x <= 1, and exp(x) >= 1 + x + x^2/2 + x^3/6 for x >= 0;
  * rational bounds L of log(m/r), certified by a truncated exponential series: an upper L when
    sum_{k<n} L^k/k! >= m/r, a lower L when an upper bound of exp(L/2^s) raised to 2^s is <= m/r.
The decimal root (60 digits) is printed for orientation only; the verdict rests on Fractions.

usage:
  rate_cert.py PROFILE.json N            exit 0 iff (R) holds at z1(N) AND fails at 1 - (N+1)/10^10
  rate_cert.py PROFILE.json --scan       the largest N whose (R) the enclosure proves, and the refusal of N + 1
PROFILE.json: {"children": {"r": n_r, ...}, "W": W, "H": H, "kc": kc, "m": m}.
Prints one JSON object.
"""
import json, sys
from fractions import Fraction
from decimal import Decimal, getcontext, ROUND_CEILING, ROUND_FLOOR

getcontext().prec = 80
K, D1 = 100, 2          # z1 = 1 - (N K + D1)/(10^10 K): the Lean rate fact's exponent
DD = 12                 # decimals of the rational log bounds


def log_upper(m, r):
    """rational L >= log(m/r), with the number of exp-series terms that certify it"""
    x = (Decimal(m).ln() - Decimal(r).ln()) * Decimal(10) ** DD
    L = Fraction(int(x.to_integral_value(rounding=ROUND_CEILING)) + 1, 10 ** DD)
    target = Fraction(m, r)
    s, term, n = Fraction(0), Fraction(1), 0
    while s < target:
        s += term
        n += 1
        term = term * L / n
        assert n < 400
    return L, n


def log_lower(m, r):
    """rational 0 < L <= log(m/r): an upper bound U of exp(L/2^s) with U^(2^s) <= m/r"""
    x = (Decimal(m).ln() - Decimal(r).ln()) * Decimal(10) ** DD
    L = Fraction(int(x.to_integral_value(rounding=ROUND_FLOOR)) - 1, 10 ** DD)
    assert L > 0
    s = 0
    while L / 2 ** s > Fraction(1, 2):
        s += 1
    y = L / 2 ** s
    target = Fraction(m, r)
    n = 2
    while True:
        S, term = Fraction(0), Fraction(1)
        for i in range(n):
            S += term
            term = term * y / (i + 1)
        fact = 1
        for i in range(2, n + 1):
            fact *= i
        U = S + y ** n * (n + 1) / (fact * n)
        if U ** (2 ** s) <= target:
            return L, s, n
        n += 1
        assert n < 60


def exp_up(x):
    assert 0 <= x <= 1
    return 1 + x + x * x / 2 + x ** 3 / 6 + x ** 4 * Fraction(5, 96)


def exp_lo(x):
    assert x >= 0
    return 1 + x + x * x / 2 + x ** 3 / 6


def check_profile(p):
    ch = {int(r): int(n) for r, n in p['children'].items() if int(n)}
    W, H, kc, m = int(p['W']), int(p['H']), int(p['kc']), int(p['m'])
    assert 0 < H < W and 1 < kc < m
    assert all(1 <= r < kc for r in ch), 'not separated: every paid block must be smaller than every job'
    return ch, W, H, kc, m


def holds(p, N):
    """exact: (R) at z1 = 1 - (N K + D1)/(10^10 K); returns (bool, relative margin)"""
    ch, W, H, kc, m = check_profile(p)
    eps = Fraction(N * K + D1, 10 ** 10 * K)
    lhs = Fraction(0)
    for r, n in ch.items():
        L, _ = log_upper(m, r)
        lhs += n * Fraction(r, m) * exp_up(eps * L)
    Ll, _, _ = log_lower(m, kc)
    rhs = (W - H) + H * Fraction(kc, m) * exp_lo(eps * Ll)
    return lhs <= rhs, float((rhs - lhs) / W)


def fails(p, N1):
    """exact: (R) is false at the grid exponent 1 - N1/10^10"""
    ch, W, H, kc, m = check_profile(p)
    eps = Fraction(N1, 10 ** 10)
    lhs = Fraction(0)
    for r, n in ch.items():
        L, _, _ = log_lower(m, r)
        lhs += n * Fraction(r, m) * exp_lo(eps * L)
    Lu, _ = log_upper(m, kc)
    rhs = (W - H) + H * Fraction(kc, m) * exp_up(eps * Lu)
    return lhs > rhs, float((lhs - rhs) / W)


def root(p):
    """60-digit decimal root b of (R) (orientation only)"""
    ch, W, H, kc, m = check_profile(p)
    getcontext().prec = 60
    c = [(Decimal(n), (Decimal(r) / m).ln()) for r, n in ch.items()]
    t = [(Decimal(W - H), Decimal(0)), (Decimal(H), (Decimal(kc) / m).ln())]
    lo, hi = Decimal('0.97'), Decimal(1)
    for _ in range(200):
        z = (lo + hi) / 2
        if sum(n * (l * z).exp() for n, l in c) > sum(n * (l * z).exp() for n, l in t):
            lo = z
        else:
            hi = z
    getcontext().prec = 80
    return 1 - (lo + hi) / 2


def certify(p, N):
    ok, margin = holds(p, N)
    no, over = fails(p, N + 1)
    return dict(N=N, figure='%.7e' % (N / 1e10), holds_at_N=ok, relative_margin=margin,
                fails_at_N_plus_1=no, relative_excess_at_N_plus_1=over, root_decimal=str(root(p))[:40],
                PASS=bool(ok and no))


def scan(p):
    b = root(p)
    N = int(b * 10 ** 10)
    while not holds(p, N)[0]:
        N -= 1
    while holds(p, N + 1)[0]:
        N += 1
    return certify(p, N)


def main():
    p = json.load(open(sys.argv[1]))
    out = scan(p) if sys.argv[2] == '--scan' else certify(p, int(sys.argv[2]))
    print(json.dumps(out, indent=1))
    sys.exit(0 if out['PASS'] else 1)


if __name__ == '__main__':
    main()
