"""ex-legit: which roots of unity does upstream's `order` supply, and how much room is left under 1024 n^3?

Port of OAI/Computability/FourierTransform/{PrimeSelection,LengthSelection}.lean (r6 candidate tree =
upstream fd4aeeb): stackN, silo, harvestN, reap, bed, fill, bead, barge, orderE, Ports, express, receipt.
upstream proves `receipt_bound : 0 < receipt n ∧ receipt n < 8*n^3` (LengthSelection.lean:262) and
`transform_main_order` (Main.lean:122): run order n = receipt n.

For each n we print: Ports n, d = receipt n, v2(d) (so zeta_{2^j} is a power of the supplied root for j <= v2(d)),
the odd part's prime factors, and the headroom 1024 n^3 / d (any factor below it can be multiplied into d).
"""
import sys
from math import log2


def nat_log2(x):          # Nat.log 2 x
    return x.bit_length() - 1 if x > 0 else 0


def stackN(n):
    return nat_log2(2 * n) + 1


def silo(k):
    return (16 * (k + 4)) ** 2


def harvestN(n):
    return silo(stackN(n))


_sieve_cache = {}


def primes_below(P):
    if P in _sieve_cache:
        return _sieve_cache[P]
    s = bytearray([1]) * max(P, 2)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(P ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    _sieve_cache[P] = s
    return s


def reap(N, P):           # iterate p = 0 .. P-1: push p if Bar p (odd prime) and vol*p <= N
    s = primes_below(P)
    vol, qs = 1, []
    for p in range(P):
        if p > 2 and s[p] and vol * p <= N:
            vol *= p
            qs.append(p)
    return vol, qs


def bed(n):
    return reap(2 * n, harvestN(n) + 1)


def fill(N, p, i):
    q = 1
    for _ in range(i):
        if p * q < N:
            q *= 2
    return q


def bead(n):
    return fill(2 * n, bed(n)[0], stackN(n))


def orderE(n):
    return nat_log2(2 * harvestN(n)) + 2


def Ports(n):
    return 2 ** orderE(n) <= n


def receipt(n):
    if Ports(n):
        vol, qs = bed(n)
        return (2 * n) * (vol * bead(n)) * 2 ** orderE(n), qs
    return n, None


def v2(x):
    return (x & -x).bit_length() - 1


def main():
    ns = [2 ** k for k in range(1, 41, 3)] + [3 ** 13, 10 ** 6, 10 ** 7 + 19, 2 ** 24 + 1, 10 ** 9 + 7, 2 ** 40 - 1]
    first_ports = None
    lo, hi = 1, 2 ** 40
    while lo < hi:        # Ports is monotone in n (2^orderE n grows like log^2 n): find the first n with Ports
        mid = (lo + hi) // 2
        if Ports(mid):
            hi = mid
        else:
            lo = mid + 1
    first_ports = lo
    print(f"first n with Ports n (upstream switches from d = n to d = express n): {first_ports} "
          f"(2^{log2(first_ports):.2f}); orderE there = {orderE(first_ports)}")
    print(f"{'n':>16} {'Ports':>5} {'log2 d':>8} {'v2(d)':>6} {'1024n^3/d':>12} {'odd primes in vol':>18}")
    for n in sorted(ns):
        d, qs = receipt(n)
        assert 0 < d < 8 * n ** 3, n
        assert d % n == 0, n
        room = 1024 * n ** 3 / d
        print(f"{n:>16} {str(Ports(n)):>5} {log2(d):8.2f} {v2(d):6d} {room:12.4g} "
              f"{('%d primes, max %d' % (len(qs), max(qs))) if qs else '-':>18}")
    print("checks: 0 < d < 8 n^3 and n | d for every n above (asserted)")


if __name__ == "__main__":
    sys.exit(main())
