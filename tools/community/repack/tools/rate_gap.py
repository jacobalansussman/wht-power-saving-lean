"""x-pack: where the grid points between the exact moment root and the Lean kernel numerator go, and what a
higher-order exponential bound (Real.exp_bound' with n terms: exp x <= sum_{i<n} x^i/i! + x^n (n+1)/(n! n) on
[0,1]; the engine's term_le uses n = 3, i.e. 1 + x + x^2/2 + (2/9) x^3) would give.
Usage: rate_gap.py <schedule.json.gz> [more]"""
import sys, os, json, gzip
from fractions import Fraction as F
from math import factorial
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import foldlib as fl, gen_rates
from price_sched import invoice
from exact_moment import certify


def poly(x, n):
    tot = F(0); term = F(1)
    for i in range(n):
        tot += term; term = term * x / (i + 1)
    return tot + x ** n * F(n + 1, factorial(n) * n)


def plan_block_n(u, A, B, s, n, Dd=12):
    eps = F(A, 10 ** B); t = 2 ** s
    m, W, H = u['m'], u['W'], u['H']
    rows = []
    for r, cnt in [(r, H[r]) for r in sorted(H)] + [(m - 1, 1), (1, 1)]:
        L, N = gen_rates.log_upper(m, r, Dd)
        x = eps * L
        assert x <= 1
        rows.append(cnt * F(r, m) * poly(x, n))
    nH = len(H)
    exact = (t - 1) * sum(rows[:nH]) + W * (rows[nH] + rows[nH + 1])
    return t * W - exact > 0


def best_A(u, A0, n, s=40, B=10):
    A = A0
    while not plan_block_n(u, A, B, s, n): A -= 1
    return A


if __name__ == '__main__':
    for p in sys.argv[1:]:
        sch = json.loads(gzip.decompress(open(p, 'rb').read()))
        hist, W = invoice(sch['pads'])
        u = dict(m=73, W=W, D=2560, H=hist)
        rate = certify(hist, 73, W)
        A0 = int(F(rate['root_lower']) * 10 ** 10)
        Aex = fl.best_units(u, 10, 40)[0]            # exact (high-precision decimal) finite-fill s=40 limit
        A3 = best_A(u, A0, 3)
        assert fl.plan_block(u, A3, 10, 40) is not None and fl.plan_block(u, A3 + 1, 10, 40) is None
        A4 = best_A(u, A0, 4); A5 = best_A(u, A0, 5)
        print(f"{os.path.basename(p)}: banks={len(sch['pads'])} exact root grid={A0}  finite-fill s=40 exact limit={Aex}  "
              f"Lean term_le (n=3) A={A3}  n=4 A={A4}  n=5 A={A5}")
