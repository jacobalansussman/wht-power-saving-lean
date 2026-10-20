"""many_primes.py -- ex2-modp: a NON-CONSTRUCTIVE existence proof from F_p points alone (PROOF.md Proposition 5).

  python3 many_primes.py PREFIX [PREFIX2 ...] [--tries T] [--start P]

The joined system (integer coefficients) is rebuilt over Q.  If it had no complex zero, Krick-Pardo-Sombra Theorem 1
gives a nonzero integer a = sum g_i f_i with ln|a| <= H := 4n(n+1)d^n(h + ln s + (n+7) ln(n+1) d); every prime p at
which the system has an F_p point divides a, so at most H/ln(P) such primes are >= P.  We look for F_p points (random
alpha in F_p, beta from the linear system mod p) at consecutive primes >= P = 2^30; finding more than H/ln(P) of them
proves that the system has a complex zero, with no characteristic-0 computation at all.  No point is lifted or
reconstructed; the proof says nothing about where the zero is."""
import sys, os, math, random, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zsys import *

def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):   # deterministic below 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

ap = argparse.ArgumentParser()
ap.add_argument('prefixes', nargs='+'); ap.add_argument('--tries', type=int, default=20)
ap.add_argument('--start', type=int, default=2 ** 30); ap.add_argument('--maxprimes', type=int, default=100000)
a = ap.parse_args()
systems = [load_Q(pf) for pf in a.prefixes]
names, polys = joint(systems); n = len(names)
assert all(c.denominator == 1 for P in polys for c in P.values()), 'clear denominators first'
ai, bi = split_linear(names, polys)
d = max(degree(P) for P in polys); s = len(polys); h = max(math.log(abs(c)) for P in polys for c in P.values())
H = 4 * n * (n + 1) * d ** n * (h + math.log(s) + (n + 7) * math.log(n + 1) * d)
K = math.floor(H / math.log(a.start)) + 1
print(f"{' + '.join(os.path.basename(x) for x in a.prefixes)}: n={n} d={d} s={s} h={h:.3f}: KPS ln|a| <= {H:.4g}; "
      f"need F_p points at {K} primes >= {a.start}", flush=True)
rng = random.Random(1); p = a.start; found = 0; tested = 0
while found < K and tested < a.maxprimes:
    p += 1
    if not is_prime(p): continue
    tested += 1
    for _ in range(a.tries):
        alpha = {i: rng.randrange(p) for i in ai}
        A, r = beta_system(names, polys, ai, bi, alpha, p)
        if solve_mod(A, r, p) is not None:
            found += 1; break
    if tested in (1, 10, 100) or tested % 1000 == 0:
        print(f"   primes tested {tested}, with F_p points {found}", flush=True)
if found >= K:
    print(f"EXISTS [non-constructive: F_p points at {found} distinct primes in ({a.start}, {p}] > KPS bound {H:.4g}/ln({a.start})]")
else:
    print(f"NO CERTIFICATE: F_p points at {found} of {tested} primes (needed {K})")
