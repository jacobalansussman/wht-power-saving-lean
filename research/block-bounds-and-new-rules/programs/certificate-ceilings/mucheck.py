# b56-ceilings: the numbers mu_l = dim (S^(n-l,l))^H / dim S^(n-l,l), H = <(x x^e1)> = (Z/2)^(n/2) (the pairs of dir(e1)).
# dim (M^l)^H = number of H-orbits on l-subsets = sum_q C(P,q) C(P-q, l-2q)  (q full pairs, l-2q half pairs), P = n/2.
# M^l = M^(l-1) (+) S^(n-l,l) for l <= n/2, so dim (S^(n-l,l))^H = orb(l) - orb(l-1).
# Claim checked here EXACTLY for n = 4, 8, ..., 256: mu_1 = (n-2)/(2(n-1)) and mu_l < mu_1 for 2 <= l <= n/2
# (the spectral ceiling for all Venn certificates uses max_{l>=1} mu_l = mu_1).
from math import comb
from fractions import Fraction
def C(x, y): return comb(x, y) if 0 <= y <= x else 0
for m in range(2, 11):
    n = 2 ** m; P = n // 2
    orb = [sum(C(P, q) * C(P - q, l - 2 * q) for q in range(l // 2 + 1)) for l in range(n // 2 + 1)]
    mu = [Fraction(orb[l] - (orb[l - 1] if l else 0), C(n, l) - C(n, l - 1)) for l in range(n // 2 + 1)]
    assert mu[0] == 1 and mu[1] == Fraction(n - 2, 2 * (n - 1))
    ok = all(mu[l] < mu[1] for l in range(2, n // 2 + 1))
    dec = all(mu[l + 1] <= mu[l] for l in range(n // 2))
    print(f"m={m} n={n}: mu_1 = {mu[1]} = {float(mu[1]):.6f}; mu_2 = {float(mu[2]):.6f}; mu_last = {float(mu[-1]):.6f}; "
          f"all mu_l < mu_1 for l>=2: {ok}; mu decreasing: {dec}; 2/(1-mu_1) = {Fraction(2) / (1 - mu[1])}")
