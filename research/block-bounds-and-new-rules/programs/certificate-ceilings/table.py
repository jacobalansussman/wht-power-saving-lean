# b56-ceilings: the table of family ceilings for block sizes m = 2..10 (closed forms; exact fractions).
#   hook (V (x) Lambda^k, N = {x in T}, every f = k+1):  cap = (n/2) C(n-2,k-1)  [exact],
#        dem <= min( k C(n,k) - C(n,k-1),  k C(n,k+1) )   [psi-subspace bound, module bound],
#        so dem/cap <= min( 2(n-1)/(n-k+1), 2(n-1)/(k+1) ),  family max at k = n/2:  H = 4(n-1)/(n+2).
#   all Venn sets (any number of exterior-power factors, any S_n-invariant coordinate set):  V = 4(n-1)/n  (spectral),
#        and the trivial dem <= m cap.
#   implied: a family with ceiling C can never prove more than C moves per array, so it can never exclude a saving
#   below s* = ln(m/C)/ln m (block m).  Gap = m - C.
from fractions import Fraction as F
from math import comb, log
print(" m    n   plain | hook ceiling H=4(n-1)/(n+2)   s*(H)  gap | Venn ceiling min(m,4(n-1)/n)  s*(V)  gap")
for m in range(2, 11):
    n = 2 ** m
    # hook per-k ceilings, check the closed form of the max
    best = max(F(min(k * comb(n, k) - comb(n, k - 1), k * comb(n, k + 1)), (n // 2) * comb(n - 2, k - 1)) for k in range(1, n - 1))
    H = F(4 * (n - 1), n + 2)
    assert best == H, (m, best, H)
    V = min(F(m), F(4 * (n - 1), n))
    sH = log(m / H) / log(m) if H < m else 0.0
    sV = log(m / V) / log(m) if V < m else 0.0
    print(f"{m:2d} {n:4d} {m:5d}   | {str(H):>9s} = {float(H):.4f}        {sH:.4f} {float(m - H):5.2f} | {str(V):>7s} = {float(V):.4f}            {sV:.4f} {float(m - V):5.2f}")
