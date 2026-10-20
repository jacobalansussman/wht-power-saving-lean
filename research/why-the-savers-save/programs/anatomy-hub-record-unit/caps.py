"""an-hub: closed-form caps of the one-unpaid-coordinate architecture, and two lineage checks. Stdlib only.
a(m, u, cw, reff) = -ln(1 - (u - cw)/m) / ln(m / reff): u unpaid coordinates per family, cw = centre charge per
family (c/w), reff = effective child rank.  Exact form of the first-order law used by the hub (a ~ D/E)."""
import math
def a(m, u, cw, reff):
    x = (u - cw) / m
    return -math.log(1 - x) / math.log(m / reff) if x > 0 else 0.0
print('m   | rank-1, cw=0 | reff=4, cw=0.643 | reff=4, cw=0 | reff=m/sqrt2, cw=0.643 | reff=m/sqrt2, cw=0')
for m in (4, 5, 6, 8, 12, 16, 21, 28, 45, 73, 145):
    R = m / math.sqrt(2)
    print(f'{m:4d} | {a(m,1,0,1):.4f} | {a(m,1,0.643,min(4,R)):.4f} | {a(m,1,0,min(4,R)):.4f} | {a(m,1,0.643,R):.4f} | {a(m,1,0,R):.4f}')
# T74 check of the first-order law
m, w = 73, 896 / 9
print('\nT74 per call: D = w - 64 =', w - 64, ' x = D/(m w) =', (w - 64) / (m * w))
print('T74 a from closed form with reff=3.9448:', a(73, 1, 64 / w, 3.9448))
# PR1625 helper-donor kernel: H1 +1, H2 -2, H3 +1 per call, 666 calls; predicted change of a (first order)
mass, reff, calls = 4816512, 3.9448, 666
E = mass * math.log(m / reff)
dE = -calls * (0 - 2 * 2 * math.log(2) + 3 * math.log(3))
print('kernel 566: mass change per call', 1 - 4 + 3, ' dE/E =', dE / E, ' predicted da =', -0.0016796718 * dE / E,
      ' (hub kappa 1.676855229e-3 -> 1.676896669e-3, i.e. d kappa =', 1.676896669e-3 - 1.676855229e-3, ')')
# outer conversion: kappa = a/(1+a) before PR1673, kappa ~ a after
b = 1.67975495612554e-3
print('a/(1+a) at b =', b / (1 + b), ' (hub pivot-554 package kappa 1.676938111e-3); gain of PR1673 outer =', b / (b / (1 + b)) - 1)
