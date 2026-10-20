# selftest.py -- sanity checks of the Krawczyk test and the exclusion test on systems with known answers
from certnum import *
from fractions import Fraction as Fr
ctx.prec = 200
# G1: x^2 - 2 = 0 (one unknown)
G1 = PolySys(1, [{(2,): (Fr(1), Fr(0)), (0,): (Fr(-2), Fr(0))}])
print('x^2-2 at 1.41421356, rho 1e-6 :', krawczyk(G1, newton_flint(G1, [1.41421356]), 1e-6)[0], '(expect True)')
print('x^2-2 at 1.3 (no Newton), rho 1e-3:', krawczyk(G1, [acb(1.3)], 1e-3)[0], '(expect False: no root in the box)')
print('x^2-2 at 0 (singular J), rho 0.1 :', krawczyk(G1, [acb(0.0)], 0.1)[0], '(expect False)')
# G2: x^2 + y^2 - 1, x - y : roots +-(1/sqrt2, 1/sqrt2)
G2 = PolySys(2, [{(2, 0): (Fr(1), Fr(0)), (0, 2): (Fr(1), Fr(0)), (0, 0): (Fr(-1), Fr(0))}, {(1, 0): (Fr(1), Fr(0)), (0, 1): (Fr(-1), Fr(0))}])
y = newton_flint(G2, [0.7, 0.7]); ok, K = krawczyk(G2, y, 1e-10)
print('circle/line near (0.707,0.707):', ok, K[0], '(expect True, 0.7071...)')
ok, _ = krawczyk(G2, [acb(0.5), acb(0.5)], 1e-3); print('circle/line at (0.5,0.5), rho 1e-3:', ok, '(expect False)')
# a big box containing both roots: uniqueness must fail
ok, _ = krawczyk(G2, [acb(0.0), acb(0.0)], 2.0); print('circle/line, box containing two roots:', ok, '(expect False)')
# exclusion: (x^2+y^2-1, x-y, x-0.9) has no zero; box around (0.707, 0.707) radius 0.05 must be excluded
G3 = PolySys(2, G2.polys + [{(1, 0): (Fr(1), Fr(0)), (0, 0): (Fr(-9, 10), Fr(0))}])
print('exclusion of an empty overdetermined box:', exclude_box(G3, [0.7071, 0.7071], 0.05), '(expect True)')
print('exclusion must FAIL on a box containing a true zero of G2:', exclude_box(G2, [0.7071, 0.7071], 0.05), '(expect False)')
# certify_point on G3' = G2 plus a redundant row (x^2 - y^2): exists; on G3: refused
G4 = PolySys(2, G2.polys + [{(2, 0): (Fr(1), Fr(0)), (0, 2): (Fr(-1), Fr(0))}])
r = certify_point_S(G4, np.array([0.70710678 + 0j, 0.70710678 + 0j]), fix_vars=[])
print('G2 + redundant row, Route S:', r.get('verdict'), '| rational closes:', r.get('simplest_rational_closes'))
r = certify_point_S(G3, np.array([0.70710678 + 0j, 0.70710678 + 0j]), fix_vars=[])
print('G2 + inconsistent row, Route S:', r.get('verdict'))
