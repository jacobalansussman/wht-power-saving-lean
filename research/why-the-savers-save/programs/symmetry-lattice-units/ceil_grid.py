"""Float version of the hub's skeleton ceiling X(s,h,k=1,x) (format-limits-corrected skel.py, theorem variant),
for exploring which (h, x) a new unit needs to beat the E8 units.  Certified values: hubcopy/skel.py."""
import math, sys
sys.path.insert(0, 'hubcopy')
from skel import profile
from fractions import Fraction as Fr

def Xf(s, h, x, k=1):
    prof, m, D = profile(s, h, k, Fr(x).limit_denominator(10**9) if not isinstance(x, Fr) else x)
    if D <= 0: return 0.0
    P = [(float(r), float(z)) for r, z in prof.items()]
    D = float(D); m = float(m)
    f = lambda q: sum(z*r*((m/r)**q-1) for r, z in P) - D
    lo, hi = 0.0, 0.5
    for _ in range(100):
        mid = (lo+hi)/2
        if f(mid) < 0: lo = mid
        else: hi = mid
    return lo

def best(h, x, S=range(3, 41, 2)):
    return max((Xf(s, h, x), s) for s in S)

def x_needed(h, target, s=None):
    """largest x with X(h, x) >= target (at best s, or at given s)."""
    g = (lambda x: best(h, x)[0]) if s is None else (lambda x: Xf(s, h, x))
    if g(0) < target: return None
    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = (lo+hi)/2
        if g(mid) >= target: lo = mid
        else: hi = mid
    return lo

if __name__ == '__main__':
    E8cx = best(9, Fr(64,112)); E8bit = best(8, Fr(63,111))
    print('E8 complex P112 (h=9, x=64/112): X* = %.4e at s=%d; s=9: %.4e' % (E8cx[0], E8cx[1], Xf(9, 9, Fr(64,112))))
    print('E8 bit P111 (h=8, x=63/111): X* = %.4e at s=%d; s=9: %.4e' % (E8bit[0], E8bit[1], Xf(9, 8, Fr(63,111))))
    for h in (8, 9, 10, 12, 16, 17, 18, 20, 24, 32):
        b0 = best(h, 0)
        xs = x_needed(h, E8cx[0]); x9 = x_needed(h, Xf(9, 9, Fr(64,112)), s=9)
        xb = x_needed(h, E8bit[0])
        print('h=%2d  X*(x=0)=%.4e (s=%d)  x needed to match E8cx X*: %s ; at s=9: %s ; to match E8bit X*: %s' % (
            h, b0[0], b0[1], 'impossible' if xs is None else '%.4f' % xs, 'impossible' if x9 is None else '%.4f' % x9,
            'impossible' if xb is None else '%.4f' % xb))
