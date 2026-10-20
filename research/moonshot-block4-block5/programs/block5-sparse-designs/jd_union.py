"""jd_union.py (w8-b5sparse): a FINITE NON-ORBIT design family for the 2-in-8: two sigma-cycles through one base stream.

Streams of input l: {u1^i beta_0 : 0 <= i < o1} U {u2^j beta_0 : 1 <= j < o2}, with u1, u2 two DIFFERENT accepted
elements (same acceptance as the orbit modes: every own-branch Lemma-Q row type has an ALIVE non-class producer inside
each cycle; swap or cycle type), beta_0 the base stream (z-producer of the other array's leftover in the *Z modes).
Each stream's own-branch rows are produced inside its own cycle (Lemma E1-sigma); the cross rows between the two cycles
are left to the exact peel.  This is not an orbit of any cyclic group (u1, u2 generate a large group); it is the
smallest frame-closed-per-cycle design that is not an orbit.  Decided by peel2 + modp + LM (jd2) and, where the d = 4
closure is binomial, alpha3's exact V4 test.
usage: as jdesign.py with mode orbitSNZ / orbitCYZ / orbitSN / orbitCY (the mode picks the sigma type)."""
import os, sys
import jd2
import alpha3
import jdesign as J
from endsolve import lift_P
from sp2 import mul as smul, order as sorder

U2 = {}
orig_sigma = J.orbit_sigma
orig_inj = J.orbit_injections


def sigma2(rng, cf, L):
    u1, o1, fix1, st1 = orig_sigma(rng, cf, L)
    if u1 is None:
        return u1, o1, fix1, st1
    for t in range(20):
        u2, o2, fix2, st2 = orig_sigma(rng, cf, L)
        if u2 is None:
            break
        if u2 != u1 and sorder(smul(u2, u1)) != 1 and all(u2 != x for x in powers(u1, o1)):
            U2['u'] = u2; U2['o'] = o2
            print("  union: second sigma order %d (first %d); alive %s" % (o2, o1, st2.get('alive')), flush=True)
            return u1, o1, fix1, st1
    print("  union: no second sigma", flush=True)
    return None, None, None, st1


def powers(u, o):
    out = []; x = u
    for _ in range(o):
        out.append(x); x = smul(u, x)
    return out


def inj2(rng, cf, L, l, u, n=0):
    js = orig_inj(rng, cf, L, l, u)
    if 'u' not in U2:
        return js
    b0 = js[0]['letter']; M0m = js[0]['mat']
    u2 = U2['u']; o2 = U2['o']
    Um = lift_P(J.EXACT[0], u2)
    b = smul(u2, b0); Mi = Um @ M0m
    for j in range(1, o2):
        js.append(dict(site='t1', l=l, where='union%d/%d' % (o2, j), Y=None, letter=b, tau_order=o2, lam=0,
                       mat=J.EXACT[0].snap(Mi)))
        b = smul(u2, b); Mi = Um @ Mi
    assert b == b0
    return js


J.orbit_sigma = sigma2
J.orbit_injections = inj2

if __name__ == "__main__":
    print("jd_union config %s argv %s" % (jd2.CFG, sys.argv[1:]), flush=True)
    J.main()
