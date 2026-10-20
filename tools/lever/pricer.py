#!/usr/bin/env python3
"""m-anatomy pricer: exact re-implementation of the hub's complex-side moment price and kappa assembly.

Re-implements (does not import) CrocSwap/integer-mult-bounds
  research/native-complex-moment552/engines/moment.py   (moment, certify; Mercator logs 224 terms, 24-term exp, 224-bit dyadic outward rounding)
  research/native-complex-moment552/engines/outer.py    (47 strict constraints, 7 margins)
  research/gauge44-frames-x/code/assemble.py            (outer_witness: s=min(c,(1-beta)b)-1e-27, 7-step bootstrap chain from 384599e-10,
                                                         kappa = floor-strict((1-eta) q/(1+q)) on the 1e-18 grid, q = chain[-1](1-2eta))
b = largest b on the 10^-24 grid with  sum_r (r n_r)/(m W) (m/r)^b < 1   (certified by exact rational enclosures at b and b+1e-24)
kappa = assembly(b, c) where c is the bit-side saving (default: PR589's 1733427039589844234749/10^24).

Usage:
  python3 -I pricer.py                                   # #799 control (profile, m=73, W=7168, deficit 2560)
  python3 -I pricer.py --profile '{"1":183457,...}' --m 73 --W 7168 [--deficit 2560] [--c 0.0017334...] [--bad]
  python3 -I pricer.py --json FILE                        # FILE has keys profile, m, W, optional deficit, c
  python3 -I pricer.py --local-delta '{"1":8,"2":-4,"3":4,"4":-4}'   # a per-call change of #799's local word (x72 calls);
        the change in residual width (= change in local rank mass, since copies and data climbs are fixed) is moved to or from
        the idle remainder at rank 48 (1.5 rank-48 children per width unit per call; a 24 left over is paid as one rank-24 child).
        This decomposition (72 x local {1:2548,2:338,3:430,4:194,5:56,8:8} + idle {1:1,3:1,44:1,48:2396}) reproduces #799 and,
        with #685's promoted local word {1:2572,2:326,...}, #685's kappa 0.001660339298410333 exactly (check_models.py).
  add --fast to skip the exact rational certification (Decimal bisection only; same 24-digit b in practice)
  add --nokappa to print b only
As a module: from pricer import price_b, kappa_from, b_fast, CONTROL
"""
import argparse, json, sys
from decimal import Decimal, localcontext, ROUND_FLOOR
from fractions import Fraction as F
from functools import lru_cache

DYAD = 1 << 224
def _floor(x): return F(x.numerator*DYAD//x.denominator, DYAD)
def _ceil(x): return -_floor(-x)

@lru_cache(None)
def log_bounds(x):
    """Exact rational enclosure of ln(x), x >= 1 (Mercator series, binary reduction; same scheme as the hub engine)."""
    x = F(x); assert x >= 1
    k = 0
    while x >= 2: x /= 2; k += 1
    j = int(x >= F(3, 2))
    if j: x /= F(3, 2)
    def mercator(y):
        z = y-1; assert 0 <= z <= F(1, 2)
        t = z; s = F(0)
        for n in range(1, 225):
            s += (t if n & 1 else -t)/n
            t *= z
        return s, s+t/225
    a, b = mercator(x); c, d = mercator(F(3, 2)); e, f = mercator(F(4, 3))
    return _floor(a+(k+j)*c+k*e), _ceil(b+(k+j)*d+k*f)

def exp_bounds(l, u):
    assert 0 <= l <= u < F(1, 10)
    def poly(x):
        s = t = F(1)
        for j in range(1, 25): t = t*x/j; s += t
        return s, t
    lo, _ = poly(l); hi, last = poly(u)
    return _floor(lo), _ceil(hi+last*u/25/(1-u/26))

def moment(H, m, W, a, bad=False):
    """Exact rational enclosure [lo, hi] of sum_r (r n_r)/(m W) (m/r)^a (+ the inherited bad-residue surcharge if bad)."""
    lo = hi = F(0)
    for r, n in H.items():
        assert 0 < r < m and n > 0
        l, u = log_bounds(F(m, r)); e, f = exp_bounds(a*l, a*u)
        p = F(r*n, m*W); lo += p*e; hi += p*f
    if bad:
        l, u = log_bounds(F(m)); e, f = exp_bounds(a*l, a*u)
        p = F(32*m*sum(H.values()), 10**16*W); lo += p*e; hi += p*f
    return _floor(lo), _ceil(hi)

def b_fast(H, m, W, bad=False, prec=70, iters=210, top='.003'):
    """Decimal bisection; returns the 24-digit grid candidate as Fraction (what the hub's certify proposes).
    The hub brackets b in [0, 0.003]; price_b retries with a wider bracket [0, 0.02] if the root is above 0.003."""
    with localcontext() as c:
        c.prec = prec
        terms = [(Decimal(r*n)/Decimal(m*W), (Decimal(m)/Decimal(r)).ln()) for r, n in H.items()]
        fb = Decimal(32*m*sum(H.values()))/Decimal(10**16*W) if bad else Decimal(0)
        lm = Decimal(m).ln(); lo = Decimal(0); hi = Decimal(top)
        for _ in range(iters):
            a = (lo+hi)/2
            if sum(p*(a*l).exp() for p, l in terms)+fb*(a*lm).exp() < 1: lo = a
            else: hi = a
        q = int(((lo+hi)/2*10**24).to_integral_value(rounding=ROUND_FLOOR))
    return F(q, 10**24)

def b_float(H, m, W):
    """Unrounded b as Decimal (prec 60), for exchange rates finer than the 1e-24 grid."""
    with localcontext() as c:
        c.prec = 60
        terms = [(Decimal(r*n)/Decimal(m*W), (Decimal(m)/Decimal(r)).ln()) for r, n in H.items()]
        lo = Decimal(0); hi = Decimal('.01')
        for _ in range(200):
            a = (lo+hi)/2
            if sum(p*(a*l).exp() for p, l in terms) < 1: lo = a
            else: hi = a
        return (lo+hi)/2

def price_b(H, m, W, bad=False, certify=True):
    H = {int(r): int(n) for r, n in H.items() if int(n)}
    assert m*m > 2*max(H)**2, 'halving constraint m^2 > 2 max(r)^2 fails'
    b = b_fast(H, m, W, bad)
    if b >= F(2999999, 10**9): b = b_fast(H, m, W, bad, iters=240, top='.02')
    if certify:
        lv = moment(H, m, W, b, bad); hv = moment(H, m, W, b+F(1, 10**24), bad)
        assert lv[1] < 1 < hv[0], ('not certified', b, lv, hv)
    return b

# ---- outer assembly (re-implementation of engines/outer.py and assemble.py:outer_witness) ----
def outer_ok(a, b, kappa, eta, beta, degree_gap=F(10**6)-F(51*20161, 25), lit_gap=1):
    tau, sigma = 1-a, 1-b
    q = a*(1-2*eta); lp = 1-q; lam = (tau+lp)/2
    c = q+eta/4; eps = (1-eta)/(1+q)
    minimum = eps*q; r = (minimum+1-eps)/2; delta = eta/8
    internal = tau+(1-beta)*max(sigma-tau, F(0)); leaf = sigma+beta*(1-sigma)
    margins = dict(balanced_prefix=1-eps, coordinate_movement=a, compact_phase_layer=minimum, bulk_exposure=a,
                   Gaussian_arithmetic=min(1-eps-delta, r-delta), scalar_work=1-eps-delta, dimension=eps)
    slacks = dict(bit_positive=a, complex_above_bit=b-a, complex_below_one_over32=F(1, 32)-b, beta_positive=beta,
        beta_below_one=1-beta, leaf_saving_above_bit=(1-beta)*b-a, q_positive=q, q_below_internal=1-internal-q,
        q_below_leaf=1-leaf-q, c_positive=c, c_below_one=1-c, q_below_reservations=c-q, lambda_above_tau=lam-tau,
        lambda_above_sigma=lam-sigma, lambda_above_internal=lam-internal, lambda_prime_above_lambda=lp-lam,
        compact_leaf=lp-leaf, compact_reservations=lp-(1-c), lambda_prime_below_one=q, epsilon_positive=eps,
        epsilon_below_one=1-eps, guard_width=1-eps, K_geometry=1-eps*(1+c), K_dominates_log=eps*c, record_suffix=1-eps,
        phase_local=1-eps-delta, phase_boundary=r-delta, gamma_sublinear=1-eps-r, cell_above_band=eps-(1-r)/2,
        prime_interval_packing=1-eps, alpha_positive=r, alpha_below_one=1-r, alpha_below_one_fourth=F(1, 4)-r,
        delta_positive=delta, delta_below_one_eighth=F(1, 8)-delta, short_record_fallback=eps-a,
        small_field_exposure=1-eps-minimum, artificial_boundary=8-eps+r-delta-minimum, literal_scalar_guard=F(lit_gap),
        row_product_gap=degree_gap)
    slacks.update({k+'_above_kappa': v-kappa for k, v in margins.items()})
    assert len(slacks) == 47
    bad = [k for k, v in slacks.items() if v <= 0]
    return not bad, bad

def kappa_from(b, c, eta=F(1, 10**24), beta=F(1, 10**24), check=True):
    """kappa from complex side b and bit side c, exactly as assemble.py:outer_witness."""
    s = min(c, (1-beta)*b)-F(1, 10**27); assert 0 < s < c and s < (1-beta)*b
    chain = [F(384599, 10**10)]
    for _ in range(7):
        old = chain[-1]; a = (1-s)*s+s*old
        assert old < a < s < 1-a
        chain.append(a)
    q = chain[-1]*(1-2*eta); margin = (1-eta)*q/(1+q); ticks = margin*10**18
    k = F((ticks.numerator-1)//ticks.denominator, 10**18)
    if check:
        ok, bad = outer_ok(chain[-1], b, k, eta, beta); assert ok, bad
        ok2, _ = outer_ok(chain[-1], b, k+F(1, 10**18), eta, beta); assert not ok2, 'next grid point should fail'
    return k, s, chain[-1]

C_PR589 = F(1733427039589844234749, 10**24)
CONTROL = dict(profile={1: 183457, 2: 24336, 3: 30961, 4: 13968, 5: 4032, 8: 576, 44: 1, 48: 2396}, m=73, W=7168, deficit=2560)
B799 = F(1664400278850482327913, 10**24); K799 = F(1661634653669567, 10**18)

def dec(x, p=40):
    with localcontext() as c:
        c.prec = p; return str(Decimal(x.numerator)/Decimal(x.denominator))

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--json'); ap.add_argument('--profile'); ap.add_argument('--m', type=int); ap.add_argument('--W', type=int)
    ap.add_argument('--deficit', type=int); ap.add_argument('--c', help='bit-side saving (decimal or p/q); default PR589')
    ap.add_argument('--bad', action='store_true', help='include the inherited bad-residue surcharge (conservative model)')
    ap.add_argument('--local-delta', help='per-call change of #799 local paid profile, JSON {rank: count}')
    ap.add_argument('--fast', action='store_true'); ap.add_argument('--nokappa', action='store_true')
    a = ap.parse_args()
    d = dict(CONTROL)
    if a.json: d.update(json.load(open(a.json)))
    if a.profile: d['profile'] = json.loads(a.profile)
    if a.local_delta:
        LOC = {1: 2548, 2: 338, 3: 430, 4: 194, 5: 56, 8: 8}; IDLE = {1: 1, 3: 1, 44: 1, 48: 2396}
        dl = {int(r): int(n) for r, n in json.loads(a.local_delta).items()}
        G = {}
        for r in set(LOC) | set(dl): G[r] = 72*(LOC.get(r, 0)+dl.get(r, 0))
        dw = sum(r*n for r, n in dl.items())          # width change per call
        moved = -72*dw                               # idle mass returning to (or taken from) the paid remainder
        q, rem = divmod(moved, 48)
        for r, n in IDLE.items(): G[r] = G.get(r, 0)+n
        G[48] += q
        if rem: G[rem] = G.get(rem, 0)+1
        assert all(n >= 0 for n in G.values()), G
        d['profile'] = {r: n for r, n in G.items() if n}
    if a.m: d['m'] = a.m
    if a.W: d['W'] = a.W
    if a.deficit is not None: d['deficit'] = a.deficit
    elif a.profile or a.m or a.W: d.pop('deficit', None)
    H = {int(r): int(n) for r, n in d['profile'].items() if int(n)}; m, W = int(d['m']), int(d['W'])
    mass = sum(r*n for r, n in H.items()); D = m*W-mass
    if d.get('deficit') is not None: assert D == int(d['deficit']), f'deficit mismatch: m*W-mass = {D}, given {d["deficit"]}'
    c = F(a.c) if a.c else (F(d['c']) if 'c' in d else C_PR589)
    b = price_b(H, m, W, a.bad, certify=not a.fast)
    out = dict(profile=dict(sorted(H.items())), m=m, W=W, rank_mass=mass, deficit=D, calls=sum(H.values()), maxrank=max(H),
               b=f'{b.numerator}/{b.denominator}', b_decimal=dec(b), certified=not a.fast, surcharge=a.bad)
    if not a.nokappa:
        k, s, chain = kappa_from(b, c)
        out.update(c=dec(c), binding='complex' if b*(1-F(1, 10**24)) < c else 'bit', kappa=f'{k.numerator}/{k.denominator}', kappa_decimal=dec(k, 20))
    if H == CONTROL['profile'] and m == 73 and W == 7168 and not a.bad and not a.c:
        out['control_b_matches_799'] = (b == B799)
        if not a.nokappa: out['control_kappa_matches_799'] = (k == K799)
    print(json.dumps(out, indent=1))

if __name__ == '__main__':
    main()
