#!/usr/bin/env python3
"""an2-word: (a) the quotient identity, checked exactly on random frames; (b) the replica bound; (c) its price.
Own code, stdlib only.
(a) For a frame change f -> g (a climb f < g or a drop g < f, subspaces of F_2^m) and a fixed nonzero vector s, the cost
    |dim g - dim f| exceeds the cost of the projected change (f+<s>)/<s> -> (g+<s>)/<s> by exactly 1 when the change
    crosses s (s in exactly one of f, g) and by 0 otherwise.  Equal frames stay equal under the projection.
(b) Per replica of the template (9 calls: 5 GF, 4 GB; W = 896 data arrays + B banks):
    D_two(m=74) = D_one(quotient, m=73) + W - (crossings),  crossings >= 9 * 112 + B   (the lemma, per call)
    => D_two <= D_one(quotient) - 112.  With the one-sided D = 320: D_two <= 208.
(c) Price of D_two = 208 per replica at m = 74 with the U8 word's own profile (exact fill, every idle remainder at rank
    48, the 112 crossings of every call as rank-1 pieces): an optimistic upper bound for any two-sided U8 word in this
    template.  Root of the Lean fold hypothesis (Work/BlockAccounting/Moment.lean; as in v-an-build/calc/vcalib.py),
    my own Decimal bisection.
usage: python3 -I quotient.py OUT.json"""
import json, os, random, sys
from decimal import Decimal as D, getcontext
getcontext().prec = 60

# ---------- (a) ----------
def reduce(rows):
    """row-reduced basis (ints as bit vectors) of the span of rows."""
    basis = []
    for r in rows:
        for b in basis:
            if r ^ b < r: r ^= b
        if r:
            basis = [b ^ r if b ^ r < b else b for b in basis] + [r]
            basis.sort(reverse=True)
    return basis
def span_dim(rows): return len(reduce(rows))
def contains(rows, v): return span_dim(rows + [v]) == span_dim(rows)
def proj_dim(rows, s): return span_dim(rows + [s]) - 1          # dim of (f + <s>)/<s>
rnd = random.Random(7)
m = 7
bad = 0; crossings = 0; trials = 20000
for _ in range(trials):
    s = rnd.randrange(1, 1 << m)
    f = [rnd.randrange(0, 1 << m) for _ in range(rnd.randrange(0, m))]
    g = f + [rnd.randrange(0, 1 << m) for _ in range(rnd.randrange(1, 4))]     # g contains f: a climb f -> g
    for a, b in ((f, g), (g, f)):                                              # climb and drop
        cost = abs(span_dim(b) - span_dim(a))
        pcost = abs(proj_dim(b, s) - proj_dim(a, s))
        cross = contains(a, s) != contains(b, s)
        crossings += cross
        if cost - pcost != (1 if cross else 0): bad += 1
print('(a) random climbs and drops in F_2^%d: %d checks, %d crossings, identity cost - projected cost = [crossing] fails %d times'
      % (m, 2 * trials, crossings, bad))
# equal frames stay equal: trivially (the projection is a function of the frame)

# ---------- (b) ----------
Wd, calls, GF, GB, ch, cst = 896, 9, 5, 4, 112, 64
D_one = Wd - calls * cst                       # 320: the identity of an-ours / v-an-build
gain = Wd                                      # one more free line per data array
cross_min = calls * ch                         # per call: GF needs ch climbs, GB ch drops (the lemma); banks add B and cancel
D_two_max = D_one + gain - cross_min
print('(b) one-sided D per replica %d; two-sided gain +%d; crossings >= %d (+B, cancelled by the B banks); D_two <= %d'
      % (D_one, gain, cross_min, D_two_max))
print('    an-build ceiling (GB only, 4 x 112):', D_one + gain - GB * ch, ' (assumed the 5 GF calls free: they are not)')

# ---------- (c) ----------
def lhs(H, W, m, a, s=40):
    e = 1 - a; lm = D(m)
    tot = sum(D(n) * ((D(r) / lm).ln() * e).exp() for r, n in H.items())
    fill = ((D(m - 1) / lm).ln() * e).exp() + ((D(1) / lm).ln() * e).exp()
    two = D(2) ** s
    return (two - 1) * tot + D(W) * fill - two * D(W)
def root(H, W, m):
    lo, hi = D(0), D('0.05')
    assert lhs(H, W, m, lo) < 0 and lhs(H, W, m, hi) > 0
    for _ in range(110):
        mid = (lo + hi) / 2
        if lhs(H, W, m, mid) < 0: lo = mid
        else: hi = mid
    return float(lo)
BC = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs') + '/')
wj = json.load(open(BC + 'words.json'))
U8E1 = {int(k): n for k, n in wj['869/870/890_U8_exit_frames(E1)']['profile'].items()}
T = 8; W = T * Wd
def fill(mm, Drep, extra1):
    H = {r: 72 * n for r, n in U8E1.items()}
    H[1] = H.get(1, 0) + extra1
    rem = mm * W - T * Drep - sum(r * n for r, n in H.items())
    q, rr = divmod(rem, 48)
    H[48] = H.get(48, 0) + q
    if rr: H[rr] = H.get(rr, 0) + 1
    assert mm * W - sum(r * n for r, n in H.items()) == T * Drep
    return H, q, rr
res = {}
for label, mm, Drep, extra in (('one-sided U8, exact fill (calibration, 1.6647327e-3)', 73, 320, 0),
                               ('two-sided, the lemma floor: D 208, 9 x 112 rank-1 crossings', 74, 208, 72 * ch),
                               ('two-sided, an-build fixed-order tie: D 320 (GB crossings only)', 74, 320, 32 * 224),
                               ('two-sided, an-build ceiling D 768 (GF free: impossible)', 74, 768, 32 * ch)):
    H, q, rr = fill(mm, Drep, extra)
    a = root(H, W, mm)
    res[label] = dict(m=mm, D_per_replica=Drep, rank1_crossing_pieces_per_word=extra, rank48_fill=q, remainder=rr, root=a)
    print('(c) %-66s m %d  D/replica %4d  root %.7e' % (label, mm, Drep, a))
cal = res['one-sided U8, exact fill (calibration, 1.6647327e-3)']['root']
assert abs(cal - 1.6647327e-3) < 1e-10, cal
json.dump(dict(a_bad=bad, a_checks=2 * trials, a_crossings=crossings, D_one=D_one, D_two_max=D_two_max,
               prices=res), open(sys.argv[1], 'w'), indent=1)
