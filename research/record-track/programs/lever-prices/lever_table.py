#!/usr/bin/env python3
"""The recursion lever on every circuit of the record track, priced with certified brackets.

For each circuit (a paid profile, the number W of arrays, the idle rank-48 blocks and how many of them are heads)
this prints the saving b = 1 - tau, where tau is the root of the balance

    sum over paid blocks of rank^tau  =  sum over arrays of (job size)^tau,

for four rules:
  uniform  - the community's own recursion (every job has 73 axes);
  half cut - each head cut by 24 of its 48 axes: the head array's job has 49 axes, the head keeps 24 paid axes,
             the other rank-48 blocks stay as they are (largest paid rank 48 < smallest job 49: separated);
  full cut - each head removed: the head array's job has 25 axes, and every other rank-48 block is split into two
             blocks of rank 24 (largest paid rank 24 < 25: separated);
  not separated - heads removed (job 25) and the other rank-48 blocks kept (48 >= 25: no theorem covers it).

Each root is bracketed with exact rational interval arithmetic (iprice.py in this folder: logarithms by atanh
series, exp by Taylor series, both with tail bounds), to a width of 1e-16. It also finds the full cut's break-even
(the least number of heads for which the full cut beats the uniform figure) by bisection on certified comparisons.

The paid profiles are the community's published global profiles (PROFILE.json of each pull request, quoted here
as numbers); the head counts are the counts of the head test (README, section 4.5).

Run: python3 -I -B lever_table.py [out.json]      (about a minute, under 100 MB, standard library only)
"""
import sys, os, json
from fractions import Fraction as Q
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from iprice import root, sep, add, F

CIRCUITS = [
    # name, W, paid profile (rank: count), [(label, heads), ...]
    ('T74 (#1839)', 66304, {1: 1560432, 2: 336336, 3: 278724, 4: 84912, 5: 54612, 8: 5328, 48: 22748},
     [('published fill', 11010)]),
    ('#2979', 66304, {1: 1645680, 2: 470196, 3: 159468, 4: 71928, 5: 13716, 8: 5328, 48: 28190},
     [('published fill', 16086), ('every head bin empty (design A)', 16576)]),
    ('#3945', 322560, {1: 8041680, 2: 2306880, 3: 729648, 4: 353160, 5: 69120, 8: 25920, 48: 137952},
     [('published fill', 57312), ('every head bin empty (design B)', 80640)]),
    ('#3815', 322560, {1: 7989840, 2: 2349000, 3: 749088, 4: 330480, 5: 69120, 8: 25920, 48: 137952},
     [('published fill', 57312)]),
]

# Values stated elsewhere in the record track; the program checks that its brackets contain them (to 1e-10).
EXPECT = {
    ('T74 (#1839)', 11010, 'uniform'): '1.6797549561e-3', ('T74 (#1839)', 11010, 'half'): '1.6968621073e-3',
    ('T74 (#1839)', 11010, 'full'): '1.6948357905e-3', ('T74 (#1839)', 11010, 'unsep'): '1.7434773229e-3',
    ('#2979', 16086, 'uniform'): '1.7290339826e-3', ('#2979', 16086, 'half'): '1.7556473153e-3',
    ('#2979', 16086, 'full'): '1.7744822149e-3', ('#2979', 16576, 'half'): '1.7564708007e-3',
    ('#3945', 57312, 'uniform'): '1.7306632365e-3', ('#3945', 57312, 'half'): '1.7501109246e-3',
    ('#3945', 57312, 'full'): '1.7307679642e-3', ('#3945', 57312, 'unsep'): '1.8033130e-3',
    ('#3945', 80640, 'half'): '1.7581520168e-3', ('#3945', 80640, 'full'): '1.7806935406e-3',
    ('#3815', 57312, 'half'): '1.7500242253e-3',
}
BREAKEVEN = {'T74 (#1839)': 9443, '#2979': 11702, '#3945': 57262}

def rules(prof, W, h):
    I = prof[48]; base = dict(prof); del base[48]
    return {
        'uniform': (add(base, 48, I), {73: W}),
        'half': (add(add(base, 24, h), 48, I - h), {73: W - h, 49: h}),
        'full': (add(base, 24, 2 * (I - h)) if I - h else dict(base), {73: W - h, 25: h}),
        'unsep': (add(base, 48, I - h) if I - h else dict(base), {73: W - h, 25: h}),
    }

def bracket(ch, tg, width=Q(1, 10**16)):
    ch = {r: n for r, n in ch.items() if n}
    return root(ch, tg, width=width)

def beats(ch, tg, lo_u, hi_u):
    """True if the rule's root is certainly above the bracket [lo_u, hi_u], False if certainly below."""
    ch = {r: n for r, n in ch.items() if n}
    # F(b) < 0 means b is below the root
    if F(ch, tg, hi_u).hi < 0: return True
    if F(ch, tg, lo_u).lo > 0: return False
    raise RuntimeError('undecided comparison')

def main():
    out = []; bad = []
    for name, W, prof, fills in CIRCUITS:
        u_lo, u_hi = bracket(*rules(prof, W, 0)['uniform'])
        for label, h in fills:
            row = dict(circuit=name, fill=label, W=W, idle48=prof[48], heads=h)
            for key, (ch, tg) in rules(prof, W, h).items():
                lo, hi = bracket(ch, tg)
                row[key] = dict(lo=str(lo.numerator * 10**16 // lo.denominator), hi=str(-((-hi.numerator * 10**16) // hi.denominator)),
                                b='%.10e' % float(lo), separated=sep({r: n for r, n in ch.items() if n}, tg))
                e = EXPECT.get((name, h, key))
                if e is not None:
                    ev = Q(e)
                    ok = lo - Q(1, 10**10) <= ev <= hi + Q(1, 10**10)
                    row[key]['expected'] = e; row[key]['matches_expected'] = ok
                    if not ok: bad.append((name, h, key, e, float(lo)))
            out.append(row)
        # break-even of the full cut against the uniform figure, at fixed idle rank-48 count
        I = prof[48]; lo_t, hi_t = 0, I
        while hi_t - lo_t > 1:
            mid = (lo_t + hi_t) // 2
            ch, tg = rules(prof, W, mid)['full']
            if beats(ch, tg, u_lo, u_hi): hi_t = mid
            else: lo_t = mid
        be = dict(circuit=name, full_cut_break_even_heads=hi_t, of_idle48=I, expected=BREAKEVEN.get(name),
                  matches_expected=(BREAKEVEN.get(name) in (None, hi_t)))
        if not be['matches_expected']: bad.append((name, 'break-even', hi_t))
        out.append(be)
    print(f"{'circuit':12s} {'fill':34s} {'heads':>7s} {'uniform':>16s} {'half cut':>16s} {'full cut':>16s} {'not separated':>16s}")
    for r in out:
        if 'fill' in r:
            print(f"{r['circuit']:12s} {r['fill']:34s} {r['heads']:7d} {r['uniform']['b']:>16s} {r['half']['b']:>16s} "
                  f"{r['full']['b']:>16s} {r['unsep']['b']:>16s}")
        else:
            print(f"{r['circuit']:12s} full-cut break-even: {r['full_cut_break_even_heads']} heads of {r['of_idle48']} idle rank-48 blocks")
    print('every bracket has width at most 1e-16; separated: half and full cut yes, last column no' if all(
        r['half']['separated'] and r['full']['separated'] and not r['unsep']['separated'] for r in out if 'fill' in r) else 'SEPARATION CHECK FAILED')
    print('all stated values reproduced' if not bad else f'MISMATCH: {bad}')
    if len(sys.argv) > 1:
        json.dump(out, open(sys.argv[1], 'w'), indent=1)
    return 0 if not bad else 1

if __name__ == '__main__':
    sys.exit(main())
