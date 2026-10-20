"""crosscheck.py -- ex2-modp, independent re-check of every EXISTS point printed in out/*.txt against the RAW mod-p
exports (the .ms files, not the system rebuilt over Q): each exact coordinate is reduced mod p1 and mod p2 and every
raw polynomial must vanish.  Points over Q(t) are checked at every root of t's minimal polynomial mod p that exists
(the point must vanish at one of them).  Catches any error in zsys.to_Q / joint."""
import sys, os, re
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zsys import parse_ms, ev

def red(q, p):
    q = Fraction(q); return q.numerator * pow(q.denominator, -1, p) % p

def parse_val(s):
    s = s.strip()
    m = re.fullmatch(r'\(([-\d/]+)([+-])([\d/]+)\*t\)', s)
    if m:
        a = Fraction(m.group(1)); b = Fraction(m.group(3)) * (1 if m.group(2) == '+' else -1); return (a, b)
    return (Fraction(s), Fraction(0))

def check(prefixes, vals, field):
    ok_all = True
    for pr in (1, 2):
        ts = [None]
        if field:
            e, f = field
            p = parse_ms(prefixes[0] + f'.p{pr}.ms')[1]
            ts = [t for t in range(p) if (t * t + red(e, p) * t + red(f, p)) % p == 0] if p < 50 else None
            if ts is None:   # large p: find roots by Tonelli-Shanks on the discriminant
                D = (red(e, p) ** 2 - 4 * red(f, p)) % p
                if pow(D, (p - 1) // 2, p) != 1:
                    print(f'   p{pr}: t not in F_p (skipped)'); continue
                r = tonelli(D, p); inv2 = pow(2, -1, p)
                ts = [(-red(e, p) + r) * inv2 % p, (-red(e, p) - r) * inv2 % p]
        good_t = False
        for t in ts:
            allzero = True
            for l, pf in enumerate(prefixes):
                if open(pf + f'.p{pr}.ms').read().startswith('#'): continue
                names, p, polys = parse_ms(pf + f'.p{pr}.ms')
                x = []
                for v in names:
                    key = v if v[0] == 'a' or len(prefixes) == 1 else f'{v}_{l}'
                    a, b = vals.get(key, (Fraction(0), Fraction(0)))
                    x.append((red(a, p) + red(b, p) * (t or 0)) % p)
                if any(ev(P, x, p) for P in polys): allzero = False
            if allzero: good_t = True
        print(f'   p{pr}: raw export vanishes at the point: {good_t}'); ok_all &= good_t
    return ok_all

def tonelli(n, p):
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p); m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r

if __name__ == '__main__':
    pass
