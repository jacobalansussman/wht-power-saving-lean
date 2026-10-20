"""sing2ms.py (w8-b5sparse, resumed): convert a v4gb.py ball script (ideal I, char p) to msolve input, for msolve's
solver: output [-1] = no solution over the algebraic closure of F_p (EXCLUDED), [1, n, -1, []] = positive dimensional.
Coefficients are reduced mod p; variable a(k) -> ak, x1(c) -> x1c<c>, b(k) -> bk.  No powers occur (products of
distinct variables only), so msolve 0.6.5's x*x bug cannot bite.  usage: sing2ms.py in.sing out.ms"""
import re, sys
from fractions import Fraction

src = open(sys.argv[1]).read()
p = int(re.search(r'ring R = (\d+),', src).group(1))
assert p > 0
vdecl = re.search(r'ring R = \d+, \((.*)\), dp;', src).group(1)
body = re.search(r'ideal I =\n(.*?);\n', src, re.S).group(1)


def ren(m):
    name, idx = m.group(1), m.group(2)
    return name + ('c' if name[-1].isdigit() else '') + idx


def frac(m):
    f = Fraction(m.group(1))
    return '(%d)' % (f.numerator * pow(f.denominator, p - 2, p) % p)


body = re.sub(r'\((-?\d+/\d+)\)', frac, body)
body = re.sub(r'\((-?\d+)\)', lambda m: '(%d)' % (int(m.group(1)) % p), body)
body = re.sub(r'([a-z]\d?)\((\d+)\)', ren, body)
body = body.replace('(', '').replace(')', '')
polys = [q.strip() for q in body.split(',\n')]
polys = [q for q in polys if q and q != 't']
vars_ = []
for part in [v.strip() for v in vdecl.split(',')]:
    m = re.match(r'([a-z]\d?)\((\d+)\.\.(\d+)\)', part)
    if m:
        vars_ += [m.group(1) + ('c' if m.group(1)[-1].isdigit() else '') + str(i) for i in range(int(m.group(2)), int(m.group(3)) + 1)]
    elif part != 't':
        vars_.append(part)
used = set(re.findall(r'[a-z]\d?c?\d+|\bt\b', ','.join(polys)))
vars_ = [v for v in vars_ if v in used]
open(sys.argv[2], 'w').write(','.join(vars_) + '\n' + str(p) + '\n' + ',\n'.join(polys) + '\n')
print('wrote %s: %d variables, %d polynomials, char %d' % (sys.argv[2], len(vars_), len(polys), p))
