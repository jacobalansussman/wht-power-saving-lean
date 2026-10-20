"""vpoints.py -- v-ex2-modp: independent exact check of ex2-modp's EXISTS points on the RAW exports.
For every point: (1) exact evaluation over Q (or over Q(t), t^2 + c1 t + c0 = 0) on the system rebuilt by MY
reconstruction; (2) evaluation mod p1 and mod p2 directly on the raw files (rational points only).
Also re-reads the points from ex2-modp's text outputs (so the values checked are the ones claimed)."""
import os, re, sys
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlib import *
M = os.environ.get('MODP_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

def parse_point(s):
    pt = {}
    for m in re.finditer(r'(\w+)=(\([^)]*\)|-?\d+(?:/\d+)?)', s):
        pt[m.group(1)] = m.group(2)
    return pt

class QF:
    """a + b t in Q(t), t^2 = -c1 t - c0"""
    c0 = c1 = None
    def __init__(s, a, b=0): s.a, s.b = Fraction(a), Fraction(b)
    def __add__(s, o): o = o if isinstance(o, QF) else QF(o); return QF(s.a + o.a, s.b + o.b)
    __radd__ = __add__
    def __mul__(s, o):
        o = o if isinstance(o, QF) else QF(o)
        # (a + b t)(c + d t) = ac + (ad + bc) t + bd t^2,  t^2 = -c1 t - c0
        bd = s.b * o.b
        return QF(s.a * o.a - bd * QF.c0, s.a * o.b + s.b * o.a - bd * QF.c1)
    __rmul__ = __mul__
    def __pow__(s, e):
        r = QF(1)
        for _ in range(e): r = r * s
        return r
    def iszero(s): return s.a == 0 and s.b == 0
    def __repr__(s): return f'({s.a}+{s.b}t)'

def val(v, quad):
    if v.startswith('('):
        m = re.fullmatch(r'\((-?\d+(?:/\d+)?)([+-]\d+(?:/\d+)?)\*t\)', v)
        assert m, v
        return QF(Fraction(m.group(1)), Fraction(m.group(2)))
    return QF(Fraction(v)) if quad else Fraction(v)

def check(pre_list, pt, quad=None, joint=False, label=''):
    """pre_list: list of export prefixes (one per array).  joint: alphas shared, betas suffixed _l."""
    if quad: QF.c1, QF.c0 = quad
    allok = True; neq = 0; modok = True
    for l, pre in enumerate(pre_list):
        A = read_ms(pre + '.p1.ms'); B = read_ms(pre + '.p2.ms')
        names, F, _ = to_Q(A, B)
        def name_in_pt(v):
            return v if (not joint or v.startswith('a')) else f'{v}_{l}'
        x = []
        for v in names:
            w = name_in_pt(v)
            x.append(val(pt[w], quad) if w in pt else (QF(0) if quad else Fraction(0)))
        unknown = [k for k in pt if not any(name_in_pt(v) == k for v in names)]
        for P in F:
            neq += 1
            r = ev(P, x)
            z = r.iszero() if quad else (r == 0)
            if not z: allok = False
        if not quad:
            for (nn, p, PP) in (A, B):
                xm = [(xx.numerator * pow(xx.denominator, -1, p)) % p for xx in x]
                if any(ev_mod(P, xm, p) for P in PP): modok = False
    print(f"{'OK  ' if allok and modok else 'FAIL'} {label}: {neq} equations, exact={allok}" + ('' if quad else f", raw mod p1,p2={modok}")
          + (f", field t^2+{quad[0]}t+{quad[1]}" if quad else ''))
    return allok and modok

# 1. relaxation hits (rational)
txt = open(f'{M}/out/relax_hits.txt').read()
n = 0
for name, line in re.findall(r'^(\S+\.u): [^\n]*\nEXISTS \[exact rational point[^\]]*\] try \d+: ([^\n]*?)  \(all other', txt, re.M):
    n += 1
    check([f'{M}/data/ms/{name}'], parse_point(line), label=name)
print('rational relaxation points checked:', n)
# 2. the two quadratic-field hits
txt = open(f'{M}/out/padic_two_hits.txt').read()
for name, c1, c0, line in re.findall(r'^(\S+\.u): .*?EXISTS \[exact point over Q\(t\), t\^2 \+ \((-?\d+)\) t \+ \((-?\d+)\) = 0.*?\n\s+(.*?)\s+\(others 0\)', txt, re.M | re.S):
    check([f'{M}/data/ms/{name}'], parse_point(line), quad=(Fraction(c1), Fraction(c0)), label=name)
# 3. the plain unit, joint over both arrays
txt = open(f'{M}/out/cal_plain_unit.txt').read()
for chart, line in re.findall(r'^chart (\S+)\n.*?\nEXISTS \[exact rational point[^\]]*\] try \d+: (.*?)$', txt, re.M):
    pre = [f'{M}/data/dec/cal_plain4F_s1/full_c_{chart}.a{l}' for l in (0, 1)]
    check(pre, parse_point(line), joint=True, label=f'cal_plain4F_s1 chart {chart} (joint)')
# 4. special-position control
txt = open(f'{M}/out/cal_special_position.txt').read()
for name, line in re.findall(r'^(\S+\.u): .*?\nEXISTS \[exact rational point[^\]]*\] try \d+: (.*?)  \(all other', txt, re.M | re.S):
    pre = f'{M}/data/msctrl/{name}'
    if os.path.exists(pre + '.p1.ms'):
        check([pre], parse_point(line), label='msctrl ' + name)
    else:
        print('no export for', name, os.listdir(f'{M}/data/msctrl')[:5])
