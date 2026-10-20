#!/usr/bin/env python3
"""myfo.py -- v-ml-proofcx's own builder (independent of bcurve.py / nsdeg.py) for
 (1) the relaxation itself at a new prime, and
 (2) the first-order border system  alpha = A + eps D,  beta = B/eps + C :
     eps^-1 coefficient of f(alpha,beta) = 0  and  eps^0 coefficient = 0 (f includes the -t constant).
Input: msolve-format file (vars line, char line, polys); coefficients lifted to (-p/2, p/2).
Beta variables: names starting with 'b'.  Powers written x^k.
Usage: python3 -I myfo.py IN OUTPREFIX prime1,prime2,...
"""
import sys, re
from collections import defaultdict

def parse(path):
    lines = open(path).read().split('\n')
    vs = [v.strip() for v in lines[0].split(',') if v.strip()]
    p = int(lines[1].strip())
    body = ' '.join(lines[2:])
    polys = []
    for ptxt in body.split(','):
        ptxt = ptxt.strip()
        if not ptxt: continue
        P = defaultdict(int)
        for term in ptxt.replace(' ', '').split('+'):
            if not term: continue
            fac = term.split('*')
            c = 1; mon = defaultdict(int)
            for f in fac:
                if re.fullmatch(r'-?\d+', f):
                    c *= int(f)
                else:
                    if '^' in f:
                        v, e = f.split('^'); mon[v] += int(e)
                    else:
                        mon[f] += 1
            c %= p
            if c > p // 2: c -= p
            P[tuple(sorted(mon.items()))] += c
        polys.append({m: c for m, c in P.items() if c})
    return vs, p, polys

def fo_system(vs, polys):
    """returns list of polys (dict mon->int) in new vars A_x, D_x, B_y, C_y"""
    out = []
    for P in polys:
        em1 = defaultdict(int); e0 = defaultdict(int)
        for mon, c in P.items():
            al = [(v, e) for v, e in mon if not v.startswith('b')]
            be = [(v, e) for v, e in mon if v.startswith('b')]
            assert sum(e for _, e in be) <= 1, 'not linear in beta'
            # alpha product expanded to first order: prod A^e  +  eps * sum_v e*A_v^(e-1)*D_v * prod_{others}
            base = {('A_' + v): e for v, e in al}
            first = []   # list of monomial dicts
            for v, e in al:
                m = dict(base); m['A_' + v] -= 1
                if m['A_' + v] == 0: del m['A_' + v]
                m['D_' + v] = m.get('D_' + v, 0) + 1
                first.append((e, m))
            def key(m): return tuple(sorted((k, x) for k, x in m.items() if x))
            if be:
                bv = be[0][0]
                # (base + eps*first) * (B/eps + C): eps^-1: base*B ; eps^0: base*C + first*B
                m = dict(base); m['B_' + bv] = 1; em1[key(m)] += c
                m = dict(base); m['C_' + bv] = 1; e0[key(m)] += c
                for e, mf in first:
                    m = dict(mf); m['B_' + bv] = m.get('B_' + bv, 0) + 1; e0[key(m)] += c * e
            else:
                e0[key(base)] += c
        for E in (em1, e0):
            E = {m: c for m, c in E.items() if c}
            if E: out.append(E)
    return out

def write_ms(polys, p, path):
    used = sorted({v for P in polys for m in P for v, _ in m})
    body = []
    for P in polys:
        ts = []
        for m, c in sorted(P.items()):
            c %= p
            if not c: continue
            mon = '*'.join(v + (f'^{e}' if e > 1 else '') for v, e in m)
            ts.append(f'{c}*{mon}' if mon else f'{c}')
        if ts: body.append('+'.join(ts))
    open(path, 'w').write(','.join(used) + '\n' + str(p) + '\n' + ',\n'.join(body) + '\n')
    return len(used), len(body)

if __name__ == '__main__':
    inp, outp, primes = sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(',')]
    vs, p0, polys = parse(inp)
    fo = fo_system(vs, polys)
    for p in primes:
        n1 = write_ms(polys, p, f'{outp}.relax.p{p}.ms')
        n2 = write_ms(fo, p, f'{outp}.fo.p{p}.ms')
        print(inp.split('/')[-1], 'p', p, 'relax vars/eqs', n1, 'first-order vars/eqs', n2, file=sys.stderr)
