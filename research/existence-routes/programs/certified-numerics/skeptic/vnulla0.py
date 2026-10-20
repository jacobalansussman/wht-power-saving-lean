# v-ex2-certnum: independent exact degree-0 Nullstellensatz check (sum_k c_k f_k = 1 with constant c_k) and exact
# evaluation of recognised points, on ex2-certnum's exported systems.  Own code (Fractions).
import sys, json
from fractions import Fraction as Fr
from collections import Counter
def load(f):
    d = json.load(open(f)); eqs = []
    for e in d['eqs']:
        p = Counter()
        for mon, c in e:
            key = tuple(sorted(mon)); p[key] += complex(c[0], c[1])
        eqs.append({k: v for k, v in p.items() if v != 0})
    return d, eqs
def nulla0(eqs):
    mons = sorted({m for e in eqs for m in e}); assert all(v.imag == 0 for e in eqs for v in e.values())
    # unknowns c_k; equations: for each monomial, sum_k c_k coef_k(m) = [m == ()]
    rows = [[Fr(e.get(m, 0).real) for e in eqs] + [Fr(1 if m == () else 0)] for m in mons]
    if () not in mons: rows.append([Fr(0)] * len(eqs) + [Fr(1)])   # 0 = 1: inconsistent, no certificate
    M = rows; rk = 0; ncol = len(eqs)
    for c in range(ncol):
        p = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]; iv = 1 / M[rk][c]; M[rk] = [x * iv for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c]; M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        rk += 1
    consistent = all(M[i][-1] == 0 for i in range(rk, len(M)))
    if not consistent: return False, None
    c = [Fr(0)] * ncol
    piv = []
    for i in range(rk):
        j = next(jj for jj in range(ncol) if M[i][jj] != 0); c[j] = M[i][-1]
    # re-verify sum c_k f_k == 1 identically
    tot = Counter()
    for ck, e in zip(c, eqs):
        if ck:
            for m, v in e.items(): tot[m] += ck * Fr(v.real)
    tot = {m: v for m, v in tot.items() if v != 0}
    return tot == {(): Fr(1)}, sum(1 for x in c if x)
def evaluate(eqs, pt):
    bad = 0
    for e in eqs:
        s = (Fr(0), Fr(0))
        for m, v in e.items():
            t = (Fr(v.real), Fr(v.imag))
            for x in m:
                a = pt[x]; t = (t[0] * a[0] - t[1] * a[1], t[0] * a[1] + t[1] * a[0])
            s = (s[0] + t[0], s[1] + t[1])
        bad += s != (0, 0)
    return bad
if __name__ == '__main__':
    for f in sys.argv[1:]:
        d, eqs = load(f); ok, sup = nulla0(eqs)
        print(f.split('/')[-1], 'nvar', d['nvar'], 'neq', len(eqs), 'degree-0 certificate sum c_k f_k = 1:', ok, '(support %s)' % sup)
