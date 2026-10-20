# v-ex2-certnum: degree-d Nullstellensatz certificate with ALL-variable multiplier monomials of degree <= d, exact
# (sparse Fraction elimination), then the identity sum_k m_k f_k = 1 re-checked by expansion.  Own code.
import sys, itertools
from fractions import Fraction as Fr
from collections import Counter
from vnulla0 import load
def mons_upto(n, d):
    out = [()]
    for k in range(1, d + 1): out += list(itertools.combinations_with_replacement(range(n), k))
    return out
def cert(eqs, n, d):
    cols = []   # each unknown: (k, multiplier monomial) -> polynomial column {monomial: coef}
    for k, e in enumerate(eqs):
        for mu in mons_upto(n, d):
            col = Counter()
            for m, v in e.items(): col[tuple(sorted(m + mu))] += Fr(v.real)
            cols.append(((k, mu), {m: c for m, c in col.items() if c}))
    # solve sum_j a_j col_j = 1: eliminate on rows = monomials. Build row dicts
    rows = {}
    for j, (_, col) in enumerate(cols):
        for m, c in col.items(): rows.setdefault(m, {})[j] = c
    rhs = {m: Fr(1 if m == () else 0) for m in rows}
    if () not in rows: return None
    # Gaussian elimination on sparse rows (list), columns chosen greedily
    R = [(dict(r), rhs[m]) for m, r in rows.items()]
    pivots = []   # (col, rowdict, rhs)
    for r, b in R:
        r = dict(r)
        for pc, pr, pb in pivots:
            if pc in r:
                f = r[pc]
                for jj, v in pr.items():
                    nv = r.get(jj, 0) - f * v
                    if nv: r[jj] = nv
                    else: r.pop(jj, None)
                b = b - f * pb
        if not r:
            if b != 0: return None        # 0 = b != 0: inconsistent, no certificate at this degree
            continue
        pc = min(r, key=lambda jj: (len(cols[jj][1]), jj)); iv = 1 / r[pc]
        r = {jj: v * iv for jj, v in r.items()}; b = b * iv
        # keep pivots reduced against the new one (so back substitution is trivial)
        newp = []
        for qc, qr, qb in pivots:
            if pc in qr:
                f = qr[pc]
                for jj, v in r.items():
                    nv = qr.get(jj, 0) - f * v
                    if nv: qr[jj] = nv
                    else: qr.pop(jj, None)
                qb = qb - f * b
            newp.append((qc, qr, qb))
        pivots = newp + [(pc, r, b)]
    a = {pc: pb for pc, pr, pb in pivots}   # free unknowns = 0
    tot = Counter()
    for j, v in a.items():
        if v:
            for m, c in cols[j][1].items(): tot[m] += v * c
    tot = {m: c for m, c in tot.items() if c}
    return tot == {(): Fr(1)}, sum(1 for v in a.values() if v)
if __name__ == '__main__':
    d = int(sys.argv[1])
    for f in sys.argv[2:]:
        D, eqs = load(f); r = cert(eqs, D['nvar'], d)
        print(f.split('/')[-1], 'degree', d, 'certificate (all-variable multipliers):', r, flush=True)
