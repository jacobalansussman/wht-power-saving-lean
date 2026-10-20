#!/usr/bin/env python3
"""d-hpe-multigate STEP 1b: the low-end admission cones of every U8 helper, independent of the cleanup gates.

The U8 audit (vendor/audit_local_pr870.py, Kothe/icekylinx PR870 independent checker) admits a helper q ending at
frame U < FULL only if
  (inverse-row cone)  every register k in the support of row q of the inverse prefix map has endB[k] <= U, and
  (dirt-column cone)  every register k != q that holds q's initial content after the prefix has endB[k] <= U,
where prefix = precondition + compensation + A + cut + B and endB is the frame state at the end of part B.
These conditions do not mention the cleanup at all, so the least admissible end frame of q is
  Umin(q) = join of endB[q] and endB[k] over both cones,
for ANY cleanup or tail.  A helper can end below FULL only if dim Umin(q) <= 8.
Usage: python3 -I cones.py FULL-WORD.json.gz out.json"""
import sys, gzip, json
from fractions import Fraction as F
from collections import Counter, defaultdict
def basis(vs):
    piv = {}
    for a in vs:
        while a:
            bit = a.bit_length()-1
            if bit in piv: a ^= piv[bit]
            else: piv[bit] = a; break
    for k in sorted(piv):
        for j in piv:
            if j > k and piv[j] >> k & 1: piv[j] ^= piv[k]
    return tuple(piv[k] for k in sorted(piv, reverse=True))
def contained(a, b): return basis(list(a)+list(b)) == b
def registers(e): return e[1:3] if e[0] == 'add' else e[1:2]
def apply(e, rows):
    if e[0] == 'add':
        t, s = e[1:3]; q = F(e[3], e[4]); dst = rows[t]
        for k, x in list(rows[s].items()):
            a = dst.get(k, F(0))+q*x
            if a: dst[k] = a
            else: dst.pop(k, None)
    else:
        t = e[1]; q = F(e[2], e[3]); rows[t] = {k: q*x for k, x in rows[t].items()}
def invert(e):
    r = list(e)
    if r[0] == 'add': r[3] = -r[3]
    else: r[2], r[3] = r[3], r[2]
    return r
w = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v+R
ident = lambda: [{k: F(1)} for k in range(n)]
main = w['A']+w['cut']+w['B']
raw = ident()
for e in main: apply(e, raw)
dirt = [{} for _ in range(n)]
for y in range(v, 2*v): dirt[y] = {k: a for k, a in raw[y].items() if k >= 2*v}
for e in w['precondition']: apply(e, dirt)
compensation = [['add', y, k, -a.numerator, a.denominator, []] for y in range(v, 2*v) for k, a in sorted(dirt[y].items())]
prefix = w['precondition']+compensation+main
full = basis([1 << i for i in range(h)])
start = [(p,) for p in w['ports']]+[()]*(v+R)
for q, U in w.get('start_extra', {}).items(): start[int(q)] = basis(U)
cur = start.copy()
for section in ('A', 'B'):
    for e in w[section]:
        for q in registers(e): cur[q] = basis(list(cur[q])+list(e[-1]))
endB = cur.copy()
p = ident(); pinv = ident()
for e in prefix: apply(e, p)
for e in reversed(prefix): apply(invert(e), pinv)
# column supports: col[q] = {k : p[k][q] != 0}
col = defaultdict(set)
for k in range(n):
    for q in p[k]: col[q].add(k)
lowends = {int(q): basis(U) for q, U in w['helper_final'].items()}
entering = set(map(int, w['start_extra']))
out = {}; summ = Counter(); rowonly = Counter(); colonly = Counter()
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
for q in range(2*v, n):
    rowk = set(pinv[q]); colk = col[q] - {q}
    Ur = basis([x for k in rowk | {q} for x in endB[k]])
    Uc = basis([x for k in colk | {q} for x in endB[k]])
    U = basis(list(Ur)+list(Uc))
    end = lowends.get(q, full)
    out[q] = dict(entering=q in entering, start_dim=len(start[q]), endB_dim=len(endB[q]), end_dim=len(end),
                  row_support=len(rowk), row_classes=dict(Counter(cls(k) for k in rowk)),
                  col_support=len(colk), col_classes=dict(Counter(cls(k) for k in colk)),
                  Umin_row_dim=len(Ur), Umin_col_dim=len(Uc), Umin_dim=len(U), Umin=list(U),
                  low_end_ok=contained(U, end))
    summ[(len(end), len(U))] += 1; rowonly[(len(end), len(Ur))] += 1; colonly[(len(end), len(Uc))] += 1
assert all(o['low_end_ok'] for o in out.values()), 'control: every existing end must contain its cone'
print('(current end dim, least admissible end dim Umin) : helpers', dict(sorted(summ.items())))
print('row cone only', dict(sorted(rowonly.items())), ' column cone only', dict(sorted(colonly.items())))
full_ended = [q for q in out if out[q]['end_dim'] == 9]
print('FULL-ended helpers', len(full_ended), 'with Umin dim <= 8:', sum(1 for q in full_ended if out[q]['Umin_dim'] <= 8))
print('FULL-ended, by (start_dim, endB_dim, Umin_row_dim, Umin_col_dim):')
for k, x in sorted(Counter((out[q]['start_dim'], out[q]['endB_dim'], out[q]['Umin_row_dim'], out[q]['Umin_col_dim']) for q in full_ended).items()): print('  ', k, x)
json.dump(dict(summary={f'{a}->{b}': x for (a, b), x in sorted(summ.items())}, helpers=out), open(sys.argv[2], 'w'), indent=0)
