"""d-mcopy-code: small exact library over the U8 FULL-WORD (hub format, frames as F2 bases).
Program order and climb rule are the ones of the hub's pr870 independent audit (audit_local.py):
program = precondition + compensation + A + cut + B + tail(early_cleanup) + cleanup; frames climb in A, B, tail."""
import json, gzip
from fractions import Fraction as F
from collections import Counter, defaultdict

def load(p):
    return json.loads(gzip.decompress(open(p, 'rb').read()))

def basis(vs):
    piv = {}
    for a in vs:
        while a:
            bit = a.bit_length() - 1
            if bit in piv: a ^= piv[bit]
            else: piv[bit] = a; break
    for k in sorted(piv):
        for j in piv:
            if j > k and piv[j] >> k & 1: piv[j] ^= piv[k]
    return tuple(piv[k] for k in sorted(piv, reverse=True))

def contained(a, b): return basis(list(a) + list(b)) == tuple(b)
def pc(x): return bin(x).count('1')

def perp(p, h=9):
    return basis([x for x in range(1, 1 << h) if pc(x & p) % 2 == 0])

def compensation(w):
    v, R = w['report']['v'], w['report']['R']; n = 2*v + R
    rows = [{k: F(1)} for k in range(n)]
    for e in w['A'] + w['cut'] + w['B']: apply(e, rows)
    dirt = [{} for _ in range(n)]
    for y in range(v, 2*v): dirt[y] = {k: a for k, a in rows[y].items() if k >= 2*v}
    for e in w['precondition']: apply(e, dirt)
    return [['add', y, k, -a.numerator, a.denominator, []] for y in range(v, 2*v) for k, a in sorted(dirt[y].items())]

def apply(e, rows):
    if e[0] == 'add':
        t, s = e[1:3]; q = F(e[3], e[4])
        dst = rows[t]
        for k, x in list(rows[s].items()):
            a = dst.get(k, F(0)) + q*x
            if a: dst[k] = a
            else: dst.pop(k, None)
    else:
        t = e[1]; q = F(e[2], e[3]); rows[t] = {k: q*x for k, x in rows[t].items()}

def program(w, comp=None):
    """list of (phase, index, gate) in audit program order"""
    if comp is None: comp = compensation(w)
    early = w.get('early_cleanup', []); eid = {i for i, U in early}
    tail = [w['cleanup'][i][:-1] + [U] for i, U in early]
    clean = [e for i, e in enumerate(w['cleanup']) if i not in eid]
    out = []
    for ph, L in [('precondition', w['precondition']), ('compensation', comp), ('A', w['A']), ('cut', w['cut']),
                  ('B', w['B']), ('tail', tail), ('cleanup', clean)]:
        for i, e in enumerate(L): out.append((ph, i, e))
    return out

def frames(w):
    """audit climb rule: start, centres, A, (centre levels), B, tail, final.  returns start, end, chains, paid Counter"""
    v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v + R
    full = basis([1 << i for i in range(h)])
    end = [full]*v + [perp(p, h) for p in w['ports']] + [full]*R
    for q, U in w.get('helper_final', {}).items(): end[int(q)] = basis(U)
    start = [(p,) for p in w['ports']] + [()]*(v + R)
    for q, U in w.get('start_extra', {}).items(): start[int(q)] = basis(U)
    cur = list(start); chain = [[(None, start[q])] for q in range(n)]; paid = Counter(); runs = defaultdict(list)
    early = w.get('early_cleanup', []); tail = [w['cleanup'][i][:-1] + [U] for i, U in early]
    def climb(q, U, tag):
        U = basis(U)
        if not contained(cur[q], U): raise RuntimeError('not monotone %s %s %s %s' % (q, cur[q], U, tag))
        d = len(U) - len(cur[q])
        if d: paid[d] += 1; runs[q].append(d); chain[q].append((tag, U))
        cur[q] = U
    for sec in ('A', 'B'):
        if sec == 'B':
            for q, U in zip(w['centers'], w['center_frames']):
                U = basis(U)
                if cur[q] != U: raise RuntimeError('centre frame')
        for i, e in enumerate(w[sec]):
            for q in (e[1:3] if e[0] == 'add' else e[1:2]): climb(q, e[-1], (sec, i))
    for i, e in enumerate(tail):
        for q in (e[1:3] if e[0] == 'add' else e[1:2]): climb(q, e[-1], ('tail', i))
    for q in range(n): climb(q, end[q], ('final', None))
    return start, end, chain, paid, runs
