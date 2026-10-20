"""d-halfread: shared replay of a hub FULL-WORD in the audit's order (audit_local.py of PR870/#942):
program = precondition + compensation (derived) + A + cut + B + tail(early cleanup with its frames) + cleanup.
Frames: audit climb rule (A, centres, B, tail, then final ends). Values: exact rationals over symbolic columns."""
import json, gzip
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
def contained(a, b): return basis(list(a)+list(b)) == tuple(b)
def popc(x): return bin(x).count('1')
def load(p): return json.loads(gzip.decompress(open(p, 'rb').read()))
def save(w, p):
    open(p, 'wb').write(gzip.compress(json.dumps(w, separators=(',', ':')).encode(), mtime=0))
def regs(e): return e[1:3] if e[0] == 'add' else e[1:2]
def apply(e, rows):
    if e[0] == 'add':
        t, s = e[1:3]; q = F(e[3], e[4]); dst = rows[t]
        for k, x in list(rows[s].items()):
            a = dst.get(k, F(0)) + q*x
            if a: dst[k] = a
            else: dst.pop(k, None)
    else:
        t = e[1]; q = F(e[2], e[3]); rows[t] = {k: q*x for k, x in rows[t].items()}
def program(w):
    """list of (phase, gate) in the audit's order; compensation derived as in the audit."""
    v, R = w['report']['v'], w['report']['R']; n = 2*v+R
    main = w['A']+w['cut']+w['B']
    raw = [{k: F(1)} for k in range(n)]
    for e in main: apply(e, raw)
    dirt = [{} for _ in range(n)]
    for y in range(v, 2*v): dirt[y] = {k: a for k, a in raw[y].items() if k >= 2*v}
    for e in w['precondition']: apply(e, dirt)
    comp = [['add', y, k, -a.numerator, a.denominator, []] for y in range(v, 2*v) for k, a in sorted(dirt[y].items())]
    early = w.get('early_cleanup', []); ids = {i for i, U in early}
    tail = [w['cleanup'][i][:-1]+[U] for i, U in early]
    cleanup = [e for i, e in enumerate(w['cleanup']) if i not in ids]
    out = [('precondition', e) for e in w['precondition']] + [('compensation', e) for e in comp]
    out += [('A', e) for e in w['A']] + [('cut', e) for e in w['cut']] + [('B', e) for e in w['B']]
    out += [('tail', e) for e in tail] + [('cleanup', e) for e in cleanup]
    return out
def ends(w):
    v, R, h = w['report']['v'], w['report']['R'], w['report']['h']
    full = basis([1 << i for i in range(h)])
    end = [full]*v+[basis([x for x in range(1, 1 << h) if popc(x & p) % 2 == 0]) for p in w['ports']]+[full]*R
    for q, U in w.get('helper_final', {}).items(): end[int(q)] = basis(U)
    start = [(p,) for p in w['ports']]+[()]*(v+R)
    for q, U in w.get('start_extra', {}).items(): start[int(q)] = basis(U)
    return start, end
def frame_replay(w, prog=None):
    """per register: chain list of (time, frame); time index into prog. Audit climb rule (A, B, tail; centres between A and B)."""
    v, R = w['report']['v'], w['report']['R']; n = 2*v+R
    prog = prog or program(w); start, end = ends(w)
    cur = list(start); chain = [[(-1, start[q])] for q in range(n)]; bad = []
    def climb(q, U, T):
        U = basis(U)
        if not contained(cur[q], U): bad.append((T, q, cur[q], U)); return
        if U != cur[q]: cur[q] = U; chain[q].append((T, U))
    centred = False
    for T, (ph, e) in enumerate(prog):
        if ph == 'cut' and not centred:
            for q, U in zip(w['centers'], w['center_frames']): climb(q, U, T)
            centred = True
        if ph in ('A', 'B', 'tail'):
            for q in regs(e): climb(q, e[-1], T)
    for q in range(n): climb(q, end[q], len(prog))
    return chain, bad
