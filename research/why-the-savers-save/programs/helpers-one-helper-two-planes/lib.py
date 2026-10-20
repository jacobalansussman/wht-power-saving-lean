"""d-aboveplane: shared helpers for the U8 FULL-WORD (the format audit_local.py checks).
Frame replay follows audit_local.py exactly (A, centre paths, B, tail = early cleanup, end climb); the value replay
follows its program order: precondition + compensation (derived as the audit derives it) + A + cut + B + tail + cleanup.
Content rows are dicts {initial-register-index: Fraction} (the audit's identity-column form), so a row's keys < v are
sources, v..2v-1 targets' initial values, >= 2v helper dirt."""
import gzip, json
from fractions import Fraction as F
from collections import Counter, defaultdict

def load(path):
    return json.loads(gzip.decompress(open(path, 'rb').read()))

def basis(vs):
    piv = {}
    for a in vs:
        while a:
            bit = a.bit_length() - 1
            if bit in piv: a ^= piv[bit]
            else:
                piv[bit] = a; break
    for k in sorted(piv):
        for j in piv:
            if j > k and piv[j] >> k & 1: piv[j] ^= piv[k]
    return tuple(piv[k] for k in sorted(piv, reverse=True))

def contained(a, b): return basis(list(a) + list(b)) == tuple(b)
def join(a, b): return basis(list(a) + list(b))
def pc(x): return bin(x).count('1')
def dot(a, b): return pc(a & b) & 1
def perp(u, h=9): return basis([x for x in range(1, 1 << h) if dot(x, u) == 0])
def inperp(Fr, u): return all(dot(x, u) == 0 for x in Fr)

def registers(e): return e[1:3] if e[0] == 'add' else e[1:2]

def apply(e, rows):
    if e[0] == 'add':
        t, s = e[1:3]; q = F(e[3], e[4]); dst = rows[t]
        for k, x in list(rows[s].items()):
            a = dst.get(k, F(0)) + q * x
            if a: dst[k] = a
            else: dst.pop(k, None)
    else:
        t = e[1]; q = F(e[2], e[3]); rows[t] = {k: q * x for k, x in rows[t].items()}

def geometry(w):
    v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2 * v + R
    full = basis([1 << i for i in range(h)])
    end = [full] * v + [perp(p, h) for p in w['ports']] + [full] * R
    for q, U in w.get('helper_final', {}).items(): end[int(q)] = basis(U)
    start = [(p,) for p in w['ports']] + [()] * (v + R)
    for q, U in w.get('start_extra', {}).items(): start[int(q)] = basis(U)
    return v, h, R, n, full, start, end

def program(w):
    """Returns (segments, compensation) in the audit's order; segments = list of (name, gates)."""
    v, h, R, n, full, start, end = geometry(w)
    ident = lambda: [{k: F(1)} for k in range(n)]
    main = w['A'] + w['cut'] + w['B']
    raw = ident()
    for e in main: apply(e, raw)
    dirt = [{} for _ in range(n)]
    for y in range(v, 2 * v): dirt[y] = {k: a for k, a in raw[y].items() if k >= 2 * v}
    for e in w['precondition']: apply(e, dirt)
    comp = [['add', y, k, -a.numerator, a.denominator, []] for y in range(v, 2 * v) for k, a in sorted(dirt[y].items())]
    early = w.get('early_cleanup', []); earlyids = {i for i, U in early}
    tail = [w['cleanup'][i][:-1] + [U] for i, U in early]
    cleanup = [e for i, e in enumerate(w['cleanup']) if i not in earlyids]
    segs = [('precondition', w['precondition']), ('compensation', comp), ('A', w['A']), ('cut', w['cut']), ('B', w['B']),
            ('tail', tail), ('cleanup', cleanup)]
    return segs, comp

def frame_replay(w):
    """Audit climb rule. Returns paid Counter, chains {q: [(t, seg, frame)]} with t a global gate index in program order
    (centre paths get t = index of first B gate), and the cur frames at the end of B (endB)."""
    v, h, R, n, full, start, end = geometry(w)
    segs, comp = program(w)
    off = {}; t = 0
    for name, g in segs: off[name] = t; t += len(g)
    cur = list(start); paid = Counter(); chains = {q: [(-1, 'start', start[q])] for q in range(n)}
    def climb(q, U, t, seg):
        U = basis(U); assert contained(cur[q], U), ('not monotone', q, cur[q], U, seg, t)
        d = len(U) - len(cur[q])
        if d: paid[d] += 1; chains[q].append((t, seg, U))
        cur[q] = U
    for name in ('A', 'B', 'tail'):
        if name == 'B':
            for q, U in zip(w['centers'], w['center_frames']):
                U = basis(U); assert cur[q] == U, ('centre frame', q)
                paid[len(U)] += 1  # centre path to the empty target level (target_levels empty)
        for i, e in enumerate(dict(segs)[name]):
            for q in registers(e): climb(q, e[-1], off[name] + i, name)
    endB = None
    for q, U in enumerate(end): climb(q, U, t, 'end')
    return paid, chains, off, segs

def value_replay(w, watch=None):
    """Full value replay in program order from identity columns. If watch is a set of registers, records their content
    snapshots after every gate touching them: snaps[q] = [(t, seg, row-copy)]."""
    v, h, R, n, full, start, end = geometry(w)
    segs, comp = program(w)
    rows = [{k: F(1)} for k in range(n)]; snaps = defaultdict(list); t = 0
    for name, g in segs:
        for e in g:
            apply(e, rows)
            if watch is not None and e[1] in watch: snaps[e[1]].append((t, name, dict(rows[e[1]])))
            t += 1
    return rows, snaps
