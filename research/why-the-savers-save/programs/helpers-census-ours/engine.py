"""g-x4-order frame engine: replays the literal frame chains of a FULL-WORD (#799 X4 format) with monotone repair,
records every register's chain (segment, op index, frame), its runs and stops, and checks legality
(target caps, helper tails, END frames, and the audit's final-helper cone conditions on the end-of-B frames).
Same semantics as build_jx.build() (hub construction/jx/build_jx.py) and independent/audit_local.py, but also repairs
A-segment frames (A frames are literal in the original word; the repair is the identity there, checked in census.py)."""
import gzip, json
from functools import lru_cache
from fractions import Fraction as Fr
from collections import Counter

V, H = 112, 9

def _basis(vs):
    piv = {}
    for a in vs:
        while a:
            b = a.bit_length() - 1
            if b in piv: a ^= piv[b]
            else: piv[b] = a; break
    for k in sorted(piv):
        for j in sorted(piv):
            if j > k and (piv[j] >> k) & 1: piv[j] ^= piv[k]
    return tuple(piv[k] for k in sorted(piv, reverse=True))

@lru_cache(maxsize=None)
def basis_t(t): return _basis(list(t))
def basis(vs): return basis_t(tuple(vs))
@lru_cache(maxsize=None)
def join(a, b): return _basis(list(a) + list(b))
def contained(a, b): return join(a, b) == b
def registers(e): return (e[1], e[2]) if e[0] == 'add' else (e[1],)
def hyper(p): return _basis([x for x in range(1, 1 << H) if bin(x & p).count("1") % 2 == 0])
FULL = _basis([1 << i for i in range(H)])

def load(path): return json.loads(gzip.decompress(open(path, 'rb').read()))

def ends(w):
    R = w['report']['R']; n = 2 * V + R
    END = [FULL] * V + [hyper(p) for p in w['ports']] + [FULL] * R
    for q, U in w['helper_final'].items(): END[int(q)] = basis(U)
    return END

def cone_caps(w, B=None):
    """the audit's final helper cone conditions: endB[k] must lie in END[q] for k in pinv[q] and for k != q with p[k][q] != 0,
    for every low-ended helper q (prefix = precondition + compensation + A + cut + B).  Returns {k: [U, ...]}."""
    B = w['B'] if B is None else B
    R = w['report']['R']; n = 2 * V + R; v = V
    def apply(e, rows):
        if e[0] == 'add':
            t, s = e[1:3]; q = Fr(e[3], e[4]); dst = rows[t]
            for k, x in list(rows[s].items()):
                a = dst.get(k, Fr(0)) + q * x
                if a: dst[k] = a
                else: dst.pop(k, None)
        else:
            t = e[1]; q = Fr(e[2], e[3]); rows[t] = {k: q * x for k, x in rows[t].items()}
    def invert(e):
        r = list(e)
        if r[0] == 'add': r[3] = -r[3]
        else: r[2], r[3] = r[3], r[2]
        return r
    ident = lambda: [{k: Fr(1)} for k in range(n)]
    main = w['A'] + w['cut'] + B
    raw = ident()
    for e in main: apply(e, raw)
    dirt = [{} for _ in range(n)]
    for y in range(v, 2 * v): dirt[y] = {k: a for k, a in raw[y].items() if k >= 2 * v}
    for e in w['precondition']: apply(e, dirt)
    comp = [['add', y, k, -a.numerator, a.denominator, []] for y in range(v, 2 * v) for k, a in sorted(dirt[y].items())]
    prefix = w['precondition'] + comp + main
    p = ident(); pinv = ident()
    for e in prefix: apply(e, p)
    for e in reversed(prefix): apply(invert(e), pinv)
    END = ends(w); caps = {}
    for qs in w['helper_final']:
        q = int(qs); U = END[q]
        ks = set(pinv[q]) | {k for k in range(n) if k != q and p[k].get(q)}
        for k in ks: caps.setdefault(k, []).append(U)
    return caps

_ACACHE = {}

def replay(w, A=None, B=None, early=None, explicit=None, caps=None, record=False, stop_at_first=True):
    """explicit: {('A'|'B'|'E', idx): frame list} (base frame for that op; result = join(base, operands)).
    Returns dict(ok, why, paid, A, B, early, chains) where chains[q] = [(seg, idx, frame), ...] (record=True)."""
    A = w['A'] if A is None else A; B = w['B'] if B is None else B
    early = w['early_cleanup'] if early is None else early; explicit = explicit or {}
    R = w['report']['R']; n = 2 * V + R; END = ends(w)
    cur = [(p,) for p in w['ports']] + [()] * (V + R)
    for q, U in w['start_extra'].items(): cur[int(q)] = basis(U)
    chains = [[('S', -1, cur[q], 0)] for q in range(n)] if record else None
    paid = Counter(); why = None
    def seg(name, ops, frames_of, capcheck):
        nonlocal why
        out = []
        for i, e in enumerate(ops):
            rs = registers(e)
            base = explicit.get((name, i))
            F = basis(base) if base is not None else basis(frames_of(i, e))
            for q in rs: F = join(F, cur[q])
            for q in rs:
                if capcheck(q, F) is False:
                    why = why or ('cap', name, i, q)
                if F != cur[q]:
                    d = len(F) - len(cur[q]); paid[d] += 1; cur[q] = F
                    if record: chains[q].append((name, i, F, d))
            out.append(F)
        return out
    capT = lambda q, F: (contained(F, END[q]) if V <= q < 2 * V else None)
    cacheable = (A is w['A']) and not record and not any(k[0] == 'A' for k in explicit)
    ck = id(w)
    if cacheable and ck in _ACACHE and _ACACHE[ck][0] is w:
        _, c_cur, c_paid, c_why, Af = _ACACHE[ck]
        cur[:] = c_cur; paid.update(c_paid); why = c_why
    else:
        Af = seg('A', A, lambda i, e: e[-1] or (), capT)
        # centre levels (target_levels: served levels; centre keeps its frame afterwards)
        lev = {int(t): basis(U) for t, U in w.get('target_levels', {}).items()}
        for t in range(V):
            if cur[V + t] != lev.get(t, ()): why = why or ('target level', t)
        for q, U in zip(w['centers'], w['center_frames']):
            U = basis(U)
            if cur[q] != U: why = why or ('centre frame', q)
            served = {lev.get(e[1] - V, ()) for e in w['cut'] if e[2] == q}
            for T in sorted(served, key=lambda x: -len(x)):
                if T != U:
                    paid[len(U) - len(T)] += 1
                    if record: chains[q].append(('C', q, T, len(U) - len(T)))
                U = T
        if cacheable: _ACACHE[ck] = (w, list(cur), Counter(paid), why, Af)
    Bf = seg('B', B, lambda i, e: e[-1] or (), capT)
    endB = list(cur)
    if caps is not None:
        for k, Us in caps.items():
            for U in Us:
                if not contained(endB[k], U): why = why or ('cone', k)
    hf = set(int(q) for q in w['helper_final'])
    cl = w['cleanup']
    Ef = []
    for j, (i, U) in enumerate(early):
        e = cl[i]; rs = registers(e)
        base = explicit.get(('E', j))
        F = basis(base if base is not None else U)
        for q in rs: F = join(F, cur[q])
        for q in rs:
            if q in hf and not contained(F, END[q]): why = why or ('tailend', j, q)
            if F != cur[q]:
                d = len(F) - len(cur[q]); paid[d] += 1; cur[q] = F
                if record: chains[q].append(('E', j, F, d))
        Ef.append([i, list(F)])
    for q in range(n):
        if not contained(cur[q], END[q]): why = why or ('end', q)
        elif END[q] != cur[q]:
            d = len(END[q]) - len(cur[q]); paid[d] += 1
            if record: chains[q].append(('END', -1, END[q], d))
    return dict(ok=why is None, why=why, paid=paid, Af=Af, Bf=Bf, Ef=Ef, chains=chains, endB=endB, END=END)

def word_with(w, A=None, B=None, res=None):
    """a new word dict with the given ops and the repaired frames of res written in."""
    import copy
    o = copy.deepcopy(w)
    A = w['A'] if A is None else A; B = w['B'] if B is None else B
    o['A'] = [list(e[:-1]) + [list(F)] for e, F in zip(A, res['Af'])]
    o['B'] = [list(e[:-1]) + [list(F)] for e, F in zip(B, res['Bf'])]
    o['early_cleanup'] = res['Ef']
    o['report'] = dict(o['report']); o['report']['local_profile'] = {str(k): v for k, v in sorted(res['paid'].items())}
    return o

def save(o, path):
    open(path, 'wb').write(gzip.compress(json.dumps(o, sort_keys=True, separators=(',', ':')).encode(), mtime=0))
