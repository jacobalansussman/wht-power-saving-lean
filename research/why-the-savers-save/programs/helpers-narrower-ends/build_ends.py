#!/usr/bin/env python3
"""d-hpe-multigate STEP 2 (surviving variant): lower the end of chosen low-ended U8 helpers to a smaller frame U by
re-annotating frames only (scalar gates, their order, entrances, the centre interface and every other end unchanged).

Rule (maximal admissible annotation under caps).  Each register r gets an end cap: its new end for the chosen helpers,
its old end otherwise.  A backward pass gives every A gate, B gate and early-cleanup (tail) gate a cap
C_g = intersection of the caps of its operands' frames after g; a centre must equal its interface frame before B.
The forward pass sets the new frame F'_g = old_g  INTERSECT  C_g, which contains both operands' current frames
(induction), so the chains stay monotone and only chains that meet the chosen helpers' cones change.
Then helper_final[q] = U.  The local profile and width inventory are recounted with the PR870 audit's climb rule and
written into report so that the unchanged audit re-derives and compares them.
Usage: python3 -I build_ends.py FULL-WORD.json.gz lazy_cones.json OUT.json.gz q1,q2,...  (U = lazy cone join of q)"""
import sys, gzip, json, hashlib
from collections import Counter
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from cones_lib import basis, contained, registers
_perp = {}
def perp(U):
    if U not in _perp:
        _perp[U] = basis([x for x in range(1, 512) if all(bin(x & b).count('1') % 2 == 0 for b in U)])
    return _perp[U]
_cap = {}
def meet(U, V):
    k = (U, V)
    if k not in _cap: _cap[k] = perp(basis(list(perp(U))+list(perp(V))))
    return _cap[k]
def build(w, newend):
    v, h, R = w['report']['v'], w['report']['h'], w['report']['R']; n = 2*v+R
    full = basis([1 << i for i in range(h)])
    end = [full]*v+[basis([x for x in range(1, 1 << h) if bin(x & p).count('1') % 2 == 0]) for p in w['ports']]+[full]*R
    for q, U in w['helper_final'].items(): end[int(q)] = basis(U)
    for q, U in newend.items():
        assert contained(U, end[q]), 'new end must lie inside the old one'; end[q] = U
    start = [(p,) for p in w['ports']]+[()]*(v+R)
    for q, U in w['start_extra'].items(): start[int(q)] = basis(U)
    ev = [('A', i) for i in range(len(w['A']))]+[('centre', k) for k in range(len(w['centers']))]+[('B', i) for i in range(len(w['B']))]+[('tail', k) for k in range(len(w['early_cleanup']))]
    def regs_of(x):
        s, i = x
        if s == 'centre': return [w['centers'][i]]
        if s == 'tail': return list(registers(w['cleanup'][w['early_cleanup'][i][0]]))
        return list(registers(w[s][i]))
    def old_frame(x):
        s, i = x
        if s == 'centre': return basis(w['center_frames'][i])
        if s == 'tail': return basis(w['early_cleanup'][i][1])
        return basis(w[s][i][-1])
    cap = list(end); C = {}
    for x in reversed(ev):
        rs = regs_of(x)
        if x[0] == 'centre':
            Z = old_frame(x)
            if not contained(Z, cap[rs[0]]): return None, 'centre interface frame exceeds cap for register %d' % rs[0]
            C[x] = Z; cap[rs[0]] = Z; continue
        c = full
        for r in rs: c = meet(c, cap[r])
        C[x] = c
        for r in rs: cap[r] = c
    for r in range(n):
        if not contained(start[r], cap[r]): return None, 'start frame of register %d exceeds its cap' % r
    w2 = json.loads(json.dumps(w)); changed = Counter()
    cur = start.copy()
    for x in ev:
        rs = regs_of(x); F = old_frame(x)
        if x[0] == 'centre':
            assert cur[rs[0]] == F or contained(cur[rs[0]], F); continue
        Fn = meet(F, C[x])
        for r in rs: assert contained(cur[r], Fn), ('non-monotone', x)
        for r in rs: cur[r] = Fn
        if Fn != F:
            changed[x[0]] += 1
            if x[0] == 'tail': w2['early_cleanup'][x[1]][1] = list(Fn)
            else: w2[x[0]][x[1]][-1] = list(Fn)
    for q, U in newend.items():
        assert contained(cur[q], U); w2['helper_final'][str(q)] = list(U)
    # recount with the audit's rule (A, centre path, B, tail, final ends)
    cur = start.copy(); paid = Counter()
    def climb(q, U):
        U = basis(U); assert contained(cur[q], U); d = len(U)-len(cur[q])
        if d: paid[d] += 1
        cur[q] = U
    for section in ('A', 'B'):
        if section == 'B':
            lev = {int(t): basis(U) for t, U in w2.get('target_levels', {}).items()}
            for q, U in zip(w2['centers'], w2['center_frames']):
                U = basis(U); assert cur[q] == U
                served = {lev.get(e[1]-v, ()) for e in w2['cut'] if e[2] == q}
                for T in sorted(served, key=lambda x: -len(x)):
                    if T != U: paid[len(U)-len(T)] += 1
                    U = T
        for e in w2[section]:
            for q in registers(e): climb(q, e[-1])
    for i, U in w2['early_cleanup']:
        for q in registers(w2['cleanup'][i]): climb(q, U)
    for q in range(n): climb(q, end[q])
    inv = Counter(len(end[q])-len(start[q]) for q in range(2*v, n))
    w2['report']['local_profile'] = {str(k): x for k, x in sorted(paid.items())}
    w2['report']['width_inventory'] = {str(k): x for k, x in sorted(inv.items())}
    return w2, dict(changed=dict(changed), profile=dict(sorted(paid.items())), widths=dict(sorted(inv.items())),
                    width_total=sum(len(end[q])-len(start[q]) for q in range(2*v, n)))
if __name__ == '__main__':
    w = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
    lz = json.load(open(sys.argv[2]))['helpers']
    qs = [int(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 and sys.argv[4] else []
    newend = {q: tuple(lz[str(q)]['Ulazy']) for q in qs}
    w2, info = build(w, newend)
    if w2 is None: print('INFEASIBLE', info); sys.exit(1)
    raw = json.dumps(w2, separators=(',', ':')).encode()
    open(sys.argv[3], 'wb').write(gzip.compress(raw, mtime=0))
    info['content_sha256'] = hashlib.sha256(raw).hexdigest(); print(json.dumps(info))
