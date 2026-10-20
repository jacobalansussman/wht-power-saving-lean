"""x-exit-words: lower a hub canonical word with EXIT-FRAME helpers, split centres on source registers, a
compensation phase and a tail (early cleanup) into gcert/1 for OUR gx checker, through Shea's source-copy route.

Adapted from Ryan Shea's diagnostic `canonical_to_gcert_probe.py` (shea256/fourier-transform-below-nlogn,
research/complex-transfer-588-596/experiments), which handles #588 atm3.  Changes (x-exit-words):
  * the centres that live on source registers (gauge44 split centres) are remapped to the source copies, like
    every other mutable use of a source;
  * the compensation phase (target += helper/centre at FRAME 0: the hub's explicit entrance-dirt cancellation)
    is NOT emitted, as in Shea's lowering: in our engine that is the derived frame-0 dirt gate, and its
    correctness is the E7 condition (gx.q2) on the x = 0 replay;
  * cleanup AND tail (early cleanup at exit frames) write helpers only once sources are copied; they are omitted
    (invocation helpers need not be restored), exactly as Shea omits the cleanup;
  * helpers end at their EXIT frames c['final'] (dim 4/8/9).
Checks with our gx.py (unchanged): scal (E5 exact scalar identity: x restored, every target gets exactly its
source), q2 (E7 entrance dirt), and check1 with the two exit-frame relaxations a tier-2 engine needs:
  E1  a slot's final frame is its exit frame (nested over its start) instead of the full frame;
  E6  N + (entrance tails) = sum_slots dim(exit) + 2 v (h-1) + cst   (instead of R h + ...).
"""
import sys, os, json, gzip, hashlib, importlib.util, inspect
from collections import Counter
from fractions import Fraction as Q

def read(p):
    raw = open(p, 'rb').read()
    return json.loads(gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw)

def lower(c, gx):
    c = json.loads(json.dumps(c))
    v, R, h = c['v'], c['R'], c['h']
    base = 2 * v + R
    remap = lambda r: base + r if r < v else r
    for part in ('A', 'B', 'precondition', 'cleanup', 'compensation', 'cut', 'tail'):
        for e in c.get(part) or []:
            e[1] = remap(e[1])
            if e[0] == 'add': e[2] = remap(e[2])
    c['centers'] = [remap(r) for r in c['centers']]
    c['A'] = [['add', base + r, r, 1, 1, c['start'][r]] for r in range(v)] + c['A']
    c['start'] += [0] * v
    c['final'] += [1] * v
    c['R'] += v
    for rank in (1, h - 1):
        c['local_profile'][str(rank)] = c['local_profile'].get(str(rank), 0) + v
    v, h, R = c['v'], c['h'], c['R']
    n = 2 * v + R
    units = [Q(1)] * n
    A, B = [], []
    def emit(e, out):
        kind, r = e[:2]
        if kind == 'scale':
            out.append(['out', e[-1], r, [], []]); f = Q(e[2], e[3]); assert f; units[r] *= f
        else:
            assert kind == 'add'
            _, t, s, num, den, f = e
            fac = Q(num, den) * units[s] / units[t]
            out.append(['out', f, s, [[t, fac.numerator, fac.denominator]], []])
    for e in c['A']:
        regs = e[1:3] if e[0] == 'add' else e[1:2]
        assert all(not v <= r < 2 * v for r in regs), 'A touches a target'
        emit(e, A)
    ret = [[k, r, c['center_frames'][k]] for k, r in enumerate(c['centers'])]
    for _, r, f in ret: A.append(['out', f, r, [], []])
    centres = {r: k for k, r, _ in ret}
    S = [[Q(0)] * len(ret) for _ in range(v)]
    for _, t, s, num, den, _ in c['cut']:
        assert v <= t < 2 * v and s in centres, ('cut', t, s)
        S[t - v][centres[s]] += Q(num, den)
    for e in reversed(c['precondition']):
        regs = e[1:3] if e[0] == 'add' else e[1:2]
        assert all(v <= r < 2 * v for r in regs) and e[-1] == 0
        if e[0] == 'scale':
            S[e[1] - v] = [x / Q(e[2], e[3]) for x in S[e[1] - v]]
        else:
            _, t, s, num, den, _ = e
            S[t - v] = [x - Q(num, den) * y for x, y in zip(S[t - v], S[s - v])]
    scatter = []
    for row in S:
        out = []
        for k, val in enumerate(row):
            val *= units[ret[k][1]]
            if val: out.append([k, val.numerator, val.denominator])
        scatter.append(out)
    for e in c['precondition'] + c['B']:
        emit(e, B)
    for r, f in enumerate(c['final']):
        B.append(['out', f, r, [], []])
    for part in ('cleanup', 'tail'):
        assert all(e[1] >= 2 * v for e in c.get(part) or []), part + ' writes a data register'
    bad = [[i, str(q)] for i, q in enumerate(units[:2 * v]) if q != 1]
    assert not bad, ('nonunit data output', bad[:5])
    dims = [len(f) for f in c['frames']]
    cur = c['start'].copy()
    hist = {k: Counter() for k in ('x', 'y', 's', 'c')}
    cls = lambda r: 'x' if r < v else 'y' if r < 2 * v else 's'
    def climb(r, f):
        if cur[r] != f:
            o = cur[r]
            assert dims[f] > dims[o] and gx.inside(c['frames'][o], c['frames'][f]), (r, o, f)
            hist[cls(r)][dims[f] - dims[o]] += 1; cur[r] = f
    for ph, gates in (('A', A), ('B', B)):
        if ph == 'B':
            for _, r, f in ret:
                assert cur[r] == f; hist['c'][dims[f]] += 1
        for g in gates:
            for r in gx.gate_regs(g)[0]: climb(r, g[1])
    for r, f in enumerate(c['final']): climb(r, f)
    total = Counter()
    for part in hist.values(): total.update(part)
    prof_ok = dict(total) == {int(r): n for r, n in c['local_profile'].items()}
    cert = dict(format='gcert/1', h=h, v=v, R=R, N=sum(r * n for r, n in total.items()),
                cst=sum(dims[f] for _, _, f in ret), ports=c['ports'], frames=c['frames'],
                start=c['start'], final=c['final'], ext=[[r, f] for r, f in enumerate(c['start']) if r >= 2 * v and f],
                ret=ret, scat={'table': scatter}, A=A, B=B,
                blocks={k: {str(r): n for r, n in sorted(p.items())} for k, p in hist.items()})
    return cert, prof_ok, dict(total)

def check1x_from(gx, drop_e4=False):
    """our gx.check1 with the two exit-frame relaxations (E1 slot final, E6 N), by source substitution"""
    src = inspect.getsource(gx.check1)
    a = '''        need(start[r] == ext.get(r, 0) and final[r] == 1, "E1: slot start / final")'''
    b = '''        need(start[r] == ext.get(r, 0) and inside(fr[start[r]], fr[final[r]]), "E1x: slot start / exit")'''
    a2 = '''    need(N == c["N"] and N + tails == R * h + 2 * v * (h - 1) + c["cst"], "E6: N")'''
    b2 = '''    need(N == c["N"] and N + tails == sum(dim[final[r]] for r in range(2 * v, nreg)) + 2 * v * (h - 1) + c["cst"], "E6x: N")'''
    assert src.count(a) == 1 and src.count(a2) == 1
    src = src.replace('def check1(', 'def check1x(').replace(a, b).replace(a2, b2)
    if drop_e4:   # E4 = the acyclic-pivot SUFFICIENT condition; Shea replaces it by an exact replay of every target basis vector
        e4 = '''    need(not (piv & ytg), "E4: a pivot is the target of a gate reading y")'''
        assert src.count(e4) == 1; src = src.replace(e4, '    pass')
    ns = dict(vars(gx)); exec(src, ns); return ns['check1x']

if __name__ == '__main__':
    wp, gxp, out = sys.argv[1], sys.argv[2], sys.argv[3]
    spec = importlib.util.spec_from_file_location('gx', gxp); gx = importlib.util.module_from_spec(spec); spec.loader.exec_module(gx)
    res = dict(word_sha256=hashlib.sha256(open(wp, 'rb').read()).hexdigest(), gx_sha256=hashlib.sha256(open(gxp, 'rb').read()).hexdigest())
    try:
        cert, prof_ok, total = lower(read(wp), gx)
    except AssertionError as e:
        res.update(lowering_failed=True, error=repr(e)); print(json.dumps(res, indent=1)); json.dump(res, open(out, 'w'), indent=1); sys.exit(1)
    res.update(local_histogram_preserved=prof_ok, histogram=total, N=cert['N'])
    open(out[:-5] + '.gcert.json.gz', 'wb').write(gzip.compress(json.dumps(cert, sort_keys=True).encode(), mtime=0))
    chk = {}
    for lab, f in (('scal_E5_exact_scalar_identity', lambda: gx.scal(cert)), ('q2_E7_entrance_dirt', lambda: gx.q2(cert)),
                   ('check1_unchanged', lambda: gx.check1(cert, scalar=False)),
                   ('check1x_exit_frames', lambda: check1x_from(gx)(cert, scalar=False)),
                   ('check1x_exit_frames_without_E4', lambda: check1x_from(gx, drop_e4=True)(cert, scalar=False))):
        try:
            f(); chk[lab] = dict(passed=True)
        except AssertionError as e:
            chk[lab] = dict(passed=False, error=str(e)[:300])
    res['checks'] = chk
    print(json.dumps(res, indent=1)); json.dump(res, open(out, 'w'), indent=1)
