#!/usr/bin/env python3
"""p-lower (round wht45/lever, 19 Oct 2026): lower a staged544 BORROW canonical word (CANONICAL-WORD, as written by
our tools/transport.py or the hub's native chain) into gcert/1w (spec: p-lower/FORMAT.md).

Extends bh-gap/k-lower's lower_x.py (itself adapted from Ryan Shea's canonical_to_gcert_probe.py): units absorb the
hub's scales, centres become retained totals, the cut a scatter table, exit frames are kept, the hub's tail is the
last part of phase B followed by one final climb per register.  New here:
  * warps: gate ['warp', r, op, f0, f1] (the hub's event shape) in A (24 on x, the lifted-target port) and B (72,
    dirty borrowing), and the operator table `warps`;
  * no copies of the sources (mutable x, as gcert/1b and the U8X0 route); the zero-frame y-only precondition P is
    the first nP gates of A (P touches the 12 lifted targets, so lower_x's P^-1-into-the-scatter is impossible) and
    the scatter is the hub's cut with units only;
  * the hub's explicit dirt phase, kept as data: D0 = the 36 twisted adds t += c W(s) expanded as warp / add /
    unwarp at frame 0; C0 = the compensation (units after P);
  * the hub's explicit full-frame cleanup, kept as data: U = zalloc of the 91 private zero buffers (registers
    632..722), the cleanup (72 warps at FULL, Z reads), zfree.
Exact checks (asserted, results in the report): frames nested at every climb, every warp carries its frame by G,
D0/C0 on frame-0 registers, U on FULL registers and zero buffers, block histogram of the main phase (warps free) ==
the hub's local profile, invoice identity, units of data registers 1 after B and of every register 1 after U,
the full program in hub order (D0, P, C0, A-rest, scatter, B, U) replayed in the free warp algebra on every input
column ends at the hub endpoint (identity, y_t += x_t, zero buffers 0 at zfree), residuals per helper == the hub
word's widths, and mutants (wrong G column, wrong warp frame, dropped unwarp, tadd coefficient, zero buffer left
dirty) are rejected.

usage: lower_w.py CANONICAL-WORD.json.gz OUT.gcert1w.json.gz REPORT.json --provider NAME [--word BORROW-WORD]
                  [--diff-1b GCERT1B.json.gz] [--mutants]
       lower_w.py --to-1b IN.gcert1w.json.gz OUT.gcert1b.json.gz      (permute warp fields, drop the 1w-only keys)
"""
import sys, json, gzip, hashlib, copy, time
from collections import Counter
from fractions import Fraction as Q


def need(ok, msg):
    if not ok:
        raise AssertionError(msg)


def read(p):
    raw = open(p, 'rb').read()
    return json.loads(gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw)


def rref(vs):
    rows = []
    for x in vs:
        for b in rows:
            if x >> (b.bit_length() - 1) & 1:
                x ^= b
        if x:
            p = x.bit_length() - 1
            rows = [b ^ x if b >> p & 1 else b for b in rows]
            rows.append(x)
            rows.sort(reverse=True)
    return rows


def inside(a, b):
    for x in a:
        for y in b:
            if x >> (y.bit_length() - 1) & 1:
                x ^= y
        if x:
            return False
    return True


def apply_G(cols, x):
    out, i = 0, 0
    while x:
        if x & 1:
            out ^= cols[i]
        x >>= 1
        i += 1
    return out


# ------------------------------------------------------------------ the lowering
def lower(c0, provider):
    c = copy.deepcopy(c0)
    v, R, h = c['v'], c['R'], c['h']
    n = 2 * v + R
    zb = list(c['zero_buffers'])
    need(zb == sorted(zb) and all(2 * v <= r < n for r in zb), 'zero buffers')
    zreg = {s: n + i for i, s in enumerate(zb)}
    nZ = len(zb)

    def reg(r):
        if isinstance(r, list):
            need(r[0] == 'Z' and r[1] in zreg, ('zero buffer', r))
            return zreg[r[1]]
        need(0 <= r < n, ('register', r))
        return r

    pre = c['precondition']
    k = 0
    while k < len(pre) and pre[k][0] == 'tadd':
        k += 1
    tadds, P = pre[:k], pre[k:]
    need(all(e[0] != 'tadd' for e in P), 'twisted add after the precondition proper')
    need(all(e[-1] == 0 and all(v <= r < 2 * v for r in ([e[1], e[2]] if e[0] == 'add' else [e[1]])) for e in P),
         'P is not zero-frame y-only')
    for part in ('A', 'cut', 'B', 'tail', 'compensation'):
        need(all(e[0] != 'tadd' for e in c[part]), 'twisted add outside the precondition')
        need(all(not isinstance(x, list) for e in c[part] for x in e[1:3]), 'zero buffer outside the cleanup')
    units = [Q(1)] * (n + nZ)

    def emit(e, out):
        kind = e[0]
        if kind == 'warp':
            _, r, op, f0, f1 = e
            out.append(['warp', reg(r), op, f0, f1])
        elif kind == 'scale':
            _, r, a, b, f = e
            r = reg(r)
            q = Q(a, b)
            need(q != 0, 'zero scale')
            out.append(['out', f, r, [], []])
            units[r] *= q
        else:
            need(kind == 'add', ('gate kind', kind))
            _, t, s, a, b, f = e
            t, s = reg(t), reg(s)
            need(t != s, 'add onto itself')
            q = Q(a, b) * units[s] / units[t]
            out.append(['out', f, s, [[t, q.numerator, q.denominator]], []])

    D0 = []
    for e in tadds:
        _, t, s, ops, a, b, f = e
        need(f == 0 and v <= t < 2 * v and s >= 2 * v and all(ops), ('twisted add', e))
        for op in reversed(ops):
            emit(['warp', s, op, 0, 0], D0)
        emit(['add', t, s, a, b, 0], D0)
        for op in ops:
            emit(['warp', s, -op, 0, 0], D0)
    A = []
    for e in P:
        emit(e, A)
    nP = len(A)
    C0 = []
    for e in c['compensation']:
        need(e[0] == 'add' and e[-1] == 0, ('compensation gate', e))
        emit(e, C0)
    levels = {int(t): f for t, f in c['target_levels'].items()}
    lifted = {v + t for t in levels}
    yA = Counter()
    for e in c['A']:
        regs = [e[1]] + ([e[2]] if e[0] == 'add' else [])
        ys = [r for r in regs if v <= r < 2 * v]
        if ys:
            need(e[0] != 'warp' and all(r in lifted for r in ys) and e[-1] in (0, levels[ys[0] - v]),
                 ('phase A names a y register that is not a lifted target at its level', e))
            yA[e[0]] += 1
        emit(e, A)
    ret = [[k_, r, c['center_frames'][k_]] for k_, r in enumerate(c['centers'])]
    for _, r, f in ret:
        A.append(['out', f, r, [], []])
    cidx = {r: k_ for k_, r, _ in ret}
    S = [[Q(0)] * len(ret) for _ in range(v)]
    served = Counter()
    for e in c['cut']:
        _, t, s, a, b, f = e
        need(v <= t < 2 * v and s in cidx and f is None, ('cut', e))
        S[t - v][cidx[s]] += Q(a, b)
        served[levels.get(t - v, 0)] += 1
    need(set(served) == {0}, ('a centre serves a lifted target', dict(served)))
    table = []
    for t, row in enumerate(S):
        out = []
        for k_, val in enumerate(row):
            val = val * units[ret[k_][1]] / units[v + t]
            if val:
                out.append([k_, val.numerator, val.denominator])
        table.append(out)
    B = []
    for e in c['B']:
        emit(e, B)
    nB0 = len(B)
    for e in c['tail']:
        need(all(not v <= r < 2 * v for r in ([e[1], e[2]] if e[0] == 'add' else [e[1]])), 'tail touches a target')
        emit(e, B)
    for r, f in enumerate(c['final']):
        B.append(['out', f, r, [], []])
    bad = [[i, str(u)] for i, u in enumerate(units[v:2 * v], v) if u != 1]
    need(not bad, ('nonunit target after B', bad[:5]))
    slot_units_after_B = sum(1 for u in units[2 * v:n] if u != 1)
    x_units_after_B = sum(1 for u in units[:v] if u != 1)
    U = [['zalloc', z] for z in range(n, n + nZ)]
    for e in c['cleanup']:
        if e[0] == 'warp':
            _, r, op, f0, f1 = e
            need(f0 == f1 == 1, 'cleanup warp off FULL')
        else:
            need(e[-1] == 1, ('cleanup gate off FULL', e))
        emit(e, U)
    U += [['zfree', z] for z in range(n, n + nZ)]
    need(all(u == 1 for u in units[v:2 * v]), 'nonunit target after U')
    units_end = [[i, u.numerator, u.denominator] for i, u in enumerate(units[:n]) if u != 1]
    z_units = sum(1 for u in units[n:] if u != 1)        # irrelevant: a zero buffer is 0 at zfree in either scale
    fr = c['frames']
    dim = [len(f) for f in fr]
    E = [r for r in range(2 * v, n) if c['final'][r] != 1]
    need(E == zb, 'exits are not the zero-buffer slots')
    ext = [[r, f] for r, f in enumerate(c['start']) if r >= 2 * v and f]
    residuals = [dim[c['final'][r]] - dim[c['start'][r]] for r in range(2 * v, n)]
    cert = dict(format='gcert/1w', provider=provider, h=h, v=v, R=R, nZ=nZ, copies=0, ports=c['ports'], frames=fr,
                start=c['start'], final=c['final'], ext=ext, exits=E, residuals=residuals,
                levels={str(t): f for t, f in sorted(levels.items())}, nP=nP, A=A, ret=ret, scat={'table': table}, B=B,
                D0=D0, C0=C0, U=U, order=['D0', 'A[:nP]', 'C0', 'A[nP:]', 'scat', 'B', 'U'],
                tadd=[[e[1], e[2], list(e[3]), e[4], e[5]] for e in tadds],
                zero_buffers=[[zreg[s], s] for s in zb],
                warps=[dict(op=w['op'], G_columns=w['G_columns'], from_frame=w['from_frame'], to_frame=w['to_frame'],
                            lender=w['lender'], B_pos=w['B_pos']) for w in c['warps']],
                helper_origin=c['helper_origin'], removed=c['removed'], borrow=c['borrow'],
                units_end=units_end,
                hub_local_profile={str(k_): x for k_, x in sorted((int(a), b) for a, b in c['local_profile'].items())})
    info = dict(nA=len(A), nP=nP, nB_without_tail=nB0, tail=len(c['tail']), nB=len(B), nD0=len(D0), nC0=len(C0),
                nU=len(U), phase_A_y_gates=dict(yA), slot_units_not_one_after_B=slot_units_after_B,
                x_units_not_one_after_B=x_units_after_B,
                units_end_not_one=dict(x=sum(1 for r, _, _ in units_end if r < v), slots=sum(1 for r, _, _ in units_end if r >= 2 * v),
                                       exits=sum(1 for r, _, _ in units_end if r in set(E)), zero_buffers=z_units))
    return cert, info


# ------------------------------------------------------------------ checks on a gcert/1w file
def regs_of(g):
    if g[0] == 'warp':
        return [g[1]]
    if g[0] == 'out':
        return [g[2]] + [t for t, _, _ in g[3]] + list(g[4])
    if g[0] == 'in':
        return [g[2]] + [s for s, _, _ in g[3]]
    if g[0] == 'neg':
        return [g[2]]
    if g[0] in ('zalloc', 'zfree'):
        return [g[1]]
    raise AssertionError(('gate kind', g[0]))


def check_frames(c, table=True):
    """frame chains of the whole program (section 5.1 of FORMAT.md) and the price data of the main phase"""
    v, R, h, nZ = c['v'], c['R'], c['h'], c['nZ']
    n = 2 * v + R
    fr = [tuple(f) for f in c['frames']]
    dim = [len(f) for f in fr]
    need(fr[0] == () and fr[1] == tuple(1 << i for i in range(h - 1, -1, -1)), 'frames 0 and 1')
    need(all(list(f) == rref(list(f)) for f in fr) and len(set(fr)) == len(fr), 'frame table not reduced echelon')
    W = {}
    for i, w in enumerate(c['warps']):
        need(w['op'] == i + 1 and len(rref(w['G_columns'])) == h, ('warp table row', i))
        W[w['op']] = w['G_columns']
        a, b = fr[w['from_frame']], fr[w['to_frame']]
        need(not table or tuple(rref([apply_G(w['G_columns'], x) for x in a])) == b, ('table: G does not carry from -> to', w['op']))
    levels = {int(t): f for t, f in c['levels'].items()}
    cls = lambda r: 'x' if r < v else 'y' if r < 2 * v else 's'
    cur = list(c['start'])
    H = {k: Counter() for k in 'xysc'}
    steps, wsteps = set(), set()
    exts = {r for r, _ in c['ext']} | set(c['exits'])

    def climb(r, f, where):
        o = cur[r]
        if o != f:
            if (o, f) not in steps:
                need(dim[o] < dim[f] and inside(fr[o], fr[f]), ('climb not nested', where, r, o, f))
                steps.add((o, f))
            H[cls(r)][dim[f] - dim[o]] += 1
            cur[r] = f

    def warp(g, where):
        _, r, op, f0, f1 = g
        need(abs(op) in W, ('unknown warp op', where))
        if r >= n:
            need(f0 == f1 == 1, ('zero-buffer warp off FULL', where))
            return
        need(cur[r] == f0, ('warp frame mismatch', where, g, cur[r]))
        if (op, f0, f1) not in wsteps:
            a, b = (f0, f1) if op > 0 else (f1, f0)
            need(dim[a] == dim[b] and tuple(rref([apply_G(W[abs(op)], x) for x in fr[a]])) == fr[b],
                 ('warp G does not carry the frame', where, g))
            wsteps.add((op, f0, f1))
        cur[r] = f1

    def zero_phase(gates, name):
        for i, g in enumerate(gates):
            if g[0] == 'warp':
                need(g[3] == g[4] == 0 and g[1] < n, (name, 'warp not at frame 0', i))
                warp(g, (name, i))
                continue
            need(g[0] == 'out' and g[1] == 0, (name, 'gate not at frame 0', i))
            rr = regs_of(g)
            need(all(r < n and cur[r] == 0 and r not in exts for r in rr), (name, 'register off frame 0', i))
    zero_phase(c['D0'], 'D0')
    nP = c['nP']
    for i, g in enumerate(c['A'][:nP]):
        need(g[0] == 'out' and g[1] == 0 and all(v <= r < 2 * v for r in regs_of(g)), ('P', i))
    zero_phase(c['C0'], 'C0')
    for i, g in enumerate(c['A']):
        rr = regs_of(g)
        need(all(0 <= r < n for r in rr) and len(set(rr)) == len(rr), ('A registers', i))
        if g[0] == 'warp':
            need(not v <= g[1] < 2 * v, ('A: warp on y', i))
            warp(g, ('A', i))
            continue
        for r in rr:
            if v <= r < 2 * v and i >= nP:
                need(r - v in levels and g[1] in (0, levels[r - v]), ('A: y named off its level', i))
            climb(r, g[1], ('A', i))
    for k, s, f in c['ret']:
        need(cls(s) != 'y' and cur[s] == f, 'retained total not at its frame')
        H['c'][dim[f]] += 1
    need(all(cur[v + t] == levels.get(t, 0) for t in range(v)), 'a y role is not at its level at the scatter')
    for i, g in enumerate(c['B']):
        rr = regs_of(g)
        need(all(0 <= r < n for r in rr) and len(set(rr)) == len(rr), ('B registers', i))
        if g[0] == 'warp':
            warp(g, ('B', i))
        else:
            for r in rr:
                climb(r, g[1], ('B', i))
    need(cur == c['final'], 'B does not end at the final frames')
    live = set()
    for i, g in enumerate(c['U']):
        if g[0] == 'zalloc':
            need(n <= g[1] < n + nZ and g[1] not in live, ('U zalloc', i)); live.add(g[1]); continue
        if g[0] == 'zfree':
            need(g[1] in live, ('U zfree', i)); live.discard(g[1]); continue
        if g[0] == 'warp':
            need(g[3] == g[4] == 1 and (g[1] in live or (g[1] < n and cur[g[1]] == 1)), ('U warp', i))
            warp(g, ('U', i))
            continue
        need(g[0] == 'out' and g[1] == 1, ('U gate off FULL', i))
        for r in regs_of(g):
            need(r in live if r >= n else cur[r] == 1, ('U register not FULL / not a live zero buffer', i, r))
    need(not live, 'zero buffer not freed')
    blocks = {k: {str(r): x for r, x in sorted(w.items())} for k, w in H.items()}
    tot = Counter()
    for w in H.values():
        tot.update(w)
    N = sum(r * x for r, x in tot.items())
    cst = sum(dim[f] for _, _, f in c['ret'])
    lhs = N + sum(dim[f] for _, f in c['ext'])
    rhs = sum(dim[c['final'][r]] for r in range(2 * v, n)) + sum(dim[c['final'][r]] - dim[c['start'][r]] for r in range(v)) \
        + sum(dim[c['final'][v + t]] for t in range(v)) + cst
    need(lhs == rhs, ('invoice identity', lhs, rhs))
    return dict(blocks=blocks, hist={str(r): x for r, x in sorted(tot.items())}, N=N, cst=cst, invoice=[lhs, rhs],
                frame_steps=len(steps), warp_steps=len(wsteps), max_block_rank=max(tot))


def warp_row(row, op):
    out = {}
    for (col, word), val in row.items():
        out[col, word[1:] if word and word[0] == -op else (op,) + word] = val
    return out


def axpy(d, src, q):
    for key, u in src.items():
        w = d.get(key, 0) + q * u
        if w:
            d[key] = w
        else:
            d.pop(key, None)


def run(reg, gates, zcheck=None):
    for i, g in enumerate(gates):
        k = g[0]
        if k == 'warp':
            reg[g[1]] = warp_row(reg[g[1]], g[2])
        elif k == 'out':
            src = dict(reg[g[2]])
            for t, a, b in g[3]:
                axpy(reg[t], src, Q(a, b))
        elif k == 'in':
            for s, a, b in g[3]:
                axpy(reg[g[2]], reg[s], Q(a, b))
        elif k == 'neg':
            reg[g[2]] = {key: -u for key, u in reg[g[2]].items()}
        elif k == 'zalloc':
            reg[g[1]] = {}
        elif k == 'zfree':
            need(not reg[g[1]], ('zero buffer not zero at zfree', g[1], i))
            reg[g[1]] = None
        else:
            raise AssertionError(k)


def full_replay(c):
    """section 5.3 of FORMAT.md: the whole program in hub order on every input column (free warp algebra)"""
    v, R, nZ = c['v'], c['R'], c['nZ']
    n = 2 * v + R
    reg = [{(r, ()): Q(1)} for r in range(n)] + [None] * nZ
    nP = c['nP']
    run(reg, c['D0'])
    run(reg, c['A'][:nP])
    run(reg, c['C0'])
    run(reg, c['A'][nP:])
    tot = {k: dict(reg[s]) for k, s, _ in c['ret']}
    for t in range(v):
        for k, a, b in c['scat']['table'][t]:
            axpy(reg[v + t], tot[k], Q(a, b))
    run(reg, c['B'])
    run(reg, c['U'])
    bad = []
    ue = {r: Q(a, b) for r, a, b in c['units_end']}
    for r in range(n):
        want = {(r, ()): 1 / ue.get(r, Q(1))}
        if v <= r < 2 * v:
            want[(r - v, ())] = Q(1)
        if reg[r] != want:
            bad.append(r)
    need(not bad, ('hub endpoint fails on rows', bad[:10]))
    need(all(z is None for z in reg[n:]), 'zero buffer not freed')
    return dict(endpoint='exact on all %d columns: identity, y_t += x_t, %d zero buffers 0 at zfree' % (n, nZ))


def main_phase_words(c):
    """the clean main phase (A, scatter, B) on every input column: warp-worded entries left in x / y rows"""
    v, R = c['v'], c['R']
    n = 2 * v + R
    reg = [{(r, ()): Q(1)} for r in range(n)]
    run(reg, c['A'])
    tot = {k: dict(reg[s]) for k, s, _ in c['ret']}
    for t in range(v):
        for k, a, b in c['scat']['table'][t]:
            axpy(reg[v + t], tot[k], Q(a, b))
    run(reg, c['B'])
    hid = all({(j, w): u for (j, w), u in reg[v + t].items() if j < v} == {(t, ()): 1} for t in range(v))
    return dict(x_row_entries_with_words=sum(1 for t in range(v) for (_, w) in reg[t] if w),
                y_row_entries_with_words=sum(1 for t in range(v) for (_, w) in reg[v + t] if w),
                hid_y_from_x_identity=hid)


# ------------------------------------------------------------------ gcert/1b view and diff
def to_1b(c):
    d = {k: copy.deepcopy(c[k]) for k in ('h', 'v', 'R', 'ports', 'frames', 'start', 'final', 'ext', 'ret', 'scat',
                                          'exits', 'levels', 'nP', 'helper_origin', 'borrow', 'copies', 'residuals')}
    conv = lambda g: ['warp', g[4], g[1], g[2], g[3]] if g[0] == 'warp' else g
    d['A'] = [conv(g) for g in c['A']]
    d['B'] = [conv(g) for g in c['B']]
    d['warps'] = [dict(op=w['op'], G_columns=w['G_columns'], from_frame=w['from_frame'], to_frame=w['to_frame'])
                  for w in c['warps']]
    d['format'] = 'gcert/1b'
    return d


def diff(a, b):
    out = {}
    for k in sorted(set(a) | set(b)):
        if k not in a or k not in b:
            out[k] = 'only in ' + ('ours' if k in a else 'theirs')
        elif a[k] != b[k]:
            if isinstance(a[k], list) and isinstance(b[k], list):
                idx = [i for i in range(min(len(a[k]), len(b[k]))) if a[k][i] != b[k][i]]
                out[k] = dict(len=[len(a[k]), len(b[k])], first_diffs=[[i, a[k][i], b[k][i]] for i in idx[:3]], ndiff=len(idx))
            else:
                out[k] = 'differs'
    return out


# ------------------------------------------------------------------ mutants
def mutants(c):
    res = {}

    def expect_fail(name, mut, fn):
        d = copy.deepcopy(c)
        mut(d)
        try:
            fn(d)
        except AssertionError as e:
            res[name] = 'rejected: ' + str(e)[:160]
            return
        raise SystemExit('MUTANT ACCEPTED: ' + name)

    # (1) a wrong G column: take the first B warp; change one column so G stays invertible but no longer carries the frame
    wi = next(i for i, g in enumerate(c['B']) if g[0] == 'warp')
    op = abs(c['B'][wi][2])
    w = c['warps'][op - 1]
    fr = c['frames']
    cols = w['G_columns']
    choice = None
    for j in range(9):
        for x in range(1, 512):
            G2 = list(cols)
            G2[j] = x
            if len(rref(G2)) == 9 and rref([apply_G(G2, u) for u in fr[w['from_frame']]]) != fr[w['to_frame']]:
                choice = (j, x)
                break
        if choice:
            break
    j, x = choice

    def m1(d):
        d['warps'][op - 1]['G_columns'][j] = x
    expect_fail('wrong G column (warp %d, column %d := %d, still invertible)' % (op, j, x), m1, check_frames)

    expect_fail('wrong G column, warp-table check skipped (the gate-level frame check alone)', m1,
                lambda d: check_frames(d, table=False))
    # (2) a warp gate with a wrong end frame
    def m2(d):
        g = d['B'][wi]
        g[4] = next(f for f in range(2, len(fr)) if f != g[4] and len(fr[f]) == len(fr[g[4]]))
    expect_fail('B warp with a wrong end frame', m2, check_frames)
    # (3) a dropped unwarp in U: the scalar endpoint must fail
    ui = next(i for i, g in enumerate(c['U']) if g[0] == 'warp')
    expect_fail('U unwarp dropped', lambda d: d['U'].pop(ui), full_replay)
    # (4) a twisted-add coefficient doubled
    di = next(i for i, g in enumerate(c['D0']) if g[0] == 'out')
    expect_fail('D0 twisted-add coefficient doubled', lambda d: d['D0'][di][3][0].__setitem__(1, 2 * d['D0'][di][3][0][1]), full_replay)
    # (5) a zero buffer left dirty: drop the last gate of U that writes a zero buffer
    n = 2 * c['v'] + c['R']
    zi = max(i for i, g in enumerate(c['U']) if g[0] == 'out' and any(t >= n for t, _, _ in g[3]))
    expect_fail('U: last write to a zero buffer dropped', lambda d: d['U'].pop(zi), full_replay)
    # (6) a twisted add with its warp sign flipped (unwarp before warp)
    wi0 = next(i for i, g in enumerate(c['D0']) if g[0] == 'warp')
    expect_fail('D0 warp sign flipped', lambda d: d['D0'][wi0].__setitem__(2, -d['D0'][wi0][2]), full_replay)
    return res


def width_inventory(word, origin):
    full = [1 << i for i in range(9)]
    w = [len(rref(word['helper_final'].get(str(r), full))) - len(rref(word['start_extra'].get(str(r), []))) for r in origin]
    return w


def main(a):
    if a[0] == '--to-1b':
        c = read(a[1])
        payload = json.dumps(to_1b(c), sort_keys=True, separators=(',', ':')).encode()
        open(a[2], 'wb').write(gzip.compress(payload, mtime=0))
        print(hashlib.sha256(payload).hexdigest())
        return
    t0 = time.monotonic()
    wp, outp, repp = a[0], a[1], a[2]
    opt = lambda k: a[a.index(k) + 1] if k in a else None
    raw = open(wp, 'rb').read()
    c0 = read(wp)
    rep = dict(canonical=wp, canonical_file_sha256=hashlib.sha256(raw).hexdigest(),
               canonical_content_sha256=hashlib.sha256(json.dumps(c0, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
    cert, info = lower(c0, opt('--provider'))
    rep['lowering'] = info
    chk = check_frames(cert)
    cert['N'], cert['cst'], cert['blocks'] = chk['N'], chk['cst'], chk['blocks']
    chk['hist_equals_hub_local_profile'] = chk['hist'] == cert['hub_local_profile']
    need(chk['hist_equals_hub_local_profile'], ('histogram differs from the hub profile', chk['hist'], cert['hub_local_profile']))
    rep['frames_and_price'] = chk
    rep['full_replay'] = full_replay(cert)
    rep['main_phase_words'] = main_phase_words(cert)
    cnt = lambda gs: sum(g[0] == 'warp' for g in gs)
    rep['counts'] = dict(R=cert['R'], nZ=cert['nZ'], exits=len(cert['exits']), ext=len(cert['ext']), tadds=len(cert['tadd']),
                         warp_events=dict(A=cnt(cert['A']), B=cnt(cert['B']), U=cnt(cert['U']), D0=cnt(cert['D0']),
                                          main_and_cleanup=cnt(cert['A']) + cnt(cert['B']) + cnt(cert['U'])),
                         zero_buffer_incidences=sum(1 for g in cert['U'] if g[0] == 'out' for r in regs_of(g) if r >= 2 * cert['v'] + cert['R']),
                         residuals_hist={str(k): x for k, x in sorted(Counter(cert['residuals']).items())})
    if opt('--word'):
        word = read(opt('--word'))
        hw = width_inventory(word, cert['helper_origin'])
        rep['width_inventory'] = dict(word=opt('--word'), word_file_sha256=hashlib.sha256(open(opt('--word'), 'rb').read()).hexdigest(),
                                      hist={str(k): x for k, x in sorted(Counter(hw).items())},
                                      equals_residuals_per_helper=hw == cert['residuals'])
        need(hw == cert['residuals'], 'residuals differ from the hub widths')
    if '--mutants' in a:
        rep['mutants'] = mutants(cert)
    payload = json.dumps(cert, sort_keys=True, separators=(',', ':')).encode()
    open(outp, 'wb').write(gzip.compress(payload, mtime=0))
    rep['out'] = outp
    rep['out_file_sha256'] = hashlib.sha256(open(outp, 'rb').read()).hexdigest()
    rep['out_content_sha256'] = hashlib.sha256(payload).hexdigest()
    if opt('--diff-1b'):
        theirs = read(opt('--diff-1b'))
        mine = to_1b(cert)
        dd = diff(mine, theirs)
        rep['diff_vs_g_heads'] = dict(file=opt('--diff-1b'), equal=not dd, differences=dd,
                                      our_1b_content_sha256=hashlib.sha256(json.dumps(mine, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
                                      their_content_sha256=hashlib.sha256(json.dumps(theirs, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
    rep['seconds'] = round(time.monotonic() - t0, 2)
    json.dump(rep, open(repp, 'w'), indent=1)
    print(json.dumps({k: x for k, x in rep.items() if k != 'mutants'}, indent=1))


if __name__ == '__main__':
    main(sys.argv[1:])
