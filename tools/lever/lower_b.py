#!/usr/bin/env python3
"""g-heads (wht45/lever, 19 Oct 2026): lower a staged544 BORROW canonical word (the hub's native transport output,
bw_transport.py's CANONICAL-WORD.json.gz) into OUR certificate format, extended: gcert/1b.

OWN code. It follows bh-gap/k-lower's lower_x.py (the T74 lowering: copies of the 112 sources, units, exit frames,
the hub's tail lowered as the last gates of phase B, centres as retained totals, the cut as a scatter table) and adds
what the staged544 words of #3734 -> #3815 -> #3945 use and T74 did not:
  (w) WARPS: monomial Clifford reference changes K_G (PR2226 dirty borrowing, PR2436 lifted-target port). A warp is a
      free relabelling of one register's frame f0 -> f1 = G(f0); gcert/1b gate ['warp', f1, r, op, f0] (op < 0 = the
      inverse of warp |op|); the operator table `warps` (G by the images of the basis vectors) is carried along.
  (l) LIFTED TARGETS: 12 y targets are raised in phase A to a level frame (`levels`, here dim 1) by the 84-event prelude
      (12 zero-frame markers, then 12 groups warp x / 4 adds into lifted y / unwarp). gcert/1 forbids y in phase A
      (E2/E3); gcert/1b allows y in A for (a) the zero-frame y-only precondition P and (b) gates on lifted targets at
      their level. At the scatter every y sits at levels.get(t, 0).
  (p) PRECONDITION P stays where the hub has it (zero-frame phase, before the lifted prelude): it is emitted at the
      start of phase A (after the copies), and the scatter table is the hub's cut divided by the y units (no P^-1
      compensation as in lower_x, because P now precedes the cut). Reason: P touches all 12 lifted targets, so P
      after the cut (lower_x's place) would put frame-0 gates on y registers already at their level.
  Dropped exactly as lower_x drops them (the engine derives dirt and undo gates from the clean product): the
  compensation (zero-frame y += c s), the 36 initial corrections (zero-frame twisted additions y += c W(s) on dirty
  helper content, PR2226 order: they belong to the compensation) and the full cleanup with the private zero buffers
  ['Z', c] (the undo). Kept: removed-helper compaction (R = 408) as the hub transports it.
Then (as k-combine's partial_uncopy with keep = {}): the copies are folded back into the sources (mutable sources,
U8X0 route, no copies). Output: the no-copies gcert/1b (default) and, with --with-copies, also the copied one.
Exact checks (all asserted, results in the JSON report):
  frames nested at every climb, warp frames carried by G, block histogram (warps free) == the hub's native
  local_profile, invoice identity, residual per helper on the exit route == the hub's width (dim helper_final -
  dim start_extra of the BORROW-WORD / h - dim start), the total matrix of the clean main phase in the FREE warp
  algebra (rows {(column, warp word): Q}) and lower_x's exit conditions (R), (C) and Shea's scalar facts on it.
usage: lower_b.py CANONICAL-WORD.json.gz OUT.gcert1b.json.gz REPORT.json [--word BORROW-WORD.json.gz] [--with-copies P]
"""
import sys, json, gzip, hashlib, copy
from collections import Counter
from fractions import Fraction as Q


def read(p):
    raw = open(p, 'rb').read()
    return json.loads(gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw)


def need(c, msg):
    if not c:
        raise AssertionError(msg)


# ---------------------------------------------------------------- GF(2) frames (own code)
def red(u, basis):                       # basis reduced echelon, pivots descending
    for b in basis:
        if u >> (b.bit_length() - 1) & 1:
            u ^= b
    return u


def inside(a, b):
    return all(red(u, b) == 0 for u in a)


def echelon(vs):
    piv = {}
    for x in vs:
        for k in sorted(piv, reverse=True):
            if x >> k & 1:
                x ^= piv[k]
        if x:
            k = x.bit_length() - 1
            for j in piv:
                if piv[j] >> k & 1:
                    piv[j] ^= x
            piv[k] = x
    return [piv[k] for k in sorted(piv, reverse=True)]


def Gimg(cols, x):
    out = 0
    for i, c in enumerate(cols):
        if x >> i & 1:
            out ^= c
    return out


# ---------------------------------------------------------------- the lowering
def lower(c0):
    c = copy.deepcopy(c0)
    v, R, h = c['v'], c['R'], c['h']
    n0 = 2 * v + R
    rep = {}
    # 0. parts the engine derives (dirt / undo), dropped as lower_x drops compensation and cleanup
    pre = c['precondition']
    ntadd = 0
    while ntadd < len(pre) and pre[ntadd][0] == 'tadd':
        ntadd += 1
    need(all(e[0] != 'tadd' for e in pre[ntadd:]), 'initial corrections not first')
    for e in pre[:ntadd]:
        need(e[-1] == 0 and c['start'][e[1]] == 0 and c['start'][e[2]] == 0 and v <= e[1] < 2 * v and e[2] >= 2 * v,
             'initial correction is not a zero-frame y += c W(helper)')
    P = pre[ntadd:]
    need(all(e[-1] == 0 and all(v <= r < 2 * v for r in (e[1:3] if e[0] == 'add' else e[1:2])) for e in P),
         'precondition not zero-frame y-only')
    need(all(e[0] != 'tadd' for part in ('A', 'cut', 'B', 'tail') for e in c[part]), 'twisted add outside the zero-frame phase')
    need(all(not isinstance(r, list) for part in ('A', 'cut', 'B', 'tail') for e in c[part] for r in e[1:3] if e[0] == 'add')
         and all(not isinstance(e[1], list) for part in ('A', 'cut', 'B', 'tail') for e in c[part]),
         'a zero buffer outside the cleanup')
    rep['dropped'] = dict(initial_corrections=ntadd, compensation=len(c['compensation']), cleanup=len(c['cleanup']),
                          cleanup_warps=sum(e[0] == 'warp' for e in c['cleanup']), zero_buffers=len(c['zero_buffers']))
    # 1. copies of the sources (lower_x): every x reference -> copy register base + r
    base = n0
    remap = lambda r: base + r if r < v else r
    for part in ('A', 'B', 'cut', 'tail'):
        for e in c[part]:
            e[1] = remap(e[1])
            if e[0] == 'add':
                e[2] = remap(e[2])
    P = [list(e) for e in P]                       # y only: unchanged by the remap
    centres0 = list(c['centers'])
    c['centers'] = [remap(r) for r in c['centers']]
    copies = [['add', base + r, r, 1, 1, c['start'][r]] for r in range(v)]
    c['start'] += [0] * v
    c['final'] += [1] * v
    R += v
    n = 2 * v + R
    frames = c['frames']
    dims = [len(f) for f in frames]
    WT = {w['op']: w for w in c['warps']}
    levels = {int(t): f for t, f in c['target_levels'].items()}
    lifted = {v + t for t in levels}
    units = [Q(1)] * n
    A, B = [], []

    def emit(e, out):
        kind, r = e[:2]
        if kind == 'scale':
            out.append(['out', e[-1], r, [], []]); f = Q(e[2], e[3]); need(f, 'zero scale'); units[r] *= f
        elif kind == 'warp':
            _, r, op, f0, f1 = e
            need(abs(op) in WT, 'unknown warp op')
            out.append(['warp', f1, r, op, f0])
        else:
            need(kind == 'add', kind)
            _, t, s, num, den, f = e
            fac = Q(num, den) * units[s] / units[t]
            out.append(['out', f, s, [[t, fac.numerator, fac.denominator]], []])

    for e in copies:
        emit(e, A)
    for e in P:
        emit(e, A)
    ny_A = Counter()
    for e in c['A']:
        regs = [e[1]] + ([e[2]] if e[0] == 'add' else [])
        ys = [r for r in regs if v <= r < 2 * v]
        if ys:
            need(e[0] != 'warp' and all(r in lifted for r in ys) and e[-1] in (0, levels[ys[0] - v]),
                 ('phase A names a y register that is not a lifted target at its level', e))
            ny_A[e[0]] += 1
        emit(e, A)
    rep['phase_A_y_gates'] = dict(ny_A)
    ret = [[k, r, c['center_frames'][k]] for k, r in enumerate(c['centers'])]
    for _, r, f in ret:
        A.append(['out', f, r, [], []])
    centres = {r: k for k, r, _ in ret}
    S = [[Q(0)] * len(ret) for _ in range(v)]
    served = Counter()
    for e in c['cut']:
        _, t, s, num, den, f = e
        need(v <= t < 2 * v and s in centres and f is None, ('cut', e))
        S[t - v][centres[s]] += Q(num, den)
        served[levels.get(t - v, 0)] += 1
    need(set(served) == {0}, ('a centre serves a lifted target: the c-block accounting would change', dict(served)))
    scatter = []
    for t, row in enumerate(S):
        out = []
        for k, val in enumerate(row):
            val = val * units[ret[k][1]] / units[v + t]
            if val:
                out.append([k, val.numerator, val.denominator])
        scatter.append(out)
    for e in c['B']:
        emit(e, B)
    nB0 = len(B)
    for e in c['tail']:
        regs = e[1:3] if e[0] == 'add' else e[1:2]
        need(all(r >= 2 * v for r in regs), 'tail touches a data register')
        emit(e, B)
    for r, f in enumerate(c['final']):
        B.append(['out', f, r, [], []])
    bad = [[i, str(q)] for i, q in enumerate(units[:2 * v]) if q != 1]
    need(not bad, ('nonunit data output', bad[:5]))
    E = [r for r in range(2 * v, n) if c['final'][r] != 1]
    ext = [[r, f] for r, f in enumerate(c['start']) if r >= 2 * v and f]
    cert = dict(format='gcert/1b', h=h, v=v, R=R, ports=c['ports'], frames=frames, start=c['start'], final=c['final'],
                ext=ext, ret=ret, scat={'table': scatter}, A=A, B=B, exits=E,
                levels={str(t): f for t, f in sorted(levels.items())},
                warps=[dict(op=w['op'], G_columns=w['G_columns'], from_frame=w['from_frame'], to_frame=w['to_frame'])
                       for w in c['warps']],
                nP=len(P), helper_origin=c['helper_origin'], borrow=c['borrow'])
    rep.update(nA=len(A), nB_without_tail=nB0, nB=len(B), tail=len(c['tail']), nP=len(P),
               warp_gates=dict(A=sum(g[0] == 'warp' for g in A), B=sum(g[0] == 'warp' for g in B)))
    return cert, rep


def regs_of(g):
    if g[0] == 'warp':
        return [g[2]]
    if g[0] == 'out':
        return [g[2]] + [t for t, _, _ in g[3]] + list(g[4])
    raise AssertionError(g[0])


def uncopy(c):
    """k-combine partial_uncopy with keep = {} (verbatim logic), plus the warp gate: fold the copies back into the
    sources; source r then carries its copy's content, start frame start[r] and final frame final[copy]"""
    v, R = c['v'], c['R']; base = 2 * v + R - v
    cp = c['A'][:v]
    for r, g in enumerate(cp):
        need(g[0] == 'out' and g[2] == r and g[3] == [[base + r, 1, 1]] and g[4] == [] and g[1] == c['start'][r], ('copy', g))
    m = lambda q: q - base if q >= base else q

    def mg(g):
        g = copy.deepcopy(g)
        if g[0] == 'warp':
            g[2] = m(g[2])
        else:
            g[2] = m(g[2]); g[3] = [[m(t), a, b] for t, a, b in g[3]]; g[4] = [m(t) for t in g[4]]
        return g
    d = copy.deepcopy(c)
    d['A'] = [mg(g) for g in c['A'][v:]]
    d['B'] = [mg(g) for g in c['B']]
    d['ret'] = [[k, m(s), f] for k, s, f in c['ret']]
    d['ext'] = [[m(r), f] for r, f in c['ext']]
    d['exits'] = [m(r) for r in c['exits']]
    start = [None] * base; final = [None] * base
    for q in range(base):
        start[q] = c['start'][q]; final[q] = c['final'][q]
    for r in range(v):
        final[r] = c['final'][base + r]
    d['start'], d['final'] = start, final
    d['R'] = R - v
    d['format'] = 'gcert/1b'
    d['copies'] = 0
    return d


def check(c, hub_profile=None):
    """frames, warps, block histogram, scatter cut with levels, invoice identity; returns (hist by class, total, N)"""
    v, R, h = c['v'], c['R'], c['h']
    n = 2 * v + R
    fr = c['frames']; dim = [len(f) for f in fr]
    WT = {w['op']: w for w in c['warps']}
    for w in c['warps']:
        need(len(echelon(w['G_columns'])) == h, 'non-invertible G')
    levels = {int(t): f for t, f in c['levels'].items()}
    cur = list(c['start'])
    cls = lambda r: 'x' if r < v else 'y' if r < 2 * v else 's'
    H = {k: Counter() for k in 'xysc'}
    steps = set(); wsteps = set()

    def climb(r, f):
        o = cur[r]
        if o != f:
            if (o, f) not in steps:
                need(dim[o] < dim[f] and inside(fr[o], fr[f]), ('step not nested', r, o, f))
                steps.add((o, f))
            H[cls(r)][dim[f] - dim[o]] += 1
            cur[r] = f

    def warp(g):
        _, f1, r, op, f0 = g
        need(cur[r] == f0, ('warp frame mismatch', g, cur[r]))
        if (op, f0, f1) not in wsteps:
            G = WT[abs(op)]['G_columns']
            a, b = (f0, f1) if op > 0 else (f1, f0)
            need(echelon([Gimg(G, x) for x in fr[a]]) == list(fr[b]), ('G does not carry the frame', g))
            wsteps.add((op, f0, f1))
        cur[r] = f1
    for ph in ('A', 'B'):
        if ph == 'B':
            for k, s, f in c['ret']:
                need(cls(s) != 'y' and cur[s] == f, 'retained total not at its frame')
                H['c'][dim[f]] += 1
            need(all(cur[v + t] == levels.get(t, 0) for t in range(v)), 'a y role is not at its level at the scatter')
        for i, g in enumerate(c[ph]):
            rr = regs_of(g)
            need(all(0 <= r < n for r in rr) and len(set(rr)) == len(rr), ('registers', ph, i))
            if g[0] == 'warp':
                warp(g)
            else:
                for r in rr:
                    climb(r, g[1])
    for r in range(n):
        climb(r, c['final'][r])
    tot = Counter()
    for w in H.values():
        tot.update(w)
    N = sum(r * k for r, k in tot.items())
    lhs = N + sum(dim[f] for _, f in c['ext'])
    cst = sum(dim[f] for _, _, f in c['ret'])
    rhs = sum(dim[c['final'][r]] for r in range(2 * v, n)) + sum(dim[c['final'][r]] - dim[c['start'][r]] for r in range(v)) \
        + sum(dim[c['final'][v + t]] for t in range(v)) + cst
    # (x rows climb from their port frame (dim 1) to their final frame; y from 0; slots from their entrance)
    need(lhs == rhs, ('invoice identity', lhs, rhs))
    out = dict(hist_by_class={k: {str(r): x for r, x in sorted(w.items())} for k, w in H.items()},
               hist={str(r): x for r, x in sorted(tot.items())}, N=N, cst=cst, invoice=[lhs, rhs],
               frame_steps=len(steps), warp_steps=len(wsteps), max_block_rank=max(tot))
    if hub_profile is not None:
        out['hist_equals_hub_native_local_profile'] = out['hist'] == {str(k): x for k, x in sorted((int(a), b) for a, b in hub_profile.items())}
    return out


# ---------------------------------------------------------------- total matrix in the free warp algebra
def warp_row(row, op):
    out = {}
    for (var, word), val in row.items():
        word = word[1:] if word and word[0] == -op else (op,) + word
        out[var, word] = val
    return out


def total_matrix(c):
    v, R = c['v'], c['R']
    n = 2 * v + R
    reg = [{(r, ()): Q(1)} for r in range(n)]

    def axpy(t, src, co):
        d = reg[t]
        for k, u in src.items():
            w = d.get(k, 0) + co * u
            if w: d[k] = w
            else: d.pop(k, None)

    def run(gates):
        for g in gates:
            if g[0] == 'warp':
                reg[g[2]] = warp_row(reg[g[2]], g[3])
                continue
            src = dict(reg[g[2]])
            for t, a, b in g[3]:
                axpy(t, src, Q(a, b))
    run(c['A'])
    tot = {k: dict(reg[s]) for k, s, _ in c['ret']}
    for t in range(v):
        for k, a, b in c['scat']['table'][t]:
            axpy(v + t, tot[k], Q(a, b))
    run(c['B'])
    return reg


def exit_conditions(c, reg, copies):
    v, R = c['v'], c['R']
    n = 2 * v + R
    E = set(c['exits'])
    lo = 2 * v if copies else 0                        # slot rows: helpers (+ copies); no copies: x rows are slots too
    slot = lambda i: i >= 2 * v or (not copies and i < v)
    badR = [e for e in sorted(E) if set(reg[e]) - {(e, ())}]
    colsE = {}
    for i in range(n):
        for (j, w) in reg[i]:
            if j in E and (j != i or w):
                colsE.setdefault(j, []).append(i)
    badC = {e: [i for i in rows if slot(i)] for e, rows in colsE.items()}
    badC = {e: r for e, r in badC.items() if r}
    words = Counter(len(w) for row in reg for (_, w) in row)
    xrow_words = sum(1 for t in range(v) for (_, w) in reg[t] if w)
    yrow_words = sum(1 for t in range(v) for (_, w) in reg[v + t] if w)
    out = dict(R_rows_offdiag_zero=not badR, R_bad_count=len(badR), R_bad=badR[:10],
               C_cols_slot_offdiag_zero=not badC, C_bad_count=len(badC), C_bad=[[e, r[:6]] for e, r in list(badC.items())[:10]],
               exit_diagonals=dict(Counter(str(reg[e].get((e, ()), 0)) for e in E)),
               entries_by_warp_word_length=dict(sorted(words.items())), x_row_entries_with_words=xrow_words,
               y_row_entries_with_words=yrow_words)
    if copies:
        out['hx_x_rows_unit'] = all(reg[t] == {(t, ()): 1} for t in range(v))
    out['hid_y_from_x_identity'] = all({(j, w): u for (j, w), u in reg[v + t].items() if j < v} == {(t, ()): 1} for t in range(v))
    out['hy_y_columns_unit'] = all(not any(v <= j < 2 * v and (j != i or w) for (j, w) in reg[i]) for i in range(n)) and \
        all(reg[v + t].get((v + t, ())) == 1 for t in range(v))
    return out


def dirt_check(c0, d):
    """the hub's zero-frame phase in its own order: initial corrections (expanded as bw_transport's expand_tadd:
    warps, add, unwarps), the precondition P (already the first nP gates of our phase A), the compensation (its
    coefficients divided by the y units P leaves, our gates being written relative to the units)"""
    v = d['v']
    z = []
    pre = c0['precondition']
    k = 0
    while pre[k][0] == 'tadd':
        _, t, s, ops, a, b, f = pre[k]
        z += [['warp', 0, s, op, 0] for op in reversed(ops)] + [['out', 0, s, [[t, a, b]], []]] + [['warp', 0, s, -op, 0] for op in ops]
        k += 1
    u = Counter()
    unit = lambda r: u.get(r, Q(1))
    for e in pre[k:]:
        if e[0] == 'scale':
            u[e[1]] = unit(e[1]) * Q(e[2], e[3])
    comp = []
    for e in c0['compensation']:
        need(e[0] == 'add' and e[-1] == 0 and e[2] >= 2 * v, 'compensation')
        q = Q(e[3], e[4]) / unit(e[1])
        comp.append(['out', 0, e[2], [[e[1], q.numerator, q.denominator]], []])
    nP = d['nP']
    dd = dict(d); dd['A'] = z + d['A'][:nP] + comp + d['A'][nP:]
    reg = total_matrix(dd)
    bad = [t for t in range(v) if reg[v + t] != {(v + t, ()): 1, (t, ()): 1}]
    worded_y = sum(1 for t in range(v) for (_, w) in reg[v + t] if w)
    return dict(initial_correction_events=len(z), compensation=len(comp), y_rows_exactly_y_plus_x=not bad,
                bad_rows=bad[:10], worded_y_entries_left=worded_y)


def residuals(c):
    """exit route (k-combine make_sched_x / Lean res_ok): res[q] = dim(final) - dim(start) for an exit slot,
    h - dim(start) otherwise; slot q = helper q (register 2v+q) in helper_origin order"""
    v, R, h = c['v'], c['R'], c['h']
    dim = [len(f) for f in c['frames']]
    E = set(c['exits'])
    return [(dim[c['final'][2 * v + q]] if 2 * v + q in E else h) - dim[c['start'][2 * v + q]] for q in range(R)]


def hub_widths(word, origin):
    """the hub's width per helper (heads_fill.py / owner verify.py): dim helper_final - dim start_extra (GF(2) rank)"""
    full = [1 << i for i in range(9)]
    return [len(echelon(word['helper_final'].get(str(r), full))) - len(echelon(word['start_extra'].get(str(r), []))) for r in origin]


if __name__ == '__main__':
    a = sys.argv[1:]
    wp, outp, repp = a[0], a[1], a[2]
    wordp = a[a.index('--word') + 1] if '--word' in a else None
    cpp = a[a.index('--with-copies') + 1] if '--with-copies' in a else None
    raw = open(wp, 'rb').read()
    c0 = read(wp)
    rep = dict(canonical_word=wp, canonical_word_sha256=hashlib.sha256(raw).hexdigest(),
               canonical_content_sha256=hashlib.sha256(json.dumps(c0, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
    cert, info = lower(c0)
    rep['lowering'] = info
    rep['with_copies'] = check(cert, c0['local_profile'])
    regc = total_matrix(cert)
    rep['with_copies']['exit_conditions'] = exit_conditions(cert, regc, True)
    del regc
    d = uncopy(cert)
    rep['no_copies'] = check(d, c0['local_profile'])
    regd = total_matrix(d)
    rep['no_copies']['exit_conditions'] = exit_conditions(d, regd, False)
    del regd
    # the dirt the engine must derive: prepend the hub's dropped zero-frame phase (initial corrections, expanded as
    # warp / add / unwarp exactly as bw_transport's expand_tadd, then the compensation) to the no-copies word; every
    # y row must then be exactly y_t + x_t (the hub's endpoint; slots are left to the engine's end gate)
    rep['no_copies']['dirt_check'] = dirt_check(c0, d)
    res = residuals(d)
    rep['residuals'] = dict(R=d['R'], exits=len(d['exits']), hist={str(k): x for k, x in sorted(Counter(res).items())}, sum=sum(res),
                            exit_route_equals_net=res == [len(d['frames'][d['final'][2 * d['v'] + q]]) - len(d['frames'][d['start'][2 * d['v'] + q]])
                                                          for q in range(d['R'])])
    d['residuals'] = res
    if wordp:
        w = read(wp if False else wordp)
        hw = hub_widths(w, c0['helper_origin'])
        rep['residuals']['hub_word'] = wordp
        rep['residuals']['hub_word_sha256'] = hashlib.sha256(open(wordp, 'rb').read()).hexdigest()
        rep['residuals']['equal_hub_widths_per_helper'] = res == hw
        rep['residuals']['mismatches'] = [[q, res[q], hw[q]] for q in range(len(res)) if res[q] != hw[q]][:20]
    payload = json.dumps(d, sort_keys=True, separators=(',', ':')).encode()
    open(outp, 'wb').write(gzip.compress(payload, mtime=0))
    rep['out'] = outp
    rep['out_sha256'] = hashlib.sha256(open(outp, 'rb').read()).hexdigest()
    rep['out_content_sha256'] = hashlib.sha256(payload).hexdigest()
    if cpp:
        pc = json.dumps(cert, sort_keys=True, separators=(',', ':')).encode()
        open(cpp, 'wb').write(gzip.compress(pc, mtime=0))
        rep['with_copies_out'] = cpp
        rep['with_copies_out_sha256'] = hashlib.sha256(open(cpp, 'rb').read()).hexdigest()
    json.dump(rep, open(repp, 'w'), indent=1)
    print(json.dumps(rep, indent=1))
