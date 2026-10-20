#!/usr/bin/env python3
"""p-lower (round wht45/lever, 19 Oct 2026): OWN transport of a staged544 BORROW-WORD to the hub's native
canonical form (CANONICAL-WORD), written from the documented semantics of the hub's bw_transport.py (read, never
run).  The construction is PR2226's anticipatory dirty borrowing (Cinna Davis, @1337isnot1337); the transport stage
of the native chain is the hub's (GK / #433 lineage).  This file is our own implementation of the same map:

  * frame ids: frames are reduced-echelon GF(2)^9 bases (pivots descending); ids are given in order of first
    appearance, in the hub's order of appearance (ports, start_extra, y hyperplanes, helper_final, the warp table,
    precondition, compensation, A, cut, B, restore, cleanup, centre frames, target levels); 0 = bottom, 1 = FULL;
  * register compaction: the k removed helpers are dropped, registers renumbered in order (R' = R - k);
  * zero buffers: in the cleanup a register n0 + c (n0 = 2v + R) is the private zero buffer of restored helper c,
    written ['Z', new id of c];
  * warps: operator table (op = row + 1), warp events inserted into B before the gate at position B_pos (ascending
    op among equal positions); A-phase warps kept in place; cleanup warps become ['warp', r, op, 1, 1];
  * initial corrections -> twisted adds ['tadd', t, s, word, num, den, 0] in front of the precondition;
  * checks (all asserted): exact replay in the free warp algebra over every input column, both orientations
    (forward, inverse-transpose); dyadic primitives; frame chains (nested climbs, warps carry frames by G, zero-frame
    phase on frame-0 registers, centres, cleanup at FULL); the paid profile equals the word's report.

usage: transport.py BORROW-WORD.json[.gz] OUT-CANONICAL.json.gz OUT-AUDIT.json [--no-replay]
"""
import sys, json, gzip, hashlib, time
from collections import Counter
from fractions import Fraction as Q


def need(ok, msg):
    if not ok:
        raise AssertionError(msg)


def read_json(path):
    raw = open(path, 'rb').read()
    return json.loads(gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw), raw


# ------------------------------------------------------------------ GF(2)^9 subspaces (own code)
def rref(vectors):
    """reduced row echelon basis of the span, as a tuple with pivots (leading bits) descending"""
    rows = []                                    # kept reduced against each other
    for x in vectors:
        for b in rows:
            if x >> (b.bit_length() - 1) & 1:
                x ^= b
        if x:
            p = x.bit_length() - 1
            rows = [b ^ x if b >> p & 1 else b for b in rows]
            rows.append(x)
            rows.sort(reverse=True)
    return tuple(rows)


def inside(a, b):
    """span(a) <= span(b); b reduced echelon"""
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


# ------------------------------------------------------------------ free warp algebra
def axpy(dst, src, q):
    for key, val in src.items():
        w = dst.get(key, 0) + q * val
        if w:
            dst[key] = w
        else:
            dst.pop(key, None)


def warp_row(row, op):
    out = {}
    for (col, word), val in row.items():
        out[col, word[1:] if word and word[0] == -op else (op,) + word] = val
    return out


def dyadic(q):
    d = q.denominator
    return d & (d - 1) == 0


def transport(word):
    need(word.get('kind') == 'BORROW-WORD', 'not a BORROW-WORD')
    rep = word['report']
    v, h, R0 = rep['v'], rep['h'], rep['R']
    need(h == 9 and v == 112, 'sizes')
    n0 = 2 * v + R0
    removed = sorted(word['removed'])
    rem = set(removed)
    need(len(rem) == len(removed) and all(2 * v <= r < n0 for r in removed), 'removed helpers')
    kept = [r for r in range(n0) if r not in rem]
    new = {r: i for i, r in enumerate(kept)}
    n = len(kept)
    R = n - 2 * v
    need(R == R0 - len(removed), 'compaction')
    ports = word['ports']
    need(len(ports) == v == len(set(ports)), 'ports')
    FULL = rref([1 << i for i in range(h)])
    frames, fids = [(), FULL], {(): 0, FULL: 1}

    def fid(vs):
        b = rref(vs)
        if b not in fids:
            fids[b] = len(frames)
            frames.append(b)
        return fids[b]

    def reg(r):                                   # BORROW numbering -> canonical register (int) or ['Z', c]
        if isinstance(r, list):
            return list(r)
        if r >= n0:
            return ['Z', new[r - n0]]
        return new[r]

    # entrance frames
    start = [fid([p]) for p in ports] + [0] * (v + R)
    for r, f in word['start_extra'].items():
        r = int(r)
        need(2 * v <= r < n0 and r not in rem, 'start_extra register')
        start[new[r]] = fid(f)
    # final frames: y_t ends at the hyperplane orthogonal to its port
    final = [1] * v + [fid([x for x in range(1, 1 << h) if bin(x & p).count('1') % 2 == 0]) for p in ports] + [1] * R
    hfinal = {}
    for r, f in word['helper_final'].items():
        r = int(r)
        need(2 * v <= r < n0 and r not in rem, 'helper_final register')
        hfinal[new[r]] = fid(f)
        final[new[r]] = hfinal[new[r]]
    zero_buffers = sorted(hfinal)

    def conv(e):
        e = [list(x) if isinstance(x, list) else x for x in e]
        k = e[0]
        if k == 'warp':
            need(len(e) == 5, 'warp event shape')
            return ['warp', reg(e[1]), e[2], fid(e[3]), fid(e[4])]
        need(k in ('add', 'scale'), ('gate kind', k))
        e[1] = reg(e[1])
        if k == 'add':
            e[2] = reg(e[2])
        e[-1] = None if e[-1] is None else fid(e[-1])
        return e

    # warp operator table
    warps = []
    for i, (pos, lender, G, f_from, f_to) in enumerate(word['warps']):
        ff = fid(f_from)
        ft = fid(f_to)
        warps.append({'op': i + 1, 'B_pos': pos, 'lender': new[lender], 'G_columns': list(G),
                      'from_frame': ff, 'to_frame': ft})
    tadds = []
    for t in word['initial_corrections']:
        q = Q(t['coefficient'])
        tadds.append(['tadd', new[t['target']], new[t['source']], list(t['operator_word']), q.numerator, q.denominator, 0])
    precondition = tadds + [conv(e) for e in word['precondition']]
    compensation = [conv(e) for e in word['compensation']]
    A = [conv(e) for e in word['A']]
    cut = [conv(e) for e in word['cut']]
    at = {}
    for w in warps:
        at.setdefault(w['B_pos'], []).append(w)
    B = []
    for p in range(len(word['B']) + 1):
        for w in sorted(at.get(p, []), key=lambda w: w['op']):
            B.append(['warp', w['lender'], w['op'], w['from_frame'], w['to_frame']])
        if p < len(word['B']):
            B.append(conv(word['B'][p]))
    need(sum(e[0] == 'warp' for e in B) == len(warps), 'every warp placed in B')
    tail = [conv(e) for e in word['restore']]
    cleanup = []
    full_ends = 0
    for e in word['cleanup']:
        if e[0] == 'warp':
            need(len(e) == 4 and rref(e[3]) == FULL, 'cleanup warp not [warp, r, op, FULL]')
            full_ends += 1
            cleanup.append(['warp', reg(e[1]), e[2], 1, 1])
        else:
            cleanup.append(conv(e))
    centers = [new[r] for r in word['centers']]
    center_frames = [fid(f) for f in word['center_frames']]
    levels = {str(t): fid(f) for t, f in word['target_levels'].items()}
    c = dict(h=h, v=v, R=R, ports=ports, frames=[list(f) for f in frames], start=start, final=final,
             precondition=precondition, compensation=compensation, A=A, cut=cut, B=B, tail=tail, cleanup=cleanup,
             centers=centers, center_frames=center_frames, target_levels=levels, warps=warps,
             zero_buffers=zero_buffers, removed=removed, helper_origin=kept[2 * v:],
             borrow=dict(k=len(removed), base=word.get('base'), original_helpers=R0))
    info = dict(cleanup_warps_at_full=full_ends)
    return c, info


def expand_tadd(e):
    _, t, s, ops, a, b, f = e
    return ([['warp', s, op, 0, 0] for op in reversed(ops)] + [['add', t, s, a, b, 0]] +
            [['warp', s, -op, 0, 0] for op in ops])


def program(c):
    pre = []
    for e in c['precondition'] + c['compensation']:
        pre += expand_tadd(e) if e[0] == 'tadd' else [e]
    return [('pre', e) for e in pre] + [('A', e) for e in c['A']] + [('cut', e) for e in c['cut']] + \
        [('B', e) for e in c['B']] + [('tail', e) for e in c['tail']] + [('cleanup', e) for e in c['cleanup']]


def key_of(r):
    return ('Z', r[1]) if isinstance(r, list) else ('L', r)


def replay(c, prog, back):
    """exact replay in the free warp algebra; rows over (canonical column, warp word); returns the bad keys"""
    v, R = c['v'], c['R']
    n = 2 * v + R
    rows = {('L', r): {(r, ()): Q(1)} for r in range(n)}
    for z in c['zero_buffers']:
        rows[('Z', z)] = {}
    for part, e in prog:
        if e[0] == 'warp':
            k = key_of(e[1])
            rows[k] = warp_row(rows[k], e[2])
        elif e[0] == 'add':
            dst, src, a = key_of(e[1]), key_of(e[2]), e[3]
            if back:
                dst, src, a = src, dst, -a
            need(dst != src, 'add onto itself')
            axpy(rows[dst], rows[src], Q(a, e[4]))
        else:
            need(e[0] == 'scale', e[0])
            k = key_of(e[1])
            q = Q(e[3], e[2]) if back else Q(e[2], e[3])
            need(q != 0, 'zero scale')
            rows[k] = {j: q * x for j, x in rows[k].items()}
    want = {('L', r): {(r, ()): Q(1)} for r in range(n)}
    for t in range(v):
        if back:
            want[('L', t)][(v + t, ())] = Q(-1)
        else:
            want[('L', v + t)][(t, ())] = Q(1)
    for z in c['zero_buffers']:
        want[('Z', z)] = {}
    return [k for k in want if rows[k] != want[k]]


def frame_audit(c, word_profile):
    v, R, h = c['v'], c['R'], c['h']
    n = 2 * v + R
    fr = [tuple(f) for f in c['frames']]
    start, final = c['start'], c['final']
    levels = c['target_levels']
    lv = lambda t: levels.get(str(t), 0)
    zb = set(c['zero_buffers'])
    wrapped = {r for r in range(2 * v, n) if start[r] != 0 or final[r] != 1}
    G = {}
    for w in c['warps']:
        cols = w['G_columns']
        need(len(rref(cols)) == h, ('G not invertible', w['op']))
        a, b = fr[w['from_frame']], fr[w['to_frame']]
        need(len(a) == len(b) and rref([apply_G(cols, x) for x in a]) == b, ('G does not carry from -> to', w['op']))
        G[w['op']] = cols
    cur = list(start)
    paid, edges, wedges = Counter(), set(), set()

    def climb(r, f):
        o = cur[r]
        need(inside(fr[o], fr[f]), ('climb not nested', r, o, f))
        if o != f:
            paid[len(fr[f]) - len(fr[o])] += 1
            edges.add((o, f))
            cur[r] = f

    def warp(e):
        _, r, op, f0, f1 = e
        if isinstance(r, list):
            need(r[0] == 'Z' and f0 == f1 == 1, 'zero-buffer warp off FULL')
            return
        need(cur[r] == f0, ('warp frame mismatch', e, cur[r]))
        cols = G[abs(op)]
        a, b = (f0, f1) if op > 0 else (f1, f0)
        need(rref([apply_G(cols, x) for x in fr[a]]) == fr[b], ('warp G does not carry the frame', e))
        cur[r] = f1
        wedges.add((1 if op > 0 else -1, abs(op), f0, f1))

    lregs = lambda e: [r for r in ([e[1], e[2]] if e[0] == 'add' else [e[1]]) if not isinstance(r, list)]
    for e in c['precondition'] + c['compensation']:
        if e[0] == 'tadd':
            need(start[e[1]] == 0 and start[e[2]] == 0 and e[-1] == 0 and e[2] not in wrapped and all(e[3]),
                 ('twisted add', e))
            for x in expand_tadd(e):
                if x[0] == 'warp':
                    warp(x)
            continue
        need(e[-1] == 0 and all(start[r] == 0 for r in lregs(e)), ('zero-frame phase on a nonzero entrance', e))
    for e in c['A']:
        if e[0] == 'warp':
            warp(e)
        else:
            for r in lregs(e):
                climb(r, e[-1])
    need(all(cur[r] == f for r, f in zip(c['centers'], c['center_frames'])), 'centres not at their frames')
    need(all(cur[v + t] == lv(t) for t in range(v)), 'targets not at their levels')
    for r, f in zip(c['centers'], c['center_frames']):
        served = {lv(e[1] - v) for e in c['cut'] if e[2] == r}
        need(len({len(fr[s]) for s in served}) == len(served), 'tie in served levels (order would be ambiguous)')
        o = f
        for s in sorted(served, key=lambda s: -len(fr[s])):
            need(inside(fr[s], fr[o]), 'served level not inside the centre frame')
            if s != o:
                paid[len(fr[o]) - len(fr[s])] += 1
                edges.add((o, s))
                edges.add((s, o))
                o = s
    for e in c['B']:
        if e[0] == 'warp':
            warp(e)
        else:
            for r in lregs(e):
                climb(r, e[-1])
    for e in c['tail']:
        need(e[0] != 'warp', 'warp in the tail')
        for r in lregs(e):
            climb(r, e[-1])
    for r in range(n):
        climb(r, final[r])
    zreads = 0
    for e in c['cleanup']:
        if e[0] == 'warp':
            warp(e)
            continue
        for r in ([e[1], e[2]] if e[0] == 'add' else [e[1]]):
            if isinstance(r, list):
                need(r[0] == 'Z' and r[1] in zb and e[-1] == 1, ('zero buffer', e))
                zreads += 1
            else:
                need(final[r] == 1 and e[-1] == 1, ('cleanup on a non-FULL end', e))
    need(all(len(fr[final[r]]) < h for r in zb) and zb == {r for r in range(2 * v, n) if final[r] != 1}, 'zero buffers')
    prof = dict(sorted(paid.items()))
    need(prof == {int(k): x for k, x in word_profile.items()}, ('paid profile differs from the report', prof))
    return prof, sorted(edges), sorted(wedges), zreads


def main(argv):
    wp, outp, audp = argv[0], argv[1], argv[2]
    t0 = time.monotonic()
    word, raw = read_json(wp)
    c, info = transport(word)
    prog = program(c)
    aud = dict(word=wp, word_file_sha256=hashlib.sha256(raw).hexdigest(),
               word_content_sha256=hashlib.sha256(gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw).hexdigest())
    if '--no-replay' not in argv:
        for back in (False, True):
            bad = replay(c, prog, back)
            need(not bad, ('free-algebra endpoint', 'inverse-transpose' if back else 'forward', bad[:8]))
        aud['free_algebra_replay'] = 'both orientations exact on all %d columns' % (2 * c['v'] + c['R'])
    for part, e in prog:
        if e[0] == 'add':
            need(dyadic(Q(e[3], e[4])), ('non-dyadic add', e))
        elif e[0] == 'scale':
            q = Q(e[2], e[3])
            a = abs(q.numerator)
            need(dyadic(q) and a & (a - 1) == 0, ('scale not a signed power of two', e))
    prof, edges, wedges, zreads = frame_audit(c, word['report']['local_profile'])
    c['local_profile'] = prof
    c['edges'] = [list(x) for x in edges]
    c['warp_edges'] = [list(x) for x in wedges]
    payload = json.dumps(c, sort_keys=True, separators=(',', ':')).encode()
    open(outp, 'wb').write(gzip.compress(payload, mtime=0))
    fr = c['frames']
    aud.update(info, canonical=outp, canonical_file_sha256=hashlib.sha256(open(outp, 'rb').read()).hexdigest(),
               canonical_content_sha256=hashlib.sha256(payload).hexdigest(),
               v=c['v'], R=c['R'], original_helpers=c['borrow']['original_helpers'], removed=c['borrow']['k'],
               zero_buffers=len(c['zero_buffers']), zero_buffer_reads=zreads, warps=len(c['warps']),
               warp_events=dict(A=sum(e[0] == 'warp' for e in c['A']), B=sum(e[0] == 'warp' for e in c['B']),
                                cleanup=sum(e[0] == 'warp' for e in c['cleanup']),
                                tadd_expansion=sum(2 * len(e[3]) for e in c['precondition'] if e[0] == 'tadd')),
               tadds=sum(e[0] == 'tadd' for e in c['precondition']),
               scalar_operations=sum(e[0] != 'warp' for _, e in prog), local_profile={str(k): x for k, x in prof.items()},
               local_rank_mass=sum(k * x for k, x in prof.items()), frames=len(fr), directed_edges=len(edges),
               warp_edge_classes=len(wedges),
               widths_hist={str(k): x for k, x in sorted(Counter(len(fr[c['final'][r]]) - len(fr[c['start'][r]])
                                                                  for r in range(2 * c['v'], 2 * c['v'] + c['R'])).items())},
               seconds=round(time.monotonic() - t0, 2))
    json.dump(aud, open(audp, 'w'), indent=1)
    print(json.dumps(aud, indent=1))


if __name__ == '__main__':
    main(sys.argv[1:])
