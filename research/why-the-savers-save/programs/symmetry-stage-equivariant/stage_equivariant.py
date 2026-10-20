#!/usr/bin/env python3
"""stage_equivariant.py (d-I3-stage-equivariant, 13 Oct 2026).  Standard library only; run with python3 -I.

Question: can a total (scratch sum of several ports' sources) built in stage s of the hub's complex nine-stage
geometry be handed over to stage s+1 without a copy, for the hub's P112 decoder or any decoder in its symmetry orbit?

What this script does, all exact over F_2:
 (1) Rebuilds the hub's per-port data routes rho_{s,q} in O(73,2) from the hub's own specification
     (split-centres-global/proofs/INHERITED-COMPLEX-GEOMETRY.txt, "Compressed window and data routes") and the chart code
     the hub's checkers use (punctured112-global/code/audit_112_independent.py; #1042's t16/indep16.py, same lines).
     Cross-checks every port against the hub's DATA-PHASES.json: the idle-child orthogonal axes and the 3,840 stage
     translations.  Then A_q^(s) = rho_{s,q} rho_{s+1,q}: the stage-(s+1) invocation that continues the data piece of
     port q processed at stage s by invocation d is d * A_q^(s) (route semantics of the template:
     "data: physical cover class d*rho(stage,port)"; "work: private family at class d, erased before reuse").
 (2) Checks the two boundary identities (target end = next zero frame, full frame = next source line) for every port at
     all 8 boundaries, and the scalar hand-over identity at the 5 same-pair boundaries.
 (3) Hub decoder (8 totals R_a = live ports F2-orthogonal to the punctured label q_a): the continuation invocations of
     the 56 ports of each total, at each of the 5 hand-over boundaries; the frame tests of the stage map; the pre-load
     escape.
 (4) The general statement: for ANY route family satisfying the target-end identity, A_T(E_{s+1}) = q_T^perp + E_s,
     so distinct labels never share a continuation.  Checked exactly here for the hub routes, and by exhaustive
     enumeration over all routes in O(5,2) and O(7,2) on the small analogues (h = 3, 4) in small_case.py.
usage: python3 -I stage_equivariant.py [--hub-root DIR] [--k8-word FILE]
"""
import gzip, hashlib, itertools, json, os, sys
from collections import Counter

HUBDEF = __import__('os').environ.get('HUB', 'hub')  # a clone of CrocSwap/integer-mult-bounds
K8DEF = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/hub-inputs/K8-CANONICAL-WORD.json.gz'))
H, M, BL = 9, 73, 8          # local dimension, global dimension, exterior block size

def arg(name, default):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default

# ------------------------------------------------------------------ F_2 linear algebra on int bit-vectors
def dot(x, y): return bin(x & y).count('1') & 1
def reduce_basis(vs):
    basis = {}
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h in basis: v ^= basis[h]
            else: basis[h] = v; break
    return basis
def rank(vs): return len(reduce_basis(vs))
def span_eq(A, B): return rank(A) == rank(B) == rank(list(A) + list(B))
def contains(big, small): return rank(list(big) + list(small)) == rank(big)
def apply(cols, x):            # matrix given by images of e_0..e_{n-1}
    y = 0; i = 0
    while x:
        if x & 1: y ^= cols[i]
        x >>= 1; i += 1
    return y
def compose(A, B):             # (A o B)(x) = A(B(x))
    return [apply(A, c) for c in B]
def is_orthogonal(cols): return all(dot(cols[i], cols[j]) == int(i == j) for i in range(len(cols)) for j in range(i, len(cols)))
def ident(n): return [1 << i for i in range(n)]
def perp_local(q):             # basis of q^perp inside F_2^9
    return list(reduce_basis([x for x in range(1, 1 << H) if dot(x, q) == 0]).values())
def block(s): return [1 << (H + BL * (s - 1) + k) for k in range(BL)]       # exterior block s (1-based)
def E(s): return [v for j in range(1, s + 1) for v in block(j)]           # window base E_s
U = [1 << i for i in range(H)]

# ------------------------------------------------------------------ the hub's chart U_q and routes rho_{s,q}
def chart(q):
    """hub audit_112_independent.py / indep16.py: swap e0 with the first zero coordinate z, then the even transvection
    x -> x + (a0.x) a0 with a0 = q + e_z.  Column 0 is q; columns 1..8 are an orthonormal basis of q^perp (local)."""
    z = next(i for i in range(H) if not (q >> i) & 1)
    a0 = q ^ (1 << z); assert bin(a0).count('1') % 2 == 0
    order = list(range(H)); order[0], order[z] = order[z], order[0]
    B = [(1 << i) ^ (a0 if (a0 >> i) & 1 else 0) for i in order]
    assert B[0] == q and is_orthogonal(B)
    return B, z

def rho(s, q):
    """INHERITED-COMPLEX-GEOMETRY: rho_{0,q} = I; for s > 0, rho_{s,q} exchanges the complement (columns 1..8 of U_q)
    with exterior block s and fixes q and all other exterior coordinates.  Pairing column k <-> block coordinate k-1
    (the order of the hub's idle-child axes; pinned below by the hub's stage translations)."""
    if s == 0: return ident(M)
    B, _ = chart(q); bs = block(s)
    img = {}                                   # image of each basis vector of the adapted basis
    img[B[0]] = B[0]
    for k in range(1, H): img[B[k]] = bs[k - 1]; img[bs[k - 1]] = B[k]
    cols = []
    for i in range(M):
        e = 1 << i
        if i < H:                              # expand e_i in the orthonormal basis B: coefficient of B[k] = e_i . B[k]
            y = 0
            for k in range(H):
                if dot(e, B[k]): y ^= img[B[k]]
            cols.append(y)
        elif e in img: cols.append(img[e])
        else: cols.append(e)
    return cols

def is_involution(R): return compose(R, R) == ident(M)
def is_symmetric(R): return all(((R[j] >> i) & 1) == ((R[i] >> j) & 1) for i in range(M) for j in range(M))
def local_block(A): return [[(A[j] >> i) & 1 for j in range(H)] for i in range(H)]
def img_space(A, vs): return [apply(A, v) for v in vs]

# ------------------------------------------------------------------ main
def main():
    out = {}
    hub = arg('--hub-root', HUBDEF)
    php = os.path.join(hub, 'research/punctured112-global/vendor/pr514/inputs/work/complex-nine-native/DATA-PHASES.json')
    tpp = os.path.join(hub, 'research/punctured112-global/vendor/pr514/inputs/work/complex-nine-native/GLOBAL-TEMPLATE.json.gz')
    k8p = arg('--k8-word', K8DEF)
    ph = json.load(open(php)); tp = json.loads(gzip.decompress(open(tpp, 'rb').read()))
    k8 = json.loads(gzip.decompress(open(k8p, 'rb').read()))
    out['inputs'] = {p: hashlib.sha256(open(p, 'rb').read()).hexdigest() for p in (php, tpp, k8p)}
    print('inputs:'); [print('   %s  %s' % (h[:16], p)) for p, h in out['inputs'].items()]
    assert ph['m'] == M and tp['m'] == M
    print('template route semantics:', json.dumps(tp['route_semantics']))
    Q120 = [r['q'] for r in sorted(ph['recipes'], key=lambda r: r['port'])]
    live = list(k8['ports']); assert k8['v'] == 112 and k8['h'] == H and len(set(live)) == 112
    assert set(live) <= set(Q120)
    punct = [q for q in Q120 if q not in live]
    print('labels: DATA-PHASES 120 ports; #1042 k8 word 112 live ports (all distinct); punctured', len(punct),
          'at canonical indices', [Q120.index(q) for q in punct])
    assert [Q120.index(q) for q in punct] == [0, 13, 22, 35, 40, 51, 54, 119]
    assert all(bin(q).count('1') % 2 == 1 and q != (1 << H) - 1 for q in Q120)

    # ---- (1a) cross-check the charts against the hub's idle-child axes, and the routes against its translations
    ax_ok = tr_ok = 0; tr_fail = []
    for r in ph['recipes']:
        q = r['q']; B, z = chart(q)
        ch0 = r['idle_children'][0]
        if ch0['rank'] == 16 and ch0['orthogonal_axes'] == B[1:] + block(1): ax_ok += 1
        for s, al in enumerate(r['stage_alignment']):
            if s == 0:
                ok = al == {'entry_X': 0, 'entry_Y': 0, 'exit_X': 0, 'exit_Y': 0}
            else:
                ez = 1 << z; want_exit = ez ^ apply(rho(s, q), ez)
                ok = (al['entry_X'] == al['entry_Y'] == q ^ ez) and (al['exit_X'] == al['exit_Y'] == want_exit)
            if ok: tr_ok += 1
            else: tr_fail.append((r['port'], s))
    print('(1a) chart columns 1..8 + block 1 = hub idle-child-0 axes: %d / 120 ports' % ax_ok)
    print('     hub stage translations reproduced (s=0: zero; s>0: entry = q+e_z, exit = e_z + rho_s(e_z)): %d / 1080'
          % tr_ok, '' if not tr_fail else 'failures %s' % tr_fail[:5])
    out['chart_axes_match'] = ax_ok; out['translations_match'] = tr_ok
    assert ax_ok == 120 and tr_ok == 1080

    # ---- (1b) rho: orthogonal symmetric involutions; A_q^(s) = rho_s rho_{s+1}
    R = {(s, q): rho(s, q) for s in range(9) for q in live}
    assert all(is_orthogonal(R[k]) and is_involution(R[k]) and is_symmetric(R[k]) for k in R)
    print('(1b) rho_{s,q}: 9 x 112 matrices in O(73,2), all symmetric involutions')

    # ---- (2) boundary identities at all 8 boundaries, every port (frames in the piece's own class coordinates)
    bad = []
    for s in range(8):
        for q in live:
            qp = perp_local(q)
            H_s = img_space(R[s, q], U + E(s));          P_n = img_space(R[s + 1, q], [q] + E(s + 1))
            K_s = img_space(R[s, q], qp + E(s));         Z_n = img_space(R[s + 1, q], E(s + 1))
            if not (span_eq(H_s, P_n) and span_eq(K_s, Z_n) and span_eq(K_s, qp + E(s)) and span_eq(H_s, U + E(s))):
                bad.append((s, q))
    print('(2) boundary identities H_s = P_(s+1) and K_s = 0_(s+1) [= data levels X_s = U+E_s, Y_s = q^perp+E_s]: '
          '%d of %d (boundary, port) pairs fail' % (len(bad), 8 * 112))
    assert not bad
    # scalar hand-over identity: the template's recipe equals stagemap.py's B_K recipe; replay it
    sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-stage-hand-over')))
    import stagemap
    assert tp['recipe'][-1] == {'kind': 'terminal', 'pairs': 4}          # the template's name for the signed exchange
    rec = [dict(kind=e['kind'], **({'stage': e['stage'], 'pair': e['pair'], 'reverse': e['reverse']} if e['kind'] == 'local'
                                 else {'dst': e['dst'], 'src': e['src'], 'c': e['coefficient']}))
           for e in tp['recipe'][:-1]] + [{'kind': 'exchange'}]
    assert rec == stagemap.recipe(4)
    log, idw = stagemap.scalar(4)
    HB = [(si, sj) for si, sj, p, *_ in log]
    print('    scalar replay of the template recipe: word = identity: %s; same-bank hand-over boundaries %s' % (idw, HB))
    print('    (content exact at %s; a bridge term in between at %s)' % ([(a, b) for a, b, p, oi, oj, sb, ex, xt in log if ex],
                                                                         [(a, b) for a, b, p, oi, oj, sb, ex, xt in log if not ex]))
    assert HB == [(0, 1), (1, 2), (3, 4), (5, 6), (7, 8)]
    out['handover_boundaries'] = HB

    # ---- (1c) the table of A_q^(s) at the 5 hand-over boundaries
    A = {(s, q): compose(R[s, q], R[s + 1, q]) for s, _ in HB for q in live}
    table = []
    for s, t in HB:
        mats = [tuple(A[s, q]) for q in live]
        distinct = len(set(mats))
        loc_ok = all(local_block(A[s, q]) == [[(q >> i) & (q >> j) & 1 for j in range(H)] for i in range(H)] for q in live)
        base_ok = all(span_eq(img_space(A[s, q], E(s + 1)), perp_local(q) + E(s)) for q in live)
        fix_ok = all(apply(A[s, q], q) == q for q in live)
        order = Counter()
        for q in live:
            X, k = A[s, q], 1
            while X != ident(M): X = compose(X, A[s, q]); k += 1
            order[k] += 1
        moved = Counter(sum(1 for i in range(M) if A[s, q][i] != 1 << i) for q in live)
        row = dict(boundary='%d->%d' % (s, t), distinct_A=distinct, local_block_is_qqT=loc_ok, A_maps_E_next_to_qperp_plus_E=base_ok,
                   A_fixes_q=fix_ok, orders=dict(order), coordinates_moved=dict(moved))
        table.append(row)
        print('(1c) boundary %d->%d: A_q = rho_%d,q rho_%d,q: %d distinct over 112 ports; 9x9 local block = q q^T: %s; '
              'A_q(E_%d) = q^perp + E_%d: %s; A_q q = q: %s; orders %s; coordinates moved %s'
              % (s, t, s, t, distinct, loc_ok, t, s, base_ok, fix_ok, dict(order), dict(moved)))
        assert distinct == 112 and loc_ok and base_ok and fix_ok
    out['A_table_summary'] = table
    # write every A_q^(s) explicitly (columns = images of e_0..e_72, as hex) for the record
    full = {'%d->%d' % (s, t): {str(q): [format(c, 'x') for c in A[s, q]] for q in live} for s, t in HB}
    json.dump(full, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'A_table.json'), 'w'))
    print('     all 5 x 112 matrices written to out/A_table.json (column j = A(e_j), hex)')

    # ---- (3) the hub decoder: 8 totals R_a
    totals = []
    for qa in punct:
        Ra = [q for q in live if dot(q, qa) == 0]
        Pa = list(reduce_basis(Ra).values())
        assert span_eq(Pa, perp_local(qa))
        totals.append((qa, Ra, Pa))
    print('(3) hub P112 decoder: %d totals, read sizes %s, widths %s (star cost %d); P_a = q_a^perp for the punctured q_a'
          % (len(totals), [len(t[1]) for t in totals], [rank(t[2]) for t in totals], sum(rank(t[2]) for t in totals)))
    print('    punctured labels pairwise dot:', sorted(Counter(dot(a, b) for a, b in itertools.combinations(punct, 2)).items()))
    res = []
    for s, t in HB:
        transportable = 0; conts = []
        frame_in_window = preload_legal = 0
        for qa, Ra, Pa in totals:
            cont = len({tuple(A[s, q]) for q in Ra})        # distinct continuation invocations d*A_q (right action is free)
            conts.append(cont)
            if cont == 1: transportable += 1
            for q in Ra:
                # stage-map frame test: in the continuation invocation of q, is the total's frame (P_a + E_s, invocation-d
                # coordinates) inside the stage-(s+1) window [E_(s+1), U + E_(s+1)]?  In d*A_q coordinates it is A_q^-1(P_a+E_s).
                F = Pa + E(s)
                # A_q^-1 (F) contains E_(s+1)  <=>  F contains A_q(E_(s+1)) = q^perp + E_s
                if contains(F, perp_local(q) + E(s)): frame_in_window += 1
            for q in live:
                # pre-load escape: add the total into target y_q, which must end inside q^perp + E_s: needs P_a inside q^perp
                if contains(perp_local(q), Pa): preload_legal += 1
        # data test of the pre-load escape: host q_S, read port q_T: is q_T's piece in the host's continuation invocation?
        cop = max(sum(1 for q in Ra if A[s, q] == A[s, qs]) for qa, Ra, Pa in totals for qs in live)
        res.append(dict(boundary='%d->%d' % (s, t), transportable=transportable, continuation_invocations_per_total=conts,
                        total_frame_contains_new_base_pairs=frame_in_window, preload_frame_legal_pairs=preload_legal,
                        max_pieces_of_a_total_copresent_in_a_host_continuation=cop))
        print('    boundary %d->%d: transportable totals %d / 8; distinct continuation invocations per total %s; '
              '(total, port) pairs whose frame reaches the new window %d / 448; frame-legal pre-loads %d / 896; '
              'most pieces of one total co-present in any host continuation %d'
              % (s, t, transportable, conts, frame_in_window, preload_legal, cop))
        assert transportable == 0 and frame_in_window == 0 and preload_legal == 0 and cop == 1
    out['hub_decoder'] = res

    # ---- (4) the general statement, on the hub's labels: for ANY routes g_s, g_{s+1} with
    #      g_s^-1(q^perp + E_s) = g_{s+1}^-1(E_{s+1}), A = g_s g_{s+1}^-1 satisfies A(E_{s+1}) = q^perp + E_s.
    #      The map q -> q^perp + E_s is injective (intersect with U, take perp), checked here on all 120 labels:
    keys = {frozenset(reduce_basis(perp_local(q) + E(3)).items()) for q in Q120}
    print('(4) q -> q^perp + E_s injective on all 120 labels (s = 3 shown): %d distinct subspaces' % len(keys))
    assert len(keys) == 120
    # orbit statement: any port relabelling (octad stabiliser, W(E8), free relabellings of the decoder) permutes labels;
    # a total is transportable iff all its ports share one label, i.e. iff it reads one port.
    out['orbit_best'] = 0
    print('    => transportable totals of ANY decoder (any orbit element, any replica, any of the 5 boundaries) = totals that')
    print('       read a single port.  Hub decoder 0 / 8; best over its orbit 0 / 8; a decoder with only such totals costs >= 112.')
    # ---- (5) paid re-routing: co-locate a total's ports by re-framing their pieces at the boundary.  Generous cost
    #      model (copyfloor H9 relaxed: a non-nested step f -> g costs the subspace distance dim f + dim g - 2 dim(f cap g);
    #      the hub's strict rule forbids it outright).  One common continuation A needs every Y piece at the one frame
    #      A(E_(s+1)) and the X pieces at pairwise distinct frames A(<q_T> + E_(s+1)) (distinct labels), all of the same
    #      dimension as the frames they start at; two distinct same-dimension subspaces are at distance >= 2.
    lb = []
    for qa, Ra, Pa in totals:
        Ys = [frozenset(reduce_basis(perp_local(q)).items()) for q in Ra]          # Y frames relative to the common E_s
        best_Y = min(sum(0 if span_eq(perp_local(q), perp_local(c)) else 2 for q in Ra) for c in range(1, 1 << H))
        best_X = 2 * (len(Ra) - 1)          # X pieces all start at U + E_s and must end at |Ra| distinct frames
        lb.append(best_Y + best_X)
    print('(5) paid co-location of a total\'s ports (generous non-nested cost): >= %s moves per total per boundary, '
          'against the 8 its star copy costs' % sorted(set(lb)))
    out['paid_colocation_lower_bound_per_total'] = sorted(set(lb))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', 'summary.json'), 'w'), indent=1)
    print('PASS  (out/summary.json, out/A_table.json)')

if __name__ == '__main__':
    main()
