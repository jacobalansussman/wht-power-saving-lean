#!/usr/bin/env python3
"""cx-medium / e7call.py (own code; numpy from our pylib; python3 -I).  An explicit crossing call (GF: y_a += x_a on
every port) for a root design, built here and checked EXACTLY on raw complex arrays at the window's block (m = k).
Parts: source lines (X_a starts at <u_a>), sink complements (Y_a ends at u_a^perp), the centre (the natural
coordinate decoder: totals T_q = sum_b (f_q . r_b) x_b gathered at the full frame by free gates, one fresh scratch
copy per total descending full -> 0 and scattered into the targets at frame 0, then discarded), junk cleaning by side
information:
  * FREE meetings: X_b climbs one level at a time along an orthonormal flag; at a frame G of its flag a fresh copy of
    x_b is taken (free); after the scatter a target Y_a whose own flag passes through the same G subtracts j_ab times
    it (free); the copy is discarded (zero moves: a zero array changes frame for free, a copy never moves);
  * PAID carriers for the rest: a copy of x_b taken at a frame F of X_b's flag climbs (paid, while holding content)
    to a frame H of Y_a's flag with F inside H, and every target whose flag passes H subtracts it.
The flags are chosen greedily, dimension by dimension, to maximise free meetings (bounded above by
v * sum_i min(2^(i-1), 2^(k-i-1)): a frame of dimension i holds at most 2^(i-1) odd labels).
Raw-array semantics (an2-smallnet frames.py): Phi_U = J diag(i^{w(proj_U x) mod 4}) J for a frame with an orthonormal
basis; every gate asserts one common frame; every step asserts nesting; a step of an array holding content costs its
dimension change.  Output: exact Y_a = Phi_{u_a^perp}(y_a + Phi_{u_a}^{-1} X_a(in)), X_a = Phi_full Phi_{u_a}^{-1} X_a(in)."""
import sys, os, json, random, itertools, time
pass  # numpy/scipy: install with pip (the original run used a local copy)
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from designs_cx import design, labels, ip, rank2, dot

def pc(x): return bin(x).count('1')

def span_key(basis):
    """canonical key of a subspace: its reduced row echelon basis (pivot = highest bit), sorted"""
    rows = []
    for v in basis:
        for r in rows:
            if (v >> (r.bit_length() - 1)) & 1: v ^= r
        if v:
            for i, r in enumerate(rows):
                if (r >> (v.bit_length() - 1)) & 1: rows[i] = r ^ v
            rows.append(v)
    return tuple(sorted(rows))

def in_span(v, basis):
    for x in sorted(span_basis(basis), reverse=True):
        v = min(v, v ^ x)
    return v == 0

def span_basis(vs):
    b = []
    for v in vs:
        for x in b: v = min(v, v ^ x)
        if v: b.append(v)
    return b

def build(name, seed=0, verbose=True):
    rng = random.Random(seed)
    L = design(name); labs, k, B = labels(L); n = len(L)
    full = (1 << k) - 1
    odd = [x for x in range(1, 1 << k) if pc(x) % 2]
    junk = {(a, b) for a in range(n) for b in range(n) if B[a][b] == 0}
    # natural decoder (coordinates as in natural.py)
    if name == 'E7':
        F = [tuple([1, -1] + [0] * 6)] + [tuple(1 if j == i else 0 for j in range(8)) for i in range(2, 8)]
    else:
        F = [tuple(1 if j == i else 0 for j in range(len(L[0]))) for i in range(len(L[0]))]
    G = np.array([[ip(f, l) for l in L] for f in F], dtype=float)
    Mmat = sum(np.outer(G[q], G[q]) / ip(F[q], F[q]) for q in range(len(F)))
    c0 = 1.0 / Mmat[0, 0]; Mmat = c0 * Mmat
    for a in range(n):
        for b in range(n):
            if B[a][b] == 1: assert abs(Mmat[a, b] - (1 if a == b else 0)) < 1e-12
    J = {(a, b): Mmat[a, b] for (a, b) in junk if abs(Mmat[a, b]) > 1e-12}
    # ---- flags: X_b basis starts with u_b (dims 1..k), Y_a basis inside u_a^perp (dims 0..k-1)
    Xb = [[labs[b]] for b in range(n)]
    Ya = [[] for _ in range(n)]
    Xfr = [{1: span_key(Xb[b])} for b in range(n)]
    Yfr = [{0: ()} for _ in range(n)]
    covered = set()
    for i in range(1, k):
        # candidate extensions: Y_a at dim i (needs a new odd c orthogonal to basis and to u_a), X_b at dim i (i >= 2)
        def completable(basis, extra_perp):
            # the orthogonal complement (inside the target space) of span(basis) must be 0 or contain an odd vector
            rest = [x for x in range(1, 1 << k) if all(dot(x, y) == 0 for y in basis + extra_perp)]
            return (not rest) or any(pc(x) % 2 for x in rest)
        def yext(a): return [c for c in odd if dot(c, labs[a]) == 0 and all(dot(c, x) == 0 for x in Ya[a]) and not in_span(c, Ya[a])
                             and completable(Ya[a] + [c], [labs[a]])]
        def xext(b): return [c for c in odd if all(dot(c, x) == 0 for x in Xb[b]) and not in_span(c, Xb[b]) and completable(Xb[b] + [c], [])]
        yopt = {a: {} for a in range(n)}
        for a in range(n):
            for c in yext(a): yopt[a].setdefault(span_key(Ya[a] + [c]), c)
        xopt = {b: {} for b in range(n)}
        if i == 1:
            for b in range(n): xopt[b] = {span_key(Xb[b]): None}
        else:
            for b in range(n):
                for c in xext(b): xopt[b].setdefault(span_key(Xb[b] + [c]), c)
        ydone, xdone = set(), set()
        while True:
            score = {}
            for a in range(n):
                if a in ydone: continue
                for g in yopt[a]:
                    score.setdefault(g, [set(), set()])[1].add(a)
            for b in range(n):
                if b in xdone: continue
                for g in xopt[b]:
                    if g in score: score[g][0].add(b)
            best, bs = None, 0
            for g, (xs, ys) in score.items():
                s = sum(1 for a in ys for b in xs if (a, b) in J and (a, b) not in covered)
                if s > bs or (s == bs and s > 0 and rng.random() < 0.3): best, bs = g, s
            if best is None or bs == 0: break
            xs, ys = score[best]
            for a in ys:
                Ya[a].append(yopt[a][best]); Yfr[a][i] = best; ydone.add(a)
                for b in xs:
                    if (a, b) in J: covered.add((a, b))
            for b in xs:
                if i >= 2 and b not in xdone:
                    Xb[b].append(xopt[b][best]); Xfr[b][i] = best; xdone.add(b)
        for a in range(n):
            if a not in ydone:
                opts = list(yopt[a].items()); g, c = rng.choice(opts)
                Ya[a].append(c); Yfr[a][i] = g
        if i >= 2:
            for b in range(n):
                if b not in xdone:
                    opts = list(xopt[b].items()); g, c = rng.choice(opts)
                    Xb[b].append(c); Xfr[b][i] = g
    # complete: X to full (last vector), Y's dim k-1 frame is u_a^perp automatically
    for b in range(n):
        while len(Xb[b]) < k:
            c = next(c for c in odd if all(dot(c, x) == 0 for x in Xb[b]) and not in_span(c, Xb[b]))
            Xb[b].append(c)
        Xfr[b][k] = span_key(Xb[b])
    for a in range(n):
        assert len(Ya[a]) == k - 1 and all(dot(c, labs[a]) == 0 for c in Ya[a])
    # free meetings (recount exactly from the flags)
    free = {}
    for (a, b) in J:
        for i in range(1, k):
            if Xfr[b].get(i) == Yfr[a][i] and Yfr[a][i] is not None:
                free[(a, b)] = i; break
    resid = [p for p in J if p not in free]
    # paid carriers: copy of x_b at X_b's frame of dim i climbs to Y's frame H of dim j >= i containing it.
    # greedy cover: for each b, candidate (i, H) -> targets served; cost j - i
    carriers = []
    left = set(resid)
    Ykeys = {}
    for a in range(n):
        for j in range(1, k):
            Ykeys.setdefault((j, Yfr[a][j]), set()).add(a)
    while left:
        bestc, bestv = None, -1
        bs = {b for (_, b) in left}
        for b in bs:
            for (j, H), tg in Ykeys.items():
                serve = {a for a in tg if (a, b) in left}
                if not serve: continue
                # cheapest source frame of X_b inside H
                Hb = list(H)
                best_i = max((i for i in range(1, min(j, k) + 1) if all(in_span(v, Hb) for v in Xfr[b][i])), default=None)
                if best_i is None: continue
                cost = j - best_i
                val = len(serve) / max(cost, 1e-9) if cost > 0 else 1e9
                if val > bestv: bestv, bestc = val, (b, best_i, j, H, serve, cost)
        b, i, j, H, serve, cost = bestc
        carriers.append(dict(b=b, i=i, j=j, H=H, targets=sorted(serve), cost=cost))
        left -= {(a, b) for a in serve}
    out = dict(design=name, n=n, k=k, junk_entries=len(J), free=len(free), resid=len(resid), carriers=len(carriers),
               carrier_levels=sum(c['cost'] for c in carriers),
               bound_free=n * sum(min(2 ** (i - 1), 2 ** (k - i - 1)) for i in range(1, k)))
    if verbose:
        print(f"{name}: junk {len(J)}; free meetings {len(free)} (bound {out['bound_free']}); residual {len(resid)} -> "
              f"{len(carriers)} paid carriers, {out['carrier_levels']} paid levels; free carriers zero-cost", flush=True)
    return dict(L=L, labs=labs, k=k, B=B, F=F, G=G, c0=c0, J=J, Xb=Xb, Ya=Ya, Xfr=Xfr, Yfr=Yfr, free=free,
                carriers=carriers, summary=out)

# ------------------------------------------------------------------ raw arrays
def walsh(m):
    n = 1 << m
    H = np.array([[(-1) ** dot(x, y) for y in range(n)] for x in range(n)], dtype=complex)
    return H / np.sqrt(n)

class Raw:
    def __init__(s, k, f=2, seed=0):
        s.k, s.N, s.J = k, 1 << k, walsh(k)
        s.rng = np.random.default_rng(seed); s.f = f
        s.A = {}; s.paid = 0; s.cache = {}; s.children = []
    def P(s, basis):
        key = span_key(basis)
        if key not in s.cache:
            bs = list(basis)
            for i, x in enumerate(bs):
                for j, y in enumerate(bs): assert dot(x, y) == (1 if i == j else 0)
            d = np.empty(s.N, dtype=complex)
            for x in range(s.N):
                p = 0
                for v in bs:
                    if dot(x, v): p ^= v
                d[x] = 1j ** (pc(p) % 4)
            s.cache[key] = s.J @ np.diag(d) @ s.J
        return s.cache[key]
    def new(s, r, basis, content=None):
        s.A[r] = [list(basis), np.zeros((s.N, s.f), complex) if content is None else content.copy()]
    def step(s, r, basis, who):
        old, c = s.A[r]
        lo, hi = (old, basis) if len(span_basis(old)) <= len(span_basis(basis)) else (basis, old)
        assert all(in_span(v, hi) for v in lo), 'not nested'
        d = abs(len(span_basis(basis)) - len(span_basis(old)))
        if np.abs(c).max() > 1e-12 and d:
            s.paid += d; s.children.append((who, d))
            c = s.P(basis) @ np.linalg.inv(s.P(old)) @ c
        s.A[r] = [list(basis), c]
    def gate(s, tgt, terms):
        F = span_key(s.A[tgt][0])
        acc = s.A[tgt][1].copy()
        for coef, src in terms:
            assert span_key(s.A[src][0]) == F, f'gate without a common frame {tgt} <- {src}'
            acc = acc + coef * s.A[src][1]
        s.A[tgt][1] = acc

def run_raw(D, seed=1):
    L, labs, k, n = D['L'], D['labs'], D['k'], len(D['L'])
    R = Raw(k, seed=seed)
    xin = [R.rng.normal(size=(R.N, 2)) + 1j * R.rng.normal(size=(R.N, 2)) for _ in range(n)]
    yin = [R.rng.normal(size=(R.N, 2)) + 1j * R.rng.normal(size=(R.N, 2)) for _ in range(n)]
    for b in range(n): R.new(('X', b), [labs[b]], xin[b])
    for a in range(n): R.new(('Y', a), [], yin[a])
    # copies needed: free meetings use copy (b, i) at X_b's frame of dim i; paid carriers copy (b, i) too
    need = {}
    for (a, b), i in D['free'].items(): need.setdefault((b, i), None)
    for ci, c in enumerate(D['carriers']): need.setdefault((c['b'], c['i']), None)
    # phase A: X_b climbs its flag one level at a time where a copy is needed, else in one block to the next need
    for b in range(n):
        stops = sorted(i for (bb, i) in need if bb == b)
        for i in stops + [k]:
            if i > len(span_basis(R.A[('X', b)][0])):
                R.step(('X', b), D['Xb'][b][:i], 'x')
            if i < k or (b, k) in need:
                R.new(('K', b, i), D['Xb'][b][:i]); R.gate(('K', b, i), [(1.0, ('X', b))])
    # totals at the full frame, copies descend to 0, scatter
    fullb = [1 << j for j in range(k)]
    for b in range(n): assert span_key(R.A[('X', b)][0]) == span_key(fullb)
    for q in range(len(D['F'])):
        R.new(('C', q), D['Xb'][0])          # the full frame (any orthonormal basis of it)
        R.gate(('C', q), [(D['G'][q][b], ('X', b)) for b in range(n) if D['G'][q][b] != 0])
        R.step(('C', q), [], 'centre')
        w = D['c0'] / ip(D['F'][q], D['F'][q])
        for a in range(n):
            if D['G'][q][a] != 0: R.gate(('Y', a), [(w * D['G'][q][a], ('C', q))])
    # paid carriers climb now (phase B start), each from its copy
    for ci, c in enumerate(D['carriers']):
        src = R.A[('K', c['b'], c['i'])]
        R.new(('P', ci), src[0], src[1]); R.step(('P', ci), yflag_basis(D, c), 'carrier')
    # phase B: Y_a climbs level by level along its flag; at each level subtract free copies and paid carriers
    car_at = {}
    for ci, c in enumerate(D['carriers']):
        for a in c['targets']: car_at.setdefault((a, c['j']), []).append((ci, c['b']))
    for a in range(n):
        for j in range(1, k):
            terms = []
            for (aa, b), i in D['free'].items():
                if aa == a and i == j: terms.append((-D['J'][(a, b)], ('K', b, i)))
            for ci, b in car_at.get((a, j), []): terms.append((-D['J'][(a, b)], ('P', ci)))
            if terms or j == k - 1:
                R.step(('Y', a), D['Ya'][a][:j], 'y')
                if terms: R.gate(('Y', a), terms)
    # checks
    ok = True
    Pfull = R.P(fullb)
    for a in range(n):
        Pa = R.P([labs[a]]); Pp = R.P(D['Ya'][a])
        want_y = Pp @ (yin[a] + np.linalg.inv(Pa) @ xin[a])
        want_x = Pfull @ np.linalg.inv(Pa) @ xin[a]
        ok &= np.allclose(R.A[('Y', a)][1], want_y) and np.allclose(R.A[('X', a)][1], want_x)
        ok &= span_key(R.A[('Y', a)][0]) == span_key([v for v in range(1, 1 << k) if dot(v, labs[a]) == 0])
    return ok, R

def yflag_basis(D, c):
    a = c['targets'][0]
    return D['Ya'][a][:c['j']]

if __name__ == '__main__':
    name = sys.argv[1] if len(sys.argv) > 1 else 'E7'
    best = None
    for seed in range(int(os.environ.get('SEEDS', '6'))):
        D = build(name, seed, verbose=False)
        s = D['summary']
        print(f"seed {seed}: free {s['free']} resid {s['resid']} carriers {s['carriers']} paid carrier levels {s['carrier_levels']}", flush=True)
        if best is None or s['carrier_levels'] < best['summary']['carrier_levels']: best = D
    s = best['summary']
    print(f"{name}: junk entries {s['junk_entries']}, free meetings {s['free']} (bound {s['bound_free']}), residual {s['resid']},"
          f" {s['carriers']} paid carriers, {s['carrier_levels']} paid carrier levels per call")
    t0 = time.time()
    ok, R = run_raw(best)
    from collections import Counter
    byw = Counter(); hist = Counter()
    for who, d in R.children: byw[who] += d; hist[(who, d)] += 1
    rec = dict(summary=s, raw_check_exact=bool(ok), paid_levels=R.paid, by_class=dict(byw),
               children={f'{w}:{d}': c for (w, d), c in sorted(hist.items())}, seconds=round(time.time() - t0, 1))
    print('raw-array check, block', best['k'], ': exact', ok, '| paid levels in the call', R.paid, dict(byw))
    print('children (class:rank -> count):', rec['children'])
    json.dump(rec, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'call_{name}.json'), 'w'), indent=1, default=str)
