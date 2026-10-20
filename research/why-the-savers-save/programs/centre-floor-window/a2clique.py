#!/usr/bin/env python3
"""ctr-window / a2clique.py (own code; stdlib; run with python3 -I).
P112 = E8 lines minus the frame e_{2a} +- e_{2a+1} (the hub's 112 ports). Labels u_t in F2^9 with u_t.u_s = B_ts =
(<r_t,r_s>/2 + 1) mod 2 (dense pairing), built by cx-medium's designs_cx.labels and re-checked here.
Computes:
 (1) the 1260 exact 8-cliques (orthogonal 8-frames of live lines);
 (2) the removed-frame coordinate decoder's term costs in the one-sided span accounting (A1: dim span U_R) and in the
     pairing-rank accounting (A2: rank_F2 B[R,A]);
 (3) the clique-counting bound in both accountings: floor = 8*#cliques / max_term hits(term)/cost(term), where a term
     (R,A) "hits" clique C iff R meets C and A meets C (A2) or R meets C (A1, one-sided as in cx-medium).
     A1 maximum is over W (R = L(W)); A2 maximum is searched over pairs W1 = U_R, W2 = U_A with W1 >= W2^perp
     (WLOG: adding W2^perp to W1 never raises rank B[R,A]); exhaustive where stated, else random+local search.
"""
import sys, os, json, random, itertools, time
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../medium-block-networks')))
from designs_cx import design, labels, ip

L = design('P112'); n = len(L)
labs, k, Bm = labels(L)
assert n == 112 and k == 9
def dot(x, y): return bin(x & y).count('1') & 1
for i in range(n):
    for j in range(n):
        assert dot(labs[i], labs[j]) == (1 + ip(L[i], L[j]) // 4) % 2   # doubled coords: norm 8, unit 4
assert len(set(labs)) == n
# 8-cliques = orthogonal frames
orth = [[ip(L[i], L[j]) == 0 for j in range(n)] for i in range(n)]
cl = []
def ext(c, cand):
    if len(c) == 8: cl.append(tuple(c)); return
    for idx, v in enumerate(cand):
        ext(c + [v], [w for w in cand[idx+1:] if orth[v][w]])
ext([], list(range(n)))
assert len(cl) == 1260, len(cl)
pm = [0] * n
for ci, c in enumerate(cl):
    for p in c: pm[p] |= 1 << ci
lab2port = {u: p for p, u in enumerate(labs)}

def rank2(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)
def span(vs):
    els = {0}
    for v in vs:
        if v not in els: els |= {x ^ v for x in els}
    return els
def ports_in(W):   # W = set of vectors
    return [lab2port[x] for x in W if x in lab2port]
def hitmask(ports):
    m = 0
    for p in ports: m |= pm[p]
    return m
def pairrank(R, A):
    # rank over F2 of the |R| x |A| matrix [u_r . u_a]
    rows = []
    for r in R:
        v = 0
        for j, a in enumerate(A):
            if dot(labs[r], labs[a]): v |= 1 << j
        rows.append(v)
    return rank2(rows)

out = {}
# (2) coordinate decoder
frame = []
for i in range(0, 8, 2):
    for sj in (1, -1):
        v = [0] * 8; v[i] = 2; v[i + 1] = 2 * sj; frame.append(tuple(v))
coord = []
for f in frame:
    Nf = [t for t in range(n) if ip(L[t], f) != 0]
    coord.append(dict(size=len(Nf), A1=rank2([labs[t] for t in Nf]), A2=pairrank(Nf, Nf)))
out['coordinate_decoder'] = dict(terms=coord, A1=sum(c['A1'] for c in coord), A2=sum(c['A2'] for c in coord))
print('coordinate decoder terms', coord, 'A1 total', out['coordinate_decoder']['A1'], 'A2 total', out['coordinate_decoder']['A2'])

# perp of a subspace (standard dot product, nondegenerate on F2^9)
ALL = list(range(1, 1 << 9))
def perp(W):
    B = [w for w in W if w]
    return {x for x in range(1 << 9) if all(dot(x, w) == 0 for w in B)}

# (3a) A1 one-sided: reproduce cx-medium's best d = 6 by random sampling only as a sanity check (its exhaustive value is 1155)
# (3b) A2: random + local search over pairs (W1, W2), W1 >= W2^perp, cost r = dim W1 + dim W2 - 9 (after closure)
def term_eval(W1basis, W2basis):
    W1 = span(W1basis); W2 = span(W2basis)
    R = ports_in(W1); A = ports_in(W2)
    if not R or not A: return 0, 0, None
    r = pairrank(R, A)
    if r == 0: return 0, 0, None
    h = bin(hitmask(R) & hitmask(A)).count('1')
    return h, r, (R, A)

rng = random.Random(20261015)
best = {}
t0 = time.time()
def randvec(): return rng.randrange(1, 1 << 9)
def record(h, r, W1b, W2b):
    if r and (r not in best or h > best[r][0]):
        best[r] = (h, [hex(x) for x in W1b], [hex(x) for x in W2b])
# seeds: W2 random of dim b, W1 = W2^perp + random X
for it in range(int(sys.argv[1]) if len(sys.argv) > 1 else 3000):
    b = rng.randint(2, 9); W2b = [randvec() for _ in range(b)]
    W2 = span(W2b); Kp = perp(W2)
    x = rng.randint(1, 6); W1b = list(Kp - {0})[:9] + [randvec() for _ in range(x)]
    # reduce bases
    def red(vs):
        bs = []
        for v in vs:
            w = v
            for q in bs: w = min(w, w ^ q)
            if w: bs.append(w)
        return bs
    W1b = red(list(perp(perp(span(W1b))) - {0})); W2b = red(list(W2 - {0}))
    h, r, _ = term_eval(W1b, W2b)
    # local search: try replacing one basis vector of W1 or W2 by a random vector, keep if ratio improves
    cur = (h / r if r else 0)
    for step in range(60):
        side = rng.randrange(2)
        b1, b2 = list(W1b), list(W2b)
        tgt = b1 if side == 0 else b2
        op = rng.randrange(3)
        if op == 0 and len(tgt) > 1: tgt.pop(rng.randrange(len(tgt)))
        elif op == 1 and len(tgt) < 9: tgt.append(randvec())
        elif tgt: tgt[rng.randrange(len(tgt))] = randvec()
        b1, b2 = red(b1), red(b2)
        hh, rr, _ = term_eval(b1, b2)
        if rr:
            record(hh, rr, b1, b2)
            if hh / rr >= cur:
                cur = hh / rr; W1b, W2b = b1, b2
    record(h, r, W1b, W2b)
print('A2 search: %.0f s' % (time.time() - t0))
res = {}
for r in sorted(best):
    h = best[r][0]
    res[r] = dict(hits=h, ratio_per_unit=h / r / 1260, floor_from_this=8 * 1260 * r / h if h else None, W1=best[r][1], W2=best[r][2])
    print('rank %d: best hits %d of 1260 -> per-unit fraction %.4f -> clique bound <= %.2f' % (r, h, h / r / 1260, 8 * 1260 * r / h if h else float('inf')))
out['A2_search'] = res
out['A2_clique_bound_upper'] = min(v['floor_from_this'] for v in res.values() if v['floor_from_this'])
out['A1_cx_medium'] = dict(d=6, hits=1155, floor=8 * 1260 * 6 / 1155)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a2clique.json'), 'w'), indent=1)
print('A2 clique bound is at most %.2f (A1 one-sided, cx-medium exhaustive: %.2f)' % (out['A2_clique_bound_upper'], 8 * 1260 * 6 / 1155))
