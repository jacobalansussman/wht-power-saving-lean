#!/usr/bin/env python3
"""d-halfread STEP 1: exact small case for pre-merge half reads with r = 1..4 consumers (generalises g-small's
halfread_toy.py). Ambient F2^6, packet plane P = <e1,e2,e3> with labels u1=e1, u2=e2, u3=e3, u4=e1+e2+e3.
Packet p = 1*x1 + 2*x2 + 3*x3 + 5*x4.  Halves: alpha = x1 + 2*x2 on L_a = <e1,e2> (carries p's own coefficients
1:2 on that line), beta = 3*x3 + 5*x4 on L_b = <e3,u4>.  Merge a += b at P (k = 1).
Two shapes of consumers:
  I (incomparable): consumer caps G_i = P + <g_i> pairwise meeting exactly in P (g = e4, e5, e6, e4+e5+e6);
  N (U8 copy shape): one consumer on a's branch (cap P+e4) and a nested pair (P+e5 inside P+<e5,e6>) that the
    two-copy design serves from a second copy h3.
Consumers are targets (start 0, end at their cap) or relaying helpers (start 0, restored).  All helpers and targets
carry symbolic dirt; helper dirt is cancelled by differential reads at frame 0.  Checks: exact endpoint over Q
(symbolic sources and dirt), every gate at a frame containing both registers' frames, monotone chains, targets end
inside their caps, helpers and sources restored, and which registers ever sit at P.
Also an exhaustive frame-level check of the branch lemma: in shape I no two consumers can read one post-merge
register (their read frames would have to nest inside G_i and G_j, whose meet is P).
Usage: python3 -I halfread_r.py [out.json]"""
import sys, json, itertools
from fractions import Fraction as Fr
E = [1 << i for i in range(6)]
def span(vs):
    bb = []
    for y in vs:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return bb
def inside(H, G):
    for x in H:
        y = x
        for z in span(G): y = min(y, y ^ z)
        if y: return False
    return True
def same(G, H): return inside(G, H) and inside(H, G)
P = [E[0], E[1], E[2]]; La = [E[0], E[1]]; Lb = [E[2], E[0] ^ E[1] ^ E[2]]
LAB = {'x1': E[0], 'x2': E[1], 'x3': E[2], 'x4': E[0] ^ E[1] ^ E[2]}
COEF = {'x1': 1, 'x2': 2, 'x3': 3, 'x4': 5}
class Prog:
    def __init__(s, regs):
        s.regs = regs; s.val = {}; s.frame = {}; s.chain = {}
        for r in regs:
            s.val[r] = {('x' if r in LAB else 'd', r): Fr(1)}
            s.frame[r] = [LAB[r]] if r in LAB else []; s.chain[r] = [list(s.frame[r])]
        s.ok = True; s.gates = 0
    def add(s, a, b, k, G):
        assert inside(s.frame[a], G) and inside(s.frame[b], G), ('frame', a, b, G)
        for r in (a, b):
            if not same(G, s.frame[r]): s.frame[r] = list(G); s.chain[r].append(list(G))
        for key, x in s.val[b].items():
            y = s.val[a].get(key, 0) + k*x
            if y: s.val[a][key] = y
            else: s.val[a].pop(key, None)
        s.gates += 1
    def monotone(s):
        return all(inside(c[i], c[i+1]) and not same(c[i], c[i+1]) for c in s.chain.values() for i in range(len(c)-1))
    def at_P(s): return sorted(r for r in s.regs if any(same(f, P) for f in s.chain[r]))
    def dims(s, r): return [len(span(f)) for f in s.chain[r]]
FULL = E[:]
def build(consumers, halfread_ids, m_extra):
    """consumers: list of (name, cap, branch) where branch is 'a', 'b' (needs a half read) or 'c<j>' (j-th copy).
    m_extra: number of copies at P.  Helper dirt reaching a consumer is cancelled by compensation reads at frame 0,
    derived from a first pass (as the hub's audit derives its compensation). Returns the checked program summary."""
    copies = ['c%d' % j for j in range(m_extra)]
    regs = list(LAB) + [n for n, _, _ in consumers] + ['a', 'b'] + copies
    helpers = ['a', 'b'] + copies
    gl = []
    gl += [('a', 'x1', 1, [E[0]]), ('a', 'x2', 2, La), ('b', 'x3', 3, [E[2]]), ('b', 'x4', 5, Lb)]   # halves
    gl += [(n, 'a', 1, La) for n, cap, br in consumers if br == 'b']                                   # pre-merge half reads at L_a
    gl += [('a', 'b', 1, P)]                                                                           # the merge
    gl += [(c, 'a', 1, P) for c in copies]                                                             # second copies at P
    gl += [(n, br, 1, cap) for n, cap, br in consumers]                                                # post-merge reads
    gl += [(c, 'a', -1, FULL) for c in copies]                                                         # restore
    gl += [('a', 'b', -1, FULL), ('a', 'x1', -1, FULL), ('a', 'x2', -2, FULL), ('b', 'x3', -3, FULL), ('b', 'x4', -5, FULL)]
    g0 = Prog(regs)
    for x in gl: g0.add(*x)
    comp = [(n, key[1], -c, []) for n, _, _ in consumers for key, c in g0.val[n].items() if key[0] == 'd' and key[1] != n]
    g = Prog(regs)
    for x in comp + gl: g.add(*x)
    want = {('x', s): Fr(k) for s, k in COEF.items()}
    res = {}
    for n, cap, br in consumers:
        exp = dict(want); exp[('d', n)] = Fr(1)
        res[n] = g.val[n] == exp
    restored = all(g.val[h] == {('d', h): Fr(1)} for h in helpers) and all(g.val[x] == {('x', x): Fr(1)} for x in LAB)
    caps_ok = all(inside(g.frame[n], cap) for n, cap, _ in consumers)
    return dict(consumers_get_p=all(res.values()), helpers_sources_restored=restored, monotone=g.monotone(),
                final_frames_in_caps=caps_ok, registers_at_P=g.at_P(), n_at_P=len(g.at_P()),
                consumer_chain_dims={n: g.dims(n) for n, _, _ in consumers}, gates=g.gates, compensation_reads=len(comp),
                reads_per_consumer={n: sum(1 for x in gl if x[0] == n) for n, _, _ in consumers},
                half_reads=sum(1 for _, _, br in consumers if br == 'b'))
CAPS_I = [P + [E[3]], P + [E[4]], P + [E[5]], P + [E[3] ^ E[4] ^ E[5]]]
out = {'shape_I': {}, 'shape_N': {}, 'branch_lemma': {}}
for r in range(1, 5):
    caps = CAPS_I[:r]
    # with half reads: consumer 1 on a, the last on b (half read), the rest on copies
    if r == 1: plan = [('t1', caps[0], 'a')]; mex = 0
    else:
        plan = [('t1', caps[0], 'a')] + [('t%d' % (i+1), caps[i], 'c%d' % (i-1)) for i in range(1, r-1)] + [('t%d' % r, caps[r-1], 'b')]
        mex = r - 2
    withh = build(plan, None, mex)
    # two-copy interface (no half reads): b idle after the merge, consumers 2..r on copies
    plan2 = [('t1', caps[0], 'a')] + [('t%d' % (i+1), caps[i], 'c%d' % (i-1)) for i in range(1, r)]
    without = build(plan2, None, r - 1)
    out['shape_I'][r] = dict(with_half_reads=withh, two_copy_interface=without)
    print('shape I r=%d: with half reads -> %d registers at P %s, checks %s; without -> %d registers at P, checks %s' % (
        r, withh['n_at_P'], withh['registers_at_P'], all(withh[k] for k in ('consumers_get_p', 'helpers_sources_restored', 'monotone', 'final_frames_in_caps')),
        without['n_at_P'], all(without[k] for k in ('consumers_get_p', 'helpers_sources_restored', 'monotone', 'final_frames_in_caps'))))
    print('      consumer chain dims (with half reads):', withh['consumer_chain_dims'])
# shape N: t_A on a (cap P+e4); t_H1 (cap P+<e5,e6>, reads at P+e5 first... ) and t_H2 (cap P+<e5,e6>+... nested)
GA = P + [E[3]]; GH1 = P + [E[4]]; GH2 = P + [E[4], E[5]]
planN_two_copy = [('tA', GA, 'a'), ('tH1', GH1, 'c0'), ('tH2', GH2, 'c0')]
planN_half = [('tA', GA, 'a'), ('tH1', GH1, 'b'), ('tH2', GH2, 'b')]
nt = build(planN_two_copy, None, 1); nh = build(planN_half, None, 0)
out['shape_N'] = dict(two_copy_interface=nt, with_half_reads=nh)
for nm, x in (('two-copy', nt), ('half reads', nh)):
    print('shape N (%s): registers at P %s; checks %s; consumer dims %s' % (nm, x['registers_at_P'],
          all(x[k] for k in ('consumers_get_p', 'helpers_sources_restored', 'monotone', 'final_frames_in_caps')), x['consumer_chain_dims']))
# exhaustive branch lemma on shape I (all subspaces of F2^6 containing P inside each cap)
def subspaces_between(lo, hi):
    extra = [v for v in range(64) if inside([v], hi) and not inside([v], lo)]
    seen = set(); res = []
    for k in range(0, 4):
        for comb in itertools.combinations(extra, k):
            b = span(lo + list(comb)); key = frozenset(_closure(b))
            if key not in seen and inside(b, hi): seen.add(key); res.append(b)
    return res
def _closure(b):
    s = {0}
    for x in b: s |= {y ^ x for y in s}
    return s
viol = 0; pairs = 0
for i, j in itertools.combinations(range(4), 2):
    Fi = [f for f in subspaces_between(P, CAPS_I[i]) if not same(f, P)]
    Fj = [f for f in subspaces_between(P, CAPS_I[j]) if not same(f, P)]
    for f in Fi:
        for h in Fj:
            pairs += 1
            if inside(f, h) or inside(h, f): viol += 1
out['branch_lemma'] = dict(pairs_checked=pairs, nested_pairs=viol)
print('branch lemma, shape I: read-frame pairs checked %d, nested pairs (would let one register serve two) %d' % (pairs, viol))
if len(sys.argv) > 1: json.dump(out, open(sys.argv[1], 'w'), indent=1, default=str)
