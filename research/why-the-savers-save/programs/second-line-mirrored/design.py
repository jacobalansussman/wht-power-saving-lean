#!/usr/bin/env python3
"""cx-mirror / design.py (own code, stdlib).  What a mirrored word must do that U8 does not, made exact.

Frame rules (gcert/1 semantics, two-sided contract per call): x_b from <u_b> inside v_b^perp, y_b from <v_b> inside
u_b^perp; every gate at ONE common frame; frames only grow except where a register drops (paid).  Hence, for a
register between two of its drops, with N[a][b] = g_a . u_b:
  joins x_b1 then x_b2 :  its frame holds u_b1 and lies in v_b2^perp  ->  N[b2][b1] = 0
  joins y_t1 then y_t2 :  holds v_t1, lies in u_t2^perp               ->  N[t1][t2] = 0
  joins x_b  then y_t  :  holds u_b,  lies in u_t^perp                ->  u_b . u_t = 0
  joins y_t  then x_b  :  holds v_t,  lies in v_b^perp                ->  g_t . g_b = 1
(the same conditions bind any register whose frame was made equal to one that met those arrays: frames propagate
through gates).
Computes, for the natural mirror g_b = j + u_b (N[a][b] = 1 exactly on junk pairs):
 (1) the 8 totals of U8: support S_k and readers R_k; the least number of drop-free chains that can gather / clean
     S_k (x side) and serve R_k (y side) = clique cover numbers in the non-orthogonality graph;
 (2) the span of a largest clique (a clean chain's frame);
 (3) for ANY g (even, u_b . g_b = 0, g_a . g_b = 1 on junk): the bound L <= P7 + 1 + omega on a drop-free x chain,
     P7 = most ports in a 7-dimensional subspace of F_2^9 (exhaustive over the 2-dimensional annihilators)."""
import gzip, json, os, itertools, random
from fractions import Fraction as F
P = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz'))
d = json.load(gzip.open(P)); v = d['v']; R = d['R']; n = 2 * v + R; U = d['ports']
dot = lambda a, b: bin(a & b).count('1') & 1
out = {}
# supports and readers of the totals (x inputs only)
cont = [({r: F(1)} if r < v else {}) for r in range(n)]
for g in d['A']:
    if g[0] == 'out':
        src = dict(cont[g[2]])
        for t, a, b in g[3]:
            for k, val in src.items():
                nv = cont[t].get(k, 0) + F(a, b) * val
                if nv: cont[t][k] = nv
                else: cont[t].pop(k, None)
    else:
        for s_, a, b in g[3]:
            for k, val in list(cont[s_].items()):
                nv = cont[g[2]].get(k, 0) + F(a, b) * val
                if nv: cont[g[2]][k] = nv
                else: cont[g[2]].pop(k, None)
ret = sorted(d['ret'])
Sk = [sorted(cont[s_]) for _, s_, _ in ret]
Rk = [sorted(t for t, row in enumerate(d['scat']['table']) if any(k == kk for kk, _, _ in row)) for k in range(len(ret))]
print('supports = readers for every total:', all(Sk[k] == Rk[k] for k in range(8)), '| sizes', [len(x) for x in Sk])
junk_in = [sum(1 for a, b in itertools.combinations(Sk[k], 2) if dot(U[a], U[b]) == 0) for k in range(8)]
print('junk (orthogonal) pairs inside each support:', junk_in, 'of', 56 * 55 // 2)
out['junk_pairs_in_support'] = junk_in
adj = {a: {b for b in range(v) if b != a and dot(U[a], U[b]) == 1} for a in range(v)}

def clique_cover_lb_ub(S, tries=4000, seed=1):
    """lower bound: ceil(|S| / omega(S)); upper bound: best randomised greedy cover"""
    rnd = random.Random(seed); best = None
    for _ in range(tries):
        rest = list(S); rnd.shuffle(rest); parts = []
        rest = set(rest)
        while rest:
            order = sorted(rest, key=lambda x: rnd.random())
            c = []
            for x in order:
                if all(x in adj[y] for y in c): c.append(x)
            parts.append(c); rest -= set(c)
        if best is None or len(parts) < len(best): best = parts
    return best

def omega(nodes):
    best = []
    def bb(Rr, Pp):
        nonlocal best
        if len(Rr) + len(Pp) <= len(best): return
        if not Pp: best = list(Rr); return
        for x in list(Pp):
            if len(Rr) + len(Pp) <= len(best): return
            bb(Rr + [x], [y for y in Pp if y in adj[x] and y > x])
            Pp = [y for y in Pp if y != x]
    bb([], sorted(nodes)); return best

def rank(vs):
    basis = []
    for x in vs:
        for bb_ in basis: x = min(x, x ^ bb_)
        if x: basis.append(x)
    return len(basis)

covers = []
for k in range(8):
    w = omega(Sk[k]); cov = clique_cover_lb_ub(Sk[k])
    covers.append(dict(total=k, omega=len(w), lb=-(-56 // len(w)), greedy=len(cov), span_of_parts=sorted(rank([U[x] for x in c]) for c in cov)))
    print('total %d: omega in its support %d -> at least %d drop-free clean chains (greedy cover %d; spans of the parts %s)'
          % (k, len(w), -(-56 // len(w)), len(cov), covers[-1]['span_of_parts']))
out['covers'] = covers
# span of maximal cliques overall
allmax = []
def all_max_cliques(limit=2000):
    res = []
    def bk(Rr, Pp, Xx):
        if len(res) >= limit: return
        if not Pp and not Xx:
            res.append(Rr); return
        for x in list(Pp):
            bk(Rr + [x], Pp & adj[x], Xx & adj[x]); Pp = Pp - {x}; Xx = Xx | {x}
    bk([], set(range(v)), set()); return res
mc = all_max_cliques()
sizes = {}
for c in mc:
    key = (len(c), rank([U[x] for x in c])); sizes[key] = sizes.get(key, 0) + 1
print('maximal cliques of the non-orthogonality graph (size, span): counts', dict(sorted(sizes.items())), '(first', len(mc), ')')
out['max_clique_size_span'] = {str(k): c for k, c in sorted(sizes.items())}
# P7: most ports in a 7-dim subspace = ports orthogonal to a 2-dim subspace W = <w1, w2>
best7 = 0
for w1 in range(1, 512):
    for w2 in range(w1 + 1, 512):
        if (w1 ^ w2) < w2: continue          # count each 2-dim subspace once: w2 > w1 ^ w2 > ... (canonical pair)
        c = sum(1 for u in U if dot(u, w1) == 0 and dot(u, w2) == 0)
        best7 = max(best7, c)
best8 = max(sum(1 for u in U if dot(u, w) == 0) for w in range(1, 512))
print('most ports in a 7-dim subspace P7 =', best7, '| in an 8-dim subspace P8 =', best8)
print('any g: a drop-free x chain has at most P7 + 1 + omega =', best7 + 1 + 8, 'members (a total has 56)')
out['P7'] = best7; out['P8'] = best8; out['chain_bound_any_g'] = best7 + 9
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'design.json'), 'w'), indent=1)
# (4) the V-drop of a total's merged scatter copy: its frame holds the lines of all 56 readers
S10 = 1 << 9
Vm = [S10 | (((1 << 9) - 1) ^ u) for u in U]
vr = [rank([Vm[t] for t in Rk[k]]) for k in range(8)]
print('natural mirror: rank of the 56 readers\' lines per total (the merged copy\'s V-drop):', vr, '-> sum', sum(vr), 'per call')
out['V_drop_per_total'] = vr
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'design.json'), 'w'), indent=1)
