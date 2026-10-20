#!/usr/bin/env python3
"""cx-mirror / lines.py (own code, stdlib).  The port geometry of the U8X word and the second lines v_b = s + g_b.

Reads the U8X lowering (gcert/1, our own data file) as data.  Computes:
 (1) the delivered map of the centre: after phase A the 8 retained totals hold combinations of the x's; the scatter
     adds sum_k c_tk tot_k to y_t.  Delivered D[t][b]; junk pairs = off-diagonal nonzero entries;
 (2) the binary Gram of the ports on junk pairs and on the other pairs;
 (3) the natural mirror g_b = j + u_b (j = all-ones of F_2^9): parity of g_b, u_b . v_b, v_a . v_b against u_a . u_b,
     and the round-2 orthogonality test g_a . g_b = 1 on every junk pair;
 (4) R2: for each total, the meet of u_t^perp over its readers (in F_2^10 = window + s), and over all 112 y's;
 (5) cliques of the binary non-orthogonality graph (omega), exact by branch and bound.
Local space F_2^10: bits 0..8 = the window, bit 9 = s.  v . w = parity(v & w)."""
import gzip, json, os, sys
from fractions import Fraction as F
from itertools import combinations
P = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz'))
d = json.load(gzip.open(P))
h, v, R = d['h'], d['v'], d['R']; n = 2 * v + R
U = d['ports']; S = 1 << 9; J9 = (1 << 9) - 1
dot = lambda a, b: bin(a & b).count('1') & 1
out = {}

# (1) totals after phase A
cont = [({r: F(1)} if r < v else {}) for r in range(n)]   # x inputs only (slots start clean; exterior dirt never reaches y, E7)
def addto(t, s, c):
    for k, val in list(cont[s].items()):
        nv = cont[t].get(k, 0) + c * val
        if nv: cont[t][k] = nv
        else: cont[t].pop(k, None)
for g in d['A']:
    if g[0] == 'out':
        src = dict(cont[g[2]])
        for t, a, b in g[3]:
            for k, val in src.items():
                nv = cont[t].get(k, 0) + F(a, b) * val
                if nv: cont[t][k] = nv
                else: cont[t].pop(k, None)
    else:
        for s_, a, b in g[3]: addto(g[2], s_, F(a, b))
ret = sorted(d['ret']); tot = [cont[s_] for _, s_, _ in ret]
assert all(all(k < v for k in T) for T in tot), 'totals hold only x data'
Dl = [[F(0)] * v for _ in range(v)]
readers = {k: [] for k in range(len(ret))}
for t, row in enumerate(d['scat']['table']):
    for k, a, b in row:
        readers[k].append(t)
        for bb, val in tot[k].items(): Dl[t][bb] += F(a, b) * val
diag = sorted({Dl[t][t] for t in range(v)})
junk = [(t, b) for t in range(v) for b in range(v) if b != t and Dl[t][b] != 0]
jp = {t: sorted(b for tt, b in junk if tt == t) for t in range(v)}
print('totals: support sizes', [len(T) for T in tot], '| diagonal of the delivered map', [str(x) for x in diag])
print('junk entries', len(junk), '| per target', sorted({len(x) for x in jp.values()}),
      '| symmetric', all((b, t) in set(junk) for t, b in junk))
out['total_support'] = [len(T) for T in tot]; out['junk_per_target'] = sorted({len(x) for x in jp.values()})
# (2) binary Gram on junk pairs / others
jset = set(junk)
gj = {dot(U[t], U[b]) for t, b in junk}
non = [(a, b) for a in range(v) for b in range(v) if a < b and (a, b) not in jset]
gn = {}
for a, b in non: gn[dot(U[a], U[b])] = gn.get(dot(U[a], U[b]), 0) + 1
print('binary u.u on junk pairs:', gj, '| on the other', len(non), 'pairs: counts', gn)
out['junk_binary_gram'] = sorted(gj); out['nonjunk_binary_gram_counts'] = gn
# (3) natural mirror
g = [J9 ^ u for u in U]; V = [S | x for x in g]
print('mirror: g even', all(bin(x).count('1') % 2 == 0 for x in g), '| u_b.v_b = 0 for all b', all(dot(U[b], V[b]) == 0 for b in range(v)),
      '| v.v = 1', all(dot(x, x) == 1 for x in V), '| v-Gram = u-Gram', all(dot(V[a], V[b]) == dot(U[a], U[b]) for a in range(v) for b in range(v)),
      '| g_a.g_b = 1 on junk', all(dot(g[a], g[b]) == 1 for a, b in junk), '| lines distinct', len(set(V)) == v)
print('mirror cross terms u_a.v_b = 1 + u_a.u_b:', all(dot(U[a], V[b]) == 1 ^ dot(U[a], U[b]) for a in range(v) for b in range(v)))
out['mirror'] = dict(g_even=True, passes_junk_orthogonality=all(dot(g[a], g[b]) == 1 for a, b in junk))
# (4) R2: meet of u_t^perp over readers, in F_2^10
def rank(vs):
    basis = []
    for x in vs:
        for bb in basis:
            x = min(x, x ^ bb)
        if x: basis.append(x)
    return len(basis)
def perp_meet(ts):
    """dimension of {w in F_2^10 : w . u_t = 0 for all t} = 10 - rank(u_t)"""
    return 10 - rank([U[t] for t in ts])
r2 = {k: perp_meet(ts) for k, ts in readers.items()}
print('R2: dim of the meet of u_t^perp over the readers of each total (F_2^10):', r2, '| over all 112:', perp_meet(range(v)))
out['R2_meet_dims'] = r2
# the largest set of y's that could share one frame holding all their lines:
# need v_a in u_b^perp for all a, b in the set  (frame contains every v, lies in every u^perp)
def compat(a, b): return dot(U[b], V[a]) == 0 and dot(U[a], V[b]) == 0
# (5) cliques
def omega(nodes, adj):
    best = []
    def bb(R_, P_):
        nonlocal best
        if len(R_) + len(P_) <= len(best): return
        if not P_:
            best = list(R_); return
        for x in list(P_):
            if len(R_) + len(P_) <= len(best): return
            bb(R_ + [x], [y for y in P_ if y in adj[x] and y > x])
            P_ = [y for y in P_ if y != x]
    bb([], sorted(nodes))
    return best
adj_non = {a: {b for b in range(v) if b != a and dot(U[a], U[b]) == 1} for a in range(v)}
w = omega(range(v), adj_non)
print('omega of the binary NON-orthogonality graph of the ports:', len(w), w)
adj_c = {a: {b for b in range(v) if b != a and compat(a, b)} for a in range(v)}
wc = omega(range(v), adj_c)
print('largest set of y\'s (mirror lines) that can share one frame containing all their lines:', len(wc))
out['omega_u'] = len(w); out['shared_frame_max'] = len(wc)
deg = sorted({len(adj_non[a]) for a in range(v)})
print('non-orthogonal partners per port:', deg)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lines.json'), 'w'), indent=1, default=str)
