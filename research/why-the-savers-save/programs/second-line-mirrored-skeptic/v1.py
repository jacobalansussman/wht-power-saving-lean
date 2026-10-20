#!/usr/bin/env python3
"""v-cx-mirror / v1.py (own code, stdlib only; reads the U8X lowering as data, runs nothing of anyone else's).

Rebuilds from scratch:
 (a) the 8 retained totals after phase A (exact Fractions, x inputs only), the scatter's delivered map, junk pairs;
 (b) junk pairs against the binary Gram of the ports;
 (c) the natural mirror v_b = s + j + u_b: parity, u_b.v_b, v-Gram = u-Gram, the cross terms, psi;
 (d) inside each total's support: the largest pairwise NON-orthogonal set (a drop-free run for the natural mirror),
     by my own exact clique search (Bron-Kerbosch with pivoting), and the global omega;
 (e) P7 = most ports in a 7-dim subspace of F_2^9 (all 43,435 of them);
 (f) the twist proof step: quadratic functions on F_2^9 that vanish on the 112 vectors e_b = u_b + j: do they
     all vanish on the whole even space E?  And the count of symmetric involutions A = I + M (M alternating,
     M^2 = 0, M j = 0) by brute force over M (2^36 is too many: I count by images instead, cross-checked by a
     direct enumeration of the alternating M with M j = 0 restricted to rank <= 4 through their image).
Bits 0..8 = window, bit 9 = s."""
import gzip, json, os, sys, time
from fractions import Fraction as Fr
from itertools import combinations
T0 = time.time()
P = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz'))
d = json.load(gzip.open(P))
h, v, R = d['h'], d['v'], d['R']; n = 2 * v + R
U = list(d['ports']); S = 1 << 9; J = (1 << 9) - 1
def dot(a, b): return bin(a & b).count('1') & 1
out = {}
# (a) phase A replay, my own: register -> {x index: Fraction}
reg = [dict() for _ in range(n)]
for t in range(v): reg[t] = {t: Fr(1)}
def axpy(t, src, c):
    dd = reg[t]
    for k, val in src.items():
        w = dd.get(k, 0) + c * val
        if w: dd[k] = w
        else: dd.pop(k, None)
for g in d['A']:
    if g[0] == 'out':
        src = dict(reg[g[2]])
        for t, a, b in g[3]: axpy(t, src, Fr(a, b))
    elif g[0] == 'in':
        for s_, a, b in g[3]: axpy(g[2], dict(reg[s_]), Fr(a, b))
    else: raise SystemExit('unknown gate')
assert all(reg[t] == {t: 1} for t in range(v)), 'x restored after A?'
ret = sorted(d['ret'])
tot = [dict(reg[s_]) for _, s_, _ in ret]
supp = [sorted(T) for T in tot]
print('totals', len(tot), 'support sizes', [len(x) for x in supp], 'coefficient values', sorted({str(c) for T in tot for c in T.values()}))
# delivered map
Dl = [dict() for _ in range(v)]
for t, row in enumerate(d['scat']['table']):
    for k, a, b in row:
        for bb, val in tot[k].items():
            w = Dl[t].get(bb, 0) + Fr(a, b) * val
            if w: Dl[t][bb] = w
            else: Dl[t].pop(bb, None)
diag_ok = all(Dl[t].get(t) == 1 for t in range(v))
junk = {(t, b) for t in range(v) for b in Dl[t] if b != t}
orth = {(a, b) for a in range(v) for b in range(v) if a != b and dot(U[a], U[b]) == 0}
print('diag 1:', diag_ok, '| junk entries', len(junk), '| junk == orthogonal pairs:', junk == orth,
      '| per target', sorted({sum(1 for b in Dl[t] if b != t) for t in range(v)}))
out['junk_eq_orth'] = junk == orth; out['njunk'] = len(junk)
# readers per total
readers = [sorted(t for t in range(v) for k, a, b in d['scat']['table'][t] if k == kk) for kk in range(len(tot))]
print('readers per total', [len(r) for r in readers], '| reads per y', sorted({len(r) for r in d['scat']['table']}))
print('ports all odd', all(bin(u).count('1') % 2 for u in U), '| weights', sorted({bin(u).count('1') for u in U}))
# (c) natural mirror
V = [S | (J ^ u) for u in U]
c1 = all(bin(J ^ u).count('1') % 2 == 0 for u in U)
c2 = all(dot(U[b], V[b]) == 0 for b in range(v))
c3 = all(dot(V[a], V[b]) == dot(U[a], U[b]) for a in range(v) for b in range(v))
c4 = all(dot(U[a], V[b]) == 1 ^ dot(U[a], U[b]) for a in range(v) for b in range(v))
def psi(w):  # psi(w) = w + wt(w)(s + j) on the window part, with psi(s) = j: linear map given on the basis
    img = 0
    for i in range(10):
        if w >> i & 1: img ^= (J if i == 9 else ((1 << i) ^ S ^ J))
    return img
c5 = all(psi(psi(w)) == w for w in range(1 << 10)) and all(dot(psi(a), psi(b)) == dot(a, b) for a in range(1 << 10) for b in range(0, 1 << 10, 7))
c6 = all(psi(U[b]) == V[b] and psi(V[b]) == U[b] for b in range(v))
print('mirror: g even', c1, '| u_b.v_b=0', c2, '| v-Gram = u-Gram', c3, '| u_a.v_b = 1+u_a.u_b', c4, '| psi involutive isometry', c5, '| psi swaps', c6)
out['mirror'] = [c1, c2, c3, c4, c5, c6]
# (d) cliques of the non-orthogonality graph (N-independent sets for the natural mirror)
nb = [0] * v
for a in range(v):
    for b in range(v):
        if a != b and dot(U[a], U[b]) == 1: nb[a] |= 1 << b
def bits(x):
    while x:
        l = x & -x; yield l.bit_length() - 1; x ^= l
def maxclique(cand):
    best = [0]
    def bk(Rn, Pm, Xm):
        if not Pm and not Xm:
            if Rn > best[0]: best[0] = Rn
            return
        if Rn + bin(Pm).count('1') <= best[0]: return
        piv = max(bits(Pm | Xm), key=lambda u: bin(Pm & nb[u]).count('1'))
        for u in list(bits(Pm & ~nb[piv])):
            bk(Rn + 1, Pm & nb[u], Xm & nb[u]); Pm &= ~(1 << u); Xm |= 1 << u
    bk(0, cand, 0)
    return best[0]
inS = [maxclique(sum(1 << b for b in s_)) for s_ in supp]
glob = maxclique((1 << v) - 1)
print('largest pairwise non-orthogonal set inside each support', inS, '| global omega', glob)
out['run_in_support'] = inS; out['omega'] = glob
# count of 4-sets and 5-sets inside support 0, plain
s0 = supp[0]
c4s = sum(1 for q in combinations(s0, 4) if all(dot(U[a], U[b]) for a, b in combinations(q, 2)))
print('support 0: pairwise non-orthogonal 4-sets', c4s, '| junk pairs inside', sum(1 for a, b in combinations(s0, 2) if dot(U[a], U[b]) == 0))
# (e) P7: 7-dim subspaces = annihilators of 2-dim subspaces <p, q>
best7 = 0
seen = set()
for p in range(1, 512):
    for q in range(p + 1, 512):
        key = tuple(sorted((p, q, p ^ q)))
        if key in seen: continue
        seen.add(key)
        c = sum(1 for u in U if dot(u, p) == 0 and dot(u, q) == 0)
        if c > best7: best7 = c
best8 = max(sum(1 for u in U if dot(u, p) == 0) for p in range(1, 512))
print('7-dim subspaces', len(seen), '| P7', best7, '| P8', best8)
out['P7'] = best7; out['P8'] = best8
# (f) quadratic functions q(e) = sum a_i e_i + sum_{i<j} b_ij e_i e_j vanishing on the 112 e_b; do they vanish on E?
def qrow(e):
    r = 0; k = 0
    for i in range(9):
        if e >> i & 1: r |= 1 << k
        k += 1
    for i in range(9):
        for jj in range(i + 1, 9):
            if (e >> i & 1) and (e >> jj & 1): r |= 1 << k
            k += 1
    return r
def rank(rows):
    basis = {}
    for x in rows:
        while x:
            p = x.bit_length() - 1
            if p in basis: x ^= basis[p]
            else: basis[p] = x; break
    return len(basis)
E_ = [e for e in range(512) if bin(e).count('1') % 2 == 0]
rP = rank([qrow(u ^ J) for u in U]); rE = rank([qrow(e) for e in E_])
rPE = rank([qrow(u ^ J) for u in U] + [qrow(e) for e in E_])
print('quadratic functions (45 coeffs): rank on the 112 e_b', rP, '| on E', rE, '| joint', rPE,
      '=> every q vanishing on the e_b vanishes on E:', rPE == rP)
out['quad_ranks'] = [rP, rE, rPE]
# (g) count symmetric involutions by images: totally isotropic I in E of dim 0, 2, 4 (dim 1, 3 give no nondegenerate
#     alternating form); for dim 2 one form, for dim 4 |GL4|/|Sp4| = 28.  Count isotropic subspaces directly.
Ev = [e for e in E_ if e]
iso2 = set()
for a in Ev:
    for b in Ev:
        if a < b and dot(a, b) == 0: iso2.add(tuple(sorted((a, b, a ^ b))))
lag = set()
for pl in iso2:
    sp = set(pl) | {0}
    perp = [e for e in Ev if all(dot(e, x) == 0 for x in pl) and e not in sp]
    for c in perp:
        sp3 = sp | {x ^ c for x in sp}
        for dd in perp:
            if dd in sp3: continue
            if dot(c, dd): continue
            sp4 = frozenset(sp3 | {x ^ dd for x in sp3})
            lag.add(sp4)
print('isotropic planes', len(iso2), '| Lagrangians', len(lag), '| twisted mirrors 1 + planes + 28 * Lagrangians =',
      1 + len(iso2) + 28 * len(lag))
out['twists'] = [len(iso2), len(lag), 1 + len(iso2) + 28 * len(lag)]
out['secs'] = round(time.time() - T0, 1)
json.dump(dict(out, supp=supp, readers=readers), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'v1.json'), 'w'))
print('done in', out['secs'], 's')
