"""t-egraph extraction step: exact MILP over ALL commuting reorders + frames (label space), with the word's own frame
vocabulary as the candidate pool.

Variables x[g,U] (gate g at frame U), u[r,U] (register r uses frame U), z[r,a,b] (r's chain climbs from dim a to dim b).
Constraints: one frame per gate; u >= x; per register at most one frame per dimension; cross-dimension comparability
(big-M per frame); f(g) <= f(h) on every edge of the transitively reduced dependency DAG (non-commuting pairs on a common
register, section order A < B < tail); centres reach their centre frame in A; per register a dim path start -> used dims
-> end.  Objective: sum d ln(73/d) over climbs.  Realisability (lsearch.py): order = sort by dim, ties topological.
"""
import sys, os, time, json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
from collections import defaultdict
sys.path.insert(0, '.')
import fw, search
from search import le, join, meet, FC
from lsearch import LS

W = search.Word(fw.load(os.environ.get('WORD', 'K8-FULL-WORD.json.gz')))
L = LS(W); N = L.N
TL = float(sys.argv[1]) if len(sys.argv) > 1 else 600
# DAG-closed bounds
lo = list(L.lo); hi = list(L.hi)
for _ in range(60):
    ch = False
    for g in range(N):
        for h in L.before[g]:
            nl = join(lo[g], lo[h])
            if nl != lo[g]: lo[g] = nl; ch = True
        for h in L.after[g]:
            nh = meet(hi[g], hi[h])
            if nh != hi[g]: hi[g] = nh; ch = True
    if not ch: break
assert all(le(lo[g], L.f[g]) and le(L.f[g], hi[g]) for g in range(N))
vocab = set(L.f) | set(W.start) | set(W.end) | set(W.cframe.values())
WIDE = int(os.environ.get('WIDE', '0'))
from ilp_wide import all_between
cand = []
for g in range(N):
    base = {U for U in vocab if le(lo[g], U) and le(U, hi[g])}
    gap = len(hi[g]) - len(lo[g])
    if WIDE and gap <= WIDE: base |= all_between(lo[g], hi[g])
    elif WIDE:
        nb = set()
        for r in L.regs[g]:
            for h in L.G[r]: nb.add(L.f[h])
        for U in nb:
            for X in (join(lo[g], U), meet(hi[g], join(lo[g], U))):
                if le(lo[g], X) and le(X, hi[g]): base.add(X)
    c = sorted(base, key=lambda U: (len(U), U))
    assert L.f[g] in c
    cand.append(c)
print('vocab', len(vocab), 'candidates per gate: mean', round(sum(map(len, cand)) / N, 2), 'max', max(map(len, cand)),
      'forced', sum(1 for c in cand if len(c) == 1), flush=True)
# transitive reduction of the DAG (original order is topological)
pos = {g: i for i, g in enumerate(range(N))}
reach = [0] * N; red = defaultdict(list); nedges = 0
for g in reversed(range(N)):
    covered = 0
    for h in sorted(L.after[g]):
        nedges += 1
        if not (covered >> h) & 1: red[g].append(h)
        covered |= reach[h] | (1 << h)
    reach[g] = covered
print('DAG edges', nedges, 'reduced', sum(map(len, red.values())), flush=True)
# variables
col = 0; xid = {}
for g in range(N):
    for U in cand[g]: xid[(g, U)] = col; col += 1
Cr = defaultdict(set)
for g in range(N):
    for r in L.regs[g]: Cr[r] |= set(cand[g])
uid = {}
for r in range(W.n):
    for U in Cr[r]: uid[(r, U)] = col; col += 1
zid = {}; zcost = []
for r in range(W.n):
    s, e = len(W.start[r]), len(W.end[r])
    for a in range(s, e + 1):
        for b in range(a + 1, e + 1):
            zid[(r, a, b)] = col; col += 1
rows, cols, vals, lo_b, hi_b = [], [], [], [], []
def add(ent, l, h):
    k = len(lo_b)
    for c_, v_ in ent: rows.append(k); cols.append(c_); vals.append(v_)
    lo_b.append(l); hi_b.append(h)
INF = np.inf
for g in range(N): add([(xid[(g, U)], 1) for U in cand[g]], 1, 1)
for g in range(N):
    for r in L.regs[g]:
        for U in cand[g]: add([(uid[(r, U)], 1), (xid[(g, U)], -1)], 0, INF)
for r in range(W.n):
    s, e = len(W.start[r]), len(W.end[r])
    byd = defaultdict(list)
    for U in Cr[r]:
        if U in (W.start[r], W.end[r]): continue
        byd[len(U)].append(U)
    for d, Us in byd.items():
        if len(Us) > 1: add([(uid[(r, U)], 1) for U in Us], 0, 1)
    for V in Cr[r]:
        bad = [U for U in Cr[r] if len(U) < len(V) and not le(U, V)]
        if bad:
            M = len(bad)
            add([(uid[(r, U)], 1) for U in bad] + [(uid[(r, V)], M)], -INF, M)
    # flow over dims
    if s == e: continue
    for d in range(s, e + 1):
        inn = [(zid[(r, a, d)], 1) for a in range(s, d)]
        out = [(zid[(r, d, b)], 1) for b in range(d + 1, e + 1)]
        yd = [(uid[(r, U)], -1) for U in byd.get(d, [])]
        if d == s: add(out, 1, 1)
        elif d == e: add(inn, 1, 1)
        else:
            add(inn + yd, 0, 0); add(out + yd, 0, 0)
for g, hs in red.items():
    for h in hs:
        for U in cand[g]:
            ent = [(xid[(g, U)], 1)] + [(xid[(h, V)], -1) for V in cand[h] if le(U, V)]
            add(ent, -INF, 0)
for q, A in L.centreA.items():
    if A: add([(xid[(g, W.cframe[q])], 1) for g in A if (g, W.cframe[q]) in xid], 1, INF)
c = np.zeros(col)
for (r, a, b), k in zid.items(): c[k] = FC[b - a]
A_ = coo_matrix((vals, (rows, cols)), shape=(len(lo_b), col)).tocsr()
integ = np.zeros(col); integ[:len(xid)] = 1; integ[len(xid):len(xid) + len(uid)] = 1
print('MILP vars', col, '(x', len(xid), 'u', len(uid), 'z', len(zid), ') rows', len(lo_b), flush=True)
# sanity: the current labelling must be feasible with its own cost
x0 = np.zeros(col)
for g in range(N): x0[xid[(g, L.f[g])]] = 1
for r in range(W.n):
    used = {L.f[g] for g in L.G.get(r, ())}
    for U in used: x0[uid[(r, U)]] = 1
    S = sorted(used | {W.start[r], W.end[r]}, key=len)
    for a_, b_ in zip(S, S[1:]): x0[zid[(r, len(a_), len(b_))]] = 1
viol = A_ @ x0
bad = [i for i in range(len(lo_b)) if viol[i] < lo_b[i] - 1e-9 or viol[i] > hi_b[i] + 1e-9]
print('current labelling feasible:', not bad, 'violations', len(bad), 'cost', round(float(c @ x0), 4), flush=True)
t = time.time()
res = milp(c, constraints=LinearConstraint(A_, lo_b, hi_b), integrality=integ, bounds=Bounds(0, 1),
           options=dict(time_limit=TL, disp=True, mip_rel_gap=1e-9))
print('status', res.status, res.message, 'obj', res.fun, 'bound', getattr(res, 'mip_dual_bound', None), 'secs', round(time.time() - t, 1), flush=True)
out = dict(status=int(res.status), obj=res.fun, dual_bound=getattr(res, 'mip_dual_bound', None), current=float(c @ x0))
if res.x is not None:
    newf = [max(cand[g], key=lambda U: res.x[xid[(g, U)]]) for g in range(N)]
    L.f = newf
    out['relabelled_cost'] = L.total(); out['changed_gates'] = sum(1 for g in range(N) if newf[g] != W.f[g])
    print('re-costed', out['relabelled_cost'], 'changed gates', out['changed_gates'], flush=True)
    if out['relabelled_cost'] < out['current'] - 1e-6:
        w2 = L.to_word(); fw.save(w2, 'W-lmilp-w%d.json.gz' % WIDE); print('saved W-lmilp.json.gz', w2['report']['local_profile'], flush=True)
json.dump(out, open('lmilp_%s_w%d.json' % (os.path.basename(os.environ.get('WORD', 'k8')).split('.')[0], WIDE), 'w'), indent=1)
