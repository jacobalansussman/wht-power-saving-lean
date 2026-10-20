#!/usr/bin/env python3
"""fc-reopt: chronology-free bound on the lever price at fixed words, for every separated cut depth (o-mip's cf_lp.py
generalised from the half cut to the full cut, partial depths and optional rank-16/32 first-touch heads). No hub code.

Separated rule with paid-rank cap C (24 <= C <= 48; k_min = C + 1; PROOF.md hypothesis (S)):
  - a literal head of rank r (empty first-touch bin) is cut by d_r = min(r, 72 - C) axes: its paid block r -> r - d_r,
    its array's target 73 -> 73 - d_r (>= C + 1);  C = 24: full cut (rank-48 heads removed, target 25);
    C = 48: half cut (rank-48 heads 48 -> 24, target 49).
  - every other paid block of rank l > C is split into floor(l/C) blocks of C and one of l mod C (the cheapest split
    into pieces <= C, by concavity of l^tau); split_block in Work.Mixed.HeadM1 allows any a + b split.
Variables: x[c, q] >= 0 bins of capacity c with pattern q (non-head); x_head[r] >= 0 empty bins of capacity r that are
heads (r = 48: children 8/10; with --ext also r = 32: children 4/6 and r = 16: children 0/2), each <= 224 T.
Rows: bins per capacity; roles per width (need_w, the whole fill).
Objective (Dinkelbach, exact for this fractional problem): G(tau) = sum cost x + sum_local n_r r^tau - W 73^tau with
  cost(x[c, q]) = phi_C(l_q),  cost(x_head[r]) = (r - d_r)^tau + 73^tau - (73 - d_r)^tau.
Windows, the one-role-per-bin-per-call rule and the compact prefix are dropped, so the LP value bounds every valid fill
of these words (MILP: the best integer mix ignoring chronology).
usage: fc_lp.py INV.json --C 24 [--ext] [--milp] [--b0 B] [--no-heads] [--out R.json]
"""
import sys, json, argparse
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

ap = argparse.ArgumentParser()
ap.add_argument('inv'); ap.add_argument('--C', type=int, default=24); ap.add_argument('--ext', action='store_true')
ap.add_argument('--milp', action='store_true'); ap.add_argument('--b0', type=float, default=1.78e-3)
ap.add_argument('--no-heads', action='store_true'); ap.add_argument('--rounds', type=int, default=10)
ap.add_argument('--out')
a = ap.parse_args()
inv = json.load(open(a.inv))
T = inv['T']; W = 896 * T; C = a.C
assert 24 <= C <= 48
need = {int(k): v for k, v in inv['need'].items() if v}
local = {int(k): v for k, v in inv['local'].items()}
bins = {int(k): v for k, v in inv['bins_by_cap'].items()}
widths = sorted(need); caps = sorted(bins)
headranks = [] if a.no_heads else ([48, 32, 16] if a.ext else [48])
Hmax = 224 * T


def patterns(cap):
    out = []
    def rec(i, rest, cur):
        if i == len(widths): out.append(tuple(cur)); return
        for k in range(rest // widths[i] + 1): rec(i + 1, rest - k * widths[i], cur + [k])
    rec(0, cap, [])
    return out


def pieces(l):
    if l <= C: return {l: 1} if l else {}
    q, m = divmod(l, C); p = {C: q}
    if m: p[m] = p.get(m, 0) + 1
    return p


def dcut(r): return min(r, 72 - C)

cols = []   # (kind, cap, pattern, left)
for c in caps:
    for q in patterns(c):
        cols.append(('bin', c, q, c - sum(k * w for k, w in zip(q, widths))))
for r in headranks:
    cols.append(('head', r, tuple(0 for _ in widths), r))
nv = len(cols)
A = np.zeros((len(caps) + len(widths), nv)); rhs = np.zeros(len(caps) + len(widths))
for j, (kind, c, q, l) in enumerate(cols):
    A[caps.index(c), j] = 1
    for i, k in enumerate(q):
        if k: A[len(caps) + i, j] = k
for i, c in enumerate(caps): rhs[i] = bins[c]
for i, w in enumerate(widths): rhs[len(caps) + i] = need[w]
ub = np.array([Hmax if kind == 'head' else np.inf for (kind, c, q, l) in cols])


def costs(tau):
    out = []
    for kind, c, q, l in cols:
        if kind == 'bin':
            out.append(sum(n * p ** tau for p, n in pieces(l).items()))
        else:
            d = dcut(c)
            out.append(((c - d) ** tau if c - d else 0.0) + 73 ** tau - (73 - d) ** tau)
    return np.array(out)


def model(x):
    """children and targets of the solution x."""
    ch = dict((r, float(n)) for r, n in local.items()); tg = {73: float(W)}
    for j, (kind, c, q, l) in enumerate(cols):
        if x[j] <= 1e-9: continue
        if kind == 'bin':
            for p, n in pieces(l).items(): ch[p] = ch.get(p, 0) + n * x[j]
        else:
            d = dcut(c)
            if c - d: ch[c - d] = ch.get(c - d, 0) + x[j]
            tg[73] -= x[j]; tg[73 - d] = tg.get(73 - d, 0) + x[j]
    return ch, tg


def root(ch, tg):
    def G(t): return sum(n * r ** t for r, n in ch.items() if n) - sum(n * k ** t for k, n in tg.items() if n)
    lo, hi = 0.99, 0.99999
    for _ in range(200):
        m = (lo + hi) / 2
        if G(m) > 0: lo = m
        else: hi = m
    return 1 - (lo + hi) / 2


def solve(b):
    tau = 1 - b; cst = costs(tau)
    if a.milp:
        res = milp(cst, constraints=LinearConstraint(A, rhs, rhs), integrality=np.ones(nv),
                   bounds=Bounds(np.zeros(nv), ub), options={'disp': False, 'time_limit': 900, 'mip_rel_gap': 1e-13})
        assert res.status == 0, res.message
        return res.x, float(res.fun), None
    res = linprog(cst, A_eq=A, b_eq=rhs, bounds=list(zip(np.zeros(nv), ub)), method='highs')
    assert res.status == 0, res.message
    duals = dict(eq=[float(v) for v in res.eqlin.marginals], upper=[float(v) for v in res.upper.marginals])
    return res.x, float(res.fun), duals


def label(kind, c, q):
    s = ','.join('%d' % w for w, k in zip(widths, q) for _ in range(k)) or 'empty'
    return ('head%d' % c) if kind == 'head' else '%d:%s' % (c, s)

b = a.b0; hist = []
for it in range(a.rounds):
    x, obj, duals = solve(b)
    ch, tg = model(x); nb = root(ch, tg)
    heads = {c: float(x[j]) for j, (kind, c, q, l) in enumerate(cols) if kind == 'head'}
    idle = {}
    for j, (kind, c, q, l) in enumerate(cols):
        if kind == 'bin' and l and x[j] > 1e-9: idle[l] = idle.get(l, 0) + float(x[j])
    hist.append(dict(b_in=b, b_out=nb, heads=heads, idle={str(k): round(v, 6) for k, v in sorted(idle.items())}))
    print('C=%d ext=%s milp=%s round %d  b_in %.13e -> %.13e  heads %s  idle-remainders %s' % (
        C, a.ext, a.milp, it, b, nb, heads, {k: round(v, 3) for k, v in sorted(idle.items())}), flush=True)
    if abs(nb - b) < 1e-15: break
    b = nb
mix = {label(kind, c, q): float(x[j]) for j, (kind, c, q, l) in enumerate(cols) if x[j] > 1e-9}
out = dict(inv=a.inv, T=T, C=C, ext=a.ext, milp=a.milp, headranks=headranks, Hmax=Hmax, rounds=hist, b=nb,
           tau=1 - nb, heads=heads, mix=mix, children={str(k): v for k, v in sorted(ch.items())},
           targets={str(k): v for k, v in sorted(tg.items())}, ncols=nv, widths=widths, caps=caps)
if duals is not None:
    out['duals_at_b_in'] = dict(b_in=hist[-1]['b_in'], **duals)
if a.out: json.dump(out, open(a.out, 'w'), indent=1)
