#!/usr/bin/env python3
"""fc-reopt: exact-arithmetic dual certificate for the chronology-free bound of fc_lp.py (standard library only).
Generalises o-mip's certify_cf.py from the half cut to every separated cut cap C and optional rank-16/32 heads.

Claim certified at b* (tau* = 1 - b*): for EVERY assignment of the words' helper roles to the idle bins (any windows,
any chronology, any per-call rule), with heads only on empty bins of the head classes (at most 224 T per rank) and the
cut of cap C (fc_lp.py's model),
    G(tau*) = sum_bins cost(bin) + sum_local n_r r^tau* - W 73^tau*  >=  G_lb,
where G_lb is the Lagrangian dual value of a dual vector y (one entry per capacity, one per width) and s_r >= 0 (the
head caps), minus bins[c] x (most negative reduced cost of capacity c) if any reduced cost is negative.
If G_lb > 0, every such fill has its cut root tau > tau*, i.e. its price b < b* (separation makes G decreasing).
y comes from complementary slackness on the tight columns of the LP optimum (fc_lp.py's mix + HiGHS duals), solved
in 90-digit Decimal; every pattern of every capacity is then checked.
usage: fc_cert.py INV.json LP.json --bstar B [--out R.json]      (C, ext read from LP.json)
"""
import sys, json, argparse
from decimal import Decimal as D, getcontext
getcontext().prec = 90
ap = argparse.ArgumentParser(); ap.add_argument('inv'); ap.add_argument('lp'); ap.add_argument('--bstar', required=True)
ap.add_argument('--out')
a = ap.parse_args()
inv = json.load(open(a.inv)); lp = json.load(open(a.lp))
T = inv['T']; W = 896 * T; C = lp['C']; headranks = lp['headranks']; Hmax = 224 * T
need = {int(k): v for k, v in inv['need'].items() if v}; widths = sorted(need)
bins = {int(k): v for k, v in inv['bins_by_cap'].items()}; caps = sorted(bins)
local = {int(k): D(v) for k, v in inv['local'].items()}
bstar = D(a.bstar); tau = 1 - bstar
_pw = {}
def pw(r):
    if r not in _pw: _pw[r] = (D(r).ln() * tau).exp() if r else D(0)
    return _pw[r]
def phi(l):
    if l <= 0: return D(0)
    q, m = divmod(l, C)
    return q * pw(C) + (pw(m) if m else D(0))
def hcost(r):
    d = min(r, 72 - C)
    return (pw(r - d) if r - d else D(0)) + pw(73) - pw(73 - d)

def patterns(cap):
    out = []
    def rec(i, rest, cur):
        if i == len(widths): out.append(tuple(cur)); return
        for k in range(rest // widths[i] + 1): rec(i + 1, rest - k * widths[i], cur + [k])
    rec(0, cap, [])
    return out

def parse(key):
    if key.startswith('head'): return ('head', int(key[4:]), tuple(0 for _ in widths))
    c, p = key.split(':'); ws = [] if p == 'empty' else [int(x) for x in p.split(',')]
    return ('bin', int(c), tuple(ws.count(w) for w in widths))

allcols = []
for c in caps:
    for q in patterns(c): allcols.append(('bin', c, q))
for r in headranks: allcols.append(('head', r, tuple(0 for _ in widths)))
def cost(col):
    kind, c, q = col
    if kind == 'head': return hcost(c)
    return phi(c - sum(k * w for k, w in zip(q, widths)))

# tight columns: the LP's support (x > 0), plus columns with zero float reduced cost under HiGHS's duals
supp = [parse(k) for k, v in lp['mix'].items() if v > 1e-7]
heads_at_cap = {parse(k)[1] for k, v in lp['mix'].items() if k.startswith('head') and abs(v - Hmax) < 1e-6}
eq = lp['duals_at_b_in']['eq']
yf = {c: eq[i] for i, c in enumerate(caps)}; ywf = {w: eq[len(caps) + i] for i, w in enumerate(widths)}
tau_f = 1 - lp['duals_at_b_in']['b_in']
def rcf(col):
    kind, c, q = col
    if kind == 'head':
        d = min(c, 72 - C); v = ((c - d) ** tau_f if c - d else 0) + 73 ** tau_f - (73 - d) ** tau_f
    else:
        l = c - sum(k * w for k, w in zip(q, widths))
        if l <= 0: v = 0.0
        else:
            qq, m = divmod(l, C); v = qq * C ** tau_f + (m ** tau_f if m else 0)
    return v - yf[c] - sum(k * ywf[w] for k, w in zip(q, widths))
eqcols = [col for col in supp if not (col[0] == 'head' and col[1] in heads_at_cap)]
extra = [col for col in allcols if col not in supp and col[0] == 'bin' and abs(rcf(col)) < 1e-7]
nun = len(caps) + len(widths)

def solve(cols):
    rows = []
    for col in cols:
        kind, c, q = col
        coef = [D(0)] * nun; coef[caps.index(c)] = D(1)
        for i, k in enumerate(q): coef[len(caps) + i] = D(k)
        rows.append(coef + [cost(col)])
    A = [r[:] for r in rows]; m = len(A); r0 = 0; pivcols = []
    for col in range(nun):
        if r0 >= m: break
        best = max(range(r0, m), key=lambda i: abs(A[i][col]))
        if abs(A[best][col]) < D('1e-40'): continue
        A[r0], A[best] = A[best], A[r0]
        for i in range(m):
            if i != r0 and A[i][col] != 0:
                f = A[i][col] / A[r0][col]; A[i] = [x - f * y for x, y in zip(A[i], A[r0])]
        pivcols.append(col); r0 += 1
    if len(pivcols) < nun: return None, len(pivcols)
    y = [D(0)] * nun
    for i, col in enumerate(pivcols): y[col] = A[i][nun] / A[i][col]
    resid = max((abs(A[i][nun]) for i in range(nun, m)), default=D(0))
    return (y, resid), nun

sol, rk = solve(eqcols)
used_extra = 0
if sol is None:
    sol, rk = solve(eqcols + extra); used_extra = len(extra)
assert sol is not None, ('rank-deficient tight system', rk, nun)
y, resid = sol
ycap = {c: y[i] for i, c in enumerate(caps)}; yw = {w: y[len(caps) + i] for i, w in enumerate(widths)}
minrc = {c: None for c in caps}; argmin = {}; zero = []; npat = 0
for col in allcols:
    kind, c, q = col
    rc = cost(col) - ycap[c] - sum(k * yw[w] for k, w in zip(q, widths))
    if kind == 'head':
        continue
    npat += 1
    if minrc[c] is None or rc < minrc[c]: minrc[c], argmin[c] = rc, q
    if abs(rc) < D('1e-60'): zero.append('%d:%s' % (c, ','.join(str(w) for w, k in zip(widths, q) for _ in range(k)) or 'empty'))
s = {}
for r in headranks:
    rc_h = hcost(r) - ycap[r]
    s[r] = -rc_h if rc_h < 0 else D(0)       # multiplier of x_head_r <= Hmax; the head column's rc becomes >= 0
Dval = sum(ycap[c] * bins[c] for c in caps) + sum(yw[w] * need[w] for w in widths) - sum(s.values()) * Hmax
const = sum(n * pw(r) for r, n in local.items()) - W * pw(73)
neg = sum(bins[c] * minrc[c] for c in caps if minrc[c] < 0)
G_lb = Dval + const + neg
out = dict(inv=a.inv, lp=a.lp, C=C, headranks=headranks, Hmax=Hmax, b_star=str(bstar), tau_star=str(tau),
           patterns_checked=npat, tight_support=len(eqcols), extra_tight_used=used_extra, duals=nun,
           consistency_residual=str(resid), min_reduced_cost={str(c): str(v) for c, v in minrc.items()},
           negative_rc_correction=str(neg), head_multipliers={str(r): str(v) for r, v in s.items()},
           zero_rc_patterns=zero, y_cap={str(c): str(v) for c, v in ycap.items()}, y_w={str(w): str(v) for w, v in yw.items()},
           G_lower_bound_all_fills=str(G_lb), certified=G_lb > 0,
           claim=('every fill of these words, any chronology, priced with this cut, has b < %s' % bstar) if G_lb > 0 else 'NOT certified at this b*')
print(json.dumps({k: out[k] for k in out if k not in ('y_cap', 'y_w', 'zero_rc_patterns')}, indent=1))
if a.out: json.dump(out, open(a.out, 'w'), indent=1)
