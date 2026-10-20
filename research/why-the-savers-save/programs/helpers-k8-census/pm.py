#!/usr/bin/env python3
"""d-pair-market pm: ONE Hamming block's helper count at the visit level (g-search m1 generalised; HiGHS via scipy).
Needs of a block: every (target, 4-subset S of the block, sign pattern) in the demand.  Under U8's own partition
(u8_blocks.json): 12 face planes with 8 readers and one pattern, 2 own-block planes with 4 readers (Hadamard).
Readers of a face plane are interchangeable, so face options are INTEGER reader counts (aggregated); own-block
needs are binaries.  Options per need:
  P    read the packet row at S's plane
  LL   two pair rows (3 splits), both reads raise the copy
  LF   two pair rows, the read on one side is the target's FIRST read (frame 0 -> that line: its low frame), so
       that copy is not raised (a copy read only first is not terminal at its line)
  H    half read (g-small-1): pre-merge first read of the half on line L inside P from h1, then the other half
       from h2 at P after the merge h1 += h2 (h2 is a holder at P, read raised).  One merge = one full copy + one
       half holder; each half holder serves q readers.
Rows: (support, ratio) [UNSIGNED=1: support only], frame = span of labels.  Copies c_r >= 1 if used, raised copies
>= raised readers / q  [FIX=1: c_r = 2 for every used row, readers ignored = pair_toy's k = 3].
Holders: V(F) >= sum_r c_r (+ half holders at planes) + DUP_F (DUP_F = 1 if a row at F has >= 2 copies; at a plane a
used half holder plays the third register).  Creation: V(F) >= 2 if a row of >= 2 sources is used at F.
Terminal [TERMINAL=1]: raised copies of a LINE row sit in registers that visit no plane afterwards.
Registers: <= 1 line and <= 1 plane, nested; a source only at frames containing its label.  Helpers minimised.
FB = budget of low-frame uses (LF + H) in this block (one per target per call; coupled across blocks outside).
SCOPE=faces drops the own-block needs.  BLOCKS = partition file.  Env: TL (s), Q, FB, OPTS (e.g. P,LL,LF,H).
Usage: VP -I pm.py BLOCK   -> one JSON line"""
import sys, os, json, itertools
from fractions import Fraction as Fr
from collections import defaultdict, Counter
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
here = os.path.dirname(os.path.abspath(sys.argv[0]))
E = os.environ
bi = int(sys.argv[1]); q = int(E.get('Q', 4)); FIX = E.get('FIX') == '1'; UNS = E.get('UNSIGNED') == '1'
TERM = E.get('TERMINAL', '1') == '1'; SCOPE = E.get('SCOPE', 'all'); OPTS = set(E.get('OPTS', 'P').split(','))
FB = E.get('FB'); TL = float(E.get('TL', 60))
d = json.load(open(here + '/a0_demand.json')); ports = d['ports']; v = len(ports)
blocks = json.load(open(E.get('BLOCKS', here + '/u8_blocks.json')))['blocks']; S8 = blocks[bi]
dem = {int(t): {int(s): Fr(x) for s, x in dd.items()} for t, dd in d['demand'].items()}
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return tuple(sorted(bb, reverse=True))
def frame(ss): return basis([ports[s] for s in ss])
def sub(a, b): return len(basis(list(a) + list(b))) == len(b)
def rkey(sup, m):
    ss = sorted(sup)
    if UNS: return tuple(ss)
    x0 = m[ss[0]]; return tuple((s, m[s]/x0) for s in ss)
# needs grouped by (S, pattern)
grp = defaultdict(list)
for t in range(v):
    m = {s: x for s, x in dem[t].items() if s in S8}
    if not m: continue
    assert len(m) == 4, 'this model handles 4-needs only (U8 partition)'
    ss = sorted(m); x0 = m[ss[0]]; grp[(tuple(ss), tuple(m[s]/x0 for s in ss))].append(t)
items = []  # (S, pattern dict, readers list, kind)
for (ss, pat), ts in grp.items():
    kind = 'face' if len(grp_s := [k for k in grp if k[0] == ss]) == 1 else 'own'
    if SCOPE == 'faces' and kind == 'own': continue
    items.append((ss, dict(zip(ss, pat)), ts, kind))
rows = {}
def row(sup, m):
    k = rkey(sup, m)
    if k not in rows: rows[k] = dict(frame=frame(sup), size=len(sup))
    return k
names, cost, lb, ub, integ = [], [], [], [], []
def var(nm, c=0.0, lo=0, hi=1):
    names.append(nm); cost.append(c); lb.append(lo); ub.append(hi); integ.append(1); return len(names) - 1
cons = []
def add(cf, lo, hi): cons.append((cf, lo, hi))
INF = np.inf
raised = defaultdict(dict); first = defaultdict(dict)   # row -> {var: coef}
lowuse = {}                                            # var -> coef (low-frame uses)
halfs = []   # (plane frame, packet row, pre row frame, post row frame, mu var)
optvars = defaultdict(list)
def bump(dct, r, x, k=1): dct[r][x] = dct[r].get(x, 0) + k
for ii, (ss, m, ts, kind) in enumerate(items):
    n = len(ts); vs = []
    S = list(ss); splits = [((S[0], a), tuple(x for x in S if x not in (S[0], a))) for a in S[1:]]
    def ov(nm): x = var((nm, ii), 0, 0, n); vs.append(x); optvars[nm].append(x); return x
    if 'P' in OPTS:
        x = ov('P'); bump(raised, row(ss, m), x)
    if 'LL' in OPTS:
        for p1, p2 in splits:
            x = ov('LL'); bump(raised, row(p1, m), x); bump(raised, row(p2, m), x)
    lowok = not (E.get('FACELOW') == '0' and kind == 'face')
    if 'LF' in OPTS and lowok:
        for p1, p2 in splits:
            for a, b in ((p1, p2), (p2, p1)):
                x = ov('LF'); bump(first, row(a, m), x); bump(raised, row(b, m), x); lowuse[x] = 1
    if 'H' in OPTS and lowok:
        rP = row(ss, m)
        for p1, p2 in splits:
            for a, b in ((p1, p2), (p2, p1)):
                x = ov('H'); lowuse[x] = 1
                mu = var(('mu', ii, a), 0, 0, 8)
                add({mu: q, x: -1}, 0, INF)                  # each half holder serves q readers
                halfs.append((rows[rP]['frame'], rP, frame(a), frame(b), mu, x))
    add({x: 1 for x in vs}, n, n)
rk = list(rows)
U = {r: var(('U', r)) for r in rk}; C = {r: var(('C', r), 0, 0, 16) for r in rk}
CR = {r: var(('CR', r), 0, 0, 16) for r in rk}      # raised copies
for r in rk:
    rd = dict(raised[r]); fs = dict(first[r]); allr = dict(rd)
    for x, k in fs.items(): allr[x] = allr.get(x, 0) + k
    add({U[r]: 16, **{x: -k for x, k in allr.items()}}, 0, INF)       # used if read
    add({C[r]: 1, U[r]: -16}, -INF, 0); add({C[r]: 1, U[r]: -1}, 0, INF)
    if FIX:
        add({C[r]: 1, U[r]: -2}, 0, INF); add({CR[r]: 1, C[r]: -1}, 0, 0) if rows[r]['size'] == 2 else None
    else:
        add({CR[r]: q, **{x: -k for x, k in rd.items()}}, 0, INF)       # q * raised copies >= raised readers
        add({C[r]: 1, CR[r]: -1}, 0, INF)
lines = sorted({r['frame'] for r in rows.values() if len(r['frame']) == 2} | {h[2] for h in halfs} | {h[3] for h in halfs})
pls = sorted({r['frame'] for r in rows.values() if len(r['frame']) == 3})
frames = lines + pls
DUP = {f: var(('DUP', f)) for f in frames}
MUany = {}
YS = {(s, f): var(('ys', s, f)) for s in S8 for f in frames if sub(frame([s]), f)}
YT = {(s, f): var(('yt', s, f)) for s in S8 for f in lines if (s, f) in YS}
HC = {}
for L in lines:
    HC['L', L] = var(('h', 'L', L), 1.0, 0, 32); HC['Lt', L] = var(('h', 'Lt', L), 1.0, 0, 32)
    for P in pls:
        if sub(L, P): HC['LP', L, P] = var(('h', 'LP', L, P), 1.0, 0, 32)
for P in pls: HC['P', P] = var(('h', 'P', P), 1.0, 0, 32)
def V(f):
    cf = {YS[s, f]: 1 for s in S8 if (s, f) in YS}
    for k, i in HC.items():
        if (k[0] in ('L', 'Lt') and k[1] == f) or (k[0] == 'P' and k[1] == f) or (k[0] == 'LP' and (k[1] == f or k[2] == f)): cf[i] = 1
    return cf
# pseudo rows of the half-read device: pre half at line a (h1, non-terminal), post half at line b (h2, non-terminal)
extra_hold = defaultdict(dict); half_at = defaultdict(dict); created = defaultdict(list)
for (Pf, rP, fa, fb, mu, x) in halfs:
    extra_hold[fa][mu] = extra_hold[fa].get(mu, 0) + 1; extra_hold[fb][mu] = extra_hold[fb].get(mu, 0) + 1
    half_at[Pf][mu] = half_at[Pf].get(mu, 0) + 1
    add({U[rP]: 8, x: -1}, 0, INF)                                    # the packet row is created
    created[fa].append(mu); created[fb].append(mu)
    add({DUP[fa]: 8, mu: -1}, -1, INF); add({DUP[fb]: 8, mu: -1}, -1, INF)   # mu >= 2 at a line -> third register
# full copies of a plane's packet row >= merges feeding half holders there
fullneed = defaultdict(dict)
for (Pf, rP, fa, fb, mu, x) in halfs: fullneed[rP][mu] = fullneed[rP].get(mu, 0) + 1
for rP, cf in fullneed.items(): add({C[rP]: 1, **{m_: -k for m_, k in cf.items()}}, 0, INF)
for f in frames:
    rf = [r for r in rk if rows[r]['frame'] == f]
    Vf = V(f)
    for r in rf:
        if rows[r]['size'] >= 2: add({**Vf, U[r]: Vf.get(U[r], 0) - 2}, 0, INF)
        add({C[r]: 1, DUP[f]: -16}, -INF, 1)
    for mu in created[f]:
        ind = var(('cr', f, mu)); add({ind: 8, mu: -1}, 0, INF); add({**Vf, ind: -2}, 0, INF)
    hold = dict(Vf)
    for r in rf: hold[C[r]] = hold.get(C[r], 0) - 1
    for mu, k in extra_hold[f].items(): hold[mu] = hold.get(mu, 0) - k
    h1 = dict(hold); h1[DUP[f]] = h1.get(DUP[f], 0) - 1; add(h1, 0, INF)
    if half_at[f]:
        h2 = dict(hold)
        for mu, k in half_at[f].items(): h2[mu] = h2.get(mu, 0) - k
        add(h2, 0, INF)
    if len(f) == 2 and TERM:
        cf = {YT[s, f]: 1 for s in S8 if (s, f) in YT}; cf[HC['Lt', f]] = 1
        for r in rf:
            if rows[r]['size'] == 2: cf[CR[r]] = cf.get(CR[r], 0) - 1
        add(cf, 0, INF)
for (s, f), y in YT.items(): add({YS[s, f]: 1, y: -1}, 0, INF)
for s in S8:
    add({YS[s, f]: 1 for f in lines if (s, f) in YS}, 0, 1)
    add({YS[s, f]: 1 for f in pls if (s, f) in YS}, 0, 1)
    for L in lines:
        for P in pls:
            if (s, L) in YS and (s, P) in YS and not sub(L, P): add({YS[s, L]: 1, YS[s, P]: 1}, 0, 1)
            if (s, L) in YT and (s, P) in YS: add({YT[s, L]: 1, YS[s, P]: 1}, 0, 1)
if FB is not None and lowuse: add(dict(lowuse), 0, int(FB))
LAM = float(E.get('LAM', 0))
for x in lowuse: cost[x] += LAM          # Lagrangian price per low-frame use (per-target budget dualised)
Ar, Ac, Av, lo, hi = [], [], [], [], []
for k, (cf, l, h) in enumerate(cons):
    for c_, x in cf.items(): Ar.append(k); Ac.append(c_); Av.append(float(x))
    lo.append(l); hi.append(h)
A = coo_matrix((Av, (Ar, Ac)), shape=(len(cons), len(names))).tocsr()
res = milp(np.array(cost), constraints=[LinearConstraint(A, np.array(lo, float), np.array(hi, float))],
           integrality=np.array(integ), bounds=Bounds(np.array(lb, float), np.array(ub, float)),
           options=dict(time_limit=TL, disp=False))
out = dict(block=bi, q=q, FIX=FIX, UNSIGNED=UNS, TERM=TERM, SCOPE=SCOPE, OPTS=sorted(OPTS), FB=FB, status=int(res.status),
           helpers=None if res.x is None else int(round(sum(res.x[i] for i in range(len(names)) if names[i][0] == 'h'))), objective=None if res.x is None else float(res.fun), LAM=E.get('LAM'), FACELOW=E.get('FACELOW'), bound=getattr(res, 'mip_dual_bound', None),
           needs=sum(len(i[2]) for i in items), items=len(items), rows=len(rows))
if res.x is not None:
    xv = res.x
    out['options'] = {k: int(round(sum(xv[i] for i in vs))) for k, vs in optvars.items()}
    out['low_uses'] = int(round(sum(xv[i] for i in lowuse)))
    hk = Counter()
    for k, i in HC.items(): hk[k[0]] += int(round(xv[i]))
    out['helper_kinds'] = dict(hk)
    out['source_visits'] = dict(Counter(len(f) for (s, f), i in YS.items() if xv[i] > .5))
    out['half_merges'] = int(round(sum(xv[h[4]] for h in halfs)))
    if E.get('DUMP'):
        json.dump(dict(out=out, sol={str(names[i]): round(float(xv[i])) for i in range(len(names)) if abs(xv[i]) > .5}),
                  open(E['DUMP'], 'w'), indent=0, default=str)
print(json.dumps(out, default=str))
