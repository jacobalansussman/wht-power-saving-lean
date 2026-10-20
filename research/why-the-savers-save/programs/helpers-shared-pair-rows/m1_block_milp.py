#!/usr/bin/env python3
"""g-search m1: the U8 helper network of ONE Hamming block as an integer programme (our tools/e8 xmilp style:
a pool of signed rows, every target need covered exactly once by uses, rows built, copies; solved by HiGHS via
scipy.optimize.milp), at the hub's VISIT level (ingredients of a created row ride inside the same chains, as in the
hub's packet-plane bound; this is a lower-bound/what-if model, not a scalar construction).

Needs (from a0/a2): every (target, block-subset, sign pattern) of the block: 12 face planes (one pattern, 6-8
readers), 2 own-block planes (4 readers, 4 patterns), and whole-block needs (8 sources).
Options per 4-need: 'P' one packet row at its plane; 'LL' two pair rows (3 splits); 'S' four single rows.
Options per 8-need: two packet rows of a parallel split (7), or 'S' 8 singles.
Rows: support (with normalised signs) and frame = span of the support's labels (rank 1 line, 2 line, 3 plane).
Copies: c_r >= readers(r)/q.  Visitors at frame F: >= 2 if a row of >= 2 sources is created there; >= sum c_r + 1
if some row at F has >= 2 copies (the hub's third-visitor lemma), else >= sum c_r.
Registers: each source / helper visits at most one rank-2 and one rank-3 frame, nested (L inside P); a source can
visit F only if its label is in F. A register that HOLDS a retained copy at a line L for target reads is terminal
(does not go on to a plane) unless variant 'first': every reader reads that row as its FIRST read (target frame 0,
so the read does not raise the copy), and each target has at most one first read.
Singles: a rank-1 row x_s needs helper copies (c >= readers/q), each a terminal helper.
Usage: venv-python -I m1_block_milp.py BLOCK [q] [variants...]   variants: LL S first nosrcterm"""
import sys, json, itertools
from fractions import Fraction as Fr
from collections import defaultdict
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
HERE = sys.path[0] if False else '.'
bi = int(sys.argv[1]); q = int(sys.argv[2]) if len(sys.argv) > 2 else 4; var = set(sys.argv[3:])
d = json.load(open('a0_demand.json')); ports = d['ports']; v = len(ports)
blocks = json.load(open('a1_blocks.json'))['blocks']; S8 = blocks[bi]
dem = {int(t): {int(s): Fr(x) for s, x in dd.items()} for t, dd in d['demand'].items()}
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return tuple(sorted(bb, reverse=True))
def span_key(ss): return basis([ports[s] for s in ss])
def canon(sup, sg):
    ss = sorted(sup); x0 = sg[ss[0]]
    return tuple((s, int(sg[s]/x0)) for s in ss)
needs = []
for t in range(v):
    m = {s: x for s, x in dem[t].items() if s in S8}
    if m: needs.append((t, m))
planes = []
for quad in itertools.combinations(S8, 4):
    if len(span_key(quad)) == 3 and (lambda L: L[0] ^ L[1] ^ L[2] ^ L[3] == 0)([ports[s] for s in quad]): planes.append(frozenset(quad))
assert len(planes) == 14, len(planes)
rows = {}  # key -> dict(frame, size)
def row(sup, sg):
    k = canon(sup, sg)
    if k not in rows: rows[k] = dict(frame=span_key(sup), size=len(sup))
    return k
opts = []  # per need: list of (name, [rows])
for t, m in needs:
    o = []
    S = frozenset(m)
    if len(S) == 4:
        o.append(('P', [row(S, m)]))
        if 'LL' in var:
            a = sorted(S)
            for p1 in [(a[0], a[1]), (a[0], a[2]), (a[0], a[3])]:
                p2 = tuple(x for x in a if x not in p1)
                o.append(('LL', [row(p1, m), row(p2, m)]))
    else:
        for P in planes:
            if P <= S and min(P) == min(S):  # split S into P and S-P (parallel planes)
                Q = S - P
                if Q in planes: o.append(('PP', [row(P, m), row(Q, m)]))
    if 'S' in var: o.append(('S', [row([s], m) for s in sorted(S)]))
    opts.append(o)
frames = sorted({r['frame'] for r in rows.values()}, key=lambda f: (len(f), f))
rk = list(rows)
lines = [f for f in frames if len(f) == 2]; pls = [f for f in frames if len(f) == 3]
def inF(s, f): return len(basis(list(f) + [ports[s]])) == len(f)
def sub(a, b): return len(basis(list(a) + list(b))) == len(b)
# ---- variables
names = []; cost = []; integ = []; lb = []; ub = []
def var_(nm, c=0.0, lo=0, hi=1, it=1):
    names.append(nm); cost.append(c); lb.append(lo); ub.append(hi); integ.append(it); return len(names) - 1
X = {}
for i, o in enumerate(opts):
    for j in range(len(o)): X[i, j] = var_(('x', i, j))
U = {r: var_(('u', r)) for r in rk}
C = {r: var_(('c', r), 0, 0, 64) for r in rk}
DUP = {f: var_(('dup', f)) for f in frames}
FIRST = {}
if 'first' in var:
    for i, o in enumerate(opts):
        for j, (nm, rr) in enumerate(o):
            for r in rr:
                if rows[r]['size'] == 2: FIRST[i, r] = var_(('first', i, r))
    NT = {r: var_(('nonterm', r)) for r in rk if rows[r]['size'] == 2}   # row read only as first reads -> copy not terminal
# sources: visits
YS = {}
for s in S8:
    for f in lines + pls:
        if inF(s, f):
            YS[s, f] = var_(('ys', s, f))
YT = {}  # source holds terminal copy at line f
for s in S8:
    for f in lines:
        if inF(s, f): YT[s, f] = var_(('yt', s, f))
# helpers: chains (L, P) with L in P, (L,), (P,), terminal line holder (Lt,), single copy holder ('1',)
HC = {}
for L in lines:
    HC[('L', L)] = var_(('h', 'L', L), 1.0, 0, 64)
    HC[('Lt', L)] = var_(('h', 'Lt', L), 1.0, 0, 64)
    for P in pls:
        if sub(L, P): HC[('LP', L, P)] = var_(('h', 'LP', L, P), 1.0, 0, 64)
for P in pls: HC[('P', P)] = var_(('h', 'P', P), 1.0, 0, 64)
HS = var_(('h', 'single'), 1.0, 0, 512)
rowsC = []
def add(coefs, lo, hi): rowsC.append((coefs, lo, hi))
INF = np.inf
for i, o in enumerate(opts): add({X[i, j]: 1 for j in range(len(o))}, 1, 1)
readers = defaultdict(dict)
for i, o in enumerate(opts):
    for j, (nm, rr) in enumerate(o):
        for r in rr:
            add({U[r]: 1, X[i, j]: -1}, 0, INF)
            readers[r][X[i, j]] = readers[r].get(X[i, j], 0) + 1
for r in rk:
    cf = {C[r]: q}; 
    for x, k in readers[r].items(): cf[x] = cf.get(x, 0) - k
    add(cf, 0, INF)                      # q c_r >= readers
    add({C[r]: 1, U[r]: -64}, -INF, 0)   # copies only if used
    add({C[r]: 1, U[r]: -1}, 0, INF)     # used -> >= 1 copy
for f in frames:
    for r in rk:
        if rows[r]['frame'] == f: add({C[r]: 1, DUP[f]: -64}, -INF, 1)  # c_r <= 1 + 64 dup_f
# visitors at frame f
def visitors(f):
    cf = {}
    for s in S8:
        if (s, f) in YS: cf[YS[s, f]] = 1
    if len(f) == 2:
        cf[HC[('L', f)]] = 1
        cf[HC[('Lt', f)]] = 1
        for P in pls:
            if ('LP', f, P) in HC: cf[HC[('LP', f, P)]] = 1
    else:
        cf[HC[('P', f)]] = 1
        for L in lines:
            if ('LP', L, f) in HC: cf[HC[('LP', L, f)]] = 1
    return cf
for f in frames:
    rf = [r for r in rk if rows[r]['frame'] == f]
    if len(f) == 1:
        # singles: helper copies; source itself not used as a reader copy (it is restored / climbs to its ladder)
        continue
    V = visitors(f)
    for r in rf:
        cf = dict(V); cf[U[r]] = cf.get(U[r], 0) - 2; add(cf, 0, INF)     # creation: 2 visitors
    cf = dict(V)
    for r in rf: cf[C[r]] = cf.get(C[r], 0) - 1
    cf[DUP[f]] = cf.get(DUP[f], 0) - 1
    add(cf, 0, INF)                                                        # holders (+1 if duplicated)
    if len(f) == 2:
        # copy holders at a line are terminal unless the row is read only as first reads
        cf = {}
        for s in S8:
            if (s, f) in YT: cf[YT[s, f]] = 1
        cf[HC[('Lt', f)]] = 1
        for r in rf:
            if 'first' in var:
                tc = var_(('tc', r), 0, 0, 64)            # terminal copies of r: tc >= c_r - 64 NT_r (fixed per row)
                add({tc: 1, C[r]: -1, NT[r]: 64}, 0, INF)
                cf[tc] = cf.get(tc, 0) - 1
            else:
                cf[C[r]] = cf.get(C[r], 0) - 1
        add(cf, 0, INF)
# singles
cf = {HS: 1}
for r in rk:
    if rows[r]['size'] == 1: cf[C[r]] = -1
add(cf, 0, INF)
# terminal holders are also visitors at their line? they count as holders: add them to the holder constraint
# (a terminal holder visits the line).  Implemented by requiring YT <= YS and Lt helpers counted as visitors:
for (s, f), y in YT.items(): add({YS[s, f]: 1, y: -1}, 0, INF)
# NOTE: Lt helpers count as visitors at f:
for f in lines:
    pass
# source nesting: <=1 line, <=1 plane, nested; a terminal line holder visits no plane
for s in S8:
    add({YS[s, f]: 1 for f in lines if (s, f) in YS}, 0, 1)
    add({YS[s, f]: 1 for f in pls if (s, f) in YS}, 0, 1)
    for L in lines:
        for P in pls:
            if (s, L) in YS and (s, P) in YS and not sub(L, P): add({YS[s, L]: 1, YS[s, P]: 1}, 0, 1)
        if (s, L) in YT:
            for P in pls:
                if (s, P) in YS: add({YT[s, L]: 1, YS[s, P]: 1}, 0, 1)
# first reads: target has at most one first read; NT_r = 1 requires every reader of r to read it first
if 'first' in var:
    byt = defaultdict(dict)
    for (i, r), fv in FIRST.items():
        t = needs[i][0]; byt[t][fv] = 1
        # first var only if the option containing r is chosen
        add({fv: 1, **{X[i, j]: -1 for j, (nm, rr) in enumerate(opts[i]) if r in rr}}, -INF, 0)
    for t, cf in byt.items(): add(cf, 0, 1)
    FB = __import__('os').environ.get('FB')       # global budget of first reads in this block (a target has ONE first read per call)
    if FB is not None: add({fv: 1 for fv in FIRST.values()}, 0, int(FB))
    for r in NT:
        # nonterm_r <= first_(i,r) for every chosen reader i:  nonterm_r + x_ij - first_(i,r) <= 1
        for i, o in enumerate(opts):
            for j, (nm, rr) in enumerate(o):
                if r in rr: add({NT[r]: 1, X[i, j]: 1, FIRST[i, r]: -1}, -INF, 1)
# Lt helpers visit their line: include in visitors
# (patch: add Lt to every visitor constraint at that line)
for k, (cf, lo, hi) in enumerate(rowsC):
    pass
A_r, A_c, A_v, lo, hi = [], [], [], [], []
for k, (cf, l, h) in enumerate(rowsC):
    for c_, x in cf.items(): A_r.append(k); A_c.append(c_); A_v.append(float(x))
    lo.append(l); hi.append(h)
A = coo_matrix((A_v, (A_r, A_c)), shape=(len(rowsC), len(names))).tocsr()
res = milp(np.array(cost), constraints=[LinearConstraint(A, np.array(lo, float), np.array(hi, float))],
           integrality=np.array(integ), bounds=Bounds(np.array(lb, float), np.array(ub, float)), options=dict(time_limit=float(__import__("os").environ.get("TL","120")), disp=False))
xv = res.x
out = dict(block=bi, q=q, variants=sorted(var), status=int(res.status), bound=getattr(res, "mip_dual_bound", None), helpers=None if xv is None else round(res.fun),
           needs=len(needs), rows_in_pool=len(rows))
if xv is not None:
    from collections import Counter
    ch = Counter(); 
    for i, o in enumerate(opts):
        for j, (nm, rr) in enumerate(o):
            if xv[X[i, j]] > .5: ch[(nm, len(needs[i][1]))] += 1
    out['options_chosen'] = {f'{a}/{b}': x for (a, b), x in sorted(ch.items())}
    hk = Counter()
    for k_, i in HC.items():
        if xv[i] > .5: hk[k_[0]] += int(round(xv[i]))
    hk['single'] = int(round(xv[HS]))
    out['helper_kinds'] = dict(hk)
    out['source_visits'] = dict(Counter(len(f) for (s, f), i in YS.items() if xv[i] > .5))
    out['source_terminal_line_holders'] = sum(1 for k_, i in YT.items() if xv[i] > .5)
print(json.dumps(out))
