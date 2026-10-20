#!/usr/bin/env python3
"""mtlp2.py (w8-rmexact): multi-target rank-cover LP with general partitions (pieces R x C of supp F).
Same relaxation as mtlp.py (virtual slice of the H-game: virtual targets h h^T).  For each form F and each chosen
partition {R_1..R_p} of supp F:  slack(F) <= sum_{i,j} p(R_i, R_j)   (p symmetric: one variable for {R_i,R_j}).
Each piece: p(R,C) <= c(R), p(R,C) <= c(C), p(R,C) <= slack(G) for every form G with R u C in supp G (all G, or those
with |supp G| <= --gmax).  slack(G) = sum_{x in G} c_x - a_G >= 0.  Partition families: cells for every F; all set
partitions for |supp F| <= --allpart; partitions into <= --kblocks blocks for |supp F| <= --kmax.
Options as mtlp.py for H.  --uniform forces equal c_x; --orbit forces c constant on H-orbits of directions."""
import sys, argparse, itertools
sys.path.insert(0, '.')
from forms import setup, pair
ap_ = argparse.ArgumentParser()
ap_.add_argument('--m', type=int, default=4)
ap_.add_argument('--H', default='trivial')
ap_.add_argument('--poly', type=int, default=None)
ap_.add_argument('--uniform', action='store_true')
ap_.add_argument('--conj', type=int, default=-1)
ap_.add_argument('--allpart', type=int, default=4)
ap_.add_argument('--kblocks', type=int, default=2)
ap_.add_argument('--kmax', type=int, default=6)
ap_.add_argument('--gmax', type=int, default=99)
ap_.add_argument('--dump', default='')
args = ap_.parse_args()
if args.poly is None: args.poly = {3: 0b1011, 4: 0b10011, 5: 0b100101}[args.m]
m = args.m
D = setup(m); GL = D['GL']; supp = D['supp']; nf = D['nf']; dirs = D['dirs']
mm, tr = D['matmul'], D['transpose']
exec(open('hgroups.py').read())
H = make_H(args, D)
n = len(H)
targets = [D['enc_sym'](mm(h, tr(h))) for h in H]
from fractions import Fraction as Fr
aF = [Fr(sum(pair(F, S) for S in targets), n) for F in range(nf)]
a = [float(x) for x in aF]
nd = len(dirs); di = {x: i for i, x in enumerate(dirs)}
S_of = [frozenset(di[x] for x in supp[F]) for F in range(nf)]
def set_partitions(lst, kmax=None):
    if not lst: yield []; return
    first, rest = lst[0], lst[1:]
    for p in set_partitions(rest, kmax):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        if kmax is None or len(p) < kmax: yield [[first]] + p
pieces = {}
def piece(R, C):
    key = (R, C) if (sorted(R), len(R)) <= (sorted(C), len(C)) else (C, R)
    key = tuple(sorted([tuple(sorted(R)), tuple(sorted(C))]))
    if key not in pieces: pieces[key] = nd + len(pieces)
    return pieces[key]
rows = []
for F in range(1, nf):
    P = sorted(S_of[F]); s = len(P)
    parts = [[[x] for x in P]]
    if s <= args.allpart: parts = list(set_partitions(P))
    elif s <= args.kmax: parts += [p for p in set_partitions(P, args.kblocks) if len(p) > 1 or True]
    for part in parts:
        row = {}
        for i in P: row[i] = row.get(i, 0) + 1
        for X in part:
            for Y in part:
                k = piece(frozenset(X), frozenset(Y)); row[k] = row.get(k, 0) - 1
        rows.append((row, a[F]))
for G in range(1, nf):
    rows.append(({i: -1 for i in S_of[G]}, -a[G]))
# piece bounds
Gs = [G for G in range(1, nf) if len(S_of[G]) <= args.gmax]
bysize = sorted(Gs, key=lambda G: len(S_of[G]))
for key, k in list(pieces.items()):
    R, C = frozenset(key[0]), frozenset(key[1]); U = R | C
    rows.append(({k: 1, **{i: -1 for i in R}}, 0))
    rows.append(({k: 1, **{i: -1 for i in C}}, 0))
    for G in bysize:
        if U <= S_of[G]:
            row = {i: -1 for i in S_of[G]}; row[k] = 1
            rows.append((row, -a[G]))
nv = nd + len(pieces)
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
ri, ci, vv = [], [], []
b = np.zeros(len(rows))
for r, (row, rhs) in enumerate(rows):
    for k, v in row.items(): ri.append(r); ci.append(k); vv.append(v)
    b[r] = rhs
A = coo_matrix((vv, (ri, ci)), shape=(len(rows), nv)).tocsr()
cobj = np.zeros(nv); cobj[:nd] = 1
Aeq = beq = None
if args.uniform:
    Aeq = coo_matrix(([1] * (nd - 1) + [-1] * (nd - 1), (list(range(nd - 1)) * 2, [0] * (nd - 1) + list(range(1, nd)))), shape=(nd - 1, nv)).tocsr()
    beq = np.zeros(nd - 1)
res = linprog(cobj, A_ub=A, b_ub=b, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * nv, method='highs')
print("m=%d H=%s conj=%d |H|=%d targets=%d uniform=%s allpart=%d kblocks=%d kmax=%d gmax=%d rows=%d vars=%d: status %s, LP M/W >= %.6f" % (
    m, args.H, args.conj, n, len(set(targets)), args.uniform, args.allpart, args.kblocks, args.kmax, args.gmax, len(rows), nv, res.status,
    res.fun if res.status == 0 else float('nan')), flush=True)
