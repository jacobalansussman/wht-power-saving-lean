#!/usr/bin/env python3
"""mtlp_t.py (t2-block4, copied from w8-rmexact mtlp.py; only the target list changed: explicit --targets I,r0.r1.r2.r3,...).\nmtlp.py (w8-rmexact): the multi-target rank-cover LP.
By Lemma C of e3-rmbuild, a unit of C[Sym_m x| H] (moves z_t) gives, in the channel representation Ind chi_Q = C^H,
the order-2 system (*) of ex2-families for Wn virtual arrays and Mn virtual moves, where virtual array (l,h) has target
S_h = h h^T and virtual move (t,h) has direction h z_t (dropping the left-regular structure is a relaxation).
Theorem 10 holds verbatim with a target per array: rank E[D_F, D_F] = |D_F| - #{virtual arrays with <F,S> = 1}.
Rank cover (Theorem 17) with the cell partition of every support gives the LP below (normalised per virtual array):
  minimise sum_x c_x  s.t.  for every form F:   sum_{x in F} c_x - a_F <= sum_{x in F} p_xx + 2 sum_{x<y in F} p_xy
                            p_xy <= c_x, p_xy <= c_y, p_xy <= sum_{z in G} c_z - a_G  (every G containing x, y)
                            sum_{z in G} c_z >= a_G (deletion)
  where a_F = fraction of virtual arrays odd at F.  Options: --H trivial|singer:<power>|GL|perm|orth ; --uniform forces
  c_x equal (the GL case has uniform virtual directions anyway).
Prints the LP optimum (a lower bound on M/W for the H-game)."""
import sys, argparse, itertools
sys.path.insert(0, '.')
from forms import setup, pair
ap_ = argparse.ArgumentParser()
ap_.add_argument('--m', type=int, default=4)
ap_.add_argument('--targets', default='I')
ap_.add_argument('--poly', type=int, default=None)
ap_.add_argument('--uniform', action='store_true')
ap_.add_argument('--conj', type=int, default=-1, help='conjugate H by the conj-th element of GL')
ap_.add_argument('--evenflats', action='store_true', help='only the 4-point supports as F (Theorem 17 rows)')
ap_.add_argument('--dump', default='')
args = ap_.parse_args()
if args.poly is None: args.poly = {3: 0b1011, 4: 0b10011, 5: 0b100101}[args.m]
m = args.m
D = setup(m); GL = D['GL']; supp = D['supp']; nf = D['nf']; dirs = D['dirs']
mm, tr = D['matmul'], D['transpose']
ap2 = None
targets_spec = args.targets.split(',')
def _T(sp):
    if sp == 'I': return [1 << i for i in range(m)]
    return [int(v) for v in sp.split('.')]
targets = [D['enc_sym'](_T(sp)) for sp in targets_spec]
n = len(targets)
a = [sum(pair(F, S) for S in targets) / n for F in range(nf)]
from fractions import Fraction as Fr
aF = [Fr(sum(pair(F, S) for S in targets), n) for F in range(nf)]
nd = len(dirs); di = {x: i for i, x in enumerate(dirs)}
pairs = [(i, j) for i in range(nd) for j in range(i, nd)]; pi = {p: k for k, p in enumerate(pairs)}
nv = nd + len(pairs)
rows = []   # (dict var->coef, rhs) meaning sum coef*var <= rhs
Fs = [F for F in range(1, nf) if (len(supp[F]) == 4 or not args.evenflats)]
for F in Fs:
    P = sorted(di[x] for x in supp[F])
    row = {}
    for i in P: row[i] = row.get(i, 0) + 1
    for i in P:
        k = nd + pi[(i, i)]; row[k] = row.get(k, 0) - 1
    for i, j in itertools.combinations(P, 2):
        k = nd + pi[(i, j)]; row[k] = row.get(k, 0) - 2
    rows.append((row, a[F], ('cover', F)))
for G in range(1, nf):
    P = sorted(di[x] for x in supp[G])
    rows.append(({i: -1 for i in P}, -a[G], ('del', G)))
    for i, j in itertools.combinations_with_replacement(P, 2):
        row = {i_: -1 for i_ in P}; row[nd + pi[(i, j)]] = row.get(nd + pi[(i, j)], 0) + 1
        rows.append((row, -a[G], ('piece<=slack', i, j, G)))
for (i, j), k in pi.items():
    rows.append(({nd + k: 1, i: -1}, 0, ('piece<=c', i, j)))
    if j != i: rows.append(({nd + k: 1, j: -1}, 0, ('piece<=c', j, i)))
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
A = lil_matrix((len(rows), nv)); b = np.zeros(len(rows))
for r, (row, rhs, lab) in enumerate(rows):
    for k, v in row.items(): A[r, k] = v
    b[r] = rhs
cobj = np.zeros(nv); cobj[:nd] = 1
Aeq = None; beq = None
if args.uniform:
    Aeq = lil_matrix((nd - 1, nv))
    for i in range(1, nd): Aeq[i - 1, 0] = 1; Aeq[i - 1, i] = -1
    beq = np.zeros(nd - 1)
res = linprog(cobj, A_ub=A.tocsr(), b_ub=b, A_eq=Aeq.tocsr() if Aeq is not None else None, b_eq=beq,
              bounds=[(0, None)] * nv, method='highs')
from collections import Counter
print("m=%d targets=%s conj=%d |T|=%d distinct targets=%d uniform=%s evenflats=%s rows=%d vars=%d: status %s, LP optimum M/W >= %.6f" % (
    m, args.targets, args.conj, n, len(set(targets)), args.uniform, args.evenflats, len(rows), nv, res.status, res.fun if res.status == 0 else float('nan')))
if res.status == 0:
    cs = res.x[:nd]
    print("  c_x:", " ".join("%.4f" % v for v in cs))
    if args.dump:
        import json
        y = res.ineqlin.marginals
        json.dump({'rows': [[{str(k): v for k, v in row.items()}, rhs, list(map(str, lab))] for (row, rhs, lab) in rows],
                   'dual': list(map(float, y)), 'opt': res.fun, 'nd': nd, 'nv': nv,
                   'aF': {str(F): str(aF[F]) for F in range(nf)}}, open(args.dump, 'w'))
