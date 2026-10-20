#!/usr/bin/env python3
"""vlp3.py (v-t2-block4, block-3 copy of vlp.py): my own build of the multi-target rank-cover LP (w8 Proposition L, as described in its text),
in my own conventions (vcore.py), plus my own exact dual certificate.
Per array normalisation: c_x = (#moves in direction x)/N, p_xy (x <= y) = max(rank E[x,y], rank E[y,x])/N.
Rows (all as G v >= h):
  deletion     sum_{z in G} c_z >= a_G                                   (every nonzero form G)
  cover        sum_{x in F} p_xx + 2 sum_{x<y in F} p_xy - sum_{x in F} c_x >= -a_F    (every nonzero form F)
  piece<=c     c_x - p_xy >= 0, c_y - p_xy >= 0
  piece<=slack sum_{z in G} c_z - p_xy >= a_G        (every G with x, y in supp G)
Minimise sum_x c_x.  Certificate: y >= 0 rational with G^T y <= (1..1, 0..0) exactly; bound = h.y (weak duality).
Usage: vlp3.py pairs   -> all O_3-classes {I, T}
       vlp.py set r0.r1.r2.r3,...  (rows as bitmasks, bit j = column j; 'I' allowed)"""
import sys, itertools, time
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
sys.path.insert(0, '.')
from vcore import setup, pairing, key, targets_and_orth
D = setup(3); P = D['P']; Qs = D['Qs']; NF = len(Qs)
supp = [np.nonzero(P[q])[0] for q in range(NF)]
cells = [(x, y) for x in range(7) for y in range(x, 7)]; ci = {c: k for k, c in enumerate(cells)}
NV = 7 + len(cells)

def build(aF):
    """aF: list of Fractions (odd fraction per form). Returns sparse rows (list of (dict, rhs))."""
    rows = []
    for q in range(NF):
        S = supp[q]
        rows.append(({int(x): 1 for x in S}, aF[q]))                                    # deletion
        r = {}
        for x in S: r[int(x)] = r.get(int(x), 0) - 1
        for x, y in itertools.combinations_with_replacement(S, 2):
            k = 7 + ci[(int(x), int(y))]; r[k] = r.get(k, 0) + (1 if x == y else 2)
        rows.append((r, -aF[q]))                                                         # cover
        for x, y in itertools.combinations_with_replacement(S, 2):
            r = {int(z): 1 for z in S}; k = 7 + ci[(int(x), int(y))]; r[k] = r.get(k, 0) - 1
            rows.append((r, aF[q]))                                                      # piece <= slack(G)
    for (x, y), k in ci.items():
        rows.append(({x: 1, 7 + k: -1}, Fr(0)))
        if y != x: rows.append(({y: 1, 7 + k: -1}, Fr(0)))                            # piece <= c
    return rows

def solve(targets, certify=True, maxden=10 ** 6):
    N = len(targets)
    cnt = sum(pairing(Qs, S) for S in targets)
    aF = [Fr(int(cnt[q]), N) for q in range(NF)]
    rows = build(aF)
    I_, J_, V_ = [], [], []
    for r, (d, h) in enumerate(rows):
        for k, v in d.items(): I_.append(r); J_.append(k); V_.append(-v)          # -G v <= -h
    A = coo_matrix((V_, (I_, J_)), shape=(len(rows), NV)).tocsr()
    b = np.array([-float(h) for _, h in rows])
    cobj = np.zeros(NV); cobj[:7] = 1
    res = linprog(cobj, A_ub=A, b_ub=b, bounds=[(0, None)] * NV, method='highs')
    out = dict(status=res.status, val=res.fun, nrows=len(rows))
    if certify and res.status == 0:
        y = -res.ineqlin.marginals                       # >= 0 duals of G v >= h
        nz = [r for r in range(len(rows)) if y[r] > 1e-12]
        yr = {r: Fr(float(y[r])).limit_denominator(maxden) for r in nz}
        colsum = [Fr(0)] * NV
        for r, yv in yr.items():
            for k, v in rows[r][0].items(): colsum[k] += yv * v
        cvec = [Fr(1)] * 7 + [Fr(0)] * len(cells)
        viol = [k for k in range(NV) if colsum[k] > cvec[k]]
        bound = sum(yv * rows[r][1] for r, yv in yr.items())
        out.update(cert_rows=len(nz), viol=len(viol), bound=bound if not viol else None)
    return out

def parse(sp):
    if sp == 'I': return np.eye(3, dtype=np.int64)
    rs = [int(v) for v in sp.split('.')]
    return np.array([[(rs[i] >> j) & 1 for j in range(3)] for i in range(3)], dtype=np.int64)

def rows_str(S):
    return '.'.join(str(int(sum(int(S[i, j]) << j for j in range(3)))) for i in range(3))

if __name__ == '__main__':
    t0 = time.time()
    if sys.argv[1] == 'pairs':
        T, O = targets_and_orth(D)
        seen = set(); reps = []
        for k, S in sorted(T.items()):
            if k in seen: continue
            orb = {key((g @ S @ g.T) % 2) for g in O}; seen |= orb; reps.append((S, len(orb)))
        I = np.eye(3, dtype=np.int64)
        for S, sz in reps:
            o = solve([I, S])
            print('T = %-12s class %2d  LP %.9f  cert rows %s viol %s  certified bound %s' % (
                rows_str(S), sz, o['val'], o.get('cert_rows'), o.get('viol'), o.get('bound')), flush=True)
    else:
        tg = [parse(s) for s in sys.argv[2].split(',')]
        o = solve(tg)
        print('targets %s  LP %.9f  rows %d  cert rows %s viol %s  certified bound %s' % (
            sys.argv[2], o['val'], o['nrows'], o.get('cert_rows'), o.get('viol'), o.get('bound')))
    print('time %.1f s' % (time.time() - t0))
