"""caseA_cover.py (d-bitspan) -- the pairwise lemma for CASE A.  If x_i = r(q_i) + w and x_k = r(q_k) + w are two read
vectors of weight <= 6 (support inside 6-flats F_i, F_k), then r(q_i + q_k) = x_i + x_k is supported on
lines(F_i) u lines(F_k).  Here: for every q != 0, is the live read set r(q) a union of two sets of Q-rank <= 6?
WLOG each part is r(q) cap F for a flat F spanned by lines of r(q).  Exact test: for every flat F of rank <= 6
(from out/flats6.pkl if present, else the root-orthogonal E6/D6 flats only), rank(r(q) minus F) <= 6 ?
usage: venv/bin/python -I caseA_cover.py <s2|word> [all]"""
import sys, pickle, os, collections; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import numpy as np
from fractions import Fraction
from e8geo import *
unit = sys.argv[1]
P = {'s2': PUNCT_S2, 'word': PUNCT_WORD, '120': []}[unit]
live = [i for i in range(NL) if i not in P]
U = f2_labels()
def rq(q): return [T for T in live if sum(U[T][k] for k in range(9) if q >> k & 1) % 2]
def rk(idx): return int(np.linalg.matrix_rank(L[idx].astype(float))) if idx else 0
W = {q: rq(q) for q in range(1, 512)}
wts = collections.Counter(rk(W[q]) for q in W)
print('unit', unit, 'label-space weights (rank of r(q)) over the 511 q:', sorted(wts.items()),
      ' sizes |r(q)|:', sorted(collections.Counter(len(W[q]) for q in W).items()))
O = (G == 0)
def perp(lines): return [T for T in range(NL) if all(O[T, l] for l in lines)]
fl = {}
for l in range(NL):
    for m in range(l + 1, NL):
        f = frozenset(perp([l, m]))
        fl[f] = 1
flats6 = [set(f) for f in fl]
src = 'root-orthogonal E6/D6 6-flats (4900)'
if len(sys.argv) > 2 and sys.argv[2] == 'all' and os.path.exists('out/flats6.pkl'):
    flats6 = [set(bits(m)) for m in pickle.load(open('out/flats6.pkl', 'rb'))]
    src = 'ALL rank-6 flats (%d)' % len(flats6)
print('flats used:', src)
cover = {}
for q, R in W.items():
    Rs = set(R); best = 9; wit = None
    for F in flats6:
        rest = [T for T in R if T not in F]
        if len(rest) > 36 and False: continue
        r_ = rk(rest)
        if r_ < best: best, wit = r_, F
        if best <= 6: break
    cover[q] = best
print('min over 6-flats F of rank(r(q) minus F):', sorted(collections.Counter(cover.values()).items()))
print('q with r(q) covered by two 6-flats:', [q for q in cover if cover[q] <= 6])
# Exactness of the float rank: for an integer n x 8 matrix M of rank r (rows = E8 lines, norm^2 4), Cauchy-Binet gives
# prod of the r nonzero singular values^2 = sum of squared r x r minors >= 1, so sigma_min >= sigma_max^-(r-1) >=
# sqrt(4n)^-7 >= 5.4e-10 for n <= 111, while numpy's threshold is sigma_max * n * eps <= 5.2e-13.  So rk() is exact.
print('min |r(q)| over q != 0 = %d  (so one part has >= %d lines, hence is an E6 (36) or D6 (30) flat: every other rank-6 '
      'root system has <= 21 lines)' % (min(len(R) for R in W.values()), (min(len(R) for R in W.values()) + 1) // 2))
