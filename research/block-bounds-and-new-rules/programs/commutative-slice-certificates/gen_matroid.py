# pf-logging, 13 Oct 2026.  GENERATOR (untrusted) of exclusion certificates for g-play's commutative-slice test
# (wave16/g-play/matroid.py): a multiset z_1..z_n of nonzero directions in F_2^m can serve W arrays only if, after
# deleting ANY W-1 of the n moves, the kernel's class kappa lies in the F_2-span of the remaining direction classes
# in D/B_D.  Here the classes are read through the homomorphism Phi : D -> F_2^{m + C(m,2)} of PROOF.md
#   Phi(q) = ( q(e_j) - q(0) mod 2 )_j  ,  ( (q(0) - q(e_j) - q(e_l) + q(e_j + e_l)) / 2 mod 2 )_{j<l}
# with Phi(lift z) = (z ; z_j z_l) and Phi(kappa) = (1...1 ; 0).  B_D lies in ker Phi, so a functional y with
# y.Phi(kappa) = 1 and y.Phi(z_i) = 0 for all but at most W-1 positions i proves the multiset fails the test.
# For every multiset this writes one y (uint16, canonical combinations_with_replacement order), or 0xFFFF if none
# exists (then the multiset PASSES the Phi-relaxation and is listed in the text output).
# Usage: python3 (with numpy) gen_matroid.py m W n OUTFILE
import sys, itertools
import numpy as np
m, W, n = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); out = sys.argv[4]
pairs = list(itertools.combinations(range(m), 2)); N = m + len(pairs)
def phi_dir(z):
    v = z
    for b, (j, l) in enumerate(pairs):
        if (z >> j) & 1 and (z >> l) & 1: v |= 1 << (m + b)
    return v
kap = (1 << m) - 1
dirs = list(range(1, 1 << m))
PH = [phi_dir(z) for z in dirs]
par = lambda x: bin(x).count('1') & 1
ys = [y for y in range(1 << N) if par(y & kap)]                       # y.Phi(kappa) = 1
bad = np.array([[par(y & ph) for y in ys] for ph in PH], dtype=np.float32)   # (15, |ys|): 1 if y.Phi(z) = 1
from math import comb
total = comb(len(dirs) + n - 1, n)
print(f"m={m} W={W} n={n}: N={N}, multisets {total}, candidate functionals {len(ys)}", flush=True)
res = np.empty(total, dtype=np.uint16); passing = []
CH = 200000; it = itertools.combinations_with_replacement(range(len(dirs)), n); st = 0
while True:
    blk = list(itertools.islice(it, CH))
    if not blk: break
    cnt = np.zeros((len(blk), len(dirs)), dtype=np.float32)
    for r, ms in enumerate(blk):
        for d in ms: cnt[r, d] += 1
    viol = cnt @ bad                                                   # positions with y.Phi(z_i) = 1
    best = viol.argmin(axis=1); mn = viol[np.arange(len(blk)), best]
    for r in range(len(blk)):
        if mn[r] <= W - 1: res[st + r] = ys[best[r]]
        else: res[st + r] = 0xFFFF; passing.append(tuple(dirs[d] for d in blk[r]))
    st += len(blk)
assert st == total
res.tofile(out)
print(f"certified failing: {int((res != 0xFFFF).sum())}, passing the Phi-test: {len(passing)}")
for h in passing[:40]: print("  passes:", h)
