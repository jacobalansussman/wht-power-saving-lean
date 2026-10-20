# exactsx.py W s file.npz [maxden] -- recognise every coordinate of a numerical stand-in solution as a Gaussian rational
# (simplest fraction within 1e-9, denominator <= maxden), then evaluate the program EXACTLY over Q(i) (Python Fractions)
# and compare with the exact target rho(D1 D0) in each block.  The proof is the exact evaluation; numerics only proposes.
import sys, numpy as np
from fractions import Fraction as Fr
exec(open('sxcore.py').read())
W, s = int(sys.argv[1]), int(sys.argv[2]); z = np.load(sys.argv[3]); p = z['p']
maxden = int(sys.argv[4]) if len(sys.argv) > 4 else 1 << 12
S = Sys(W, s)
def rec(x):
    a = Fr(float(x)).limit_denominator(maxden)
    return a if abs(float(a) - x) < 1e-9 else None
vals = []; bad = 0
for x in p:
    r, i = rec(x.real), rec(x.imag)
    if r is None or i is None: bad += 1; vals.append(None)
    else: vals.append((r, i))
print('coordinates %d, nonzero %d, not recognised %d' % (len(p), sum(1 for v in vals if v and (v[0] or v[1])), bad), flush=True)
if bad: sys.exit()
dens = sorted(set(v[0].denominator for v in vals) | set(v[1].denominator for v in vals)); print('denominators', dens)
def gq(M):   # exact Gaussian-rational matrix from a numpy complex matrix with dyadic entries
    re = np.empty(M.shape, object); im = np.empty(M.shape, object)
    for ix in np.ndindex(M.shape): re[ix] = Fr(float(M[ix].real)); im[ix] = Fr(float(M[ix].imag))
    return re, im
def mul(A, B): return (A[0].dot(B[0]) - A[1].dot(B[1]), A[0].dot(B[1]) + A[1].dot(B[0]))
def add(A, B): return (A[0] + B[0], A[1] + B[1])
def zero(shape):
    z_ = np.empty(shape, object); z_.fill(Fr(0)); return z_, z_.copy()
PaQ = [gq(Pa[a]) for a in range(dimA)]; RQ = gq(R); EQ = [gq(e) for e in S.E]; TQ = [gq(t) for t in S.T]
k = 0; coefs = []
for n_ in S.slots:
    row = []
    for i in range(n_):
        Mq = zero((d, d))
        for a in range(dimA):
            vr, vi = vals[k]; k += 1
            if vr or vi:
                Mq = add(Mq, (PaQ[a][0] * vr - PaQ[a][1] * vi, PaQ[a][0] * vi + PaQ[a][1] * vr))
        row.append(Mq)
    coefs.append(row)
srcs = list(EQ)
for j in range(s):
    c = zero((d, d * W))
    for i in range(len(srcs)): c = add(c, mul(coefs[j][i], srcs[i]))
    srcs.append(mul(RQ, c))
ok = True
for t in range(W):
    o = zero((d, d * W))
    for i in range(len(srcs)): o = add(o, mul(coefs[s + t][i], srcs[i]))
    diff = sum(1 for ix in np.ndindex(o[0].shape) if o[0][ix] != TQ[t][0][ix] or o[1][ix] != TQ[t][1][ix])
    print('output %d: %d of %d entries differ from the exact target' % (t, diff, o[0].size), flush=True)
    ok &= diff == 0
print('EXACT VERDICT: %s' % ('a %d-in-%d scheme with Gaussian-rational coefficients EXISTS in the V(x)V stand-in (checked exactly)' % (W, s) if ok else 'NOT exact'), flush=True)
