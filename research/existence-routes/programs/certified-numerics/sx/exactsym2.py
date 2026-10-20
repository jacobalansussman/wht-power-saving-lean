# exactsym2.py W s file.npz [maxden] -- exact check of a numerical Sym^2-coordinate stand-in solution (sym2core.Sys2):
# recognise each coordinate as a Gaussian rational, evaluate exactly with Fractions (R', K' have dyadic entries, read
# exactly from their binary floats and checked to be dyadic with denominator <= 4).
import sys, numpy as np
from fractions import Fraction as Fr
exec(open('sym2core.py').read())
W, s = int(sys.argv[1]), int(sys.argv[2]); z = np.load(sys.argv[3]); p = z['p']
maxden = int(sys.argv[4]) if len(sys.argv) > 4 else 1 << 12
S = Sys2(W, s)
def rec(x):
    a = Fr(float(x)).limit_denominator(maxden)
    return a if abs(float(a) - x) < 1e-9 else None
vals = []; bad = 0
for x in p:
    r, i = rec(x.real), rec(x.imag)
    if r is None or i is None: bad += 1; vals.append(None)
    else: vals.append((r, i))
print('coordinates %d, nonzero %d, not recognised %d' % (len(p), sum(1 for v in vals if v and (v[0] or v[1])), bad), flush=True)
if bad:
    print('unrecognised values (first 12):', [p[k] for k in range(len(p)) if vals[k] is None][:12]); sys.exit()
print('denominators', sorted(set(v[0].denominator for v in vals) | set(v[1].denominator for v in vals)))
def gq(M):
    re = np.empty(M.shape, object); im = np.empty(M.shape, object)
    for ix in np.ndindex(M.shape):
        re[ix] = Fr(float(M[ix].real)); im[ix] = Fr(float(M[ix].imag))
        assert re[ix].denominator <= 4 and im[ix].denominator <= 4
    return re, im
def mul(A, B): return (A[0].dot(B[0]) - A[1].dot(B[1]), A[0].dot(B[1]) + A[1].dot(B[0]))
def add(A, B): return (A[0] + B[0], A[1] + B[1])
def zero(shape):
    z_ = np.empty(shape, object); z_.fill(Fr(0)); return z_, z_.copy()
RQ = gq(Rp); EQ = [gq(e_) for e_ in S.E]; TQ = [gq(t) for t in S.T]
k = 0; coefs = []
for n_ in S.slots:
    row = []
    for i in range(n_):
        Mq = zero((10, 10))
        for a in range(nA):
            vr, vi = vals[k]; k += 1
            Mq[0][BLK[a]] += vr; Mq[1][BLK[a]] += vi
        row.append(Mq)
    coefs.append(row)
srcs = list(EQ)
for j in range(s):
    c = zero((10, 10 * W))
    for i in range(len(srcs)): c = add(c, mul(coefs[j][i], srcs[i]))
    srcs.append(mul(RQ, c))
ok = True
for t in range(W):
    o = zero((10, 10 * W))
    for i in range(len(srcs)): o = add(o, mul(coefs[s + t][i], srcs[i]))
    diff = sum(1 for ix in np.ndindex(o[0].shape) if o[0][ix] != TQ[t][0][ix] or o[1][ix] != TQ[t][1][ix])
    print('output %d: %d of %d entries differ from the exact target' % (t, diff, o[0].size), flush=True)
    ok &= diff == 0
print('EXACT VERDICT: %s' % ('Gaussian-rational %d-in-%d scheme on Sym^2 (block-diagonal mixing) checked exactly' % (W, s) if ok else 'NOT exact'), flush=True)
