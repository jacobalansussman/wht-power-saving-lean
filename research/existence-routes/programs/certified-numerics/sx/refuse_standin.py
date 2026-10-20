# refuse_standin.py W s -- calibration (b) inside the stand-in: the Sym^2 system of a shape below ceil(3W/2) (impossible
# by Lemma 8) as an exact PolySys, then certnum's exact Nullstellensatz certificates of low degree.
import sys, itertools, numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, '..')
exec(open('sym2core.py').read())
from certnum import PolySys, nulla, GZERO, gadd, gmul
W, s = int(sys.argv[1]), int(sys.argv[2])
S = Sys2(W, s); n = S.P
# polynomials: each source entry as a dict {exponent tuple: (Fr, Fr)}; sources are polynomial matrices
def const_mat(M):
    out = np.empty(M.shape, object)
    for ix in np.ndindex(M.shape):
        v = M[ix]; out[ix] = {} if v == 0 else {(): (Fr(float(v.real)), Fr(float(v.imag)))}
    return out
def padd(p, q):
    r = dict(p)
    for e, c in q.items():
        r[e] = gadd(r.get(e, GZERO), c)
        if r[e] == GZERO: del r[e]
    return r
def pscale(p, c):   # c Gaussian rational
    return {e: gmul(v, c) for e, v in p.items() if gmul(v, c) != GZERO}
def pmulvar(p, k):  # multiply by variable k (exponents as sorted tuples of variable ids)
    return {tuple(sorted(e + (k,))): v for e, v in p.items()}
Rq = [[(Fr(float(Rp[a, b].real)), Fr(float(Rp[a, b].imag))) for b in range(10)] for a in range(10)]
def mixapply(var0, src):   # (sum_a x_{var0 + a} Pb[a]) @ src, Pb[a] = matrix unit BLK[a]
    out = np.empty(src.shape, object)
    for ix in np.ndindex(src.shape): out[ix] = {}
    for a, (i, j) in enumerate(BLK):
        for c in range(src.shape[1]):
            if src[j, c]: out[i, c] = padd(out[i, c], pmulvar(src[j, c], var0 + a))
    return out
def Rapply(src):
    out = np.empty(src.shape, object)
    for ix in np.ndindex(src.shape): out[ix] = {}
    for a in range(10):
        for b in range(10):
            if Rq[a][b] != GZERO:
                for c in range(src.shape[1]):
                    if src[b, c]: out[a, c] = padd(out[a, c], pscale(src[b, c], Rq[a][b]))
    return out
srcs = [const_mat(E_) for E_ in S.E]; k = 0
for j in range(s):
    c = np.empty((10, 10 * W), object)
    for ix in np.ndindex(c.shape): c[ix] = {}
    for i in range(len(srcs)):
        m_ = mixapply(k, srcs[i]); k += nA
        for ix in np.ndindex(c.shape): c[ix] = padd(c[ix], m_[ix])
    srcs.append(Rapply(c))
polys = []
for t in range(W):
    o = np.empty((10, 10 * W), object)
    for ix in np.ndindex(o.shape): o[ix] = {}
    for i in range(len(srcs)):
        m_ = mixapply(k, srcs[i]); k += nA
        for ix in np.ndindex(o.shape): o[ix] = padd(o[ix], m_[ix])
    Tc = const_mat(S.T[t])
    for ix in np.ndindex(o.shape):
        neg = {e: (-v[0], -v[1]) for e, v in Tc[ix].items()}
        polys.append(padd(o[ix], neg))
assert k == n
def to_exp(p):
    q = {}
    for e, v in p.items():
        x = [0] * n
        for vv in e: x[vv] += 1
        q[tuple(x)] = v
    return q
F = PolySys(n, [to_exp(p) for p in polys if p])
nmove = sum(S.slots[:s]) * nA
print('Sym2 stand-in %d in %d: %d unknowns (%d in moves), %d equations' % (W, s, n, nmove, F.N), flush=True)
for d, sub in [(0, None), (1, list(range(nmove))), (1, None), (2, list(range(nmove)))]:
    found, msg = nulla(F, d, var_subset=sub)
    print('NulLA degree %d (%s multipliers): %s  [%s]' % (d, 'move-variable' if sub else 'all-variable',
          'CERTIFICATE (no solution over C)' if found else ('no certificate' if found is False else 'not run'), msg), flush=True)
    if found: break
