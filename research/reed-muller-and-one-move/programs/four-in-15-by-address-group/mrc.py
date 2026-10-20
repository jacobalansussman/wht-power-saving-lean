#!/usr/bin/env python3
"""mrc.py (w8-rmexact): minimum-rank-completion test on the tight forms (brief step 3), Z_n Singer power.
Virtual moves (t,k), channel k sees direction sigma^-k z_t; virtual arrays (l,k), target sigma^-k sigma^-kT.
A form G is tight when |D_G| = W n_odd(G).  On a tight support (square, Theorem 10 + Lemma C):
B_D A_D = I + Gamma_D, so K = B A (rank <= Wn, block circulant: left-regular C[Z_n] coefficients) has
K[(t,k),(s,k')] = delta whenever the two virtual moves lie in a common tight support and s >= t in time (s = t: same
physical move, Gamma is 0 there).  We solve for circulant B (Mn x Wn) and A (Wn x Mn) meeting all those entries
(least squares, random restarts).  Residual 0 = the tight structure is consistent; a floor = evidence of an obstruction.
Usage: mrc.py --n 3 --power 5 --dirs 1,2,4,8,3,1,2,4,8,3,1,2,4,8,3 [--order random] --restarts 5"""
import sys, argparse, itertools, numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from forms import setup, pair
ap_ = argparse.ArgumentParser()
ap_.add_argument('--n', type=int, default=3); ap_.add_argument('--power', type=int, default=5)
ap_.add_argument('--W', type=int, default=4); ap_.add_argument('--dirs', default='1,2,4,8,3,1,2,4,8,3,1,2,4,8,3')
ap_.add_argument('--restarts', type=int, default=4); ap_.add_argument('--seed', type=int, default=1)
ap_.add_argument('--shuffle', type=int, default=0, help='random time order with this seed (0 = as given)')
args = ap_.parse_args()
m = 4; D = setup(m); supp = D['supp']; nf = D['nf']; mm, tr = D['matmul'], D['transpose']; ap = D['apply']
poly = 0b10011
cols = []
for j in range(m):
    v = 1 << (j + 1)
    if v >> m: v ^= poly
    cols.append(v)
A1 = [sum(((cols[j] >> i) & 1) << j for j in range(m)) for i in range(m)]
I = [1 << i for i in range(m)]
def mpow(X, k):
    R = I
    for _ in range(k): R = mm(X, R)
    return R
sig = mpow(A1, args.power); n = args.n
assert mpow(sig, n) == I
sinv = mpow(sig, n - 1)
chan = [mpow(sinv, k) for k in range(n)]            # sigma^-k
dirs = [int(v) for v in args.dirs.split(',')]
if args.shuffle:
    rr = np.random.default_rng(args.shuffle); dirs = list(rr.permutation(dirs))
M = len(dirs); W = args.W
vdir = [[ap(chan[k], dirs[t]) for k in range(n)] for t in range(M)]
targ = [D['enc_sym'](mm(chan[k], tr(chan[k]))) for k in range(n)]
tight = []
for G in range(1, nf):
    Dg = [(t, k) for t in range(M) for k in range(n) if vdir[t][k] in supp[G]]
    nodd = sum(pair(G, targ[k]) for k in range(n))
    if len(Dg) == W * nodd and nodd > 0: tight.append(Dg)
spec = set()
for Dg in tight:
    for (t, k), (s, k2) in itertools.product(Dg, Dg):
        if s >= t: spec.add((t, s, (k - k2) % n))
spec = sorted(spec)
print("n=%d power=%d dirs=%s: %d tight forms; %d specified circulant entries of K (rank <= %d)" % (n, args.power, dirs, len(tight), len(spec), W * n), flush=True)
ti = np.array([a for a, b, c in spec]); si = np.array([b for a, b, c in spec]); di = np.array([c for a, b, c in spec])
want = np.array([1.0 if (a == b and c == 0) else 0.0 for a, b, c in spec])
# K[(t,k),(s,k')] = sum_l sum_j beta[t,l,(k-j)] alpha[l,s,(j-k')]  ->  kappa[t,s,d] = sum_l sum_e beta[t,l,e] alpha[l,s,(d-e)]
def unpack(x):
    z = x[:len(x) // 2] + 1j * x[len(x) // 2:]
    b = z[:M * W * n].reshape(M, W, n); a = z[M * W * n:].reshape(W, M, n)
    return b, a
def kappa(b, a):
    # circular convolution over the channel index
    fb = np.fft.fft(b, axis=2); fa = np.fft.fft(a, axis=2)
    fk = np.einsum('tlj,lsj->tsj', fb, fa)
    return np.fft.ifft(fk, axis=2)
def res(x):
    b, a = unpack(x); k = kappa(b, a)
    r = k[ti, si, di] - want
    return np.concatenate([r.real, r.imag])
nB = M * W * n
def jac(x):
    b, a = unpack(x)
    J = np.zeros((len(spec), 2 * nB), dtype=complex)
    for r, (t, s_, d) in enumerate(spec):
        for l in range(W):
            for e in range(n):
                J[r, (t * W + l) * n + e] += a[l, s_, (d - e) % n]          # d kappa / d beta[t,l,e]
                f = (d - e) % n
                J[r, nB + (l * M + s_) * n + f] += b[t, l, e]                # d kappa / d alpha[l,s,f]
    return np.block([[J.real, -J.imag], [J.imag, J.real]])
rng = np.random.default_rng(args.seed)
best = None
for rs in range(args.restarts):
    x0 = rng.normal(size=2 * 2 * M * W * n) * 0.6
    sol = least_squares(res, x0, jac=jac, method='trf', max_nfev=3000, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    f = float(np.sum(sol.fun ** 2))
    b, a = unpack(sol.x); kk = kappa(b, a)
    # rank of the full 45 x 45 K and of K on each tight support
    Kfull = np.zeros((M * n, M * n), dtype=complex)
    for t in range(M):
        for s in range(M):
            for k in range(n):
                for k2 in range(n): Kfull[t * n + k, s * n + k2] = kk[t, s, (k - k2) % n]
    sv = np.linalg.svd(Kfull, compute_uv=False)
    print("restart %d: residual %.3e  (status %s, nfev %d); singular values of K: %s" % (rs, f, sol.status, sol.nfev, np.round(sv[:14], 4)), flush=True)
    if best is None or f < best: best = f
print("best residual %.3e" % best)
