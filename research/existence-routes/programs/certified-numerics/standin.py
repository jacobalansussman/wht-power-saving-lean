# standin.py -- the block-2 representation stand-in (V(x)V) of ex-dimension: W in M with s shared paid moves.
# Step 1 (numerics): does the unit variety project dominantly onto the MOVE coefficients?  For random moves, solve the
#   output mixing (which enters linearly) by least squares; residual ~ 0 means: for generic moves an exact solution
#   exists, so an exact closure over Q(i) with random Gaussian-rational moves is possible.
# Step 2 (exact): random Gaussian-rational moves, exact sources over Q(i) (fmpq_mat real/imag), exact consistency of the
#   linear system for the output mixing (realified, rank [A|b] == rank A).  Consistent => an exact unit exists in the
#   stand-in (every coefficient in Q(i)).
import sys, time, numpy as np
W, NS = int(sys.argv[1]), int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'num'
exec(open('standin_core.py').read())
from flint import fmpq_mat, fmpq
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 7)
E = []
for l in range(W):
    Xl = np.zeros((d, d * W), complex); Xl[:, l * d:(l + 1) * d] = np.eye(d); E.append(Xl)
T = []
for t in range(W):
    X = np.zeros((d, d * W), complex); X[:, t * d:(t + 1) * d] = Kr; T.append(X)
# exact spanning basis of A from group elements (rho(b) has entries 0, +-1, +-i)
flat = np.array([rho(b).flatten() for b in Bgrp]); chosen = []; Q = np.zeros((0, flat.shape[1]), complex)
for k in range(len(Bgrp)):
    v = flat[k] - (Q.conj() @ flat[k]) @ Q if len(Q) else flat[k].copy()
    if np.linalg.norm(v) > 1e-6:
        chosen.append(k); Q = np.vstack([Q, v / np.linalg.norm(v)])
    if len(chosen) == dimA: break
Pa = [rho(Bgrp[k]) for k in chosen]
assert all(np.allclose(P, np.round(P.real) + 1j * np.round(P.imag)) for P in Pa)
print('exact basis of A from %d group elements' % len(Pa), flush=True)

def num_test(trials=3, scale=0.5):
    for tr in range(trials):
        srcs = list(E)
        for j in range(NS):
            c = sum(sum((rng.normal() + 1j * rng.normal()) * scale * Pa[a] for a in range(dimA)) @ srcs[i] for i in range(len(srcs)))
            srcs.append(R @ c)
        res = []
        for t in range(W):
            cols = np.array([(Pa[a] @ srcs[i]).reshape(-1) for i in range(len(srcs)) for a in range(dimA)]).T
            x, *_ = np.linalg.lstsq(cols, T[t].reshape(-1), rcond=None)
            res.append(float(np.linalg.norm(cols @ x - T[t].reshape(-1)) / np.linalg.norm(T[t])))
            sv = np.linalg.svd(cols, compute_uv=False)
        print('random moves trial %d: relative LS residual per output %s; rank of output map %d of %d columns' %
              (tr, ['%.2e' % r for r in res], int((sv > 1e-9 * sv[0]).sum()), cols.shape[1]), flush=True)

def to_q(z):
    from fractions import Fraction as Fr
    return Fr(z.real).limit_denominator(1 << 20), Fr(z.imag).limit_denominator(1 << 20)

def exact_test(den=2, rngint=2):
    # matrices over Q(i) as pairs (Re, Im) of fmpq_mat
    def qm(M):  # numpy complex with exact dyadic entries -> (fmpq_mat, fmpq_mat)
        from fractions import Fraction as Fr
        re = fmpq_mat(M.shape[0], M.shape[1]); im = fmpq_mat(M.shape[0], M.shape[1])
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                a = Fr(M[i, j].real); b = Fr(M[i, j].imag)
                if a: re[i, j] = fmpq(a.numerator, a.denominator)
                if b: im[i, j] = fmpq(b.numerator, b.denominator)
        return re, im
    def mul(A, B): return (A[0] * B[0] - A[1] * B[1], A[0] * B[1] + A[1] * B[0])
    def add(A, B): return (A[0] + B[0], A[1] + B[1])
    def scal(c, A): return (A[0] * c[0] - A[1] * c[1], A[0] * c[1] + A[1] * c[0])
    PaQ = [qm(P) for P in Pa]; RQ = qm(R); EQ = [qm(e) for e in E]; TQ = [qm(t) for t in T]
    t0 = time.time(); srcs = list(EQ); moves = []
    for j in range(NS):
        acc = None
        for i in range(len(srcs)):
            Mji = None
            for a in range(dimA):
                c = (fmpq(int(rng.integers(-rngint, rngint + 1)), den), fmpq(int(rng.integers(-rngint, rngint + 1)), den))
                term = scal(c, PaQ[a]); Mji = term if Mji is None else add(Mji, term)
            pr = mul(Mji, srcs[i]); acc = pr if acc is None else add(acc, pr)
        srcs.append(mul(RQ, acc))
    print('exact sources built in %.1fs' % (time.time() - t0), flush=True)
    ok_all = True
    for t in range(W):
        cols = []
        for i in range(len(srcs)):
            for a in range(dimA):
                Pq = mul(PaQ[a], srcs[i]); cols.append(Pq)
        nr = d * d * W; nc = len(cols)
        A = fmpq_mat(2 * nr, 2 * nc + 1)
        for j, (re, im) in enumerate(cols):
            for r_ in range(d):
                for c_ in range(d * W):
                    k = r_ * d * W + c_; x = re[r_, c_]; y = im[r_, c_]
                    if x != 0:
                        A[k, 2 * j] = x; A[nr + k, 2 * j + 1] = x      # coefficient real part u: (u)(x+iy)
                    if y != 0:
                        A[nr + k, 2 * j] = y; A[k, 2 * j + 1] = -y     # coefficient imag part v: (iv)(x+iy)
        for r_ in range(d):
            for c_ in range(d * W):
                k = r_ * d * W + c_
                A[k, 2 * nc] = TQ[t][0][r_, c_]; A[nr + k, 2 * nc] = TQ[t][1][r_, c_]
        t1 = time.time()
        Ab = A; Aonly = fmpq_mat(2 * nr, 2 * nc)
        for r_ in range(2 * nr):
            for c_ in range(2 * nc): Aonly[r_, c_] = A[r_, c_]
        r1 = Aonly.rank(); r2 = Ab.rank()
        print('output %d: exact linear system %d x %d over Q; rank A = %d, rank [A|b] = %d -> %s  (%.1fs)' %
              (t, 2 * nr, 2 * nc, r1, r2, 'CONSISTENT' if r1 == r2 else 'INCONSISTENT', time.time() - t1), flush=True)
        ok_all &= (r1 == r2)
    print('EXACT VERDICT (W=%d, s=%d, VV stand-in, random Gaussian-rational moves): %s' %
          (W, NS, 'a unit with every coefficient in Q(i) EXISTS' if ok_all else 'not closed by these moves'), flush=True)


def rand_move(nsrc, scale=0.5):
    return [sum((rng.normal() + 1j * rng.normal()) * scale * Pa[a] for a in range(dimA)) for i in range(nsrc)]

def lin_test(trials=2):
    """fix moves 0..s-2 and each output's coefficient on the last source at random; the rest enters LINEARLY"""
    for tr in range(trials):
        srcs = list(E)
        for j in range(NS - 1):
            Mj = rand_move(len(srcs)); srcs.append(R @ sum(Mj[i] @ srcs[i] for i in range(len(srcs))))
        nb = len(srcs)                              # sources feeding the last move and the outputs (except the last)
        Mlast = [sum((rng.normal() + 1j * rng.normal()) * 0.5 * Pa[a] for a in range(dimA)) for t in range(W)]
        blocks = []                                 # columns: last move coefficients (i, a), then output t coefficients
        colsL = []
        for i in range(nb):
            for a in range(dimA):
                colsL.append(np.concatenate([(Mlast[t] @ R @ Pa[a] @ srcs[i]).reshape(-1) for t in range(W)]))
        colsO = []
        for t in range(W):
            for i in range(nb):
                for a in range(dimA):
                    v = np.zeros(W * d * d * W, complex); v[t * d * d * W:(t + 1) * d * d * W] = (Pa[a] @ srcs[i]).reshape(-1)
                    colsO.append(v)
        A_ = np.array(colsL + colsO).T; b_ = np.concatenate([T[t].reshape(-1) for t in range(W)])
        x, *_ = np.linalg.lstsq(A_, b_, rcond=None)
        sv = np.linalg.svd(A_, compute_uv=False)
        print('fix-and-linear trial %d: %d unknowns, %d equations, rank %d, relative LS residual %.2e' %
              (tr, A_.shape[1], A_.shape[0], int((sv > 1e-9 * sv[0]).sum()), np.linalg.norm(A_ @ x - b_) / np.linalg.norm(b_)), flush=True)

if mode in ('num', 'both'): num_test(1)
if mode in ('lin', 'both'): lin_test()
if mode in ('exact', 'both'): exact_test()
