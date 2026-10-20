"""verify.py -- w5-twoshared: numerical solution, formal-identity check and raw-array check of a skeleton that Singular
called CONSISTENT.  (1) LM over every program coefficient (all free) of the summed closure residual of both outputs,
many starts; (2) at the best point, for each output the readout (closure columns with least-squares coefficients) is
expanded as a FORMAL sum of terms (array, matrix); terms with free (monomial) matrices are cancelled by input readouts;
the remainder must be exactly 1 * [K] x_l0 (this is the all-f identity); (3) raw arrays: the program, with that
readout and the input corrections, is applied to fresh random inputs at f = 1..fmax (letters act as M^(x)f).
Usage: verify.py pkl fmax starts"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import least_squares
from exact import G, build, Closure
pk, fmax, starts = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
d = pickle.load(open(pk, 'rb'))
P = build(2, d['moves'])
Cs = [Closure(P, l, 20000) for l in range(2)]
k = P.np
def fun(x):
    par = x[:k] + 1j * x[k:]
    v = np.concatenate([C.resid_vec(par) for C in Cs])
    return np.concatenate([v.real, v.imag, [0.0 * x[0]]])
rng = np.random.default_rng(1)
best = (np.inf, None)
for s in range(starts):
    r = least_squares(fun, rng.standard_normal(2 * k), method='trf', max_nfev=3000)
    c = 2 * r.cost
    if c < best[0]: best = (c, r.x)
    if c < 1e-24: break
print("best closure residual", best[0], flush=True)
par = best[1][:k] + 1j * best[1][k:]
coefs = P.coefs(par)
print("program coefficients |a|:", np.round(np.abs(par), 4), flush=True)
readouts = []
for l0, C in enumerate(Cs):
    A = C.matrix(par); t = np.zeros(C.nr, dtype=complex); t[0] = 1
    beta = np.linalg.lstsq(A, t, rcond=None)[0]
    S = {}
    for (kind, i, bid), c in C.cols.items():
        if abs(beta[c]) < 1e-14: continue
        b = G.MATS[bid]
        for j, (l, g) in enumerate(P.terms[i]):
            key = (l, G.intern(G.mul(b, G.MATS[g])))
            S[key] = S.get(key, 0) + beta[c] * coefs[i][j]
    free = {kk: v for kk, v in S.items() if G.is_free(kk[1]) and abs(v) > 1e-12}
    rest = {kk: v for kk, v in S.items() if not G.is_free(kk[1]) and abs(v) > 1e-9}
    ok = (len(rest) == 1 and (l0, G.KID) in rest and abs(rest[(l0, G.KID)] - 1) < 1e-8)
    print(f"output {l0}: formal sum: {len(rest)} non-free terms left {[(kk[0], G.dist(kk[1]), np.round(v, 6)) for kk, v in list(rest.items())[:6]]}; identity {'HOLDS' if ok else 'FAILS'}", flush=True)
    readouts.append((C, beta, free))
def act(M, X, f):
    for ax in range(f):
        X = np.moveaxis(np.tensordot(M, X, axes=([1], [ax])), 0, ax)
    return X
for f in range(1, fmax + 1):
    errs = []
    xs = [rng.standard_normal((16,) * f) + 1j * rng.standard_normal((16,) * f) for _ in range(2)]
    vals = []
    for i, T in enumerate(P.terms):
        out = [np.zeros((16,) * f, dtype=complex) for _ in range(2)]
        for j, (l, g) in enumerate(T):
            out[l] += coefs[i][j] * act(G.MATS[g], xs[l], f)
        vals.append(out)
    for l0, (C, beta, free) in enumerate(readouts):
        outv = [np.zeros((16,) * f, dtype=complex) for _ in range(2)]
        for (kind, i, bid), c in C.cols.items():
            for l in range(2):
                outv[l] += beta[c] * act(G.MATS[bid], vals[i][l], f)
        for (l, gid), v in free.items():
            outv[l] -= v * act(G.MATS[gid], xs[l], f)
        tgt = [act(G.K, xs[l0], f) if l == l0 else np.zeros((16,) * f) for l in range(2)]
        errs.append(sum(np.linalg.norm(outv[l] - tgt[l]) for l in range(2)) / np.linalg.norm(tgt[l0]))
    print(f"f = {f}: relative raw-array errors (fresh inputs) {['%.2e' % e for e in errs]}", flush=True)
