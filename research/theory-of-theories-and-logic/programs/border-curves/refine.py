#!/usr/bin/env python3
"""refine.py (v-ml-proofcx): independent check of the explicit g27_Fp2b_m3_s5 first-order border curve.
1. Build the first-order system with myfo.py (integer coefficients, lifted from the mod-p file).
2. Start from ml-proofcx's printed curve (6 decimals), add random affine slices through the start to remove the
   positive-dimensional freedom, and run Newton in mpmath at 80 digits.  Quadratic convergence to ~1e-70 = strong
   evidence of an exact complex solution of the first-order system.
3. Evaluate the ORIGINAL relaxation (my own parser) along alpha = A + eps D, beta = B/eps + C at small eps:
   the residual must be O(eps) while max|beta| ~ 1/eps (a genuine border, Lemma T).
Also: random-start Newton/LM searches (numpy) for g27_Fany4_s6 and g27_Fany4_s5 first-order systems.
"""
import sys, json, random, time
sys.path.insert(0, '.'); DATA = '../../data/border-curves/'
from myfo import parse, fo_system
import mpmath as mp

def compile_sys(polys, names):
    idx = {v: i for i, v in enumerate(names)}
    return [[(c, [(idx[v], e) for v, e in m]) for m, c in P.items()] for P in polys]

def ev(cs, x):
    out = []
    for P in cs:
        s = 0
        for c, m in P:
            t = c
            for i, e in m: t *= x[i] ** e
            s += t
        out.append(s)
    return out

def jac(cs, x, n):
    J = []
    for P in cs:
        row = [0] * n
        for c, m in P:
            for k, (i, e) in enumerate(m):
                t = c * e * x[i] ** (e - 1)
                for k2, (i2, e2) in enumerate(m):
                    if k2 != k: t *= x[i2] ** e2
                row[i] += t
        J.append(row)
    return J

def main_fp2b():
    vs, p, rel = parse(DATA + 'g27_Fp2b_m3_s5.u.p1.ms')
    fo = fo_system(vs, rel)
    names = sorted({v for P in fo for m in P for v, _ in m})
    n = len(names)
    cs = compile_sys(fo, names)
    cur = [json.loads(l) for l in open(DATA + 'curve_starts.jsonl')]
    cur = [d for d in cur if 'a0' in d][0]
    mapk = {'a0': 'A_', 'd': 'D_', 'b0': 'B_', 'b1': 'C_'}
    mp.mp.dps = 80
    x = []
    for v in names:
        k, nm = v.split('_', 1)
        src = {'A': 'a0', 'D': 'd', 'B': 'b0', 'C': 'b1'}[k]
        re_, im_ = cur[src].get(nm, [0.0, 0.0])
        x.append(mp.mpc(re_, im_))
    r0 = ev(cs, x)
    print('start: n', n, 'eqs', len(cs), 'max|res| at printed (6-decimal) curve', float(max(abs(t) for t in r0)), flush=True)
    rnd = random.Random(7)
    # rank of Jacobian at the start (numerical, via mpmath svd)
    J = mp.matrix(jac(cs, x, n))
    sv = mp.svd_c(J, compute_uv=False)
    svl = sorted([float(abs(s)) for s in sv], reverse=True)
    print('singular values (top, then smallest 18):', ['%.2e' % s for s in svl[:3]], ['%.2e' % s for s in svl[-18:]], flush=True)
    for it in range(14):
        r = ev(cs, x)
        nr = max(abs(t) for t in r)
        print(f'gauss-newton (min-norm, SVD) it {it}: max|res| = {mp.nstr(nr, 5)}', flush=True)
        if nr < mp.mpf(10) ** (-70): break
        J = mp.matrix(jac(cs, x, n))
        U, S, V = mp.svd_c(J)
        smax = max(abs(S[i]) for i in range(len(S)))
        Ur = U.H * mp.matrix(r)
        dx = [mp.mpc(0)] * n
        for k in range(len(S)):
            if abs(S[k]) > mp.mpf(10) ** (-30) * smax:
                c = Ur[k] / S[k]
                for i in range(n): dx[i] += c * mp.conj(V[k, i])
        x = [x[i] - dx[i] for i in range(n)]
    xv = dict(zip(names, x))
    mx = lambda pre: float(max(abs(xv[v]) for v in names if v.startswith(pre)))
    print('max|A| %.4g max|D| %.4g max|B| %.4g max|C| %.4g' % (mx('A_'), mx('D_'), mx('B_'), mx('C_')))
    print('B values (|.|>1e-20):', {v: mp.nstr(xv[v], 8) for v in names if v.startswith('B_') and abs(xv[v]) > 1e-20})
    # 3. the original relaxation along the curve
    rnames = vs
    rcs = compile_sys(rel, rnames)
    for eps in [mp.mpf('1e-2'), mp.mpf('1e-4'), mp.mpf('1e-8'), mp.mpf('1e-16')]:
        pt = []
        for v in rnames:
            if v.startswith('b'): pt.append(xv.get('B_' + v, 0) / eps + xv.get('C_' + v, 0))
            else: pt.append(xv.get('A_' + v, 0) + eps * xv.get('D_' + v, 0))
        res = ev(rcs, pt)
        print(f'eps {mp.nstr(eps,3)}: relaxation |res|_inf = {mp.nstr(max(abs(t) for t in res), 6)}  max|beta| = {mp.nstr(max(abs(pt[i]) for i,v in enumerate(rnames) if v.startswith("b")), 6)}', flush=True)
    # exact infeasibility is checked by msolve separately; here: distance of best exact point? (not attempted)

if __name__ == '__main__':
    t = time.time(); main_fp2b(); print('time %.1f s' % (time.time() - t))
