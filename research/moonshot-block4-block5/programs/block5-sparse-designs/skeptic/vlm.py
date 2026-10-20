"""vlm.py (skeptic v-w8-b5sparse, own code): NUMERICAL re-test of the ball systems on fresh random starts.
Complex unknowns: alpha_k (shared), one w per column per block.  Residuals: every ball row (and, D4=1, every d=4 row),
plus w_seed - 1 per block (the seed's K-coefficient is 1 in these pickles; asserted).  Levenberg-Marquardt / TRF via
scipy.optimize.least_squares on the real form with an exact sparse Jacobian.  Reports the residual norm and min |alpha|.
usage: vlm.py d4.pkl d3.pkl K starts seed [D4=1] [AREG=0: weight of a penalty pulling |alpha_k| away from 0: none]"""
import sys, os, pickle, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vball import ballrows


def blocks(d4f, d3f, K, D4):
    P4 = pickle.load(open(d4f, 'rb')); P3 = pickle.load(open(d3f, 'rb'))
    out = []
    for l in (0, 1):
        if D4:
            cl, sc, cap = P4['sys'][l]; out.append((cl, sc, None))
        if K > 0:
            cl, sc, cap = P3['sys'][l]; out.append((cl, sc, ballrows(cl[2], sc, K)))
    return out, P4['nal']


def build(bl):
    # unknown layout: alphas 1..nal -> index 0..nal-1 ; then w's
    eqs = []          # list of lists of (coef, k, wi)  (k = 0: constant 1)
    seeds = []
    wid = {}
    for bi, (cl, sc, rs) in enumerate(bl):
        rows, cols, ent, (klass, kterms) = cl
        byrow = {}
        for (r, c, cf) in ent:
            if rs is not None and r not in rs:
                continue
            key = (bi, c)
            if key not in wid:
                wid[key] = len(wid)
            for k, v in cf.items():
                byrow.setdefault(r, []).append((float(v), k, wid[key]))
        eqs += list(byrow.values())
        for c in sc:
            assert kterms[c] and all(cf == {0: 1.0} for (rk, cf) in kterms[c]) and len(kterms[c]) == 1
            seeds.append(wid[(bi, c)])
    return eqs, seeds, len(wid)


def main():
    d4f, d3f, K, starts, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    D4 = int(os.environ.get('D4', '1'))
    bl, nal = blocks(d4f, d3f, K, D4)
    eqs, seeds, nw = build(bl)
    A = sorted(set(k for e in eqs for (c, k, w) in e if k))
    amap = {k: i for i, k in enumerate(A)}
    na = len(A); INV = int(os.environ.get('INV', '0')); n = na + nw + (na if INV else 0)
    # sparse structure
    R = []; Cw = []; Ca = []; V = []
    for ri, e in enumerate(eqs):
        for (c, k, w) in e:
            R.append(ri); Cw.append(na + w); Ca.append(amap[k] if k else -1); V.append(c)
    R = np.array(R); Cw = np.array(Cw); Ca = np.array(Ca); V = np.array(V)
    m = len(eqs); ms = len(seeds)
    print('vlm: %d equations, %d seeds, %d alphas, %d w unknowns (complex)' % (m, ms, na, nw), flush=True)

    def unpack(x):
        return x[:n] + 1j * x[n:]

    def F(x):
        z = unpack(x)
        av = np.where(Ca >= 0, z[np.maximum(Ca, 0)], 1.0)
        t = V * av * z[Cw]
        f = np.zeros(m, complex); np.add.at(f, R, t)
        fs = z[na + np.array(seeds)] - 1.0
        f = np.concatenate([f, fs])
        if INV:
            f = np.concatenate([f, z[:na] * z[na + nw:] - 1.0])
        return np.concatenate([f.real, f.imag])

    def Jac(x):
        z = unpack(x)
        av = np.where(Ca >= 0, z[np.maximum(Ca, 0)], 1.0)
        rows = [R, R[Ca >= 0]]
        colsl = [Cw, Ca[Ca >= 0]]
        vals = [V * av, (V * z[Cw])[Ca >= 0]]
        rr = np.concatenate(rows + [m + np.arange(ms)]); cc = np.concatenate(colsl + [na + np.array(seeds)])
        vv = np.concatenate(vals + [np.ones(ms)])
        if INV:
            ii = np.arange(na)
            rr = np.concatenate([rr, m + ms + ii, m + ms + ii]); cc = np.concatenate([cc, ii, na + nw + ii])
            vv = np.concatenate([vv, z[na + nw:], z[:na]])
        Jc = sp.csr_matrix((vv, (rr, cc)), shape=(m + ms + (na if INV else 0), n))
        Jr, Ji = Jc.real, Jc.imag
        return sp.vstack([sp.hstack([Jr, -Ji]), sp.hstack([Ji, Jr])]).tocsr()

    rng = np.random.default_rng(seed)
    for s in range(starts):
        t = time.time()
        z0 = (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
        x0 = np.concatenate([z0.real, z0.imag])
        res = least_squares(F, x0, jac=Jac, method='trf', x_scale='jac', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000,
                            tr_solver='lsmr')
        z = unpack(res.x)
        r = np.sqrt(np.sum(res.fun ** 2))
        al = np.abs(z[:na])
        wv = np.abs(z[na:na + nw])
        # polish check: fix alpha, solve the LINEAR system in w by least squares; residual and |w|
        av = np.where(Ca >= 0, z[np.maximum(Ca, 0)], 1.0)
        Aw = sp.csr_matrix((V * av, (R, Cw - na)), shape=(m, nw))
        Sd = sp.csr_matrix((np.ones(ms), (np.arange(ms), np.array(seeds))), shape=(ms, nw))
        AA = sp.vstack([Aw, Sd]).toarray(); bb = np.concatenate([np.zeros(m), np.ones(ms)])
        wl, *_ = np.linalg.lstsq(AA, bb, rcond=None)
        rl = np.linalg.norm(AA @ wl - bb)
        sv = np.linalg.svd(AA, compute_uv=False)
        print('start %d: residual %.3e, min|alpha| %.3e, max|alpha| %.3e, max|w| %.3e, nfev %d; at this alpha: linear LSQ residual %.3e, min-norm |w| %.3e, rank(tol 1e-9) %d of %d cols (%.0fs)' % (
            s, r, al.min(), al.max(), wv.max(), res.nfev, rl, np.linalg.norm(wl), int((sv > 1e-9 * sv[0]).sum()), nw, time.time() - t), flush=True)


if __name__ == '__main__':
    main()
