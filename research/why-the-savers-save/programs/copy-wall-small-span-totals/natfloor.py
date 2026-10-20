"""natfloor.py (d-R7): natural-family (linear totals) star floor of a sub-unit U of P112.
A linear total is T -> <c, w_T> on U (c in span of U's roots); its frame is the label span of its support; cost =
dim of that span.  Floor = min-weight basis of span(U's roots), by greedy over frames F of increasing dimension d:
the totals that fit inside F are the annihilator of the roots of U outside F.  Exact up to float rank (8 x 8 Grams
of small integer vectors; tolerance 1e-8)."""
import sys, time; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
OUTER = np.einsum('ti,tj->tij', LAM, LAM).reshape(NP, 64)

def natfloor(U, dmax=8, verbose=False):
    U = np.array(sorted(U)); rU = qrank(LAM[U])
    Ub = np.zeros(NP, bool); Ub[U] = True
    # basis of span(U roots)
    u_, s_, vt_ = np.linalg.svd(LAM[U]); SP = vt_[:rU]               # rU x 8
    acc = np.zeros((0, 8)); cost = 0; per = []
    for d in range(1, dmax + 1):
        ds = duals(9 - d) if d < 9 else [()]
        D = np.array(ds, dtype=np.int64).reshape(len(ds), 9 - d)
        M = np.ones((len(ds), NP), dtype=bool)
        for i in range(9 - d): M &= (PT[D[:, i]] == 0)
        C = Ub[None, :] & ~M                                       # roots of U outside F
        C = np.unique(C, axis=0)
        G = (C.astype(float) @ OUTER).reshape(-1, 8, 8)
        # restrict to span(U): Gram in SP coordinates
        Gs = np.einsum('ai,nij,bj->nab', SP, G, SP)
        ev = np.linalg.eigvalsh(Gs)
        low = ev[:, 0] < 1e-8 * max(1.0, ev[:, -1].max())
        new = []
        for n in np.nonzero(low)[0]:
            w, v = np.linalg.eigh(Gs[n]); nv = v[:, w < 1e-8 * max(1.0, w[-1])]
            new.append((nv.T @ SP))
        before = qrank(acc) if len(acc) else 0
        if new: acc = np.vstack([acc] + new)
        after = qrank(acc) if len(acc) else 0
        inc = after - before; cost += d * inc; per.append((d, inc, int(low.sum())))
        if verbose: print('   d=%d: frames whose outside roots have rank < %d: %d; rank of available totals %d -> %d' % (d, rU, int(low.sum()), before, after), flush=True)
        if after >= rU: break
    return cost if after >= rU else None, rU, per
if __name__ == '__main__':
    t0 = time.time()
    print('P112 itself:', natfloor(range(NP), verbose=True), '%.0fs' % (time.time() - t0))
