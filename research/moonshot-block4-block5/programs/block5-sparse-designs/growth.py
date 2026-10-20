"""growth.py (w8-b5sparse): what fills a cap-limited d = 3 closure?  Histograms of the live columns (type, generator,
class or not), of the rows of live columns (array, layer, through a main or a stream term, number of non-dead
producers), and of the left translations c = D b b_seed^-1 D^-1 of the class columns (their symplectic images: how
many distinct, and do they fix the frames N_a, N_b, A_1?).  usage: as jd2.py (dlist 543)."""
import sys
from collections import Counter
import numpy as np
import jd2, peel2
import jdesign as J
from sp2 import image, mul as smul, inv as sinv, D as DS, L0

orig = J.S.closure_any


def cl(P, l, d, seeds, cap=3000):
    r = orig(P, l, d, seeds, cap)
    if d != 3:
        return r
    E, Sy = P.E, P.S
    colrows, rowprod, dead, expanded, pending = peel2.DEBUG['state']
    live = [k for k in colrows if k not in dead]
    klass = set(jd2.STASH[l]['cl'][3][1].keys())
    cols = jd2.STASH[l]['cl'][1]
    iscl = {k: (cols.get(k) in klass) for k in live}
    hc = Counter((k[0], k[1], 'class' if iscl[k] else 'non') for k in live)
    print("   GROWTH out %d: live columns by (type, generator, class): %s" % (l + 1, dict(hc)), flush=True)
    # rows
    hr = Counter()
    for k in live[:4000]:
        typ, u, b = k
        pre = E.D @ E.MATS[b] if typ == 'y' else E.MATS[b]
        for (l2, t, cf) in P.gens[u]:
            g = E.intern(E.snap(pre @ E.MATS[t]))
            rr = (l2, g)
            if rr not in rowprod:
                continue
            s = rowprod[rr]
            hr[(l2 + 1, E.layer(g), 'main' if cf == {0: 1.0} else 'stream', min(len(s), 4))] += 1
    print("   GROWTH out %d: rows of the first 4000 live columns by (array, layer, term, #live producers capped at 4): %s" % (
        l + 1, dict(sorted(hr.items()))), flush=True)
    # translations of class columns relative to the seed
    seedk = ('y', seeds[0][0], seeds[0][1])
    bs = E.MATS[seedk[2]]
    Dm = E.D; Di = np.conj(Dm).T
    tr = Counter(); fixes = Counter()
    for k in live:
        if not iscl[k] or k[0] != 'y' or k[1] != seedk[1]:
            continue
        c = E.intern(E.snap(Dm @ E.MATS[k[2]] @ np.conj(bs).T @ Di))
        cs = tuple(Sy.of(c))
        tr[cs] += 1
    print("   GROWTH out %d: class y-columns on %s: %d, distinct symplectic translations %d" % (
        l + 1, seedk[1], sum(tr.values()), len(tr)), flush=True)
    return r


J.S.closure_any = cl
if __name__ == "__main__":
    J.main()
