"""frames_lim.py -- mw2-absorb4: for a limit skeleton, group every generator's terms by frame and report, per frame, the
free letters relating the terms (ratio u = h' h^-1, its order).  Used for the syzygy argument (NOTE s2.4).
Usage: frames_lim.py lim.pkl"""
import sys, pickle
from twist import *
d = pickle.load(open(sys.argv[1], 'rb'))
sk = build_T(d['m'], d['W'], d['H'], d['recipe'], d['skseed'], d.get('twist', {}))
G = sk.G
for i, T in enumerate(sk.gens):
    byf = {}
    for j, (l, g) in enumerate(T):
        byf.setdefault((G.frm(g), l), []).append((j, g))
    desc = []
    for (f, l), lst in byf.items():
        if len(lst) == 1:
            desc.append(f"a{l} d{G.dist(lst[0][1])}: 1 term")
        else:
            g0 = lst[0][1]
            us = []
            for (j, g) in lst[1:]:
                u = G.mul(G.MATS[g], np.conj(G.MATS[g0]).T)
                us.append((G.is_monomial(u), order_of(G, u)))
            desc.append(f"a{l} d{G.dist(g0)}: {len(lst)} terms, ratios (free, order) {us}")
    print(f"g{i}: " + " | ".join(desc))
