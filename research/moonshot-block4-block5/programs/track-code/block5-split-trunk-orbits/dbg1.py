import sys, random, numpy as np
import split as S, endsolve as es, designer as dz, jdesign as J
from sp2 import *
S.LIMY = 3; S.PEEL = True
rng = random.Random(31)
cf, L = S.make_config_split(rng, (4,2), (2,4), 3, 3)
E = es.Eng(M, seed=31); Sy = es.Symp(E); J.EXACT[0] = E
P = dz.build_prog2(E, Sy, L, [], 0.7); P.mats = dict(P.mats); P.mats['c3b'] = es.lift_P(E, L['c3b'])
seeds = [S.class_seeds(P, L, cf, l, random.Random(7 + l)) for l in (0, 1)]
cl = S.closure_any(P, 0, 3, seeds[0])
trig = list(S.LAST_TRIG)
print("triggers", trig)
injs = []
for row in trig:
    l, gid = row
    th = dz.theta(Sy, L, l, gid); G = image(th, L0); Pr = image(th, cf['M0'])
    anc = [cf['P1'], cf['Q']][l]
    print("row array", l+1, "layer", E.layer(gid), "G==A1", G == cf['A1'], "d(anc,G)", dist(anc, G), "d(anc,Pr)", dist(anc, Pr), "Pr==Na", Pr == cf['N2'], "Pr==Nb", Pr == cf['N2p'])
    js = J.t2_options(rng, cf, L, Sy, (0, row) if False else row, 2)
    print("  options", [(j['where'], j['l']) for j in js])
    injs += js
P = dz.build_prog2(E, Sy, L, injs, 0.7); P.mats = dict(P.mats); P.mats['c3b'] = es.lift_P(E, L['c3b'])
# terms at Pr for the trigger row
for row in trig:
    l, gid = row
    Dinv = np.conj(E.D).T
    h = E.intern(E.snap(Dinv @ E.MATS[gid])); Fh = Sy.frame(h)
    for u, T in P.gens.items():
        for (ll, g, c) in T:
            if ll == l and Sy.frame(g) == Fh:
                bm = E.snap(E.MATS[h] @ np.conj(E.MATS[g]).T)
                nz = np.abs(bm) > 1e-9
                vals = sorted(set(np.round(bm[nz], 6)))
                print("  term at Pred: gen", u, "coef", c, "layer", E.layer(g), "is_mono", E.is_mono(bm), "entries", vals[:4])
cl = S.closure_any(P, 0, 3, seeds[0])
