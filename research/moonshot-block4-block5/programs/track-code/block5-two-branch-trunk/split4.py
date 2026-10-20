"""split4.py (b5-twobranch): DESIGNED off-S_Z companions on the close branches, lam = (4,4) (PROOF section 2, Lemma Z prime).
Config (split.py): branch a = the close middle frame for array 1 (d(P1,Na)=4, d(Q,Na)=3), branch b = the close one for
array 2 (mirror), lam = (4,4).  Lemma L4' puts the whole K-part there; Lemma Z' says the leftover needs companions w
outside S_Z whose producing columns have their own-array image Y = tau^-1(A) at layer <= LIMC.  For each array: ncomp
t1-injections of the other input with tau(M0) = other anchor, tau(N_own) = N_own and d(anchor, tau^-1(A)) <= LIMC,
tau^-1(A) != A.  Then the layered exact solve (optionally with designer rounds).
usage: split4.py seed ncomp LIMC rounds maxinj dlist LIMY [mu]"""
import sys, time, random
import numpy as np
import split as sp
from split import make_config_split, class_seeds, closure_multi, diag
from endsolve import Eng, Symp, lift_P, find_P, assemble, lm, vp_residual, pi_bar, merge
from sp2 import (M, L0, Fstar, D as DS, mul as smul, inv as sinv, image, dist, order as sorder)
import designer as dz
CAP = 4000


def companions(rng, cf, L, ncomp, LIMC, maxs=3000):
    out = []
    for own in (0, 1):
        N = cf['N2'] if own == 0 else cf['N2p']
        A = cf['A1'] if own == 0 else cf['A2']
        anc = cf['P1'] if own == 0 else cf['Q']
        oth = 1 - own
        pib = pi_bar(L, oth); pinv = sinv(pib)
        X = image(pinv, N); Yt = image(DS, N)
        got = 0; stats = {}
        for s in range(maxs):
            b = find_P(rng, X, Yt)
            if b is None:
                break
            tau = smul(pib, smul(sinv(b), DS))
            assert image(tau, N) == N
            Y = image(sinv(tau), A)
            key = (dist(anc, Y), Y == A); stats[key] = stats.get(key, 0) + 1
            if Y == A or dist(anc, Y) > LIMC:
                continue
            out.append(dict(site='t1', l=oth, where='comp%d/img%d' % (own, dist(anc, Y)), Y=None, letter=b,
                            tau_order=sorder(tau), lam=0))
            got += 1
            if got >= ncomp:
                break
        print("   array %d: %d companion injection(s) of x%d; image-layer stats over %d samples %s" % (
            own + 1, got, oth + 1, s + 1, dict(sorted(stats.items()))), flush=True)
    return out


def main():
    a = sys.argv[1:]
    seed, ncomp, LIMC, rounds, maxinj = [int(x) for x in a[:5]]
    dlist = [int(c) for c in a[5]] if len(a) > 5 else [5, 4, 3, 2, 1]
    sp.LIMY = int(a[6]) if len(a) > 6 else 4
    mu = float(a[7]) if len(a) > 7 else 0.7
    global CAP
    sp.PEEL = len(a) > 8 and a[8] == 'p'
    CAP = int(a[9]) if len(a) > 9 else 4000
    rng = random.Random(seed); t0 = time.time()
    cf, L = make_config_split(rng, (4, 3), (3, 4), 4, 4)
    if cf is None:
        print("seed %d: no config" % seed); return
    print("seed %d: close/close config %s (%.0fs)" % (seed, cf['info'], time.time() - t0), flush=True)
    injs = companions(rng, cf, L, ncomp, LIMC)
    E = Eng(M, seed=seed); Sy = Symp(E)
    TRIGS = [[], []]
    for d in dlist:
        for rd in range(rounds + 1):
            P = dz.build_prog2(E, Sy, L, injs, mu)
            P.mats = dict(P.mats); P.mats['c3b'] = lift_P(E, L['c3b'])
            seeds = [class_seeds(P, L, cf, l, random.Random(7 + l)) for l in (0, 1)]
            cls = []
            for l in (0, 1):
                cl = sp.closure_any(P, l, d, seeds[l], cap=CAP)
                TRIGS[l] = list(sp.LAST_TRIG)
                if cl is None:
                    print("  d=%d round %d: closure cap" % (d, rd), flush=True); cls = None; break
                cls.append(cl)
            if cls is None:
                break
            singles = []
            for li, cl in enumerate(cls):
                singles += [(li, r) for r in diag(P, cl, li)]
                if sp.PEEL:
                    singles += [(li, r) for r in TRIGS[li]]
            Nss = [assemble(cl, P.nal, rng) for cl in cls]
            cost, al = lm(Nss, P.nal, rng, starts=4, iters=80)
            print("  d=%d round %d: injections %d, closures %s, residual %.3e (%.0fs)" % (
                d, rd, len(injs), [(len(c[0]), len(c[1]), len(c[3][0])) for c in cls], cost, time.time() - t0), flush=True)
            if cost < 1e-20 or rd == rounds or len(injs) >= maxinj:
                break
            added = 0; seen = set()
            for (li, row) in singles:
                if row in seen:
                    continue
                seen.add(row)
                j, w = dz.design_injection(rng, cf, L, Sy, row, 'y' if rng.random() < 0.5 else 'z')
                if j is not None:
                    injs.append(j); added += 1
                if len(injs) >= maxinj:
                    break
            print("     added %d designer injections" % added, flush=True)
            if added == 0:
                break
        if cls is None or cost > 1e-12:
            if cls is not None:
                r, _ = vp_residual(Nss, al if al is not None else np.zeros(0), False)
                n0 = Nss[0][0].shape[0]
                print("     per output: %.3e %.3e" % (float(np.vdot(r[:n0], r[:n0]).real), float(np.vdot(r[n0:], r[n0:]).real)), flush=True)
            print("STOP at d=%d" % d, flush=True)
            break
    else:
        print("ALL LAYERS PASS (per output): verify kappa", flush=True)
    print("done (%.0fs)" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
