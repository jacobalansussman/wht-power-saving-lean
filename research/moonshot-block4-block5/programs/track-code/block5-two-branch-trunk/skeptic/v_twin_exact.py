"""v_twin_exact.py (skeptic v-b5-twobranch): exact-element test of the SAME-FRAME TWIN door on the close branch.
Config: split.py close/close, lam = (4,4) (branch a = close frame of array 1, branch b = close frame of array 2).
Array 1 gets: (i) ncomp S_Z companions = t1-injections of x2 with tau(N_a) = N_a and tau^-1(A1) = A1 (w in S_Z);
(ii) optionally a twin trunk term c3a' = Dbar^-1 u Dbar c3a (u in Stab(L0, M0), i.e. h' = u h, same frame N_a) with
u Z != Z (Z = Theta_a^-1 A1) and both images Theta_a(u^{+-1} Z) at layer <= LIMT for BOTH anchors.
Prediction (skeptic's algebra): without twin, output 1 fails at d = 4 (Lemma Z / Z'); with twin it can pass d = 4.
Output 1 is solved ALONE (lm on output 1's system).  Uses the claimant's exact library (copied to lib/).
usage: v_twin_exact.py seed twin(0/1) ncomp muj LIMT dlist"""
import sys, time, random, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import numpy as np
import split as sp
from split import make_config_split, class_seeds
from endsolve import Eng, Symp, lift_P, find_P, assemble, lm, vp_residual, pi_bar, merge
from sp2 import M, L0, Fstar, D as DS, mul as smul, inv as sinv, image, dist, order as sorder
import designer as dz


def sz_companions(rng, cf, L, ncomp, maxs=6000):
    out = []; N = cf['N2']; A = cf['A1']; anc = cf['P1']; oth = 1
    pib = pi_bar(L, oth); pinv = sinv(pib)
    X = image(pinv, N); Yt = image(DS, N)
    for s in range(maxs):
        b = find_P(rng, X, Yt)
        tau = smul(pib, smul(sinv(b), DS))
        assert image(tau, N) == N
        if image(sinv(tau), A) != A:
            continue
        out.append(dict(site='t1', l=oth, where='SZcomp', Y=None, letter=b, tau_order=sorder(tau), lam=0))
        if len(out) >= ncomp:
            break
    print("   %d S_Z companion injection(s) of x2 after %d samples" % (len(out), s + 1), flush=True)
    return out


def twin_letter(rng, cf, L, LIMT, maxs=4000):
    M0 = image(DS, L0); Nc = cf['N2']; A = cf['A1']; P1, Q = cf['P1'], cf['Q']
    H = smul(sinv(L['c2']), smul(DS, smul(sinv(L['c3']), DS)))
    assert image(H, L0) == Nc and image(H, M0) == cf['N1'], "Theta_a(L0)=N_a, Theta_a(M0)=N1"
    Hi = sinv(H); Z = image(Hi, A)
    for s in range(maxs):
        u = find_P(rng, M0, M0)
        assert image(u, L0) == L0 and image(u, M0) == M0
        if image(u, Z) == Z:
            continue
        ok = True
        for w in (u, sinv(u)):
            X = image(H, image(w, Z))
            if X == A or max(dist(P1, X), dist(Q, X)) > LIMT:
                ok = False
        if not ok:
            continue
        c3t = smul(sinv(DS), smul(u, smul(DS, L['c3'])))
        assert image(c3t, L0) == L0
        Ht = smul(sinv(L['c2']), smul(DS, smul(sinv(c3t), DS)))
        assert image(Ht, L0) == Nc, "twin in the same frame"
        X = image(H, image(u, Z))
        print("   twin u found after %d samples: order %d, image (d(P1,X), d(Q,X)) = (%d, %d)" % (
            s + 1, sorder(u), dist(P1, X), dist(Q, X)), flush=True)
        return c3t
    return None


def main():
    a = sys.argv[1:]
    seed, twin, ncomp = int(a[0]), int(a[1]), int(a[2]); muj = float(a[3]); LIMT = int(a[4])
    dlist = [int(c) for c in a[5]] if len(a) > 5 else [5, 4, 3]
    mu = 0.7
    sp.LIMY = int(a[6]) if len(a) > 6 else 4
    rng = random.Random(seed); t0 = time.time()
    cf, L = make_config_split(rng, (4, 3), (3, 4), 4, 4)
    if cf is None:
        print("no config"); return
    print("seed %d twin %d ncomp %d muj %g: config %s" % (seed, twin, ncomp, muj, cf['info']), flush=True)
    injs = sz_companions(rng, cf, L, ncomp)
    c3t = twin_letter(rng, cf, L, LIMT) if twin else None
    if twin and c3t is None:
        print("no twin letter"); return
    E = Eng(M, seed=seed); Sy = Symp(E)
    for d in dlist:
        P = dz.build_prog2(E, Sy, L, injs, mu)
        P.mats = dict(P.mats); P.mats['c3b'] = lift_P(E, L['c3b'])
        if c3t is not None:
            Mt = lift_P(E, c3t)
            assert Sy.of(E.intern(E.snap(Mt))) == c3t
            extra = [(l, P.mulD(Mt, g), {k: muj * v for k, v in c.items()}) for (l, g, c) in P.gens['t2']]
            P.gens['t3'] = merge(P.gens['t3'] + extra)
        seeds = class_seeds(P, L, cf, 0, random.Random(7))
        cl = sp.closure_any(P, 0, d, seeds, cap=6000)
        if cl is None:
            print("  d=%d: closure cap" % d, flush=True); break
        Ns = assemble(cl, P.nal, rng)
        cost, al = lm([Ns], P.nal, rng, starts=4, iters=80)
        nk = len(cl[3][0])
        print("  d=%d: output 1 alone: closure rows %d cols %d class cols %d, residual %.3e (%.0fs)" % (
            d, len(cl[0]), len(cl[1]), nk, cost, time.time() - t0), flush=True)
        sp.diag(P, cl, 0)
        if cost > 1e-12:
            print("STOP at d=%d" % d, flush=True); break
    print("done (%.0fs)" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
