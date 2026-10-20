"""jd_zcover.py (w8-b5sparse): step 2 of w7-b5's to-do list, frame-level: orbit designs whose producers' layer-3 NEED
rows are covered by PRESCRIBED z-streams.

For the accepted sigma (any orbit mode) the in-orbit producers of the own-branch Lemma-Q rows are the twisted class
columns of direction X = Hbar_k'^-1 sigma^j Hbar_c Z_o (Lemma E1-sigma).  In the 2-in-8 each has main rows on its own
branch k' only, at frame F = Hbar_k'(X) (pred N_k'), one per array.  A main row of array l with
d(anchor_l, F) = 3 that is not N1 or a branch frame is NEED: it has no automatic second producer.  It gets a z-producer
from any input-l t1-stream with a branch term AT F (tau(N_c) = F, possible because d(anchor_l, F) = 3 = d(M0, N_c)).
This driver adds, for every such NEED frame (both inputs), a z-stream with tau(N_c) = F, together with its own
sigma-orbit (so that its own Lemma-Q rows have in-orbit producers, as for the main orbit), and runs the exact layered
test (peel2 + modp/LM + alpha3's V4 when binomial).  Frame-level statistics are printed first.
usage: as jdesign.py (orbit modes); env ZC_MAX (max NEED frames covered, default 4)."""
import os, sys, random
import jd2
import alpha3
import jdesign as J
from endsolve import lift_P, find_P, pi_bar
from sp2 import mul as smul, inv as sinv, image, dist, D as DS, order as sorder

orig_sigma = J.orbit_sigma
orig_inj = J.orbit_injections
ST = {}
ZC_MAX = int(os.environ.get('ZC_MAX', '4'))


def need_frames(cf, u, o):
    """layer-3 NEED main rows of the in-orbit producers (2-in-8): list of (array l, frame F)"""
    sig = smul(DS, smul(u, DS))
    Hs = {'a': cf['Hs']['t3'], 'b': cf['Hs']['t3b']}; Nf = {'a': cf['N2'], 'b': cf['N2p']}
    Z = {0: image(sinv(Hs['a']), cf['A1']), 1: image(sinv(Hs['b']), cf['A2'])}
    out = []
    sj = sig
    for j in range(1, o):
        for o_ in (0, 1):
            for c in ('a', 'b'):
                Nk = image(sj, Nf[c])
                if Nk not in (Nf['a'], Nf['b']):
                    continue
                kp = 'a' if Nk == Nf['a'] else 'b'
                X = image(sinv(Hs[kp]), image(sj, image(Hs[c], Z[o_])))
                if X == Z[o_]:
                    continue                                   # twin: class family
                F = image(Hs[kp], X)
                sv = J.twisted_status(cf, X)
                if 'DEAD' in sv:
                    continue
                for l in (0, 1):
                    anc = cf['P1'] if l == 0 else cf['Q']
                    if dist(anc, F) == 3 and F != cf['N1'] and F not in (Nf['a'], Nf['b']):
                        out.append((l, F))
        sj = smul(sig, sj)
    seen = []; res = []
    for x in out:
        if x not in seen:
            seen.append(x); res.append(x)
    return res


def sigma_z(rng, cf, L):
    u, o, fixA, st = orig_sigma(rng, cf, L)
    if u is not None:
        nf = need_frames(cf, u, o)
        ST['need'] = nf
        print("  zcover: %d layer-3 NEED frames of in-orbit producers: %s" % (len(nf), [(l + 1, 'F') for (l, F) in nf]), flush=True)
    return u, o, fixA, st


def inj_z(rng, cf, L, l, u, n=0):
    js = orig_inj(rng, cf, L, l, u)
    o = sorder(u)
    added = 0
    for (ll, F) in ST.get('need', [])[:ZC_MAX * 2]:
        if ll != l or added >= ZC_MAX:
            continue
        pinv = sinv(pi_bar(L, l))
        c = (cf['N2'], cf['N2p'])[added % 2]
        b0 = find_P(rng, image(pinv, F), image(DS, c))
        if b0 is None:
            continue
        tau = smul(pi_bar(L, l), smul(sinv(b0), DS))
        assert image(tau, c) == F
        Um = lift_P(J.EXACT[0], u); Mi = lift_P(J.EXACT[0], b0); b = b0
        for i in range(o):
            js.append(dict(site='t1', l=l, where='zcov%d/%d' % (o, i), Y=None, letter=b, tau_order=o, lam=0,
                           mat=J.EXACT[0].snap(Mi)))
            b = smul(u, b); Mi = Um @ Mi
        added += 1
    print("  zcover x%d: %d z-stream orbits added (%d injections in all)" % (l + 1, added, len(js)), flush=True)
    return js


J.orbit_sigma = sigma_z
J.orbit_injections = inj_z

if __name__ == "__main__":
    print("jd_zcover config %s ZC_MAX %d argv %s" % (jd2.CFG, ZC_MAX, sys.argv[1:]), flush=True)
    J.main()
