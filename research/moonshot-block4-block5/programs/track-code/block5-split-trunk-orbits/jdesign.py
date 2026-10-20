"""jdesign.py (w7-b5): JOINT injection designs for the split-trunk engine 2 in 7 at block 5 (and the decoupled 2-in-8),
decided exactly with b5-twobranch's tools (split.py's skeletons and lifts, peel.py's exact forced-zero propagation,
endsolve's LM over the injection coefficients).  All designed injections of a round are offered AT ONCE (coefficients
free), so a stream's pollution can be produced by another designed stream (E1) or by the trunk (E2-E4).

Modes (argument `mode`):
  t2     only t2-injections (their pollution in the class columns sits at layer 3, never at layer 4; the producers of a
         layer-3 row with a layer-2 predecessor include t2-streams at that predecessor: companions);
  esc    t1-injections whose transports are chosen so that EVERY layer-4 Lemma-Q row of the stream (both outputs'
         class columns, both branches) lands on a trunk escape: E3 (stream term at N1 -> t2 y-column) or E4 (row frame
         = a branch frame at layer 4 -> t3 z-column) or E2; plus t2-injections as in `t2`;
  orbit  t1-injections in closed sigma-orbits {u^i beta} (Lemma E1-sigma): sigma in Stab(M0, N_a) of small order, all
         orbit streams in one frame; plus the class family's twisted columns come for free.
usage: jdesign.py seed profA profB lam1 lam2 mode rounds nper dlist [dec8]"""
import sys, time, random
import numpy as np
import split as S
import endsolve as es
from endsolve import lift_P, find_P, pi_bar, assemble, lm, vp_residual
from sp2 import M, L0, Fstar, D as DS, mul as smul, inv as sinv, image, dist, neighbors, order as sorder
import designer as dz
EXACT = [None]
LMS = [2, 50]
ZTRIES = [1]
DEC8MODE = [False]
SIDE = []
COPY = {}
SMART = 0
STATS = {}


def stream_tau(L, site, l, letter):
    pib = pi_bar(L, l)
    return smul(pib, smul(sinv(letter), DS if site == 't1' else L['c2']))


def t2_options(rng, cf, L, Sy, row, nper):
    """t2-injections giving `row` a y-producer (stream term of a t2-injection at the row's predecessor frame, layer 2)
    or z-producer (stream term at the row's frame); every feasible (frame, route) gets up to nper letters."""
    l, gid = row
    th = dz.theta(Sy, L, l, gid)
    G = image(th, L0); Pr = image(th, cf['M0'])
    anchor = [cf['P1'], cf['Q']][l]
    pinv = sinv(pi_bar(L, l))
    out = []
    for (N, where) in ((cf['N1'], 't2'), (cf['N2'], 't3'), (cf['N2p'], 't3b')):
        dN = dist(L0, N)
        for (X, route) in ((G, 'z'), (Pr, 'y')):
            if dist(anchor, X) != dN:
                continue
            got = 0
            for _ in range(SMART if route == 'y' else nper):
                b = find_P(rng, image(pinv, X), image(L['c2'], N))
                if b is None:
                    break
                tau = stream_tau(L, 't2', l, b)
                if route == 'y' and SMART > 0:
                    # producer: y-column on the generator holding the stream term; direction theta_term^-1(G)
                    if where == 't2':
                        Th = smul(tau, smul(sinv(L['c2']), DS)); Xd = image(sinv(Th), G)
                        Hm = smul(sinv(L['c2']), DS)
                        sv = [cstat(cf, ll, image(Hm, Xd), cf['N1']) for ll in (0, 1)]
                    else:
                        Hc = cf['Hs']['t3'] if where == 't3' else cf['Hs']['t3b']
                        Th = smul(tau, Hc); Xd = image(sinv(Th), G)
                        sv = twisted_status(cf, Xd)
                    STATS['smart_try'] = STATS.get('smart_try', 0) + 1
                    if 'DEAD' in sv:
                        continue
                    STATS['smart_ok'] = STATS.get('smart_ok', 0) + 1
                out.append(dict(site='t2', l=l, where='%s/%s' % (where, route), Y=None, letter=b,
                                tau_order=sorder(tau), lam=dist(anchor, G)))
                got += 1
                if got >= nper:
                    break
    return out


def lemmaQ_rows_ok(cf, L, j, mode_ok=('E3', 'E4', 'E2')):
    """for a t1-injection j (transport tau), list its layer-4 Lemma-Q rows (class directions of both outputs, both
    branches) and whether each lands on a trunk escape.  Returns (all_ok, list)."""
    l = j['l']; anc = cf['P1'] if l == 0 else cf['Q']
    tau = stream_tau(L, 't1', l, j['letter'])
    Hs = [cf['Hs']['t3'], cf['Hs']['t3b']]; Ns = [cf['N2'], cf['N2p']]
    Zs = [image(sinv(Hs[0]), cf['A1']), image(sinv(Hs[1]), cf['A2'])]
    res = []
    for o, Z in enumerate(Zs):
        for k, (Hk, Nk) in enumerate(zip(Hs, Ns)):
            F = image(smul(tau, Hk), Z); Pr = image(tau, Nk)
            lay = dist(anc, F)
            if lay < 4:
                continue
            how = None
            if Pr == cf['N1']:
                how = 'E3'
            elif any(Pr == Ni and dist(anc, Ni) == lay - 1 for Ni in Ns):
                how = 'E2'
            elif any(F == Ni and dist(anc, Ni) == lay for Ni in Ns):
                how = 'E4'
            res.append((o, k, lay, how))
    ok = all(h is not None and h in mode_ok for (_, _, _, h) in res)
    return ok, res


def esc_injections(rng, cf, L, l, n, tries=4000):
    """t1-injections of x_l whose layer-4 Lemma-Q rows all land on trunk escapes.  Strategy: force one coincidence
    exactly (tau(N_c) = N1 for c in {a, b}, i.e. E3 for the rows with predecessor tau(N_c)), then sample the rest."""
    anc = cf['P1'] if l == 0 else cf['Q']
    pinv = sinv(pi_bar(L, l))
    out = []; stats = {}
    Ns = [cf['N2'], cf['N2p']]
    for t in range(tries):
        c = t % 2
        # tau(N_c) = N1 with tau = pibar beta^-1 Dbar: beta(pibar^-1 N1) = Dbar N_c
        b = find_P(rng, image(pinv, cf['N1']), image(DS, Ns[c]))
        if b is None:
            continue
        j = dict(site='t1', l=l, where='esc/E3b%d' % c, Y=None, letter=b, tau_order=0, lam=0)
        ok, res = lemmaQ_rows_ok(cf, L, j)
        key = tuple(sorted((o, k, h) for (o, k, _, h) in res))
        stats[key] = stats.get(key, 0) + 1
        if ok:
            j['tau_order'] = sorder(stream_tau(L, 't1', l, b))
            j['where'] += str(sorted(set(h for (_, _, _, h) in res)))
            out.append(j)
            if len(out) >= n:
                break
    return out, stats


def cstat(cf, l, F, Pr):
    """second-order status of a main row (e1census.status, split.py's cf keys)"""
    anc = cf['P1'] if l == 0 else cf['Q']
    A = cf['A1'] if l == 0 else cf['A2']
    lam = dist(anc, F); brf = (cf['N2'], cf['N2p'])
    if lam == 5:
        return 'FREE' if F == A else 'DEAD'
    if lam == 4:
        if F in brf:
            return 'AUTO'
        return 'NEED' if dist(anc, Pr) == 3 else 'DEAD'
    if lam == 3:
        return 'AUTO' if (F == cf['N1'] or F in brf) else 'NEED'
    return 'LOW'


def twisted_status(cf, X):
    out = []
    for (H, N) in ((cf['Hs']['t3'], cf['N2']), (cf['Hs']['t3b'], cf['N2p'])):
        F = image(H, X)
        for l in (0, 1):
            out.append(cstat(cf, l, F, N))
    return out


ORB = dict(maxord=24, tries=40000, both=True, twin=None, swap=False)


def orbit_sigma(rng, cf, L):
    """sigma = Dbar ubar Dbar in Stab(M0) with sigma{N_a, N_b} = {N_a, N_b} (ORB['swap']: sigma swaps them, else fixes
    both), order <= maxord.  Orbit streams tau_i = tau_0 sigma^-i then keep all their t3 terms in the two frames
    tau_0 N_a, tau_0 N_b, so EVERY Lemma-Q row of every orbit stream has E1 producers inside the orbit (Lemma E1-sigma):
    the class column of output o twisted by sigma^j, direction Hbar_k'^-1 sigma^j Hbar_c Z_o (k' = sigma^j(c)).
    Accept sigma if for every (output, branch c) at least one twist j != 0 has no DEAD main row (DEAD twisted columns
    are forced to 0 but the others can still cancel the circulant pollution system).  ORB['twin'] as before."""
    Xa = image(DS, cf['N2']); Xb = image(DS, cf['N2p'])
    Hs = {'a': cf['Hs']['t3'], 'b': cf['Hs']['t3b']}; Nf = {'a': cf['N2'], 'b': cf['N2p']}
    Z = {0: image(sinv(Hs['a']), cf['A1']), 1: image(sinv(Hs['b']), cf['A2'])}
    st = {}
    for t in range(ORB['tries']):
        if ORB.get('swap'):
            u = find_P(rng, Xa, Xb)
            if image(u, Xb) != Xa:
                st['notswap'] = st.get('notswap', 0) + 1; continue
        else:
            u = find_P(rng, Xa, Xa)
            if image(u, Xb) != Xb:
                st['notb'] = st.get('notb', 0) + 1; continue
        o = sorder(u)
        if o is None or o < 2 or o > ORB['maxord']:
            st['ord'] = st.get('ord', 0) + 1; continue
        sig = smul(DS, smul(u, DS))
        fixA = image(sig, cf['A1']) == cf['A1'] and image(sig, cf['A2']) == cf['A2']
        if ORB['twin'] is True and not fixA or ORB['twin'] is False and fixA:
            st['twin'] = st.get('twin', 0) + 1; continue
        alive = {}
        sj = sig
        for j in range(1, o):
            for o_ in (0, 1):
                for c in ('a', 'b'):
                    Nk = image(sj, Nf[c]); kp = 'a' if Nk == Nf['a'] else 'b'
                    assert Nk == Nf[kp]
                    X = image(sinv(Hs[kp]), image(sj, image(Hs[c], Z[o_])))
                    sv = twisted_status(cf, X) if not (DEC8MODE[0] and False) else None
                    if 'DEAD' not in sv:
                        alive.setdefault((o_, c), []).append((j, sv.count('NEED'), 'twin' if X == Z[o_] else 'non'))
            sj = smul(sig, sj)
        if ORB.get('nontwin_own'):
            need = [(0, 'a'), (1, 'b')] + [(o_, c) for (o_, c) in ((0, 'b'), (1, 'a'))
                                           if ORB.get('all4') or dist(cf['M0'], image(Hs[c], Z[o_])) == 4]
            okk = all(any(x[2] == 'non' for x in alive.get(key, [])) for key in need)
            if not okk:
                st['dead'] = st.get('dead', 0) + 1; continue
        elif len(alive) < 4:
            st['dead'] = st.get('dead', 0) + 1; continue
        return u, o, fixA, dict(st, alive=alive)
    return None, None, None, st


def orbit_injections(rng, cf, L, l, u, n=0):
    """the closed sigma-orbit of x_l t1-injections: letters beta_i = u^i beta_0 (exact products of lifts)."""
    from sp2 import rand_P
    o = sorder(u)
    b0 = rand_P(rng)
    if ORB.get('zl'):
        # tau_0(N_c) = A_other (a z-producer of the other array's leftover L at A_o, W3 Lemma L3's route), c = branch
        A = cf['A1'] if l == 1 else cf['A2']
        anc = cf['Q'] if l == 1 else cf['P1']
        pinv = sinv(pi_bar(L, l))
        Nown = cf['N2'] if l == 1 else cf['N2p']        # the leftover's predecessor frame (N_a for L2, N_b for L1)
        best = None
        for tz in range(ZTRIES[0]):
            c = (cf['N2'], cf['N2p'])[tz % 2]
            Hc = cf['Hs']['t3'] if c == cf['N2'] else cf['Hs']['t3b']
            if dist(anc, A) != dist(cf['M0'], c):
                continue
            bz = find_P(rng, image(pinv, A), image(DS, c))
            if bz is None:
                continue
            tau = stream_tau(L, 't1', l, bz)
            W = image(sinv(Hc), image(sinv(tau), Nown))
            sv = []
            for (Hp, Np) in ((cf['Hs']['t3'], cf['N2']), (cf['Hs']['t3b'], cf['N2p'])):
                if DEC8MODE[0] and Np != c:
                    continue
                Pr = image(Hp, W)
                for ll in (0, 1):
                    an = cf['P1'] if ll == 0 else cf['Q']
                    lam = dist(an, Np)
                    if lam <= 2:
                        sv.append('LOW')
                    elif lam == 4:
                        if Pr == cf['N1'] or (Pr in (cf['N2'], cf['N2p']) and dist(an, Pr) >= 3):
                            sv.append('AUTO')
                        elif dist(an, Pr) == 3:
                            sv.append('NEED')
                        else:
                            sv.append('DEAD')
                    else:
                        sv.append('NEED')
            score = (sv.count('DEAD'), sv.count('NEED'))
            if best is None or score < best[0]:
                best = (score, bz)
            if score == (0, 0):
                break
        if best is not None:
            b0 = best[1]
            STATS['zl_score_x%d' % (l + 1)] = best[0]
    out = []; b = b0
    Um = lift_P(EXACT[0], u); Bm = lift_P(EXACT[0], b0); Mi = Bm
    for i in range(o):
        out.append(dict(site='t1', l=l, where='orbit%d/%d' % (o, i), Y=None, letter=b, tau_order=o, lam=0,
                        mat=EXACT[0].snap(Mi)))
        b = smul(u, b); Mi = Um @ Mi
    assert b == b0
    return out[:n] if n else out


def run(seed, profA, profB, lam1, lam2, mode, rounds, nper, dlist, dec8=False, cap=4000):
    S.LIMY = 3 if (lam1 == 3 and lam2 == 3) else 4
    S.DECOUPLE = dec8; S.PEEL = True; S.ESC = False
    rng = random.Random(seed); t0 = time.time()
    cf, L = S.make_config_split(rng, profA, profB, lam1, lam2)
    if cf is None:
        print("seed %d: no config" % seed); return
    print("seed %d mode %s%s: config %s (%.0fs)" % (seed, mode, ' 2-in-8' if dec8 else '', cf['info'], time.time() - t0), flush=True)
    E = es.Eng(M, seed=seed); Sy = es.Symp(E); EXACT[0] = E
    injs = []
    if COPY.get('s2b'):
        # theta' = qbar qbar'^-1 = g2 eps g2^-1 in Stab(Q); it fixes A2 iff eps fixes O2 = g2^-1 A2 = sigma2bar F*
        g2 = smul(DS, smul(L['a2'], DS))
        O2 = image(L['s2'], Fstar)
        assert image(smul(g2, L['s2']), Fstar) == cf['A2']
        for t in range(20000):
            eps = find_P(rng, O2, O2)
            th = smul(g2, smul(eps, sinv(g2)))
            assert image(th, cf['Q']) == cf['Q'] and image(th, cf['A2']) == cf['A2']
            if COPY.get('fixNb') and image(th, cf['N2p']) != cf['N2p']:
                continue
            if image(th, L0) == L0:
                continue                       # the copy must leave the main frames somewhere
            break
        L['s2b'] = smul(sinv(eps), L['s2'])
        prof = [(nm, dist(cf['Q'], image(th, X)), image(th, X) == X) for nm, X in
                (('L0', L0), ('N1', cf['N1']), ('Na', cf['N2']), ('Nb', cf['N2p']), ('A1', cf['A1']))]
        print("  copy (two-letter start): theta' moves/fixes %s" % prof, flush=True)
        cf['theta_copy'] = th
    if mode == 'esc':
        for l in (0, 1):
            js, st = esc_injections(rng, cf, L, l, nper)
            print("  esc x%d: %d admissible of tries; pollution patterns %s" % (l + 1, len(js), dict(list(st.items())[:12])), flush=True)
            injs += js
    if mode in ('e3orb', 'e3orb1'):
        # minimal closed design: every orbit stream's R-row with predecessor tau(N_own') lands on E3 (tau(N) = N1),
        # the other own-branch row is produced inside the orbit by NON-class twisted columns.  x2 streams:
        # tau_0(N_a) = N1, twists nu_b^j Z_2 (out2's class column via b); x1 streams: tau_0(N_b) = N1, twists nu_a^j Z_1.
        # sigma in Stab(M0) fixes N_a and N_b, so tau_i(N_a), tau_i(N_b) are the same for the whole orbit.
        Xa = image(DS, cf['N2']); Xb = image(DS, cf['N2p'])
        Hs = {'a': cf['Hs']['t3'], 'b': cf['Hs']['t3b']}
        Z = {0: image(sinv(Hs['a']), cf['A1']), 1: image(sinv(Hs['b']), cf['A2'])}
        st = {}; found = None
        for t in range(ORB['tries']):
            u = find_P(rng, Xa, Xa)
            if image(u, Xb) != Xb:
                st['notb'] = st.get('notb', 0) + 1; continue
            o = sorder(u)
            if o is None or o < 2 or o > ORB['maxord']:
                st['ord'] = st.get('ord', 0) + 1; continue
            sig = smul(DS, smul(u, DS))
            ok = True; info = []
            for (o_, c) in ((1, 'b'), (0, 'a')):
                nu = smul(sinv(Hs[c]), smul(sig, Hs[c])); X = Z[o_]; nal_ = 0
                for j in range(1, o):
                    X = image(nu, X)
                    if X == Z[o_]:
                        info.append((o_, c, j, 'twin')); continue
                    sv = twisted_status(cf, X)
                    if 'DEAD' in sv:
                        info.append((o_, c, j, 'DEAD'))
                        if mode == 'e3orb':
                            ok = False; break
                        continue
                    info.append((o_, c, j, sv.count('NEED'))); nal_ += 1
                if mode == 'e3orb1' and nal_ == 0:
                    ok = False
                if not ok:
                    break
            key = tuple(sorted(set(x[3] if isinstance(x[3], str) else 'alive' for x in info)))
            st[key] = st.get(key, 0) + 1
            if not ok:
                st['dead_or_twin'] = st.get('dead_or_twin', 0) + 1; continue
            found = (u, o, info); break
        print("  e3orb: %s" % st, flush=True)
        if found is None:
            return 'nosigma', None
        u, o, info = found
        print("  e3orb sigma order %d, twists %s" % (o, info), flush=True)
        for l, Nown in ((1, cf['N2']), (0, cf['N2p'])):
            pinv = sinv(pi_bar(L, l))
            b0 = find_P(rng, image(pinv, cf['N1']), image(DS, Nown))
            Um = lift_P(E, u); Mi = lift_P(E, b0); b = b0
            n = o if mode == 'e3orb' else 2
            for i in range(o):
                injs.append(dict(site='t1', l=l, where='e3orb%d/%d' % (o, i), Y=None, letter=b, tau_order=o, lam=0,
                                 mat=E.snap(Mi)))
                b = smul(u, b); Mi = Um @ Mi
        mode = 'orbit'
    if mode in ('side', 'sideo'):
        # Lemma T's gap: a side term gamma t1 inside t3 puts BOTH arrays' t1 main terms at W = gammabar^-1 M0 (a
        # neighbour of L0), an E2-type trunk producer for Lemma-Q rows whose predecessor is W.  Design, per input:
        # a t1-stream with tau(N_a) = N1 (E3 for R1) and W := tau(N_b) ~ L0 (then gamma: W -> M0 gives R2 an E2'').
        pinvs = {l: sinv(pi_bar(L, l)) for l in (0, 1)}
        st = {}
        for l in (1, 0):
            got = 0
            for t in range(6000):
                c_first = 0 if (mode == 'side' or t % 2 == 0) else 1
                Ns = [cf['N2'], cf['N2p']]
                b = find_P(rng, image(pinvs[l], cf['N1']), image(DS, Ns[c_first]))
                if b is None:
                    continue
                tau = stream_tau(L, 't1', l, b)
                W = image(tau, Ns[1 - c_first])
                key = ('x%d' % (l + 1), 'd(L0,W)', dist(L0, W))
                st[key] = st.get(key, 0) + 1
                if dist(L0, W) != 1:
                    continue
                g = find_P(rng, W, cf['M0'])
                SIDE.append(g)
                injs.append(dict(site='t1', l=l, where='side/E3b%d' % c_first, Y=None, letter=b, tau_order=0, lam=0))
                got += 1
                if got >= max(1, nper):
                    break
        print("  side design: %d injections, %d side terms; W layers %s" % (len(injs), len(SIDE), st), flush=True)
        if not SIDE:
            return 'noside', None
    if mode in ('orbitC', 'orbitCN'):
        ORB['nontwin'] = (mode == 'orbitCN')
        # commuting construction: c3b = c3a eps^-1 with eps in P, eps(X_a) = X_b (X_c = Dbar c2bar N_c); sigma~ = eps^m
        # (m = length of X_a's eps-orbit) commutes with eps, so sigma = c2bar^-1 Dbar sigma~ Dbar c2bar has nu_a = nu_b
        # (one cyclic family of twisted columns serves both branches).  Keep sigma if it fixes M0, is not 1, and the
        # class columns' other-branch images and some twist per (output, branch) stay non-DEAD.
        c2 = L['c2']; Xa = image(DS, image(c2, cf['N2'])); Xb = image(DS, image(c2, cf['N2p']))
        found = None; st = {}
        for t in range(ORB['tries']):
            eps = find_P(rng, Xa, Xb)
            m = 1; Y = image(eps, Xa)
            while Y != Xa and m < 200:
                Y = image(eps, Y); m += 1
            if m >= 200:
                continue
            st['m%d' % m] = st.get('m%d' % m, 0) + 1
            em = list(range(0)) or None
            em = eps
            for _ in range(m - 1):
                em = smul(eps, em)
            sig = smul(sinv(c2), smul(DS, smul(em, smul(DS, c2))))
            if image(sig, cf['M0']) != cf['M0']:
                st['notM0'] = st.get('notM0', 0) + 1; continue
            u = smul(DS, smul(sig, DS))
            o = sorder(u)
            if o is None or o < 2 or o > ORB['maxord']:
                st['ord'] = st.get('ord', 0) + 1; continue
            c3b = smul(L['c3'], sinv(eps))
            Hb = smul(sinv(c2), smul(DS, smul(sinv(c3b), DS)))
            if image(Hb, L0) != cf['N2p']:
                st['frame'] = st.get('frame', 0) + 1; continue
            cf2 = dict(cf); cf2['Hs'] = {'t3': cf['Hs']['t3'], 't3b': Hb}
            Za = image(sinv(cf['Hs']['t3']), cf['A1']); Zb = image(sinv(Hb), cf['A2'])
            if 'DEAD' in twisted_status(cf2, Za) or 'DEAD' in twisted_status(cf2, Zb):
                st['classdead'] = st.get('classdead', 0) + 1; continue
            nuA = smul(sinv(cf['Hs']['t3']), smul(sig, cf['Hs']['t3'])); nuB = smul(sinv(Hb), smul(sig, Hb))
            assert nuA == nuB, "commuting construction failed"
            if ORB.get('nontwin') and (image(nuA, Za) == Za or image(nuA, Zb) == Zb):
                st['twin'] = st.get('twin', 0) + 1; continue
            alive = {}
            for (o_, Z) in ((0, Za), (1, Zb)):
                X = Z
                for j in range(1, o):
                    X = image(nuA, X)
                    sv = twisted_status(cf2, X)
                    if 'DEAD' not in sv:
                        alive.setdefault(o_, []).append((j, sv.count('NEED')))
            st['cand'] = st.get('cand', 0) + 1
            if len(alive) < 2 and not ORB.get('any'):
                st['dead'] = st.get('dead', 0) + 1; continue
            found = (u, o, c3b, Hb, alive); break
        print("  orbitC: %s" % ({k: v for k, v in st.items()}), flush=True)
        if found is None:
            return 'nosigma', None
        u, o, c3b, Hb, alive = found
        L['c3b'] = c3b; cf['Hs']['t3b'] = Hb
        print("  orbitC sigma order %d, alive twists %s" % (o, alive), flush=True)
        for l in (0, 1):
            js = orbit_injections(rng, cf, L, l, u)
            injs += js
        mode = 'orbit'
    if mode.startswith('orbit') and mode != 'orbit':
        ORB['twin'] = {'orbitT': True, 'orbitN': False, 'orbitS': None, 'orbitZ': None, 'orbitSZ': None, 'orbitF': None, 'orbitSN': None, 'orbitSNZ': None, 'orbitSN4': None, 'orbitSNZ4': None}[mode]
        ORB['swap'] = mode in ('orbitS', 'orbitSZ', 'orbitSN', 'orbitSNZ', 'orbitSN4', 'orbitSNZ4')
        ORB['nontwin_own'] = mode in ('orbitSN', 'orbitSNZ', 'orbitSN4', 'orbitSNZ4')
        ORB['all4'] = mode in ('orbitSN4', 'orbitSNZ4')
        ORB['zl'] = mode in ('orbitZ', 'orbitSZ', 'orbitSNZ', 'orbitSNZ4')
        u, o, fixA, st = orbit_sigma(rng, cf, L)
        print("  orbit sigma: order %s fixes A1,A2 %s; rejections %s" % (o, fixA, st), flush=True)
        if u is None:
            return 'nosigma', None
        for l in (0, 1):
            js = orbit_injections(rng, cf, L, l, u)
            print("  orbit x%d: %d injections" % (l + 1, len(js)), flush=True)
            injs += js
        mode = 'orbit'
    for d in dlist:
        for rd in range(rounds + 1):
            P = dz.build_prog2(E, Sy, L, injs, 0.7)
            P.mats = dict(P.mats); P.mats['c3b'] = lift_P(E, L['c3b'])
            if SIDE and not dec8:
                extra = []
                for Gs in SIDE:
                    Gm = lift_P(E, Gs)
                    extra += [(ll, P.mulD(Gm, g), {k: 0.8 * v for k, v in c.items()}) for (ll, g, c) in P.gens['t1']]
                P.gens['t3'] = es.merge(P.gens['t3'] + extra)
            if dec8:
                t2 = P.gens['t2']
                P.gens = {'t1': P.gens['t1'], 't2': t2,
                          'ua': es.merge([(ll, P.mulD(P.mats['c3'], g), c) for (ll, g, c) in t2]),
                          'ub': es.merge([(ll, P.mulD(P.mats['c3b'], g), c) for (ll, g, c) in t2])}
                if SIDE:
                    extra = []
                    for Gs in SIDE:
                        Gm = lift_P(E, Gs)
                        extra += [(ll, P.mulD(Gm, g), {k: 0.8 * v for k, v in c.items()}) for (ll, g, c) in P.gens['t1']]
                    P.gens['ua'] = es.merge(P.gens['ua'] + extra)
            seeds = [S.class_seeds(P, L, cf, l, random.Random(7 + l)) for l in (0, 1)]
            cls = []; trig = []
            for l in (0, 1):
                cl = S.closure_any(P, l, d, seeds[l], cap=cap)
                if cl is None:
                    print("  d=%d round %d: closure cap (%d injections)" % (d, rd, len(injs)), flush=True)
                    cls = None; break
                cls.append(cl); trig.append(list(S.LAST_TRIG))
            if cls is None:
                return 'cap', d
            singles = []
            for li, cl in enumerate(cls):
                singles += [(li, r) for r in S.diag(P, cl, li)] + [(li, r) for r in trig[li]]
            Nss = [assemble(cl, P.nal, rng) for cl in cls]
            cost, al = lm(Nss, P.nal, rng, starts=LMS[0], iters=LMS[1])
            print("  d=%d round %d: injections %d, closures %s, residual %.3e (%.0fs)" % (
                d, rd, len(injs), [(len(c[0]), len(c[1]), len(c[3][0])) for c in cls], cost, time.time() - t0), flush=True)
            if cost < 1e-20:
                break
            if rd == rounds:
                break
            if mode == 'orbit' and rounds == 0:
                break
            added = 0; seen = set()
            for (li, row) in singles:
                if row in seen:
                    continue
                seen.add(row)
                js = t2_options(rng, cf, L, Sy, row, nper)
                injs += js; added += len(js)
            print("     round %d: %d single/trigger rows, added %d t2-injections; smart %s" % (rd, len(seen), added, STATS), flush=True)
            if added == 0:
                break
        if cost > 1e-12:
            print("STOP at d=%d (residual %.3e)" % (d, cost), flush=True)
            return 'stop', d
    print("ALL LAYERS PASS: candidate (verify kappa + formal identity)", flush=True)
    return 'pass', None


def main():
    a = sys.argv[1:]
    seed = int(a[0])
    pp = lambda s: None if s == '0' else (int(s[0]), int(s[1]))
    profA, profB = pp(a[1]), pp(a[2])
    lam1 = None if a[3] == '0' else int(a[3]); lam2 = None if a[4] == '0' else int(a[4])
    mode = a[5]; rounds = int(a[6]); nper = int(a[7])
    dlist = [int(c) for c in a[8]]
    dec8 = len(a) > 9 and a[9] == '8'
    DEC8MODE[0] = dec8
    import os
    ZTRIES[0] = int(os.environ.get('ZTRIES', '1')); LMS[0] = int(os.environ.get('LMSTARTS', '2')); LMS[1] = int(os.environ.get('LMITERS', '50'))
    if len(a) > 11:
        for tok in a[11].split(','):
            if tok:
                COPY[tok] = True
    global SMART
    SMART = int(a[10]) if len(a) > 10 else 0
    t0 = time.time()
    r = run(seed, profA, profB, lam1, lam2, mode, rounds, nper, dlist, dec8, cap=int(os.environ.get('CAP', '4000')))
    print("result %s (%.0fs) stats %s" % (r, time.time() - t0, STATS), flush=True)


if __name__ == "__main__":
    main()
