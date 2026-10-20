"""split3.py (b5-twobranch): FRUSTRATED SAME-FRAME FACTORS on a multi-branch split trunk (engine 2 in 7, block 5).

Mechanism (RUNLOG 20:2x): every class term of array 1 pairs (cK, cL2) with the same c (Lemma Z), unless the array-2 partner
h_k q has a same-frame companion inside t3.  With branch frames N_k at distance <= 3 from the other anchor, an injection of the
other input can put its t3 term in the coset B h_k q (t1-injection if d = 3, t2-injection if d = 2): the class column's
leftover becomes (1 + lam w') L2.  Two families per array with factors whose left annihilators generate C[B] escape Lemma Z.

Program: s_l = D sig_l x_l; t1 = D(a1 s1 + a2 s2 + t1-injections); t2 = D(c2 t1 + t2-injections);
t3 = D(sum_k mu_k c3_k) t2 with 4 branches: k = 0, 1 serve array 1 (d(P1, N_k) = 4, d(Q, N_k) in {2,3}, N_k ~ A1),
k = 2, 3 serve array 2 (mirror).  For each branch one designed injection of the OTHER input lands its t3 term (own branch) in
the coset of the branch's other-array main term, with wbar fixing the class direction Z_k (so the factor stays in the leftover
coset).  Layered exact-structure solve (endsolve), K-part normalisation, then the kappa test (natural representation).
usage: split3.py seed nfac LIMY dlist [mus]     nfac = number of factor injections per branch (0 = control)"""
import sys, time, random
import numpy as np
import endsolve as es
from endsolve import (Eng, Symp, lift_P, find_P, opp_graph_T, merge, assemble, lm, vp_residual, pi_bar, build_prog, IDM)
from sp2 import (M, L0, Fstar, D as DS, mul as smul, inv as sinv, image, dist, neighbors, order as sorder)
from split import H_of, closure_multi, diag


def make_config4(rng, LIMY=3, tries=300, nlet=300, n1=2, n2=2):
    """n1 branches for array 1 (d(P1,N)=4, d(Q,N) in {2,3}, all adjacent to A1), n2 for array 2 (mirror); every class
    family's images through the other branches at layer <= LIMY for both arrays (or on a target)."""
    M0 = image(DS, L0)
    N1s = [N for N in sorted(neighbors(L0)) if dist(M0, N) == 2]
    for t in range(tries):
        N1 = rng.choice(N1s)
        brs = [N for N in sorted(neighbors(N1)) if dist(L0, N) == 2 and dist(M0, N) == 3]
        anchors = [X for X in sorted(neighbors(M0)) if dist(X, L0) == 2 and dist(X, N1) == 3]
        P1, Q = rng.sample(anchors, 2)
        if dist(P1, Q) != 2:
            continue
        fams = []
        for (an, other, nn) in ((P1, Q, n1), (Q, P1, n2)):
            if nn == 0:
                opp = [A for A in neighbors(rng.choice([N for N in brs if dist(an, N) == 4] or brs)) if dist(an, A) == 5]
                if not opp:
                    break
                fams.append((rng.choice(sorted(opp)), []))
                continue
            byA = {}
            for N in brs:
                if dist(an, N) == 4 and dist(other, N) in (2, 3):
                    for A in neighbors(N):
                        if dist(an, A) == 5 and dist(other, A) <= 4:
                            byA.setdefault(A, []).append(N)
            good = [(A, Ns) for A, Ns in byA.items() if len(Ns) >= nn]
            if not good:
                break
            A, Ns = rng.choice(sorted(good, key=lambda x: str(x[0])))
            fams.append((A, rng.sample(Ns, nn)))
        if len(fams) < 2:
            continue
        (A1, f1), (A2, f2) = fams
        Nk = f1 + f2
        own = [0] * len(f1) + [1] * len(f2)
        if len(set(Nk)) < len(Nk) or not Nk:
            continue
        c2 = find_P(rng, N1, M0)
        for _ in range(nlet):
            c3 = [find_P(rng, image(DS, image(c2, N)), M0) for N in Nk]
            Hs = [H_of(c2, c) for c in c3]
            ok = True
            for k in range(len(Nk)):
                A = A1 if own[k] == 0 else A2
                Z = image(sinv(Hs[k]), A)
                for j in range(len(Nk)):
                    if j == k or not ok:
                        continue
                    Y = image(Hs[j], Z)
                    for (anc, tgt) in ((P1, A1), (Q, A2)):
                        if dist(anc, Y) > LIMY and Y != tgt:
                            ok = False
            if ok:
                cf = dict(M0=M0, N1=N1, P1=P1, Q=Q, A1=A1, A2=A2, Nk=Nk, Hs=Hs, own=own,
                          info=dict(prof=[(dist(P1, N), dist(Q, N)) for N in Nk], own=own, lam=(dist(Q, A1), dist(P1, A2)), tries=t))
                L = dict(c2=c2, c3=c3[0], c3x=c3[1:])
                L['a1'] = find_P(rng, M0, image(DS, P1)); L['a2'] = find_P(rng, M0, image(DS, Q))
                for l, (a, A) in enumerate(((L['a1'], A1), (L['a2'], A2))):
                    O = image(DS, image(sinv(a), image(DS, A)))
                    T = opp_graph_T(O)
                    sg = es.unip([T[i] ^ (1 << i) for i in range(M)])
                    assert image(sg, Fstar) == O
                    L['s%d' % (l + 1)] = sg
                return cf, L
    return None, None


def factor_injections(rng, cf, L, nfac, maxs=4000):
    """for each branch k: nfac injections of the OTHER input whose own-branch-k t3 term lies in the coset of the branch's
    other-array main term, with wbar fixing Z_k.  Returns endsolve-style injection dicts (site, l, letter)."""
    out = []
    for k in range(len(cf['Nk'])):
        own = cf['own'][k]
        tgt_arr = 1 - own                           # the other array (0-based): its main term h_k pi_tgt gets the factor
        an_t = cf['Q'] if own == 0 else cf['P1']
        A, an = (cf['A1'], cf['P1']) if own == 0 else (cf['A2'], cf['Q'])
        Nk, Hk = cf['Nk'][k], cf['Hs'][k]
        Z = image(sinv(Hk), A)
        dd = dist(an_t, Nk)
        pib = pi_bar(L, tgt_arr); pinv = sinv(pib)
        if dd == 3:
            site = 't1'; src = image(DS, Nk)          # tau(N_k) = N_k with tau = pibar beta^-1 Dbar
        elif dd == 2:
            site = 't2'; src = image(L['c2'], Nk)     # tau'(N_k) = N_k with tau' = pibar delta^-1 c2bar
        else:
            continue
        got = 0
        for s in range(maxs):
            b = find_P(rng, image(pinv, Nk), src)
            if b is None:
                break
            tau = smul(pib, smul(sinv(b), DS if site == 't1' else L['c2']))
            assert image(tau, Nk) == Nk
            # w = (injection's branch-k t3 term)(h_k pi)^-1; its symplectic image: wbar = (D c3k D b)bar ... computed exactly
            # below on matrices; here the frame test: theta_inj = tau Theta_k, theta_main = Theta_k, so
            # wbar^-1 acts as Theta_k^-1 tau Theta_k on global->local frames: require it to fix Z.
            Tk = Hk
            m = smul(sinv(Tk), smul(tau, Tk))
            if image(m, Z) != Z:
                continue
            out.append(dict(site=site, l=tgt_arr, where='fac%d' % k, Y=None, letter=b, tau_order=sorder(tau), lam=dd, k=k))
            got += 1
            if got >= nfac:
                break
        print("   branch %d: %d factor injection(s) (%s, %d samples)" % (k, got, site, s + 1), flush=True)
    return out


def build_prog4(E, Sy, L, injs, mus):
    Lm = {k: v for k, v in L.items() if k != 'c3x'}
    P = build_prog(E, Sy, Lm, injs)
    P.mats = dict(P.mats)
    t2 = P.gens['t2']
    t3 = [(l, g, {kk: mus[0] * v for kk, v in c.items()}) for (l, g, c) in P.gens['t3']]
    for i, cx in enumerate(L['c3x']):
        Mx = lift_P(E, cx)
        assert Sy.of(E.intern(E.snap(Mx))) == cx
        P.mats['c3x%d' % i] = Mx
        t3 += [(l, P.mulD(Mx, g), {kk: mus[i + 1] * v for kk, v in c.items()}) for (l, g, c) in t2]
    P.gens['t3'] = merge(t3)
    return P


def seeds4(P, L, cf, l, rng):
    E, Sy = P.E, P.S
    seeds = []
    A = cf['A1'] if l == 0 else cf['A2']
    t2main = [g for (ll, g, c) in P.gens['t2'] if ll == l and c == {0: 1.0}][0]
    for k in range(len(cf['Nk'])):
        if dist(cf['Nk'][k], A) != 1:
            continue
        Mk = P.mats['c3'] if k == 0 else P.mats['c3x%d' % (k - 1)]
        g = P.mulD(Mk, t2main)
        if E.layer(g) != M - 1:
            continue
        X = image(Sy.of(g), Fstar)
        b = find_P(rng, X, image(DS, L0))
        Bm = lift_P(E, b)
        r = E.intern(E.snap(E.D @ Bm @ E.MATS[g]))
        assert Sy.frame(r) == Fstar
        seeds.append(('t3', E.intern(E.snap(Bm))))
    return seeds


def kappa_test(P, cls, Nss, al, nsol=3, rng=None):
    """K-parts of random solutions of each output (natural 32-dim representation): rank of the 64 x 64 block."""
    E = P.E
    Kinv = np.conj(E.K).T
    sols = []
    for li, (cl, Ns) in enumerate(zip(cls, Nss)):
        A = Ns[0].toarray()
        for k in range(1, len(Ns)):
            if Ns[k].nnz:
                A = A + al[k - 1] * Ns[k].toarray()
        t = np.zeros(A.shape[0], complex); t[0] = 1
        w0 = np.linalg.lstsq(A, t, rcond=1e-11)[0]
        U, S, Vh = np.linalg.svd(A, full_matrices=True)
        rk = int((S > 1e-9 * S[0]).sum()) if len(S) else 0
        sols.append((cl, w0, Vh[rk:].conj().T))
    out = []
    for trial in range(nsol):
        ks = []
        for (cl, w0, Nb) in sols:
            z = (np.random.randn(Nb.shape[1]) + 1j * np.random.randn(Nb.shape[1])) if Nb.shape[1] else np.zeros(0)
            w = w0 + (Nb @ z if Nb.shape[1] else 0)
            rows, cols, ent, kl = cl
            kap = [np.zeros((E.n, E.n), complex), np.zeros((E.n, E.n), complex)]
            for (typ, u, b), c in cols.items():
                if typ != 'y' or abs(w[c]) < 1e-14:
                    continue
                pre = E.D @ E.MATS[b]
                for (l2, g2, cf) in P.gens[u]:
                    if E.layer(g2) < M - 1:
                        continue
                    G = E.snap(pre @ E.MATS[g2]); gid = E.intern(G)
                    if E.layer(gid) == M and P.S.frame(gid) == Fstar:
                        coef = cf.get(0, 0) if len(cf) == 1 else sum(v * (al[k - 1] if k else 1) for k, v in cf.items())
                        kap[l2] += w[c] * coef * (G @ Kinv)
            ks.append(kap)
        if len(ks) == 2:
            Kb = np.block([[ks[0][0], ks[0][1]], [ks[1][0], ks[1][1]]])
        else:
            Kb = np.block([ks[0][0], ks[0][1]])
        sv = np.linalg.svd(Kb, compute_uv=False)
        out.append((sv[0], sv[-1], int((sv > 1e-9 * sv[0]).sum()),
                    [int((np.linalg.svd(ks[i][i], compute_uv=False) > 1e-9 * max(1e-300, np.linalg.norm(ks[i][i]))).sum()) for i in range(len(ks))]))
    return out


def main():
    seed, nfac, LIMY = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    dlist = [int(c) for c in sys.argv[4]] if len(sys.argv) > 4 else [5, 4, 3, 2, 1]
    mus = [float(x) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 and sys.argv[5] != '-' else [1.0, 0.7, 0.8, 0.6]
    cap = int(sys.argv[6]) if len(sys.argv) > 6 else 3000
    n1 = int(sys.argv[7]) if len(sys.argv) > 7 else 2
    n2 = int(sys.argv[8]) if len(sys.argv) > 8 else 2
    outs = [0, 1] if n2 > 0 else [0]
    rng = random.Random(seed); t0 = time.time()
    cf, L = make_config4(rng, LIMY, n1=n1, n2=n2)
    if cf is None:
        print("seed %d: no 4-branch configuration" % seed); return
    print("seed %d: 4-branch split config %s (%.0fs)" % (seed, cf['info'], time.time() - t0), flush=True)
    injs = factor_injections(rng, cf, L, nfac) if nfac > 0 else []
    print("  factor injections: %s" % [(j['site'], 'x%d' % (j['l'] + 1), j['where'], 'ord%d' % j['tau_order']) for j in injs], flush=True)
    E = Eng(M, seed=seed); Sy = Symp(E)
    P = build_prog4(E, Sy, L, injs, mus)
    # check the factor terms land where designed
    for j in injs:
        k = j['k']; l = j['l']
        Mk = P.mats['c3'] if k == 0 else P.mats['c3x%d' % (k - 1)]
        t2main = [g for (ll, g, c) in P.gens['t2'] if ll == l and c == {0: 1.0}][0]
        hk = P.mulD(Mk, t2main)
        inj_terms = [g for (ll, g, c) in P.gens['t3'] if ll == l and j['alpha'] in c]
        hit = [g for g in inj_terms if Sy.frame(g) == Sy.frame(hk)]
        wm = [E.snap(E.MATS[g] @ np.conj(E.MATS[hk]).T) for g in hit]
        print("   factor injection %d (branch %d, array %d): %d t3 terms in the coset of h_k pi (monomial w: %s)" % (
            j['alpha'], k, l + 1, len(hit), [E.is_mono(x) for x in wm]), flush=True)
    seeds = [seeds4(P, L, cf, l, random.Random(7 + l)) for l in (0, 1)]
    for d in dlist:
        cls = []
        for l in outs:
            cl = closure_multi(P, l, d, seeds[l], cap=cap)
            if cl is None:
                print("  d=%d: closure cap hit" % d, flush=True); cls = None; break
            cls.append(cl)
        if cls is None:
            break
        for li, cl in enumerate(cls):
            diag(P, cl, li)
        Nss = [assemble(cl, P.nal, rng) for cl in cls]
        cost, al = lm(Nss, P.nal, rng, starts=4, iters=80)
        print("  d=%d: seeds %s, closures %s, residual %.3e (%.0fs)" % (
            d, [len(s) for s in seeds], [(len(c[0]), len(c[1]), len(c[3][0])) for c in cls], cost, time.time() - t0), flush=True)
        if cost > 1e-12:
            r, _ = vp_residual(Nss, al if al is not None else np.zeros(0), False)
            n0 = Nss[0][0].shape[0]
            print("     per output: %.3e %.3e" % (float(np.vdot(r[:n0], r[:n0]).real), float(np.vdot(r[n0:], r[n0:]).real)), flush=True)
            print("STOP at d=%d" % d, flush=True)
            break
        kt = kappa_test(P, cls, Nss, al if al is not None else np.zeros(0))
        print("     kappa test (natural rep, random solutions): %s" % ["smax %.2e smin %.2e rank %d, diag ranks %s" % x for x in kt], flush=True)
    print("done (%.0fs)" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
