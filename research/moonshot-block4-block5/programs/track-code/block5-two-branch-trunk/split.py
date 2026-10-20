"""split.py (b5-twobranch): exact-structure layered solve for the SPLIT-TRUNK engine 2 in 7 at block 5.

Program: s_l = D sig_l x_l; t1 = D(a1 s1 + a2 s2 + injections); t2 = D(c2 t1 + t2-injections); t3 = D(c3a + mu c3b) t2;
two terminal ends, contents free (endsolve.py's formulation: D c_e - kappa_e K X in V, V = C[B]{x, s, t1, t2, t3}).
The joint walk is needed only up to N1; branch a (frame Na) carries array 1's class letters (d(P1, Na) = 4, A1 ~ Na),
branch b (frame Nb) carries array 2's (d(Q, Nb) = 4, A2 ~ Nb).  Letters are built EXACTLY from F_2 data (sp2.py) and
lifted to monomial Cliffords (every lift checked).  Closure seeded by EVERY class family of the target (one seed per
branch whose frame is adjacent to the target frame); K-part normalisation (endsolve.assemble).
Injection rounds (optional): designer.design_injection gives each single-producer row an exact second producer.
usage: split.py seed profA profB lam1 lam2 mu rounds maxinj dlist
   profX = two digits (d(P1,N), d(Q,N)), e.g. 42; 0 = any.  lam = 0 any."""
import sys, time, random
import numpy as np
import endsolve as es
from endsolve import Eng, Symp, lift_P, find_P, opp_graph_T, merge, assemble, lm, vp_residual, pi_bar, build_prog
from sp2 import (M, L0, Fstar, D as DS, mul as smul, inv as sinv, image, dist, neighbors, order as sorder)
import designer as dz


def H_of(c2, c3):
    return smul(sinv(c2), smul(DS, smul(sinv(c3), DS)))


LIMY = 4
DOUBLE = False        # allow double class letters (Ya = A2 / Yb = A1: kappa rank 1 unless more families)
DECOUPLE = False      # True: the two branches are SEPARATE moves ua, ub (the 2-in-8 relaxation / target)


def make_config_split(rng, profA, profB, lam1, lam2, nlet=12, tries=400):
    M0 = image(DS, L0)
    N1s = [N for N in sorted(neighbors(L0)) if dist(M0, N) == 2]
    for t in range(tries):
        N1 = rng.choice(N1s)
        brs = [N for N in sorted(neighbors(N1)) if dist(L0, N) == 2 and dist(M0, N) == 3]
        anchors = [X for X in sorted(neighbors(M0)) if dist(X, L0) == 2 and dist(X, N1) == 3]
        for _ in range(30):
            P1, Q = rng.sample(anchors, 2)
            if dist(P1, Q) != 2:
                continue
            pr = {N: (dist(P1, N), dist(Q, N)) for N in brs}
            Nas = [N for N in brs if pr[N][0] == 4 and (profA is None or pr[N] == profA)]
            Nbs = [N for N in brs if pr[N][1] == 4 and (profB is None or pr[N] == profB)]
            rng.shuffle(Nas); rng.shuffle(Nbs)
            for Na in Nas[:3]:
                A1s = [A for A in sorted(neighbors(Na)) if dist(P1, A) == 5 and dist(Q, A) <= 4 and (lam1 is None or dist(Q, A) == lam1)]
                if not A1s:
                    continue
                for Nb in Nbs[:3]:
                    if Nb == Na:
                        continue
                    A2s = [A for A in sorted(neighbors(Nb)) if dist(Q, A) == 5 and dist(P1, A) <= 4 and (lam2 is None or dist(P1, A) == lam2)]
                    if not A2s:
                        continue
                    A1 = rng.choice(A1s); A2 = rng.choice(A2s)
                    c2 = find_P(rng, N1, M0)
                    for _l in range(nlet):
                        c3a = find_P(rng, image(DS, image(c2, Na)), M0)
                        c3b = find_P(rng, image(DS, image(c2, Nb)), M0)
                        Ha, Hb = H_of(c2, c3a), H_of(c2, c3b)
                        assert image(Ha, L0) == Na and image(Hb, L0) == Nb
                        Ya = image(Hb, image(sinv(Ha), A1)); Yb = image(Ha, image(sinv(Hb), A2))
                        ok1 = (dist(P1, Ya) <= LIMY or Ya == A1) and (dist(Q, Ya) <= LIMY or Ya == A2)
                        ok2 = (dist(P1, Yb) <= LIMY or Yb == A1) and (dist(Q, Yb) <= LIMY or Yb == A2)
                        if not DOUBLE and (Ya == A2 or Yb == A1):
                            continue
                        if ok1 and ok2:
                            cf = dict(M0=M0, N1=N1, N2=Na, N2p=Nb, P1=P1, Q=Q, A1=A1, A2=A2, Hs={'t3': Ha, 't3b': Hb},
                                      nN2=sorted(neighbors(Na)),
                                      info=dict(profA=pr[Na], profB=pr[Nb], lam=(dist(Q, A1), dist(P1, A2)),
                                                Ya=(dist(P1, Ya), dist(Q, Ya), Ya == A1, Ya == A2),
                                                Yb=(dist(P1, Yb), dist(Q, Yb), Yb == A1, Yb == A2)))
                            L = dict(c2=c2, c3=c3a, c3b=c3b)
                            L['a1'] = find_P(rng, M0, image(DS, P1)); L['a2'] = find_P(rng, M0, image(DS, Q))
                            for l, (a, A) in enumerate(((L['a1'], A1), (L['a2'], A2))):
                                O = image(DS, image(sinv(a), image(DS, A)))
                                assert dist(O, L0) == M
                                T = opp_graph_T(O)
                                S = [T[i] ^ (1 << i) for i in range(M)]
                                sg = es.unip(S)
                                assert image(sg, Fstar) == O
                                L['s%d' % (l + 1)] = sg
                            return cf, L
    return None, None


def class_seeds(P, L, cf, l, rng):
    """one class letter per branch whose frame is adjacent to A_l (exact), as interned matrices"""
    E, Sy = P.E, P.S
    seeds = []
    for k, key in ((0, 'c3'), (1, 'c3b')):
        Nk = cf['N2'] if k == 0 else cf['N2p']
        A = cf['A1'] if l == 0 else cf['A2']
        if dist(Nk, A) != 1:
            continue
        # the t3 term of array l through branch k: D M_k (t2 main term of array l)
        t2main = [g for (ll, g, c) in P.gens['t2'] if ll == l and c == {0: 1.0}][0]
        g = P.mulD(P.mats[key] if key in P.mats else lift_P(E, L[key]), t2main)
        X = image(Sy.of(g), Fstar)
        b = find_P(rng, X, image(DS, L0))
        Bm = lift_P(E, b)
        r = E.intern(E.snap(E.D @ Bm @ E.MATS[g]))
        assert Sy.frame(r) == Fstar
        seeds.append((('ua' if k == 0 else 'ub') if DECOUPLE else 't3', E.intern(E.snap(Bm))))
    return seeds


def closure_multi(P, target_l, dmin, seeds, cap=3000):
    """endsolve.closure with several seed columns y_{t3, b0}"""
    E, Sy = P.E, P.S
    Dinv = np.conj(E.D).T
    idx = {}
    for u, T in P.gens.items():
        for (l, g, c) in T:
            idx.setdefault((l, Sy.frame(g)), []).append((u, g))
    Fs = [Sy.frame(P.s[0]), Sy.frame(P.s[1])]
    FK = Fstar
    Dinv_cache = {}

    def dinv_frame(gid):
        v = Dinv_cache.get(gid)
        if v is None:
            h = E.intern(E.snap(Dinv @ E.MATS[gid])); v = (h, Sy.frame(h)); Dinv_cache[gid] = v
        return v

    def free(l, g):
        if E.layer(g) == 0:
            return True
        F = Sy.frame(g)
        if F == Fs[l] or F == FK:
            return True
        h, Fh = dinv_frame(g)
        if E.layer(h) == 0 or Fh == Fs[l]:
            return True
        return False

    rows = {}; queue = []; cols = {}; ent = []; klass = []; kterms = {}

    def add_col(typ, u, b):
        key = (typ, u, b)
        if key in cols:
            return True
        c = len(cols); cols[key] = c
        if c >= cap:
            return False
        bm = E.MATS[b]
        pre = E.D @ bm if typ == 'y' else bm
        iscl = False
        for (l2, g2, cf) in P.gens[u]:
            g3 = E.intern(E.snap(pre @ E.MATS[g2]))
            if l2 == target_l and E.layer(g3) == M and Sy.frame(g3) == FK:
                iscl = True
                kterms.setdefault(c, []).append(((l2, g3), cf))
            if E.layer(g3) < dmin:
                continue
            t3 = (l2, g3)
            r = rows.get(t3)
            if r is None:
                if free(l2, g3):
                    continue
                r = len(rows); rows[t3] = r; queue.append(t3)
            ent.append((r, c, cf))
        if iscl:
            klass.append(c)
        return True

    for (u0, b0) in seeds:
        add_col('y', u0, b0)
    while queue:
        (l, g) = queue.pop()
        cand = []
        for (u, gj) in idx.get((l, Sy.frame(g)), ()):
            cand.append(('z', u, E.MATS[g] @ np.conj(E.MATS[gj]).T))
        h, Fh = dinv_frame(g)
        for (u, gj) in idx.get((l, Fh), ()):
            cand.append(('y', u, E.MATS[h] @ np.conj(E.MATS[gj]).T))
        for (typ, u, bm) in cand:
            bm = E.snap(bm)
            if not E.is_mono(bm):
                continue
            if not add_col(typ, u, E.intern(bm)):
                return None
    return rows, cols, ent, (klass, kterms)


PEEL = False
LAST_TRIG = []
ESC = False


def closure_any(P, target_l, dmin, seeds, cap=3000):
    if not PEEL:
        return closure_multi(P, target_l, dmin, seeds, cap)
    from peel import peel_closure
    global LAST_TRIG
    rows, cols, ent, kl, st = peel_closure(P, target_l, dmin, seeds, cap=cap)
    E = P.E
    LAST_TRIG = [w for (w, k) in st['triggers']]
    tr = {}
    for (l, g) in LAST_TRIG[:200]:
        h = E.intern(E.snap(np.conj(E.D).T @ E.MATS[g]))
        key = (l + 1, E.layer(g), P.S.frame(h) == P.S.frame(g))
        tr[key] = tr.get(key, 0) + 1
    print("     peel out %d: kill triggers (array, layer, -) %s" % (target_l + 1, tr), flush=True)
    print("     peel out %d: created %d, killed %d, live %d (class live %d of %d, seeds live %d), rows %d%s" % (
        target_l + 1, st['created'], st['killed'], st['live'], st['class_live'], st['class_total'], st['seeds_live'],
        st['rows'], ' CAP' if st.get('cap') else ''), flush=True)
    if st.get('cap'):
        return None
    return rows, cols, ent, kl


def design_injection_clean(rng, cf, L, Sy, row, d, tries=40):
    """designer.design_injection, but reject injections whose t3 terms pollute the SEED class columns: for every class
    direction Z_e (Theta_own^-1(A_own), both outputs) and every branch k, the row D beta_e g_inj sits at global frame
    tau Theta_k(Z_e); require its layer d(anchor_l, .) < d (so it is not a row of the layer-d test)."""
    Zs = [image(sinv(cf['Hs']['t3']), cf['A1']), image(sinv(cf['Hs']['t3b']), cf['A2'])]
    Ths = [cf['Hs']['t3'], cf['Hs']['t3b']]
    rej = 0
    for _ in range(tries):
        j, w = dz.design_injection(rng, cf, L, Sy, row, 'y' if rng.random() < 0.5 else 'z')
        if j is None:
            return None, w, rej
        l = j['l']; anc = cf['P1'] if l == 0 else cf['Q']
        pib = pi_bar(L, l)
        tau = smul(pib, smul(sinv(j['letter']), DS if j['site'] == 't1' else L['c2']))
        ok = True
        Ns = [cf['N2'], cf['N2p']]
        for Z in Zs:
            for (Th, Nj) in zip(Ths, Ns):
                Y = image(smul(tau, Th), Z); lay = dist(anc, Y)
                if lay < d:
                    continue
                if ESC:                          # Lemma Q escapes through the trunk itself (E2-E4)
                    pred = image(tau, Nj)
                    if pred == cf['N1'] and lay - 1 == dist(anc, cf['N1']):
                        continue                 # E3: t2 main term at the predecessor
                    if any(pred == Ni and dist(anc, Ni) == lay - 1 for Ni in Ns):
                        continue                 # E2: a branch main term at the predecessor
                    if any(Y == Ni and dist(anc, Ni) == lay for Ni in Ns):
                        continue                 # E4: a branch main term at the row's frame (z-producer)
                ok = False
        if ok:
            return j, None, rej
        rej += 1
    return None, 'polluting', rej


def diag(P, cl, li):
    E, Sy = P.E, P.S
    rows_, cols_, ent_, kl_ = cl
    prod = {}
    for (r, c, cf) in ent_:
        prod.setdefault(r, set()).add(c)
    inv_rows = {v: k for k, v in rows_.items()}
    single = [r for r, cs in prod.items() if len(cs) == 1]
    colkind = {c: k[0] + '_' + k[1] for k, c in cols_.items()}
    byc = {}
    for r in single:
        c = next(iter(prod[r])); byc[colkind[c]] = byc.get(colkind[c], 0) + 1
    lays = {}
    for r in single:
        (l, g) = inv_rows[r]; lays[(l + 1, E.layer(g))] = lays.get((l + 1, E.layer(g)), 0) + 1
    print("     out %d: %d rows, %d cols (%d class); single-producer rows %d by producer %s, by (array, layer) %s" % (
        li + 1, len(rows_), len(cols_), len(kl_[0]), len(single), byc, lays), flush=True)
    return [inv_rows[r] for r in single]


def main():
    a = sys.argv[1:]
    seed = int(a[0])
    profA = None if a[1] == '0' else (int(a[1][0]), int(a[1][1]))
    profB = None if a[2] == '0' else (int(a[2][0]), int(a[2][1]))
    lam1 = None if a[3] == '0' else int(a[3]); lam2 = None if a[4] == '0' else int(a[4])
    mu = float(a[5]); rounds = int(a[6]); maxinj = int(a[7])
    dlist = [int(c) for c in a[8]] if len(a) > 8 else [5, 4, 3, 2, 1]
    cap = int(a[9]) if len(a) > 9 else 2500
    global LIMY, DECOUPLE, DOUBLE, PEEL
    LIMY = int(a[10]) if len(a) > 10 else 4
    DECOUPLE = len(a) > 11 and a[11] == '8'
    DOUBLE = len(a) > 12 and a[12] == '1'
    PEEL = len(a) > 13 and a[13] in ('p', 'pe')
    global ESC
    ESC = len(a) > 13 and a[13] == 'pe'
    print("mode: %s" % ("2 in 8 (branches decoupled: ua, ub separate moves)" if DECOUPLE else "2 in 7 split trunk t3 = D(c3a + mu c3b) t2"), flush=True)
    rng = random.Random(seed); t0 = time.time()
    cf, L = make_config_split(rng, profA, profB, lam1, lam2)
    if cf is None:
        print("seed %d: no split configuration found" % seed); return
    print("seed %d: split config %s (%.0fs)" % (seed, cf['info'], time.time() - t0), flush=True)
    E = Eng(M, seed=seed); Sy = Symp(E)
    assert Sy.of(E.intern(E.D)) == DS
    injs = []
    TRIGS = [[], []]
    for d in dlist:
        for rd in range(rounds + 1):
            P = dz.build_prog2(E, Sy, L, injs, mu)
            P.mats = dict(P.mats); P.mats['c3b'] = lift_P(E, L['c3b'])
            if DECOUPLE:
                t2 = P.gens['t2']
                P.gens = {'t1': P.gens['t1'], 't2': t2,
                          'ua': merge([(l, P.mulD(P.mats['c3'], g), c) for (l, g, c) in t2]),
                          'ub': merge([(l, P.mulD(P.mats['c3b'], g), c) for (l, g, c) in t2])}
            seeds = [class_seeds(P, L, cf, l, random.Random(7 + l)) for l in (0, 1)]
            cls = []
            for l in (0, 1):
                cl = closure_any(P, l, d, seeds[l], cap=cap)
                TRIGS[l] = list(LAST_TRIG)
                if cl is None:
                    print("  d=%d round %d: closure cap hit (%d injections)" % (d, rd, len(injs)), flush=True)
                    cls = None; break
                cls.append(cl)
            if cls is None:
                break
            singles = []
            for li, cl in enumerate(cls):
                singles += [(li, r) for r in diag(P, cl, li)]
                if PEEL:
                    singles += [(li, r) for r in TRIGS[li]]
            Nss = [assemble(cl, P.nal, rng) for cl in cls]
            cost, al = lm(Nss, P.nal, rng, starts=3, iters=60)
            print("  d=%d round %d: seeds %s, injections %d, closures %s, residual %.3e (%.0fs)" % (
                d, rd, [len(s) for s in seeds], len(injs), [(len(c[0]), len(c[1]), len(c[3][0])) for c in cls], cost,
                time.time() - t0), flush=True)
            if cost < 1e-20:
                break
            if cost > 1e-12:
                r, _ = vp_residual(Nss, al if al is not None else np.zeros(0), False)
                n0 = Nss[0][0].shape[0]
                print("     per output: %.3e %.3e" % (float(np.vdot(r[:n0], r[:n0]).real), float(np.vdot(r[n0:], r[n0:]).real)), flush=True)
            if rd == rounds or len(injs) >= maxinj:
                break
            added = 0; why = {}; seen = set()
            for (li, row) in singles:
                if row in seen:
                    continue
                seen.add(row)
                if PEEL:
                    j, w, nrej = design_injection_clean(rng, cf, L, Sy, row, d)
                else:
                    j, w = dz.design_injection(rng, cf, L, Sy, row, 'y' if rng.random() < 0.5 else 'z')
                if j is None:
                    why[str(w)] = why.get(str(w), 0) + 1; continue
                injs.append(j); added += 1
                if len(injs) >= maxinj:
                    break
            print("     added %d injections %s; rows with no possible producer (layer, pred-layer): %s" % (
                added, [(j['site'], 'x%d' % (j['l'] + 1), j['where'], 'ord%d' % j['tau_order']) for j in injs[-added:]] if added else [],
                why), flush=True)
            if added == 0:
                break
        if cls is None or cost > 1e-12:
            print("STOP at d=%d" % d, flush=True)
            break
    else:
        print("ALL LAYERS PASS: candidate (verify kappa + formal identity)", flush=True)
    print("done (%.0fs)" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
