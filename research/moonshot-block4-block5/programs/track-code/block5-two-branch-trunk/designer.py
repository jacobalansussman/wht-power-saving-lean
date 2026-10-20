"""designer.py (w3-b5): EXACT greedy collision closure for the 2-in-7 butterfly at block 5 with a TWO-BRANCH trunk
t3 = D(c3 t2 + mu c3b t2) (the door left open by Theorem W3: a second layer-4 route).

1. lam = (3,3) special-position configuration (main branch N2), letters exactly from F_2 data (endsolve.py).
2. Second branch: c3b with F(D c3b D c2) = N2' (a neighbour of N1 with d(M0, N2') = 3), chosen so that the class letters'
   branch-2 images psi^-1(A_l) are at layer <= 4 for both arrays or equal A_l (otherwise the class columns are dead).
3. Rounds: closure (endsolve.closure, layer >= d) for both outputs; every row with a SINGLE producer gets a second producer
   by an exact injection of x_l at t1 (or t2): its stream term is placed in the row's coset (z-route) or at its predecessor
   frame (y-route), F_2 transport tau with tau(M0) = anchor_l (tau'(L0) = anchor_l).  Then ONE layered solve with all
   injections at once (variable projection over their coefficients).
usage: designer.py seed rounds maxinj dlist mu"""
import sys, time, random
import numpy as np
import endsolve as es
LIM = 3
LOOK = 30
SECOND = True
from endsolve import *


def build_prog2(E, Sy, L, injs, mu):
    Lm = {k: v for k, v in L.items() if k != 'c3b'}
    P = build_prog(E, Sy, Lm, injs)
    if 'c3b' in L:
        Mb = lift_P(E, L['c3b'])
        assert Sy.of(E.intern(E.snap(Mb))) == L['c3b']
        t2 = P.gens['t2']
        extra = [(l, P.mulD(Mb, g), {k: mu * v for k, v in c.items()}) for (l, g, c) in t2]
        P.gens['t3'] = merge(P.gens['t3'] + extra)
    return P


def theta(Sy, L, l, gid):
    return smul(pi_bar(L, l), sinv(Sy.of(gid)))


def pick_branch(rng, cf, L, tries=4000, per=12):
    stats = [0, 0]
    M0, N1, N2 = cf['M0'], cf['N1'], cf['N2']
    H = smul(sinv(L['c2']), smul(DS, smul(sinv(L['c3']), DS)))          # Hbar, H L0 = N2
    assert image(H, L0) == N2
    cands = [N for N in sorted(neighbors(N1)) if N not in (N2, L0) and dist(M0, N) == 3]
    rng.shuffle(cands)
    for N2p in cands:
        X = image(DS, image(L['c2'], N2p))                                 # c3b X = M0
        for _ in range(per):
            c3b = find_P(rng, X, M0)
            Hb = smul(sinv(L['c2']), smul(DS, smul(sinv(c3b), DS)))
            if image(Hb, L0) != N2p:
                continue
            psi_inv = smul(Hb, sinv(H))
            stats[0] += 1
            ok = True
            for A, own in ((cf['A1'], 0), (cf['A2'], 1)):
                Ai = image(psi_inv, A)
                if Ai == A:
                    continue
                if dist(cf['P1'], Ai) > LIM or dist(cf['Q'], Ai) > LIM:
                    ok = False
            if ok and SECOND:
                psi = smul(H, sinv(Hb)); P1, Q = cf['P1'], cf['Q']
                cnt = []
                for anc in (Q, P1):
                    g = 0
                    for Y in neighbors(N2p):
                        if dist(M0, Y) != dist(anc, N2):
                            continue
                        Z = image(psi, Y)
                        if max(dist(P1, Y), dist(Q, Y)) <= 3 and max(dist(P1, Z), dist(Q, Z)) <= 3:
                            g += 1
                    cnt.append(g)
                if not (cnt[0] and cnt[1]):
                    ok = False
            if ok:
                cf['Hs'] = {'t3': H, 't3b': Hb}
                return c3b, N2p, (dist(cf['P1'], N2p), dist(cf['Q'], N2p), 'branches tried %d' % stats[0])
    return None, None, stats


def producers(cl):
    rows, cols, ent, kl = cl
    prod = {}
    for (r, c, cf) in ent:
        prod.setdefault(r, set()).add(c)
    inv = {v: k for k, v in rows.items()}
    return [(inv[r], cs) for r, cs in prod.items()]


def design_injection(rng, cf, L, Sy, row, route_pref):
    """an exact injection whose stream gives `row` a new producer; None if impossible"""
    l, gid = row
    E = Sy.E
    th = theta(Sy, L, l, gid)
    G = image(th, L0); Pr = image(th, cf['M0'])
    anchor = [cf['P1'], cf['Q']][l]
    lay_r = dist(anchor, G)
    pib = pi_bar(L, l); pinv = sinv(pib)
    opts = []
    # t1 injection: stream at tau(L0) [t1, layer d(M0,L0)=1], tau(N1) [t2, 2], tau(N2), tau(N2') [t3, 3]
    for (N, where) in ((L0, 't1'), (cf['N1'], 't2'), (cf['N2'], 't3'), (cf.get('N2p'), 't3b')):
        if N is None:
            continue
        dN = dist(cf['M0'], N)
        for (X, route) in ((G, 'z'), (Pr, 'y')):
            if dist(anchor, X) == dN:
                opts.append(('t1', N, X, route, where))
    # t2 injection: tau'(L0) = anchor; stream at tau'(N1) [t2, 1], tau'(N2), tau'(N2') [t3, 2]
    for (N, where) in ((cf['N1'], 't2'), (cf['N2'], 't3'), (cf.get('N2p'), 't3b')):
        if N is None:
            continue
        dN = dist(L0, N)
        for (X, route) in ((G, 'z'), (Pr, 'y')):
            if dist(anchor, X) == dN:
                opts.append(('t2', N, X, route, where))
    if not opts:
        return None, (lay_r, dist(anchor, Pr))
    opts.sort(key=lambda o: (o[3] != route_pref, rng.random()))
    Hs = cf['Hs']
    best = None
    for (site, N, X, route, where) in opts[:4]:
        tgt = image(DS, N) if site == 't1' else image(L['c2'], N)
        wb = where if where in Hs else None
        for smp in range(LOOK):
            b = find_P(rng, image(pinv, X), tgt, tries=300000)
            if b is None:
                break
            tau = smul(pib, smul(sinv(b), DS)) if site == 't1' else smul(pib, smul(sinv(b), L['c2']))
            score = 0
            if wb is not None:                      # lookahead on the main-stream rows the new column brings
                ti = sinv(tau)
                base = image(ti, Pr if route == 'z' else G)
                for k, Hk in Hs.items():
                    F = image(smul(Hk, sinv(Hs[wb])), base)
                    dmax = max(dist(cf['P1'], F), dist(cf['Q'], F))
                    score += max(0, dmax - 3)
            if best is None or score < best[0]:
                best = (score, site, route, where, b, tau)
            if score == 0:
                break
        if best is not None and best[0] == 0:
            break
    if best is None:
        return None, 'nofind'
    score, site, route, where, b, tau = best
    return dict(site=site, l=l, where=where + '/' + route + '/la%d' % score, Y=None, letter=b, tau_order=sorder(tau), lam=lay_r), None


def main():
    seed, rounds, maxinj = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    dlist = [int(c) for c in sys.argv[4]] if len(sys.argv) > 4 else [4, 3, 2, 1]
    mu = float(sys.argv[5]) if len(sys.argv) > 5 else 0.7
    global LIM
    LIM = int(sys.argv[6]) if len(sys.argv) > 6 else 3
    rng = random.Random(seed); t0 = time.time()
    cf = make_config(rng, 3, 3)
    L = build_letters(rng, cf)
    c3b, N2p, dd = pick_branch(rng, cf, L)
    if c3b is None:
        print("seed %d: no admissible second branch (LIM %d), tried %s" % (seed, LIM, dd)); return
    L['c3b'] = c3b; cf['N2p'] = N2p
    print("seed %d: lam (3,3); second branch N2' with (d(P1,N2'), d(Q,N2')) = %s (%.0fs)" % (seed, dd, time.time() - t0), flush=True)
    E = Eng(M, seed=seed); Sy = Symp(E)
    injs = []
    for d in dlist:
        for rd in range(rounds + 1):
            P = build_prog2(E, Sy, L, injs, mu)
            b0s = [clean_letter(P, L, l) for l in (0, 1)]
            cls = []
            for l in (0, 1):
                cl = closure(P, l, d, b0s[l], cap=2500)
                if cl is None:
                    print("  d=%d round %d: closure cap hit with %d injections" % (d, rd, len(injs)), flush=True)
                    cls = None; break
                cls.append(cl)
            if cls is None:
                break
            singles = []
            for cl in cls:
                for (row, cs) in producers(cl):
                    if len(cs) == 1:
                        singles.append(row)
            Nss = [assemble(cl, P.nal, rng) for cl in cls]
            cost, al = lm(Nss, P.nal, rng, starts=2, iters=40)
            print("  d=%d round %d: injections %d, closures %s, single-producer rows %d, residual %.3e (%.0fs)" % (
                d, rd, len(injs), [(len(c[0]), len(c[1])) for c in cls], len(singles), cost, time.time() - t0), flush=True)
            if cost < 1e-20:
                break
            if rd == rounds or len(injs) >= maxinj:
                break
            added = 0; why = {}
            seen = set()
            for row in singles:
                if row in seen:
                    continue
                seen.add(row)
                j, w = design_injection(rng, cf, L, Sy, row, 'y' if rng.random() < 0.5 else 'z')
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
        else:
            pass
        if cls is None or cost > 1e-12:
            print("STOP at d=%d (residual %.3e)" % (d, cost if cls is not None else -1), flush=True)
            break
    print("done (%.0fs)" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
