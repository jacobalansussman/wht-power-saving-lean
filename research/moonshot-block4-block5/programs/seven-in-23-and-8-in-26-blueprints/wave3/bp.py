"""bp.py -- w3-blueprint (moonshot wave 3, 14 Oct 2026): the TWO-ROUTE palindromic blueprint shape, exact structure.

Shape F O O O F under H = Z7rot (7 in 23) or Z2cube (8 in 26), block 4, engine game (final test = bm closure, free =
Z4-monomials, i.e. the widened game, which contains the engine game):
  F1 = d_x . sum_h tau_h(b0 x_0)                      (shared x-move on the inputs; d_x a tau-fixed progressing direction)
  O1_0 = D1 . (b1 x_0 + a1 F1)                          (own input = route B; F1 = route A)
  O2_0 = D1 . (c1 O1_0 + c2 O1_0 [+ extra])           (c1 progresses route A, c2 route B)
  O3_0 = D1 . (e1 O2_0 + e2 O2_0 [+ extra])
  F2 = d_x' . sum_h tau_h(f1 O3_0 [+ f2 O3_0])         (shared move at the end: route B's x-step)
Route A: x_k -> F1 (x first) -> O1 -> O2 -> O3 reaches K.  Route B: x_k -> O1 -> O2 -> O3 -> F2 (x last) reaches K.
Every array has two routes into its kernel frame, related through tau (odd order for Z7): the sigma-type two-piece
separation that mw2-absorb4 section 5 asks for, with the shared moves at both ends (D2 and its mirror hold: O1 has
co-depth 4, O3 has depth 4).  Variants: 'A' (as above), 'Bx' (extra summands from other orbit members / F1 into O2, O3).
Usage: bp.py W H variant nseeds seed0 starts cap
"""
import sys, time
from symsolve import *

def lead_of(sk, gi, l, want_dist=None):
    G = sk.G
    out = [(G.dist(g), g) for (ll, g) in sk.gens[gi] if ll == l]
    if want_dist is not None:
        out = [t for t in out if t[0] == want_dist]
    return out

def prog(G, gid, rng, Dmat=None, want=6):
    try:
        return G.progress_letters(gid, rng, want=want, Dmat=Dmat, steps=1)
    except Exception:
        return []

def build(W, H, seed, variant='A', log=print):
    R = Run(4, W, H, seed)
    sk = Skel(4, W, H, seed, R)
    G, rng = sk.G, sk.rng
    fd = fixed_dirs(R, rng)
    dfix = [d for d in fd if prog(G, G.ID, rng, Dmat=d)]
    if not dfix:
        log("no progressing fixed direction"); return None
    dx = dfix[0]
    # F1
    pl = prog(G, G.ID, rng, Dmat=dx)
    b0 = rng.choice(pl)
    if 'Tw' in variant:
        # SELECTIVE shared start (NOTE.md section 5b): pick b0 with b0^-1 c b0 NOT in K^-1 B K, so that along every path
        # exactly one branch of the orbit sum reaches the kernel frame; keep the progress condition dist(dx b0) = 3
        c = R.Hc[1][0]
        Kc = np.conj(G.K).T
        found = None
        for _ in range(4000):
            b = rand_free_b(G, rng, 4)
            if G.dist(G.intern(G.mul(dx, b))) != 3:
                continue
            if G.is_monomial(G.mul(G.K, np.conj(b).T, c, b, Kc)):
                continue
            found = b; break
        if found is None:
            log("no twisting progress letter b0"); return None
        b0 = found
        log("   shared start: b0 with b0^-1 c b0 outside K^-1 B K (selective)")
    k1 = 'A' if 'Ach' in variant else 'S'          # 'A': chi-twisted orbit sum (eigencharacter omega)
    summ1 = [(0, b0)]
    if 'M2' in variant:                              # two raw-input terms in the shared start (the loophole)
        for _ in range(200):
            bb = rand_free_b(G, rng, 4)
            if G.dist(G.intern(G.mul(dx, bb))) == 3:
                summ1.append((1 % W, bb)); break
    sk.add_move(k1, summ1, dx); iF1 = len(sk.gens) - 1
    tA = lead_of(sk, iF1, 0)
    assert tA and tA[0][0] == 3, tA
    gA = tA[0][1]
    # O1
    b1 = rng.choice(prog(G, G.ID, rng))
    a1s = prog(G, gA, rng)
    if not a1s:
        log("no progress letter for route A at F1"); return None
    a1 = rng.choice(a1s)
    if 'Sel' in variant:
        a1s = prog(G, gA, rng, want=16)
        othersF = [g for (l, g) in sk.gens[iF1] if l != 0]
        def nprog1(c):
            return sum(1 for g in othersF if G.dist(G.intern(G.mul(G.D1, c, G.MATS[g]))) == G.dist(g) - 1)
        sc = sorted(((nprog1(c), k) for k, c in enumerate(a1s)))
        a1 = a1s[sc[0][1]]
        log(f"   layer 1: selective letter progresses {sc[0][0]} of {len(othersF)} other-array F1 terms (worst {sc[-1][0]})")
    sk.add_move('F', [(0, b1), (iF1, a1)]); iO1 = sk.moves[-1]["rep"]
    # route terms of O1_0
    gB1 = G.intern(G.mul(G.D1, b1))
    gA1 = G.intern(G.mul(G.D1, a1, G.MATS[gA]))
    assert G.dist(gB1) == 3 and G.dist(gA1) == 2, (G.dist(gB1), G.dist(gA1))
    prev, gAp, gBp = iO1, gA1, gB1
    for layer in (2, 3):
        cA = prog(G, gAp, rng); cB = prog(G, gBp, rng)
        if not cA or not cB:
            log(f"no progress letters at layer {layer}"); return None
        c1, c2 = rng.choice(cA), rng.choice(cB)
        if 'Sel' in variant:
            # array-selective route A: among progress letters for array 0's route-A term, take the one that progresses
            # the FEWEST other arrays' terms at the same distance (so route A stops delivering equal coefficients)
            cA = prog(G, gAp, rng, want=16)
            others = [g for (l, g) in sk.gens[prev] if l != 0 and G.dist(g) == G.dist(gAp)]
            def nprog(c):
                return sum(1 for g in others if G.dist(G.intern(G.mul(G.D1, c, G.MATS[g]))) == G.dist(g) - 1)
            scored = sorted(((nprog(c), k) for k, c in enumerate(cA)))
            c1 = cA[scored[0][1]]
            log(f"   layer {layer}: selective letter progresses {scored[0][0]} of {len(others)} other-array terms (worst {scored[-1][0]})")
        summ = [(prev, c1), (prev, c2)]
        if variant.startswith('Bx'):
            # an extra summand from another member of the previous orbit (its route-A term aligned onto ours: collision)
            mv = sk.moves[-1]
            T = G.intern(G.mul(G.D1, c1, G.MATS[gAp]))
            done = False
            for j in mv["images"]:
                for (l, g) in sk.gens[j]:
                    b = G.mul(np.conj(G.D1).T, G.MATS[T], np.conj(G.MATS[g]).T)
                    if G.is_monomial(b) and G.dist(g) >= 1:
                        summ.append((j, b)); done = True; break
                if done:
                    break
        if variant.startswith('Sym') and str(layer) in variant:
            # symmetrise route B: summands e2 c_h^-1 . (member h of the previous orbit), every h; member h's route-B
            # term then lands on the tau_*^h-image cell of ours, so with equal coefficients route B vanishes in every
            # nontrivial character of the frame-sum relaxation (NOTE.md, toy analysis)
            mvp = sk.moves[-1]
            members = [mvp["rep"]] + list(mvp["images"])
            summ = [(prev, c1)]
            for hj, idx in enumerate(members):
                ch = R.Hc[hj][0]
                summ.append((idx, G.mul(c2, np.conj(ch).T)))
        sk.add_move('F', summ); prev = sk.moves[-1]["rep"]
        gAp = G.intern(G.mul(G.D1, c1, G.MATS[gAp]))
        gBp = G.intern(G.mul(G.D1, c2, G.MATS[gBp]))
        assert G.dist(gAp) == 4 - layer - 1 and G.dist(gBp) == 4 - layer, (layer, G.dist(gAp), G.dist(gBp))
    # F2: a tau-fixed direction that takes route B (distance 1) to the kernel frame
    best = None
    for d in fd:
        f1s = prog(G, gBp, rng, Dmat=d)
        if f1s:
            best = (d, rng.choice(f1s)); break
    kind2 = 'S'
    if best is None or variant.endswith('T'):
        # c-twisted orbit sum (tau_h u = c_h u): an absorbed move in ANY direction (w3-blueprint, tcheck.py)
        f1s = prog(G, gBp, rng)
        if not f1s:
            log("no letter finishes route B"); return None
        best = (G.D1, rng.choice(f1s)); kind2 = 'T'
    d2, f1 = best
    sk.add_move(kind2, [(prev, f1)], d2)
    cost = sum(1 + len(mv["images"]) for mv in sk.moves)
    return sk, cost

def frame_relax(sk, alpha_free=None, starts=3):
    G = sk.G
    fidx = {}
    for i, T in enumerate(sk.gens):
        for j, (l, g) in enumerate(T):
            key = (G.frm(g), l)
            if key not in fidx:
                fidx[key] = len(fidx)
    Kkey = (G.frm(G.KID), 0)
    if Kkey not in fidx:
        return None
    nrow = len(fidx)
    rowmap = [np.array([fidx[(G.frm(g), l)] for (l, g) in T]) for T in sk.gens]
    N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
    def relres(al):
        cv = sk.coefs(al)
        A = np.zeros((nrow, len(sk.gens)), dtype=complex)
        for i in range(len(sk.gens)):
            np.add.at(A[:, i], rowmap[i], cv[i])
        t = np.zeros(nrow, dtype=complex); t[fidx[Kkey]] = 1
        beta = np.linalg.lstsq(A, t, rcond=1e-10)[0]
        return A @ beta - t
    best = np.inf
    rng = random.Random(1)
    for st in range(starts):
        x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
        def full(x):
            al = np.ones(N, dtype=complex); al[free_idx] = x[:nf] + 1j * x[nf:]; return al
        if nf:
            x, c = lm(lambda x: relres(full(x)), x0, iters=40)
        else:
            c = float(np.sum(np.abs(relres(np.ones(N, dtype=complex))) ** 2))
        best = min(best, c)
    return best

def main():
    W, H, variant, nseeds, seed0, starts, cap = int(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
    full_only_one = len(sys.argv) > 8 and sys.argv[8] == 'arr0'
    t0 = time.time()
    summary = []
    for s in range(nseeds):
        seed = seed0 + s
        try:
            r = build(W, H, seed, variant)
        except AssertionError as e:
            print(f"seed {seed}: build assertion {e}", flush=True); continue
        if r is None:
            print(f"seed {seed}: build failed", flush=True); continue
        sk, cost = r
        orbs = [1 + len(mv['images']) for mv in sk.moves]
        rel = frame_relax(sk)
        print(f"seed {seed}: cost {cost} orbits {orbs} rep terms {[mv['n'] for mv in sk.moves]} frame-sum relaxation res^2 {rel} ({time.time()-t0:.0f}s)", flush=True)
        if rel is None or rel > 1e-6:
            summary.append((rel if rel is not None else 9, None, seed)); continue
        # full exact-structure test: closure for array 0 only (the others follow by symmetry)
        Wsave = sk.W
        sk.W = 1
        ok = sk.setup_final(cap)
        sk.W = Wsave
        if not ok:
            print(f"seed {seed}: closure cap {cap}", flush=True); summary.append((rel, 'cap', seed)); continue
        prep_fast(sk)
        N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
        dims = [(nr, nc) for (nr, nc, E) in sk.final]
        best = (np.inf, None)
        for st in range(starts):
            x0 = np.array([sk.rng.gauss(0, 1) for _ in range(2 * nf)])
            x, c = vp_lm(sk, N, free_idx, x0)
            if c < best[0]:
                best = (c, x)
        print(f"seed {seed}: FULL closure dims {dims} unknowns {nf} best res^2 {best[0]:.4e} ({time.time()-t0:.0f}s)", flush=True)
        summary.append((rel, best[0], seed))
        if best[0] < 1e-14:
            np.save(f"bp_sol_{W}_{H}_{variant}_{seed}.npy", best[1])
    print("summary (relax, full, seed):", summary, flush=True)

if __name__ == "__main__":
    main()
