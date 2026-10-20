"""famT.py -- w3-asym: theory-guided injections on (2,2) cores (Lemma 5 position), W = 2.
For output l (other array j) the junk J_l = f_l K p_l^-1 p_j sits in the layer-2 frame Phi_l and can only be cancelled by
a LOWER generator's x_j-term in Phi_l (NOTE.md s2.6): v_e's D ell_e D kappa x_j (x_j injected into u with letter kappa) or
v_e's D c D b_j x_j (s_j injected into v_e with letter c).  For every (l, e) and both kinds, letters landing exactly in
Phi_l are found by sampling; all are added with free coefficients, then (a) LM on the f = 2 diagonal residual and (b) LM
on the FULL closure residual; then up to `extra` greedy closure-guided injections (famI) on top.
Usage: famT.py nskel seed nsamp extra outprefix"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famI
from famI import *
from famI22 import core22

def landing_letters(rng, target_fk, pre, post, nsamp, want=1):
    """letters b with frame(pre . D b . post) == target_fk; pre/post matrices"""
    out = []
    for _ in range(nsamp):
        b = rand_engine(rng)
        g = G.mul(pre, G.D1, b, post)
        if G.frame_key(g) == target_fk:
            out.append(b)
            if len(out) >= want: break
    return out

def guided(rng, rec, nsamp):
    W = 2
    P = build_rec(W, rec)
    u = 2 * W; v = [u + 1, u + 3]
    # main elements
    pm = {}
    for (l, g) in P.terms[u]:
        if G.dist(g) == 2: pm[l] = G.MATS[g]
    ells = [rec[W + 1][0][1], rec[W + 3][0][1]]
    bs = [rec[0][0][1], rec[1][0][1]]
    rec2 = [list(m) for m in rec]
    found = []
    for l in range(W):
        j = 1 - l
        J = G.mul(G.K, np.conj(pm[l]).T, pm[j]); fk = G.frame_key(J)
        for e in range(W):
            # x_j into u (move index W) read through v_e: D ell_e D kappa
            pre = G.mul(G.D1, ells[e])
            for kap in landing_letters(rng, fk, pre, np.eye(16, dtype=complex), nsamp):
                rec2[W].append((j, kap)); found.append(('x->u', l, e))
            # s_j into v_e (move index W+1+2e): D c D b_j
            for c in landing_letters(rng, fk, np.eye(16, dtype=complex), G.mul(G.D1, bs[j]), nsamp):
                rec2[W + 1 + 2 * e].append((W + j, c)); found.append(('s->v', l, e))
    return rec2, found

if __name__ == '__main__':
    nsk, seed, nsamp, extra, outp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    rng = random.Random(seed); nrng = np.random.default_rng(seed)
    for s in range(nsk):
        t0 = time.time()
        rec = core22(rng, 2)
        if rec is None: print("core failed"); continue
        rec2, found = guided(rng, rec, nsamp)
        P = build_rec(2, rec2)
        print(f"skel {s}: guided injections {found} params {P.np}", flush=True)
        if P.np == 0: continue
        c2, par2, per2 = solve(P, 3, nrng)
        print(f"   f=2 diagonal LM: {c2:.5f} per output {per2}", flush=True)
        fr = full_res(P, 3, nrng)
        if fr is None: print("   closure cap", flush=True); continue
        print(f"   FULL closure LM: {fr[0]:.5f} per output {np.round(fr[2], 5)} closures {[(C.nr, C.nc) for C in fr[3]]} ({time.time()-t0:.0f}s)", flush=True)
        if fr[0] < 1e-9:
            pickle.dump({'W': 2, 'rec': rec2, 'par': fr[1]}, open(f"{outp}_hit{s}.pkl", 'wb')); print("   HIT (numerical)", flush=True)
