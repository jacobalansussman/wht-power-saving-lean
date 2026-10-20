"""exact_small.py -- w3-asym: exact (Groebner, Singular) decision of SMALL skeletons around the (2,2) core:
(a) the bare core (merge coefficient free); (b) core + one landing injection of x_1 into u read through v_0 (output 0's
junk frame); (c) core + one landing injection of x_0 into u read through v_1 (output 1's).  Each with its own closures
(cap 3000) and a Singular script; both outputs' equations together.  Usage: exact_small.py nskel seed nsamp outdir"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_decide import *

if __name__ == '__main__':
    nsk, seed, nsamp, outdir = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)
    W = 2
    for s in range(nsk):
        rec = core22(rng, W)
        if rec is None: print("core failed"); continue
        P0 = build_rec(W, rec)
        u = 2 * W
        pm = {l: G.MATS[g] for (l, g) in P0.terms[u] if G.dist(g) == 2}
        ells = [rec[W + 1][0][1], rec[W + 3][0][1]]
        variants = [('bare', [list(m) for m in rec])]
        allr = [list(m) for m in rec]
        for (l, e) in [(0, 0), (1, 1), (0, 1), (1, 0)]:
            j = 1 - l
            J = G.mul(G.K, np.conj(pm[l]).T, pm[j]); fk = G.frame_key(J)
            ks = landing_letters(rng, fk, G.mul(G.D1, ells[e]), np.eye(16, dtype=complex), nsamp)
            if ks:
                r2 = [list(m) for m in rec]; r2[W].append((j, ks[0]))
                variants.append((f"x{j}->u_lands_out{l}_via_v{e}", r2))
                allr[W].append((j, ks[0]))
                variants.append((f"cumulative_{len(allr[W]) - W}", [list(m) for m in allr]))
        for name, r2 in variants:
            P = build_rec(W, r2)
            scr, sizes = singular_script(P, cap=3000)
            if scr is None:
                print(f"skel {s} [{name}]: closure cap", flush=True); continue
            fn = os.path.join(outdir, f"sk{seed}_{s}_{name.replace('->','')}.sing")
            open(fn, 'w').write(scr)
            print(f"skel {s} [{name}]: params {P.np}, closures {sizes}, script {fn}", flush=True)
