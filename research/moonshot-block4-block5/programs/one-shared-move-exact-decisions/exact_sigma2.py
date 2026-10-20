"""exact_sigma2.py -- w3-asym: sigma IN THE MERGE with BOTH routes used: u = D(a_0 s_0 + c a_0 sigma s_0 + a_1 s_1);
array 0's main stream has two terms pi = D a_0 D b_0 and pi' = D a_0 sigma D b_0, so there are two different
kernel-reaching words, w = K pi^-1 and w' = K pi'^-1 (NOTE.md s5): the genuinely new "second route".  End 0 gets
multi-term letters carrying both: v_0 = D(ell u + ell' u), t_0 = D(m v_0 + m' v_0) with D m D ell pi and D m' D ell' pi'
in the kernel frame; end 1 the single path of p_1.  Variants: bare, + landing injections (x_j -> u read through v_e
landing in output l's junk frame of the pi-route).  Writes Singular scripts (they were decided on a Linux machine).
Usage: exact_sigma2.py ncores seed nsamp cap outdir"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_decide import *
from famS import order3_letters

def core_sigma2(rng, W=2, tries=400):
    sigs = order3_letters(rng)
    for _ in range(tries):
        bs = [prog_letter(rng, [G.ID]) for l in range(W)]
        hs = [G.intern(G.mul(G.D1, b)) for b in bs]
        as_ = [prog_letter(rng, [h]) for h in hs]
        ps = [G.mul(G.D1, a, G.MATS[h]) for a, h in zip(as_, hs)]
        if lam(G.mul(G.K, np.conj(ps[0]).T, ps[1])) > 3 or lam(G.mul(G.K, np.conj(ps[1]).T, ps[0])) > 3:
            continue
        rng.shuffle(sigs)
        for sg in sigs:
            q = G.intern(G.mul(G.D1, as_[0], sg, G.MATS[hs[0]]))
            if G.dist(q) == 2 and lam(G.MATS[q]) == 2 and G.frm(q) != G.frm(G.intern(ps[0])):
                pa = midpaths(rng, G.intern(ps[0]), nsamp=300); pb = midpaths(rng, q, nsamp=300)
                p1 = midpaths(rng, G.intern(ps[1]), nsamp=300)
                if pa and pb and p1:
                    rec = [[(l, bs[l])] for l in range(W)]
                    rec.append([(W, as_[0]), (W, G.mul(as_[0], sg)), (W + 1, as_[1])])
                    u = 2 * W
                    rec.append([(u, pa[0][0]), (u, pb[0][0])])          # v_0 (gen 5)
                    rec.append([(u + 1, pa[0][1]), (u + 1, pb[0][1])])  # t_0 (gen 6)
                    rec.append([(u, p1[0][0])])                         # v_1 (gen 7)
                    rec.append([(u + 3, p1[0][1])])                     # t_1 (gen 8)
                    return rec, ps, [pa[0][0], p1[0][0]]
    return None

if __name__ == '__main__':
    ncores, seed, nsamp, cap, outdir = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed); W = 2
    for s in range(ncores):
        r = core_sigma2(rng)
        if r is None: print("core failed", flush=True); continue
        rec, ps, ells = r
        variants = [('bare', [list(m) for m in rec])]
        allr = [list(m) for m in rec]
        for (l, e) in [(0, 0), (1, 1), (0, 1), (1, 0)]:
            j = 1 - l
            J = G.mul(G.K, np.conj(ps[l]).T, ps[j]); fk = G.frame_key(J)
            ks = landing_letters(rng, fk, G.mul(G.D1, ells[e]), np.eye(16, dtype=complex), nsamp)
            if ks:
                r2 = [list(m) for m in rec]; r2[W].append((j, ks[0]))
                variants.append((f"x{j}u_out{l}_v{e}", r2))
                allr[W].append((j, ks[0])); variants.append((f"cum{len(allr[W]) - 3}", [list(m) for m in allr]))
        for name, r2 in variants:
            P = build_rec(W, r2)
            scr, sizes = singular_script(P, cap=cap)
            if scr is None:
                print(f"core {s} [{name}]: closure cap", flush=True); continue
            fn = os.path.join(outdir, f"sig2_{seed}_{s}_{name}.sing")
            open(fn, 'w').write(scr)
            print(f"core {s} [{name}]: params {P.np}, closures {sizes}, {fn}", flush=True)
