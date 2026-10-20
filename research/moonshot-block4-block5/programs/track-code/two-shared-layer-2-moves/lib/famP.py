"""famP.py -- w3-asym: family P = the REDUCED form of the BRIEF3 shape (NOTE.md: Lemmas 2, 4 and the reduction):
starts s_l = D b_l x_l; merge u = D(sum_l a_l s_l [+ input injections]) with the CLEAN SPLIT lambda(K p_l^-1 p_j) <= 3 for
all l != j enforced (rejection sampling); end l: v_l = D(ell_l u + inj), t_l = D(m_l v_l + inj) with D m_l D ell_l p_l in
the K frame (single-letter private path w_l = K p_l^-1, which the reduction shows is general up to cross terms).
Injections: menus of k random engine letters from sources: into u: inputs (if injU); into v_l: every input and every
start; into t_l: every input, every start, u.  Coefficients: LM on the f = 2 diagonal-block residual vector.
Usage: famP.py W nskel starts seed k injU outprefix"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from famM2 import *

def lam(M):
    return int(round(np.log2((np.abs(M) > 1e-9).sum(axis=1).max())))

def buildP(rng, W, k, injU, tries=200):
    for _ in range(tries):
        bs = [prog_letter(rng, [G.ID]) for l in range(W)]
        hs = [G.intern(G.mul(G.D1, b)) for b in bs]
        as_ = [prog_letter(rng, [h]) for h in hs]
        ps = [G.mul(G.D1, a, G.MATS[h]) for a, h in zip(as_, hs)]
        ok = all(lam(G.mul(G.K, np.conj(ps[l]).T, ps[j])) <= 3 for l in range(W) for j in range(W) if l != j)
        if ok: break
    else:
        return None
    P = Prog(W)
    starts = [P.add_move([(l, bs[l])]) for l in range(W)]
    summ = [(starts[l], as_[l]) for l in range(W)]
    if injU:
        summ += [(l, rand_engine(rng)) for l in range(W) for _ in range(k)]
    u = P.add_move(summ)
    pl = {}
    for (l, g) in P.terms[u]:
        if G.dist(g) == 2 and l not in pl: pl[l] = g
    vs = []
    for l in range(W):
        path = midpaths(rng, pl[l], nsamp=300)
        if not path: return None
        ell, m = path[0]
        summ = [(u, ell)] + [(src, rand_engine(rng)) for src in list(range(W)) + starts for _ in range(k)]
        v = P.add_move(summ)
        summ2 = [(v, m)] + [(src, rand_engine(rng)) for src in list(range(W)) + starts + [u] for _ in range(k)]
        P.add_move(summ2)
    return P

if __name__ == '__main__':
    W, nsk, starts, seed, k, injU, outp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7]
    rng = random.Random(seed); nrng = np.random.default_rng(seed)
    allres = []
    for s in range(nsk):
        t0 = time.time()
        P = buildP(rng, W, k, injU)
        if P is None:
            print(f"skel {s}: build failed", flush=True); continue
        c, par, per = solve(P, starts, nrng)
        line = f"skel {s}: W {W} k {k} injU {injU} params {P.np} D2 min {c:.5f} per output {per} ({time.time()-t0:.0f}s)"
        if c < 1e-8:
            D3 = DTest(P, fs=(3,)).residual(par, detail=True)
            F2 = [r[0] for r in full2(P, par)]
            line += f"\n   SURVIVOR: D3 {D3}  full f=2 {F2}"
            pickle.dump({'moves': P.moves, 'par': par}, open(f"{outp}_surv{s}.pkl", 'wb'))
        print(line, flush=True)
        allres.append(c)
    print("summary: sorted D2 minima", sorted(np.round(allres, 4)))
