"""famM.py -- w3-asym: family M ("multi-term ends") for W in 3W+1 at block 4, shape of BRIEF3:
  starts s_l = D b_l x_l (progress), ONE shared move u = D(sum_l a_l s_l) at layer 2 (single progress letters), then W ends
  of two moves.  End e: v_e = D(sum_k alpha ell_{e,k} u + injections), t_e = D(sum_k alpha' m_{e,k} v_e + injections),
  where for every array l the end word has a term D m_{e,l} D ell_{e,l} reaching the K frame on array l's main stream
  p_l (a geodesic p_l -> K through a midpoint frame; different ends use different midpoints when possible), so EVERY array
  has a piece in EVERY end (W pieces per target).  Injections: every input x_l into some end's first move (Lemma Delta),
  menu of k random engine letters each.  Optional sigma: end e's letters right-multiplied... (variant 'S': the second end's
  letters for array 1 are the first end's times an order-3 Levi letter).
Objective: f = 2 diagonal-block residual over all outputs, minimised over every coefficient (BFGS, several starts).
Usage: famM.py W nskel starts seed kinj outprefix"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fam import *

def midpaths(rng, g, nsamp=600):
    """all geodesic 2-step paths g (dist 2) -> K: list of (ell, m) with dist(D ell g) = 1 and D m D ell g in the K frame,
    one per midpoint frame"""
    mids = {}
    for _ in range(nsamp):
        l = rand_engine(rng)
        h = G.intern(G.mul(G.D1, l, G.MATS[g]))
        if G.dist(h) == G.dist(g) - 1:
            fk = G.frm(h)
            if fk not in mids:
                mids[fk] = (l, h)
    out = []
    for fk, (l, h) in mids.items():
        m = prog_letter(rng, [h], tries=4000)
        if m is not None:
            out.append((l, m))
    return out

def buildM(rng, W, kinj):
    P = Prog(W)
    starts = [P.add_move([(l, prog_letter(rng, [G.ID]))]) for l in range(W)]
    summ = [(starts[l], prog_letter(rng, main_terms(P, starts[l]))) for l in range(W)]
    sh = P.add_move(summ)
    # main stream element of each array in the shared move
    pl = {}
    for (l, g) in P.terms[sh]:
        if G.dist(g) == 2: pl[l] = g
    if len(pl) < W: return None, None
    paths = {l: midpaths(rng, pl[l]) for l in range(W)}
    if any(len(v) == 0 for v in paths.values()): return None, None
    info = {l: len(paths[l]) for l in range(W)}
    inj_of = {e: [] for e in range(W)}
    for l in range(W):                     # input x_l injected into end (l+1) mod W's first move (another array's end)
        inj_of[(l + 1) % W].append(l)
    tops = []
    for e in range(W):
        ells = [paths[l][e % len(paths[l])][0] for l in range(W)]
        ms = [paths[l][e % len(paths[l])][1] for l in range(W)]
        summ = [(sh, ells[l]) for l in range(W)]
        for l in inj_of[e]:
            summ += [(l, rand_engine(rng)) for _ in range(kinj)]
        v = P.add_move(summ)
        tsumm = [(v, ms[l]) for l in range(W)]
        tops.append(P.add_move(tsumm))
    return P, info

if __name__ == '__main__':
    W, nsk, starts, seed, kinj, outp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    rng = random.Random(seed); nrng = np.random.default_rng(seed)
    allres = []
    for s in range(nsk):
        t0 = time.time()
        P, info = buildM(rng, W, kinj)
        if P is None:
            print(f"skel {s}: build failed", flush=True); continue
        DT = DTest(P, fs=(2,))
        k = P.np
        f = lambda x: DT.residual(x[:k] + 1j * x[k:])
        best = (np.inf, None)
        for st in range(starts):
            r = minimize(f, nrng.standard_normal(2 * k), method='BFGS', options={'maxiter': 300, 'gtol': 1e-12})
            if r.fun < best[0]: best = (r.fun, r.x)
        par = best[1][:k] + 1j * best[1][k:]
        tot, det = DT.residual(par, detail=True)
        line = f"skel {s}: W {W} moves {len(P.moves)} params {k} midpoints {info} D2 min {best[0]:.5f} {det} ({time.time()-t0:.0f}s)"
        if best[0] < 1e-6:
            D3 = DTest(P, fs=(3,)).residual(par, detail=True)
            F2 = [r[0] for r in full2(P, par)]
            line += f"\n   SURVIVOR: D3 {D3}  full f=2 {F2}"
            pickle.dump({'moves': P.moves, 'par': par}, open(f"{outp}_surv{s}.pkl", 'wb'))
        print(line, flush=True)
        allres.append(best[0])
    print("summary: sorted D2 minima", sorted(np.round(allres, 4)))
