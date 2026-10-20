"""fam.py -- w3-asym: generate asymmetric block-4 skeletons ("shared move early, two private layers after it, second
routes into every target") with letters chosen for (joint) progress, then minimise the f = 2 diagonal-block residual over
the coefficients (BFGS, several starts); survivors (< 1e-8) get the f = 3 diagonal block and the full f = 2 test.
Usage: fam.py family W nskel starts seed outfile"""
import sys, os, json, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from asym import *
from scipy.optimize import minimize

def prog_letter(rng, gids, target_d=None, tries=20000, letter=rand_engine):
    """letter b with dist(D1 b g) = dist(g) - 1 for every g in gids (joint progress)"""
    for _ in range(tries):
        b = letter(rng)
        ok = True
        for g in gids:
            h = G.intern(G.mul(G.D1, b, G.MATS[g]))
            if G.dist(h) != G.dist(g) - 1:
                ok = False; break
        if ok:
            return b
    return None

def main_terms(P, i, ls=None):
    """the terms of generator i of maximal height (smallest dist) per array"""
    out = {}
    for (l, g) in P.terms[i]:
        if ls is not None and l not in ls: continue
        d = G.dist(g)
        if l not in out or d < out[l][0]:
            out[l] = (d, [g])
        elif d == out[l][0]:
            out[l][1].append(g)
    return [g for (d, gs) in out.values() for g in gs]

def build(fam, rng, W=2):
    P = Prog(W)
    starts = []
    for l in range(W):
        b = prog_letter(rng, [G.ID]); starts.append(P.add_move([(l, b)]))
    inj = lambda k: [(rng.choice(list(range(W)) + starts), rand_engine(rng)) for _ in range(k)]
    if fam == 'A':      # W=2: shared u4 at layer 2; ends u5,u6 / u7,u8 joint-progress both arrays; injections at layer 3
        b0 = prog_letter(rng, main_terms(P, starts[0])); b1 = prog_letter(rng, main_terms(P, starts[1]))
        sh = P.add_move([(starts[0], b0), (starts[1], b1)])
        ends = []
        for e in range(W):
            c = prog_letter(rng, main_terms(P, sh))
            if c is None: return None
            m1 = P.add_move([(sh, c)] + inj(2))
            c2 = prog_letter(rng, main_terms(P, m1))
            if c2 is None: return None
            m2 = P.add_move([(m1, c2)] + inj(1))
            ends.append(m2)
    elif fam == 'B':    # injections INTO the shared move (x's), ends with cross links at the top
        b0 = prog_letter(rng, main_terms(P, starts[0])); b1 = prog_letter(rng, main_terms(P, starts[1]))
        sh = P.add_move([(starts[0], b0), (starts[1], b1)] + [(l, rand_engine(rng)) for l in range(W)])
        mids = []
        for e in range(W):
            c = prog_letter(rng, main_terms(P, sh))
            if c is None: return None
            mids.append(P.add_move([(sh, c)] + inj(1)))
        for e in range(W):
            c2 = prog_letter(rng, main_terms(P, mids[e]))
            if c2 is None: return None
            other = mids[(e + 1) % W]
            P.add_move([(mids[e], c2), (other, rand_engine(rng))] + inj(1))
    elif fam == 'C':    # shared move at layer 3 fed by two layer-2 privates (late merge, then the two top moves)
        mids = []
        for l in range(W):
            c = prog_letter(rng, main_terms(P, starts[l]))
            mids.append(P.add_move([(starts[l], c)] + inj(1)))
        sh = P.add_move([(mids[0], prog_letter(rng, main_terms(P, mids[0]))), (mids[1], prog_letter(rng, main_terms(P, mids[1])))] + inj(1))
        for e in range(W):
            c2 = prog_letter(rng, main_terms(P, sh))
            if c2 is None: return None
            P.add_move([(sh, c2), (mids[e], rand_engine(rng))] + inj(1))
    return P

def d2_obj(P, fs=(2,)):
    DT = DTest(P, fs=fs)
    k = P.np
    def f(x):
        return DT.residual(x[:k] + 1j * x[k:])
    return f, DT

if __name__ == '__main__':
    fam, W, nsk, starts, seed, outf = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    rng = random.Random(seed); nrng = np.random.default_rng(seed)
    res = []
    for s in range(nsk):
        t0 = time.time()
        P = build(fam, rng, W)
        if P is None:
            print(f"skel {s}: build failed", flush=True); continue
        f, DT = d2_obj(P)
        best = (np.inf, None)
        for st in range(starts):
            x0 = nrng.standard_normal(2 * P.np)
            r = minimize(f, x0, method='BFGS', options={'maxiter': 400, 'gtol': 1e-12})
            if r.fun < best[0]: best = (r.fun, r.x)
        x = best[1]; par = x[:P.np] + 1j * x[P.np:]
        line = f"skel {s}: moves {len(P.moves)} params {P.np} D2 min {best[0]:.5f} | {P.describe()} ({time.time()-t0:.0f}s)"
        if best[0] < 1e-6:
            D3 = DTest(P, fs=(3,)).residual(par, detail=True)
            F2 = [r[0] for r in full2(P, par)]
            line += f"\n   SURVIVOR: D3 {D3}  full f=2 {F2}"
            pickle.dump({'P': P.moves, 'par': par}, open(f"{outf}_surv{s}.pkl", 'wb'))
        print(line, flush=True)
        res.append(best[0])
    print("summary: min", min(res) if res else None, "sorted", sorted(np.round(res, 4)))
