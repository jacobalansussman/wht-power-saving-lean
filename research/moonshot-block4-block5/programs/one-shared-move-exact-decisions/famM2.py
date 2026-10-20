"""famM2.py -- w3-asym: family M (multi-term ends, BRIEF3 shape) with injection-placement variants and a faster
least-squares (LM) minimisation of the f = 2 diagonal-block residual vector over all coefficients.
Shape (W in 3W+1): starts s_l = D b_l x_l; shared u = D(sum_l a_l s_l [+ inputs if inj=='shared']); end e:
v_e = D(sum_l ell_{e,l} u + injections), t_e = D(sum_l m_{e,l} v_e [+ cross v_e' if cross]) with D m D ell p_l in the K
frame (geodesic through midpoint e mod #midpoints of array l): every array has a piece in every end.
inj: 'cross' (x_l -> v_(l+1)), 'own' (x_l -> v_l), 'all' (every x -> every v), 'shared' (every x -> u, plus 'cross').
sig: if 1, end e's letters for array l are end 0's letters times a random order-3 Levi letter on the right of ell
(the relative letter between the two pieces of a target is then non-central of odd order when the frames allow it).
Usage: famM2.py W nskel starts seed kinj inj cross outprefix"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from famM import *
from scipy.optimize import least_squares

def buildM2(rng, W, kinj, inj, cross):
    P = Prog(W)
    starts = [P.add_move([(l, prog_letter(rng, [G.ID]))]) for l in range(W)]
    summ = [(starts[l], prog_letter(rng, main_terms(P, starts[l]))) for l in range(W)]
    if inj == 'shared':
        summ += [(l, rand_engine(rng)) for l in range(W) for _ in range(kinj)]
    sh = P.add_move(summ)
    pl = {}
    for (l, g) in P.terms[sh]:
        if G.dist(g) == 2 and l not in pl: pl[l] = g
    if len(pl) < W: return None, None
    paths = {l: midpaths(rng, pl[l]) for l in range(W)}
    if any(len(v) == 0 for v in paths.values()): return None, None
    info = {l: len(paths[l]) for l in range(W)}
    if inj in ('cross', 'shared'): inj_of = {e: [(e - 1) % W] for e in range(W)}
    elif inj == 'own': inj_of = {e: [e] for e in range(W)}
    elif inj == 'all': inj_of = {e: list(range(W)) for e in range(W)}
    else: inj_of = {e: [] for e in range(W)}
    vs = []
    for e in range(W):
        summ = [(sh, paths[l][e % len(paths[l])][0]) for l in range(W)]
        for l in inj_of[e]:
            summ += [(l, rand_engine(rng)) for _ in range(kinj)]
        vs.append(P.add_move(summ))
    for e in range(W):
        tsumm = [(vs[e], paths[l][e % len(paths[l])][1]) for l in range(W)]
        if cross:
            tsumm += [(vs[(e + 1) % W], rand_engine(rng))]
        P.add_move(tsumm)
    return P, info

class DVec:
    """f = 2 (and 3) diagonal-block residual as a real vector (all outputs), for least squares"""
    def __init__(self, P, f=2):
        self.P = P; self.f = f; self.dim = 16 ** f; W = P.W
        self.rows = {}
        for T in P.terms:
            for (l, g) in T:
                if g not in self.rows:
                    self.rows[g] = np.array([DTest.tens(G.MATS[g][x], f) for x in range(16)])
        TK = np.array([DTest.tens(G.K[x], f) for x in range(16)])
        self.T = []
        for l0 in range(W):
            Tm = np.zeros((16, self.dim * W), dtype=complex); Tm[:, self.dim * l0:self.dim * (l0 + 1)] = TK
            self.T.append(Tm)
        self.norm = np.sqrt(np.sum(np.abs(TK) ** 2))
    def resid(self, par):
        cv = self.P.coefs(par); d = self.dim; W = self.P.W
        S = []
        for i, T in enumerate(self.P.terms):
            X = np.zeros((16, d * W), dtype=complex)
            for j, (l, g) in enumerate(T):
                X[:, d * l:d * (l + 1)] += cv[i][j] * self.rows[g]
            S.append(X)
        S = np.vstack(S)
        Q, Rr = np.linalg.qr(S.conj().T)          # columns of Q span the row space (generic full rank)
        out = []
        for Tm in self.T:
            Rm = Tm - (Tm @ Q) @ Q.conj().T
            out.append(Rm.ravel() / self.norm)
        v = np.concatenate(out)
        return np.concatenate([v.real, v.imag])
    def per_output(self, par):
        v = self.resid(par); n = v.size // 2; v = v[:n] + 1j * v[n:]
        k = v.size // self.P.W
        return [round(float(np.sum(np.abs(v[i * k:(i + 1) * k]) ** 2)), 6) for i in range(self.P.W)]

def solve(P, starts, nrng, f=2):
    DV = DVec(P, f)
    k = P.np
    fun = lambda x: DV.resid(x[:k] + 1j * x[k:])
    best = (np.inf, None)
    for st in range(starts):
        r = least_squares(fun, nrng.standard_normal(2 * k), method='lm', max_nfev=200 * (2 * k + 1), xtol=1e-14, ftol=1e-14)
        c = 2 * r.cost
        if c < best[0]: best = (c, r.x)
    par = best[1][:k] + 1j * best[1][k:]
    return best[0], par, DV.per_output(par)

if __name__ == '__main__':
    W, nsk, starts, seed, kinj, inj, cross, outp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6], int(sys.argv[7]), sys.argv[8]
    rng = random.Random(seed); nrng = np.random.default_rng(seed)
    allres = []
    for s in range(nsk):
        t0 = time.time()
        P, info = buildM2(rng, W, kinj, inj, cross)
        if P is None:
            print(f"skel {s}: build failed", flush=True); continue
        c, par, per = solve(P, starts, nrng)
        line = f"skel {s}: W {W} inj {inj} k {kinj} cross {cross} params {P.np} mids {info} D2 min {c:.5f} per output {per} ({time.time()-t0:.0f}s)"
        if c < 1e-8:
            D3 = DTest(P, fs=(3,)).residual(par, detail=True)
            F2 = [r[0] for r in full2(P, par)]
            line += f"\n   SURVIVOR: D3 {D3}  full f=2 {F2}"
            pickle.dump({'moves': P.moves, 'par': par}, open(f"{outp}_surv{s}.pkl", 'wb'))
        print(line, flush=True)
        allres.append(c)
    print("summary: sorted D2 minima", sorted(np.round(allres, 4)))
