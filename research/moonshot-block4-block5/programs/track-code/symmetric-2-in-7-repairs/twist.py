"""twist.py -- mw2-absorb4: ABSORBED shared moves from TWISTED orbit sums.
For H = <(c, pi)> cyclic of order W acting by tau (symlib), sigma(x) = c x c^-1 on B.  A content
    y = sum_j g_j tau^j(b u),   g_0 = 1,  g_(j+1) = beta^-1 sigma(g_j)
satisfies tau(y) = beta y exactly when the norm  g_W = beta^-1 sigma(beta^-1) ... sigma^(W-1)(beta^-1)  is 1 (a scalar is
absorbed into the coefficients).  With a tau-fixed direction D_delta (c D c^-1 = beta1 D) the move u = D y has
tau(u) = beta1 (D beta D^-1) u: one paid move whose left factor beta_tot = beta1 D beta D^-1 is NON-scalar when beta is.
The frame-sum relaxation cannot see beta (augmentation sends it to 1); the full problem can.
Kind 'T' in TSkel.add_move: summ entries (src, b), plus self.twist[t] = beta for that move."""
import sys
from relax_search import *

def mono_inv(G, M):
    return np.conj(M).T

def is_scalar(M):
    return np.allclose(M, M[0, 0] * np.eye(M.shape[0])) and abs(M[0, 0]) > 1e-9

def order_of(G, M, cap=64):
    X = M.copy(); o = 1
    while not is_scalar(X) and o < cap:
        X = G.mul(X, M); o += 1
    return o

class TSkel(Skel):
    def gammas(self, beta):
        """g_0..g_(W-1) as (matrix, scalar); None if the norm is not a scalar"""
        G, R = self.G, self.R
        c = R.Hc[1][0]; ci = np.conj(c).T
        bi = np.conj(beta).T
        gs = [np.eye(G.n, dtype=complex)]
        for j in range(1, len(R.Hc) + 1):
            gs.append(G.mul(bi, c, gs[-1], ci))
        gW = gs[-1]
        if not is_scalar(gW):
            return None
        lam = gW[0, 0]
        mu = lam ** (1.0 / len(R.Hc))          # rescale: g_j -> g_j mu^-j makes the norm 1
        return [(gs[j], mu ** (-j)) for j in range(len(R.Hc))]

    def add_move(self, kind, summ, Dm_override=None, beta=None):
        if kind != 'T':
            return Skel.add_move(self, kind, summ, Dm_override)
        G, R = self.G, self.R
        gam = self.gammas(beta)
        assert gam is not None, "twist norm not scalar"
        Dm = G.D1 if Dm_override is None else Dm_override
        terms = {}; order = []; maps = []
        for k, (src, b) in enumerate(summ):
            a_idx = self.nalpha + k
            for hj, h in enumerate(R.Hc):
                gmat, gsc = gam[hj]
                si, ri = [], []
                for j, (l, g) in enumerate(self.gens[src]):
                    g1 = G.intern(G.mul(b, G.MATS[g]))
                    l2, g2 = R.tau(h, {(l, g1): 1}).popitem()[0] if hj > 0 else (l, g1)
                    g3 = G.intern(G.mul(Dm, gmat, G.MATS[g2]))
                    key = (l2, g3)
                    if key not in terms:
                        terms[key] = len(order); order.append(key)
                    si.append(j); ri.append(terms[key])
                maps.append((a_idx, src, np.array(si), np.array(ri), gsc))
        self.fixed.append(self.nalpha)
        self.nalpha += len(summ)
        self.gens.append(order); self.coef_of.append(('move', len(self.moves)))
        mv = {"kind": kind, "maps": maps, "rep": len(self.gens) - 1, "n": len(order), "images": [], "beta": beta}
        self.moves.append(mv)
        if not hasattr(self, 'testcv'):
            self.testcv = [np.ones(1, dtype=complex) for _ in range(self.W)]
            self.testal = {}
        v = np.zeros(len(order), dtype=complex)
        for (a, src, si, ri, ch) in maps:
            if a not in self.testal:
                self.testal[a] = complex(self.rng.gauss(0, 1), self.rng.gauss(0, 1))
            np.add.at(v, ri, ch * self.testal[a] * self.testcv[src][si])
        self.testcv.append(v)
        test = {t: v[j] for j, t in enumerate(order)}
        mv["absorbed"] = all(self.R.S.same_up_to_free(R.tau(h, test), test) for h in R.Hc[1:])
        return 1

def build_T(m, W, H, recipe, seed, twists, R=None):
    """recipe entries (kind, summ, Dm); twists: {move index: beta} for kind 'T'"""
    sk = TSkel(m, W, H, seed, R)
    for t, (kind, summ, Dm) in enumerate(recipe):
        if kind == 'T':
            sk.add_move('T', summ, Dm, beta=twists[t])
        else:
            sk.add_move(kind, summ, Dm)
    return sk

def beta_candidates(G, R, Dm, rng, nrand=3000):
    """free letters beta with scalar norm and D beta D^-1 free; returns list of (beta, order of beta_tot, beta_tot)"""
    m, n = G.m, G.n
    c = R.Hc[1][0]
    Dmi = np.conj(Dm).T
    beta1 = G.mul(c, Dm, np.conj(c).T, Dmi)
    cands = []
    for a in range(n):                      # Paulis X^a Z^b
        for b in range(n):
            sig = [x ^ a for x in range(n)]; q = [2 * (bin(x & b).count('1') & 1) for x in range(n)]
            cands.append(G.mono(sig, q))
    for _ in range(nrand):
        cands.append(rand_free_b(G, rng, m))
    out = []; seen = set()
    S = TSkel.__new__(TSkel); S.G = G; S.R = R
    for bt in cands:
        if is_scalar(bt):
            continue
        if S.gammas(bt) is None:
            continue
        dbd = G.mul(Dm, bt, Dmi)
        if not G.is_monomial(dbd):
            continue
        btot = G.mul(beta1, dbd)
        if is_scalar(btot):
            continue
        k = G.key(bt)
        if k in seen:
            continue
        seen.add(k)
        out.append((bt, order_of(G, btot), btot))
    return out
