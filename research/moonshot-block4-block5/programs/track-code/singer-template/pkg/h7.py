"""h7.py (mw2-singer4): Singer-separation templates at block 4, prefix of M - W moves + W tops tested by tops.py
(all separation matrices T over C[<eta>] at once; also the scalar case H = {1}).
Templates (all respect RULES D2 both ways):
  H   two x-injections u_l = D1 x_l (l = 1..W, co-top moves), then a shared conveyor of depth-2/depth-3 moves whose
      contents are Singer filters  gamma = sum_t s_t eta^{k_t} rho_t  (rho_t routes along the interval, <= 3 terms,
      s_t in {+-1, +-i}), and W tops solved.  Moves: W + (M - 2W) + W.
  C   one shared x-move u1 = D1(x_1 + h x_2 + ...) (Singer pattern h in <eta>), then injections of the other arrays
      (co-depth >= 4), conveyor, tops solved.
  PL  plain control (private chains; must pass with M = 4W).
usage: h7.py TEMPLATE W M n seed [H=eta|one] [out.json]"""
import sys, time, json, random
import numpy as np
from fractions import Fraction as Fr
from s4lib import S4, E_pool, fast_frame_key
from tops import target_space, best_tops

UNITC = [Fr(1), Fr(-1)]

class Builder:
    def __init__(self, S, rng):
        self.S = S; self.G = S.G; self.rng = rng; self.rcache = {}

    def on_geo(self, gid):
        return self.S.depth(gid) + self.S.dK(gid) == 4

    def routes(self, gid, samples=300):
        """free letters b with D1 b g one step further along a geodesic: dict frame -> b"""
        r = self.rcache.get(gid)
        if r is not None:
            return r
        S, G = self.S, self.G
        d0, dk = S.depth(gid), S.dK(gid)
        r = {}
        need = (1 << dk) - 1
        for _ in range(samples):
            if len(r) >= need:
                break
            b = S.rand_free()
            h = G.intern(G.mul(S.D1, b, G.MATS[gid]))
            if S.depth(h) == d0 + 1 and S.dK(h) == dk - 1:
                fk = S.fk(h)
                if fk not in r:
                    r[fk] = b
        self.rcache[gid] = r
        return r

    def lead(self, gen, l=None):
        """deepest geodesic term of a generator (optionally of array l)"""
        S = self.S
        ts = [(l2, g) for (l2, g) in gen if self.on_geo(g) and (l is None or l2 == l)]
        if not ts:
            return None
        ts.sort(key=lambda t: -S.depth(t[1]))
        best = S.depth(ts[0][1])
        top = [t for t in ts if S.depth(t[1]) == best]
        return self.rng.choice(top)

    def singer_filter(self, gen, nterm=None, l=None):
        """list of (coef, monomial) summands acting on gen: sum_t s_t eta^{k_t} rho_t"""
        S, rng = self.S, self.rng
        t = self.lead(gen, l)
        if t is None:
            return None
        rts = list(self.routes(t[1]).values())
        if not rts:
            return None
        k = nterm or rng.choice([1, 2, 2, 3, 3])
        out = []
        for _ in range(k):
            rho = rng.choice(rts)
            e = rng.randrange(len(S.etap))
            mode = rng.random()
            if mode < 0.5:
                b = self.G.mul(S.etap[e], rho)                       # left eta: same frame, C[<mu>] coefficient
            else:
                b = self.G.mul(S.etap[e], rho, np.conj(S.etap[e]).T)  # conjugated route: frame moved by eta^e
                # keep only if it still advances the lead term
                h = self.G.intern(self.G.mul(S.D1, b, self.G.MATS[t[1]]))
                if not (S.depth(h) == S.depth(t[1]) + 1 and S.dK(h) == S.dK(t[1]) - 1):
                    b = self.G.mul(S.etap[e], rho)
            out.append((rng.choice(UNITC), b))
        return out

    def mv(self, parts):
        """parts: list of (gen, filter) -> one move D1 * sum filter(gen)"""
        summ = []
        for gen, flt in parts:
            for (co, b) in flt:
                summ.append((co, b, gen))
        return self.S.move(summ)

def scored_routes(B, gen, samples=250):
    """routes for the geodesic terms of gen, scored by how many geodesic terms of gen stay on geodesics after D1 b
    (junk control); returns list of (score, b) best first"""
    S, G = B.S, B.G
    geo = [(l, g) for (l, g) in gen if B.on_geo(g)]
    cands = {}
    for (l, g) in geo[:6]:
        for b in B.routes(g).values():
            cands.setdefault(G.key(b), b)
    out = []
    for b in cands.values():
        sc = 0
        for (l, g) in geo:
            h = G.intern(G.mul(S.D1, b, G.MATS[g]))
            if S.depth(h) == S.depth(g) + 1 and S.dK(h) == S.dK(g) - 1:
                sc += 1
        out.append((sc, b))
    out.sort(key=lambda t: -t[0])
    return out

def build(B, tmpl, W, M):
    S, G, rng = B.S, B.G, B.rng
    ins = [S.inp(l) for l in range(W)]
    gens = []
    nprefix = M - W
    if tmpl == 'PL':
        for l in range(W):
            g = ins[l]
            for _ in range(nprefix // W):
                f = B.singer_filter(g, nterm=1)
                if f is None:
                    return None
                f = [(Fr(1), b) for (_, b) in f]
                g = B.mv([(g, f)]); gens.append(g)
        return gens
    if tmpl == 'H':
        inj = []
        for l in range(W):
            u = S.move([(1, None, ins[l])]); inj.append(u); gens.append(u)
        cur = None
        for step in range(nprefix - W):
            if step == 0:
                parts = [(u, B.singer_filter(u)) for u in inj]
            else:
                src = gens[-1] if rng.random() < 0.6 else rng.choice(gens[W:] or gens)
                parts = [(src, B.singer_filter(src))]
                if rng.random() < 0.5:                  # pattern changer: an injection routed one step
                    u = rng.choice(inj); parts.append((u, B.singer_filter(u, nterm=1)))
            if any(f is None for _, f in parts):
                return None
            gens.append(B.mv(parts))
        return gens + ins
    if tmpl == 'HW':
        # W x-injections, ONE shared spreading move a (Singer filters over the seven N<x,v>), W depth-3 moves b_j built
        # on a (junk-aware routes, eta-powers, optional pattern changer from an injection), W tops solved: 3W + 1.
        assert M >= 3 * W + 1
        nb = M - 2 * W - 1
        inj = []
        for l in range(W):
            u = S.move([(1, None, ins[l])]); inj.append(u); gens.append(u)
        parts = []
        for l in range(W):
            f = B.singer_filter(inj[l])
            if f is None:
                return None
            parts.append((inj[l], f))
        a = B.mv(parts); gens.append(a)
        sr = scored_routes(B, a)
        if not sr:
            return None
        top = [b for (sc, b) in sr[:6]]
        for j in range(nb):
            k = rng.choice([1, 1, 2, 2, 3])
            flt = []
            for _ in range(k):
                b = rng.choice(top); e = rng.randrange(len(S.etap))
                flt.append((rng.choice(UNITC), G.mul(S.etap[e], b)))
            ps = [(a, flt)]
            if rng.random() < 0.6:
                l = rng.randrange(W); f = B.singer_filter(inj[l], nterm=1)
                if f is not None:
                    ps.append((inj[l], f))
            gens.append(B.mv(ps))
        return gens + ins
    if tmpl in ('HWc', 'HWs'):
        # COHERENT variant: every path geodesic.  a = D1(sum_l gamma_l u_l) with Singer-conjugate routes
        # rho_k = eta^k rho mu^-k (array l spread over frames N<x, eta^k v0>, k in K_l, |K_l| <= 2; HWs: all arrays on the
        # SAME frame with left eta-powers, i.e. coefficients in C[<mu>]); b_j built only from monomials that advance
        # EVERY geodesic term of a (full score); W tops solved.
        nb = M - 2 * W - 1
        inj = []
        for l in range(W):
            u = S.move([(1, None, ins[l])]); inj.append(u); gens.append(u)
        g0 = next(iter(inj[0]))[1]
        rts = list(B.routes(g0).values())
        if not rts:
            return None
        rho = rng.choice(rts)
        mu = G.mul(S.D1, S.eta, np.conj(S.D1).T)
        mup = [np.eye(16, dtype=complex)]
        for _ in range(len(S.etap) - 1):
            mup.append(G.mul(mup[-1], mu))
        parts = []
        for l in range(W):
            nt = rng.choice([1, 2])
            ks = rng.sample(range(len(S.etap)), nt)
            flt = []
            for k in ks:
                if tmpl == 'HWc':
                    b = G.mul(S.etap[k], rho, np.conj(mup[k]).T)       # frame N<x, eta^k v0>
                else:
                    b = G.mul(S.etap[k], rho)                          # same frame, relative factor mu^k
                flt.append((rng.choice(UNITC), b))
            parts.append((inj[l], flt))
        a = B.mv(parts); gens.append(a)
        if not all(B.on_geo(g) for (l, g) in a):
            return 'offgeo'
        sr = scored_routes(B, a, samples=400)
        ngeo = sum(1 for (l, g) in a if B.on_geo(g))
        full = [b for (sc, b) in sr if sc == ngeo]
        if not full:
            return 'nofull'
        for j in range(nb):
            k = rng.choice([1, 1, 2])
            flt = [(rng.choice(UNITC), G.mul(S.etap[rng.randrange(len(S.etap))], rng.choice(full))) for _ in range(k)]
            gens.append(B.mv([(a, flt)]))
        return gens + ins
    if tmpl == 'C':
        # shared x-move with a Singer pattern, then the other injections (co-top), conveyor
        e = [rng.randrange(len(S.etap)) for _ in range(W)]
        sg = [rng.choice(UNITC) for _ in range(W)]
        u1 = S.move([(sg[l], S.etap[e[l]], ins[l]) for l in range(W)]); gens.append(u1)
        inj = [u1]
        for l in range(1, W):
            # injection: route u1 one step and add fresh raw x_l (this move must start a 4-chain: enforced by the conveyor)
            parts = [(u1, B.singer_filter(u1, nterm=1)), (ins[l], [(rng.choice(UNITC), S.etap[rng.randrange(len(S.etap))])])]
            if parts[0][1] is None:
                return None
            gens.append(B.mv(parts)); inj.append(gens[-1])
        for step in range(nprefix - W):
            src = gens[-1]
            parts = [(src, B.singer_filter(src))]
            if rng.random() < 0.5:
                u = rng.choice(gens[:-1]); f = B.singer_filter(u, nterm=1)
                if f is not None:
                    parts.append((u, f))
            if any(f is None for _, f in parts):
                return None
            gens.append(B.mv(parts))
        return gens + ins
    raise ValueError(tmpl)

def test_prefix(S, gens, W, Hname, nE=10, rng=None):
    Hids = S.etap_id if Hname == 'eta' else [S.G.ID]
    nh = len(Hids)
    pool = E_pool(S, gens, maxn=nE, rng=rng)
    spaces = {}; stats = {'pool': len(pool), 'cap': 0, 'dims': []}
    for idx, E in enumerate([None] + pool):
        V, st = target_space(S, gens, W, E, Hids)
        if V is None:
            stats['cap'] += 1; continue
        spaces['E%d' % idx if E is not None else 'none'] = V
        stats['dims'].append(int(V.shape[0]))
    combo, diag = best_tops(spaces, W, nh)
    stats['chardiag'] = diag
    return combo, stats

def main():
    tmpl, W, M, n, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    Hname = sys.argv[6] if len(sys.argv) > 6 else 'eta'
    out = sys.argv[7] if len(sys.argv) > 7 else None
    t0 = time.time()
    S = S4(4, seed=seed); rng = random.Random(seed); B = Builder(S, rng)
    agg = {'built': 0, 'fail_build': 0, 'pass': 0, 'maxdiag_hist': {}, 'cap': 0}
    hits = []
    for it in range(n):
        if len(S.G.MATS) > 40000:
            S = S4(4, seed=seed + 1000 + it); B = Builder(S, rng)
        gens = build(B, tmpl, W, M)
        if gens is None or isinstance(gens, str):
            agg['fail_build'] += 1
            if isinstance(gens, str):
                agg[gens] = agg.get(gens, 0) + 1
            continue
        agg['built'] += 1
        combo, st = test_prefix(S, gens, W, Hname, rng=rng)
        agg['cap'] += st['cap']
        key = str(min(st['chardiag']) if st['chardiag'] else 0)
        agg['maxdiag_hist'][key] = agg['maxdiag_hist'].get(key, 0) + 1
        agg.setdefault('dims_hist', {})
        dk = str(max(st['dims']) if st['dims'] else -1)
        agg['dims_hist'][dk] = agg['dims_hist'].get(dk, 0) + 1
        if combo is not None:
            agg['pass'] += 1; hits.append({'it': it, 'combo': combo, 'stats': st})
            print('PASS', it, combo, st, flush=True)
        if it % 10 == 9:
            print('progress', it + 1, agg, f'{time.time()-t0:.0f}s', flush=True)
    agg['seconds'] = round(time.time() - t0, 1)
    print(f'h7 {tmpl} W={W} M={M} n={n} seed={seed} H={Hname}: {agg}', flush=True)
    if out:
        json.dump({'args': sys.argv[1:], 'agg': agg, 'hits': hits}, open(out, 'w'))

if __name__ == '__main__':
    main()
