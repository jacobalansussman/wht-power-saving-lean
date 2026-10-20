"""symsearch.py -- mw-sym4: EXACT search over H-invariant programs (unions of tau-orbits of moves), engine game.
H = <(c_P, pi)> with c_P in H_max = B cap K B K^-1 (lift of P in GL(m,2)), pi a permutation of the W arrays.
Each step adds one tau-orbit of moves: representative u = D c, c = sum of 1..3 summands s_k b_k gen_k (b_k free:
progress letters of a geodesic lead term, or exact alignments onto the first summand's new term), optionally an
orbit-sum over H (so that u can be H-fixed); the orbit {tau_h u} is deduplicated up to free left factors and costs
|orbit| * rank.  End: exact span test (closure over frames + solve over Q; halving allowed) for every array.
Usage: symsearch.py m W budget Hspec n seed [rank2_frac] [orbitsum_frac]
Hspec: 'plain' (no symmetry), 'Z3rot' (coordinate 3-cycle on e1,e2,e3 + array 3-cycle), 'Z3fix' (same P, arrays fixed),
       'Z2swap' (coordinate swap e1<->e2 + array swap), 'Z7' (Singer on <e1,e2,e3>, arrays fixed),
       'Z4rot' (coordinate 4-cycle + array 4-cycle), 'V4' (double transpositions on coords+arrays),
       'Z3fixT' (Singer^? order-3 element of GL(3,2) fixing x=e4: e1->e2->e1+e2 ... arrays fixed),
       'GL32' (stabiliser of e4 acting as GL(3,2) on <e1,e2,e3>, arrays fixed; large)."""
import sys, time
from symlib import *

def perm_P(m, perm):
    """coordinate permutation e_j -> e_perm[j]"""
    return [1 << perm[j] for j in range(m)]

def hspec(S, name, W):
    m = S.m
    idW = tuple(range(W))
    cyc = lambda k: tuple((l + 1) % W for l in range(W)) if W > 1 else idW
    if name == 'plain':
        return [(np.eye(S.n, dtype=complex), idW)]
    gens = []
    if name == 'Z3rot':
        P = perm_P(m, [1, 2, 0] + list(range(3, m))); gens = [(P, cyc(1))]
    elif name == 'Z3fix':
        P = perm_P(m, [1, 2, 0] + list(range(3, m))); gens = [(P, idW)]
    elif name == 'Z2swap':
        P = perm_P(m, [1, 0] + list(range(2, m))); gens = [(P, cyc(1))]
    elif name == 'Z2swapNL':           # non-orthogonal involution e1 -> e1, e2 -> e1+e2 with array swap
        P = [1, 3] + [1 << j for j in range(2, m)]; gens = [(P, cyc(1))]
    elif name == 'Z4rot':
        P = perm_P(m, [1, 2, 3, 0] + list(range(4, m))); gens = [(P, cyc(1))]
    elif name == 'V4':
        P1 = perm_P(m, [1, 0, 3, 2]); P2 = perm_P(m, [2, 3, 0, 1])
        gens = [(P1, (1, 0, 3, 2)), (P2, (2, 3, 0, 1))]
    elif name == 'Z7':                 # Singer cycle of <e1,e2,e3>: x -> alpha x in F_8 = F_2[a]/(a^3+a+1)
        P = [2, 4, 3] + [1 << j for j in range(3, m)]; gens = [(P, idW)]
    elif name == 'Z7rot':              # Singer cycle of <e1,e2,e3> (fixing x = e4) with a 7-cycle of 7 arrays
        assert W == 7
        P = [2, 4, 3] + [1 << j for j in range(3, m)]; gens = [(P, cyc(1))]
    elif name == 'Z5rot':              # order-5 element of GL(4,2) (cube of a Singer cycle of F_16) with a 5-cycle of arrays
        P = [2, 4, 8, 15]; gens = [(P, cyc(1))]
    elif name == 'Z6rot':              # order-6 element: 3-cycle (e1 e2 e3) times e4 -> e4 + e1+e2+e3? (searched below)
        P = [2, 4, 1, 15]; gens = [(P, cyc(1))]
    elif name == 'Z3fixT':             # order-3 element of GL(3,2) acting on <e1,e2> as a Singer cycle of F_4
        P = [2, 3] + [1 << j for j in range(2, m)]; gens = [(P, idW)]
    elif name == 'Z3T_rot':            # same, with array 3-cycle (targets 3 in 10)
        P = [2, 3] + [1 << j for j in range(2, m)]; gens = [(P, cyc(1))]
    elif name == 'GL32':
        P1 = [2, 4, 3] + [1 << j for j in range(3, m)]; P2 = [2, 1, 4] + [1 << j for j in range(3, m)]
        gens = [(P1, idW), (P2, idW)]
    elif name == 'Z7Z3':               # Frobenius 7:3 on <e1,e2,e3> (Singer + Frobenius x -> x^2)
        P1 = [2, 4, 3] + [1 << j for j in range(3, m)]
        P2 = [1, 4, 6] + [1 << j for j in range(3, m)]   # e1->1, e2=a -> a^2=4, e3=a^2 -> a^4 = a^2+a = 6
        gens = [(P1, idW), (P2, idW)]
    else:
        raise ValueError(name)
    cg = []
    for (P, pi) in gens:
        r = S.hmax_elem(P)
        assert r is not None, P
        cg.append((r[0], pi))
    return S.group_closure(cg)

class Run:
    def __init__(self, m, W, Hname, seed):
        self.S = Sym(m); self.G = self.S.G; self.W = W; self.m = m
        self.H = hspec(self.S, Hname, W)
        self.rng = random.Random(seed)
        self.Hc = [(c, pi, self.G.mul(self.S.Ki, np.conj(c).T, self.S.K)) for (c, pi) in self.H]
        self._pl = {}

    def tau(self, h, gen):
        c, pi, cp = h
        G = self.G
        out = {}
        for (l, gid), co in gen.items():
            t = (pi[l], G.intern(G.mul(c, G.MATS[gid], cp)))
            out[t] = out.get(t, 0) + co
        return {t: v for t, v in out.items() if v != 0}

    def tau_term_cached(self, h, l, g):
        if not hasattr(self, '_tc'):
            self._tc = {}
        k = (id(h[0]), l, g)
        v = self._tc.get(k)
        if v is None:
            v = next(iter(self.tau(h, {(l, g): 1})))
            self._tc[k] = v
        return v

    def orbit(self, u):
        reps = [u]
        for h in self.Hc[1:]:
            v = self.tau(h, u)
            if not any(self.S.same_up_to_free(v, w) for w in reps):
                reps.append(v)
        return reps

    def depth(self, gid):
        G = self.G
        M = G.MATS[gid]
        d0 = int(round(np.log2(G.rowsupp(M).max())))
        return d0, G.dist(gid)

    def lead(self, gen):
        out = []
        for (l, g) in gen:
            d0, dk = self.depth(g)
            if d0 + dk == self.m and dk > 0:
                out.append((d0, l, g))
        out.sort(key=lambda t: -t[0])
        return out

    def plets(self, gid, steps):
        k = (gid, steps)
        if k not in self._pl:
            Dm = self.G.D1 if steps == 1 else self.G.D12
            pl = self.G.progress_letters(gid, self.rng, want=4, Dmat=Dm, steps=steps)
            if os.environ.get('MWSYM_COORD'):
                # coordinate world: keep only coordinate permutation letters (real 0/1 permutation matrices)
                pl = [b for b in pl if np.allclose(b.imag, 0) and np.allclose(np.abs(b.real) - b.real, 0)]
            self._pl[k] = pl
        return self._pl[k]

    def build(self, budget, rank2_frac=0.0, orbitsum_frac=0.3, maxsum=3, prefer_shared=0.7):
        G, rng, W = self.G, self.rng, self.W
        gens = [G.inp(l) for l in range(W)]
        # inputs must form a tau-stable set (they do: tau maps x_l to a free image of x_pi(l))
        cost = 0; prog = []
        tries = 0
        while cost < budget and tries < 60:
            tries += 1
            steps = 2 if (rng.random() < rank2_frac and budget - cost >= 2) else 1
            Dm = G.D1 if steps == 1 else G.D12
            cand = [i for i in range(len(gens)) if any(dk >= steps for (d0, l, g) in [(a, b, c) for (a, b, c) in self.lead(gens[i])] for dk in [self.depth(g)[1]])]
            if not cand:
                break
            # prefer recent generators
            i0 = rng.choices(cand, weights=[1 + 2 * k for k in range(len(cand))])[0]
            ld = [t for t in self.lead(gens[i0]) if self.depth(t[2])[1] >= steps]
            d0, l0, g0 = ld[0] if rng.random() < 0.6 else rng.choice(ld)
            pl = self.plets(g0, steps)
            if not pl:
                continue
            b0 = rng.choice(pl)
            summ = [(Fraction(1), G.lmul(b0, gens[i0]))]
            Tg = G.intern(G.mul(Dm, b0, G.MATS[g0]))
            ns = rng.randint(1, maxsum)
            if ns > 1:
                pool = []
                for i in range(len(gens)):
                    for (dd, l, g) in self.lead(gens[i]):
                        if i == i0 and g == g0:
                            continue
                        b = G.mul(np.conj(Dm).T, G.MATS[Tg], np.conj(G.MATS[g]).T)
                        if G.is_monomial(b):
                            pool.append((i, b, 2 + (1 if l != l0 else 0)))
                        elif dd == d0 and rng.random() < 0.3:
                            pls = self.plets(g, steps)
                            if pls:
                                pool.append((i, rng.choice(pls), 1 + (1 if l != l0 else 0)))
                rng.shuffle(pool)
                if rng.random() < prefer_shared:
                    pool.sort(key=lambda t: -t[2])
                for (i, b, _) in pool[:ns - 1]:
                    summ.append((Fraction(rng.choice([1, -1])), G.lmul(b, gens[i])))
            c = G.add(*summ)
            if not c:
                continue
            if len(self.Hc) > 1 and rng.random() < orbitsum_frac:
                parts = [(Fraction(1), c)] + [(Fraction(rng.choice([1, -1])) if rng.random() < 0.3 else Fraction(1), self.tau(h, c)) for h in self.Hc[1:]]
                c2 = G.add(*parts)
                if c2:
                    c = c2
            u = G.lmul(Dm, c)
            orb = self.orbit(u)
            add = len(orb) * steps
            if cost + add > budget:
                continue
            gens.extend(orb); cost += add
            prog.append((steps, len(orb), len(summ)))
        return gens, prog, cost

def main():
    m, W, budget, Hname, n, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
    r2 = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
    osf = float(sys.argv[8]) if len(sys.argv) > 8 else 0.3
    t0 = time.time()
    R = Run(m, W, Hname, seed)
    print(f"m={m} W={W} budget={budget} H={Hname} |H|={len(R.H)} n={n} seed={seed} r2={r2} osf={osf}", flush=True)
    stats = {"tested": 0, "hits": 0, "partial": 0, "cap": 0, "under": 0}
    hits = []
    for it in range(n):
        if len(R.G.MATS) > 80000:
            R = Run(m, W, Hname, seed * 1000 + it)
        gens, prog, cost = R.build(budget, r2, osf)
        if cost < budget:
            stats["under"] += 1
        res = final_all(R.S, gens, W)
        stats["tested"] += 1
        ny = sum(1 for s in res if s == "yes")
        if "cap" in res:
            stats["cap"] += 1
        if 0 < ny < W:
            stats["partial"] += 1
        if ny == W:
            stats["hits"] += 1
            hits.append({"prog": prog, "cost": cost})
            print("HIT", cost, prog, flush=True)
            if cost < W * m:
                import pickle
                pickle.dump({"gens": [{(l, R.G.key(R.G.MATS[g])): str(c) for (l, g), c in u.items()} for u in gens],
                             "prog": prog, "cost": cost}, open(f"HIT_{m}_{W}_{budget}_{Hname}_{seed}_{it}.pkl", "wb"))
        if (it + 1) % 200 == 0:
            stats["sec"] = round(time.time() - t0, 1)
            print("progress", it + 1, stats, flush=True)
    stats["sec"] = round(time.time() - t0, 1)
    print(stats, flush=True)
    if hits:
        print("hits:", hits[:5])

if __name__ == "__main__":
    main()
