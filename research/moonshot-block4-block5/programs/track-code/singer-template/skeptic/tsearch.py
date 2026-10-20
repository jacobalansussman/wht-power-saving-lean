"""tsearch.py -- ms-fano: EXACT search over Fano-gadget word templates (engine game letters, formal identities).
Register form (as tb4-beam bm.py): generators start as inputs x_l; a move makes u = D1 c, c = sum_j s_j b_j g_j with
b_j free (engine monomials: affine maps, Clifford Z4 phases, Paulis), s_j in {1,-1}; the end test asks
[K]x_l in the Q[B']-span of all generators (bm.Game.span_test: closure over frames + exact solve over Q, so the final
combination may halve), or with r spare moves the private free-last-move test (bm.Game.free_last_test).
Routes: for a term g on a geodesic L0 -> K-frame (d0 + dK = m) a ROUTE letter b sends it one step further on a
geodesic: D1 b g has d0 + 1, dK - 1.  One representative per next frame (sampled), plus EXACT alignment: to land a
second term g' on the very group element T already produced, b' = D1^-1 T g'^-1 is used when it is free.
Templates (blueprint of pf-duality section 7):
  FS  start gadget: the first k moves act on raw inputs only; each content spreads both inputs over the frames N<v>,
      v in a 3-space W3 (all seven directions used across the gadget); later moves are route sums of <= 3 summands,
      preferring Fano 3-route sums (three summands whose lead terms sit on the three points of a line, aligned onto one
      element of the line's frame) -- the halving then happens in the exact final solve.
  FN  (block 4) one shared rank-1 move u1 = D1 b_x (x1 + s x2) (both terms on N<x>), then a gadget layer whose contents
      route u1 into the seven rank-2 frames N<x,v> (several frames per content = multi-term coefficients) and inject a
      raw input term (so the x2/x1 pattern changes), then route sums.
  FE  mirrored end gadget: private plain chains to depth m-2 (or m-1), then shared route sums into the frames next to
      the kernel frame, final solve (halving) with no private last moves.
  PL  plain control: private chains only (must succeed with mW moves)."""
import sys, os, time, json, random, itertools
from fractions import Fraction
from fgeo import *

class FT:
    def __init__(self, m, seed):
        self.m = m; self.geo = Geo(m, seed); self.G = self.geo.G; self.rng = random.Random(seed)
        self._d = {}; self._routes = {}

    # ---- term geometry
    def depth(self, gid):
        v = self._d.get(gid)
        if v is None:
            M = self.G.MATS[gid]
            d0 = int(round(np.log2(self.G.rowsupp(M).max())))
            dk = self.G.dist(gid)
            v = (d0, dk); self._d[gid] = v
        return v

    def on_geo(self, gid):
        d0, dk = self.depth(gid)
        return d0 + dk == self.m

    def routes(self, gid, want=None, samples=400):
        """dict frame_key -> (b, new_gid) of next frames on geodesics"""
        r = self._routes.get(gid)
        if r is None:
            G = self.G; r = {}
            d0, dk = self.depth(gid)
            need = (1 << dk) - 1 if d0 + dk == self.m else 0
            g = G.MATS[gid]
            for _ in range(samples):
                if len(r) >= need:
                    break
                b = self.geo.rand_free()
                h = G.mul(G.D1, b, g)
                hid = G.intern(h)
                if self.depth(hid) == (d0 + 1, dk - 1):
                    fk = G.frm(hid)
                    if fk not in r:
                        r[fk] = (b, hid)
            self._routes[gid] = r
        return r

    def align(self, T_gid, gid):
        """free b with D1 b g == T exactly, or None"""
        G = self.G
        b = G.mul(np.conj(G.D1).T, G.MATS[T_gid], np.conj(G.MATS[gid]).T)
        return b if G.is_monomial(b) else None

    def lead_terms(self, gen):
        """geodesic terms of a generator, deepest first"""
        ts = [(l, g) for (l, g), c in gen.items() if self.on_geo(g)]
        ts.sort(key=lambda t: -self.depth(t[1])[0])
        return ts

    # ---- building contents
    def summand(self, gen, b, s):
        return (s, self.G.lmul(b, gen))

    def move(self, summands):
        G = self.G
        c = G.add(*summands)
        if not c:
            return None
        return G.lmul(G.D1, c)

def final_test(ft, gens, W, spare, cap=1500):
    G = ft.G
    res = []
    for l in range(W):
        if spare > 0:
            st, sol = G.free_last_test(gens, l, cap)
        else:
            st, sol = G.span_test(gens, l, cap)
            if st == "yes":
                st = "yes0"
        res.append(st)
    nfree = sum(1 for s in res if s == "yes")
    ok = all(s in ("yes", "yes0") for s in res) and nfree <= spare
    return ok, res

# ------------------------------------------------------------------ template generators
def gadget_contents(ft, W, k, W3, rng):
    """k raw-input contents spreading the inputs over frames N<v>, v in W3 (seven points); each content has 2-4 terms,
    every input appears in every content; across the gadget all seven points are used when 2-4 terms allow"""
    G = ft.G
    pts = [v for v in range(1, 1 << ft.m) if v in W3]
    reps = ft.geo_rep_letters
    used = []
    out = []
    order = pts[:]; rng.shuffle(order)
    for t in range(k):
        nt = rng.choice([2, 3, 4]) if W >= 2 else rng.choice([1, 2, 3])
        terms = []
        ls = list(range(W)); rng.shuffle(ls)
        for j in range(nt):
            v = order[(t * 3 + j) % len(order)]
            l = ls[j % W]
            b = rng.choice(reps[v])
            s = rng.choice([1, -1])
            terms.append((s, G.lmul(b, G.inp(l))))
        out.append(terms)
    return out

def route_sum(ft, gens, rng, maxs=3, prefer_fano=0.6, prefer_recent=True):
    """a content of 1..maxs summands that each advance a geodesic term of an earlier generator; with probability
    prefer_fano try a 3-route sum into one frame (exactly aligned when possible)"""
    G = ft.G
    cand = [i for i in range(len(gens)) if ft.lead_terms(gens[i])]
    if not cand:
        return None
    if prefer_recent:
        wts = [1 + i for i in cand]
    else:
        wts = [1] * len(cand)
    i0 = rng.choices(cand, weights=wts)[0]
    terms0 = ft.lead_terms(gens[i0])
    l0, g0 = terms0[0] if rng.random() < 0.7 else rng.choice(terms0)
    R = ft.routes(g0)
    if not R:
        return None
    fk, (b0, T) = rng.choice(list(R.items()))
    summ = [(rng.choice([1, -1]) if False else 1, G.lmul(b0, gens[i0]))]
    ns = rng.randint(1, maxs)
    if ns > 1:
        # other generators with a geodesic term adjacent below the frame fk: align onto T exactly if possible
        pool = []
        for i in cand:
            for (l, g) in ft.lead_terms(gens[i]):
                if ft.depth(g)[0] != ft.depth(T)[0] - 1 or (i == i0 and g == g0):
                    continue
                b = ft.align(T, g)
                if b is not None:
                    pool.append((i, b, True))
                else:
                    Rg = ft.routes(g)
                    if fk in Rg:
                        pool.append((i, Rg[fk][0], False))
        rng.shuffle(pool)
        if rng.random() < prefer_fano:
            pool.sort(key=lambda t: not t[2])
        for (i, b, ex) in pool[:ns - 1]:
            summ.append((rng.choice([1, -1]), G.lmul(b, gens[i])))
    return summ

def run_template(ft, kind, W, M, rng, spare_ok=True):
    G = ft.G
    gens = [G.inp(l) for l in range(W)]
    prog = []
    m = ft.m
    if kind == "PL":
        # private chains, each input m moves (control); M must be >= mW for success
        for l in range(W):
            cur = l
            for d in range(m):
                if len(prog) >= M:
                    break
                lt = ft.lead_terms(gens[cur])
                (ll, g) = lt[0]
                R = ft.routes(g)
                fk, (b, T) = rng.choice(list(R.items()))
                u = ft.move([(1, G.lmul(b, gens[cur]))])
                gens.append(u); prog.append(("plain", l)); cur = len(gens) - 1
        used = len(prog)
        return gens, prog, M - used
    if kind == "FS":
        W3 = ft.W3
        k = rng.choice([2, 3]) if m == 4 else rng.choice([2, 3, 4])
        k = min(k, M - 1)
        for terms in gadget_contents(ft, W, k, W3, rng):
            u = ft.move(terms)
            if u is None:
                return None
            gens.append(u); prog.append(("gadget", len(terms)))
    elif kind == "FN":
        x = ft.x_rep
        s = rng.choice([1, -1])
        u1 = ft.move([(1, G.lmul(x, G.inp(0))), (s, G.lmul(x, G.inp(1)))])
        gens.append(u1); prog.append(("shared_x", s))
        # gadget over N<x>: 2 moves routing u1 into several N<x,v>, plus one raw-input term
        for t in range(2):
            lt = ft.lead_terms(u1)
            g = lt[0][1]
            R = list(ft.routes(g).items()); rng.shuffle(R)
            nf = rng.choice([1, 2, 3])
            summ = [(rng.choice([1, -1]), G.lmul(b, u1)) for (fk, (b, T)) in R[:nf]]
            if rng.random() < 0.8:
                l = rng.randrange(W); v = rng.choice(list(ft.geo_rep_letters))
                summ.append((rng.choice([1, -1]), G.lmul(rng.choice(ft.geo_rep_letters[v]), G.inp(l))))
            u = ft.move(summ)
            if u is None:
                return None
            gens.append(u); prog.append(("gadgetNx", len(summ)))
    elif kind == "FE":
        # private chains to depth m-2 for every input, then shared route sums
        for l in range(W):
            cur = l
            for d in range(m - 2):
                lt = ft.lead_terms(gens[cur]); g = lt[0][1]
                R = ft.routes(g); fk, (b, T) = rng.choice(list(R.items()))
                u = ft.move([(1, G.lmul(b, gens[cur]))]); gens.append(u); prog.append(("chain", l)); cur = len(gens) - 1
    while len(prog) < M:
        summ = route_sum(ft, gens, rng)
        if summ is None:
            return None
        u = ft.move(summ)
        if u is None:
            return None
        gens.append(u); prog.append(("route", len(summ)))
    return gens, prog, 0

def setup(m, seed):
    ft = FT(m, seed)
    G = ft.G
    # representatives of the 2^m - 1 rank-1 interval frames over L0 (several letters each)
    reps = {}
    tries = 0
    while (len(reps) < (1 << m) - 1 or min(len(v) for v in reps.values()) < 4) and tries < 50000:
        tries += 1
        b = ft.geo.rand_free(); g = G.mul(G.D1, b)
        if ft.geo.d0(g) == 1 and ft.geo.dK(g) == m - 1:
            s = ft.geo.s_of(g)[0]
            reps.setdefault(s, [])
            if len(reps[s]) < 4:
                reps[s].append(b)
    ft.geo_rep_letters = reps
    ft.W3 = set(range(1, 8))                     # the 3-space <e1, e2, e3>
    ft.x_rep = reps[1][0]
    return ft

def main():
    m, W, M, kind, n, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
    out = sys.argv[7] if len(sys.argv) > 7 else None
    t0 = time.time()
    rng = random.Random(seed)
    ft = setup(m, seed)
    stats = {"tested": 0, "hits": 0, "one_output": 0, "cap": 0, "built_fail": 0}
    hits = []
    for it in range(n):
        if len(ft.G.MATS) > 60000:
            ft = setup(m, seed + it)            # memory: fresh interning tables
        r = run_template(ft, kind, W, M, rng)
        if r is None:
            stats["built_fail"] += 1; continue
        gens, prog, spare = r
        ok, res = final_test(ft, gens, W, spare)
        stats["tested"] += 1
        if "cap" in res:
            stats["cap"] += 1
        if sum(1 for s in res if s in ("yes", "yes0")) == 1 and W == 2:
            stats["one_output"] += 1
        if ok:
            stats["hits"] += 1
            hits.append({"prog": prog, "res": res})
    stats["seconds"] = round(time.time() - t0, 1)
    print(f"m={m} W={W} M={M} {kind} n={n} seed={seed}: {stats}", flush=True)
    if hits:
        print("first hits:", hits[:3])
    if out:
        json.dump({"args": sys.argv[1:], "stats": stats, "hits": hits[:20]}, open(out, "w"))

if __name__ == "__main__":
    main()
