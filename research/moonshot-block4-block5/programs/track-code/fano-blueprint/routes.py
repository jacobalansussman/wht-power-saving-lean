"""routes.py -- ms-fano E1: route group elements of the Fano gadget in the engine's letters, block m.
For a flag of interval frames L0 < N<v> < N<s> (s a 2-space containing v) we pick free letters b1, b2 with
D1 b1 at N<v> and D1 b2 D1 b1 at N<s>.  Three routes into N<s> (v in s): r_v = D1 b2_v D1 b1_v.  They lie in one coset
B' r; f_v = r_v r_{v0}^-1 is a monomial (free) matrix.  Question: can the free letters be chosen so that the three
routes are EQUAL group elements (f_v = 1)?  Then a content sum over the three routes adds terms of one input on one
group element (cancellation possible), which is what a halving gadget needs.  We also record the set of free factors
reachable when b1, b2 range over random choices (the coset structure)."""
import sys, itertools
from fgeo import *

def main(m, trials=4000):
    geo = Geo(m, seed=7); G = geo.G; D1 = G.D1
    # stage 1: many b1 per direction v
    b1s = {}
    for _ in range(trials):
        b = geo.rand_free(); g = G.mul(D1, b)
        if geo.d0(g) == 1 and geo.dK(g) == m - 1:
            b1s.setdefault(geo.s_of(g), []).append(b)
    print(f"m={m}: directions with b1 samples: { {k: len(v) for k, v in sorted(b1s.items())} }")
    # stage 2: for 2-spaces s, routes from each v in s
    planes = set()
    for v in range(1, 1 << m):
        for w in range(v + 1, 1 << m):
            planes.add(tuple(sorted(span_of([v, w]) - {0})))
    planes = sorted(planes)
    res = {}
    for s in planes[:7]:
        routes = {}
        for v in s:
            for b1 in b1s[(v,)][:6]:
                g1 = G.mul(D1, b1)
                for _ in range(400):
                    b2 = geo.rand_free(); g2 = G.mul(D1, b2, g1)
                    if geo.d0(g2) == 2 and geo.dK(g2) == m - 2 and geo.s_of(g2) == s:
                        routes.setdefault(v, []).append(g2)
        # compare: do routes from different v ever coincide as matrices (exactly), or up to a scalar / Pauli?
        keys = {v: set(G.key(r) for r in rs) for v, rs in routes.items()}
        allk = [k for v in keys for k in keys[v]]
        coinc = []
        vs = sorted(routes)
        for a, b in itertools.combinations(vs, 2):
            coinc.append(len(keys[a] & keys[b]))
        # up to scalars
        def skey(M):
            # normalise global phase by first nonzero entry
            idx = np.argmax(np.abs(M.ravel()) > 1e-9); ph = M.ravel()[idx] / abs(M.ravel()[idx])
            return G.key(G.snap(M / ph)) if True else None
        sk = {v: set(skey(r) for r in rs) for v, rs in routes.items()}
        coinc_s = [len(sk[a] & sk[b]) for a, b in itertools.combinations(vs, 2)]
        print(f"plane {s}: routes per v {[len(routes.get(v, [])) for v in s]}; exact coincidences between v-pairs {coinc}; up to scalar {coinc_s}")
        res[s] = routes
    return geo, b1s, res

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
