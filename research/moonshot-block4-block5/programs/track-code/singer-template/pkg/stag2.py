"""stag2.py -- ms-fano: staggered shared lineage + ONE cleanup move per output, decided EXACTLY by linear algebra.
Lineage (2 arrays, block m): G1 = D b1 x1 ; G2 = D (r2 G1 + s b2 x2) ; G_t = D r_t G_{t-1}, t = 3..m+1, every r_t a
double advance (both inputs' lead terms one step further on geodesics while possible).  So x1's kernel term appears
in G_m and x2's in G_{m+1}, each with the other input's junk riding along.
Cleanup test for output l: is K x_l in  Span_{Q[B']}(inputs, G_1..G_{m+1})  +  E . Span_{Q[B']}(same)  for some paid
element E = gJ g''^-1 at distance 1 that maps a generator term g'' exactly onto a junk term gJ of a generator holding
x_l's kernel term?  (bm.Game.span_test with E: the E-columns all share ONE extra move, whose content is solved for.)
If both outputs pass, the word has m+1+2 moves: block 3 -> 6 (plain count), block 4 -> 7 (2 in 7, the target).
Each 'yes' is an exact solve over Q of the formal identity (closure over frames).
usage: python3 stag2.py m n seed [maxE] [out.json]"""
import sys, time, json, random
from tsearch import *
from staggered import double_advance

def junk_E(ft, gens, l, maxE=12):
    G = ft.G
    out = {}
    for t, gen in enumerate(gens):
        kt = [(ll, g) for (ll, g) in gen if ll == l and ft.depth(g)[1] == 0]
        if not kt:
            continue
        for (lj, gj) in gen:
            if (lj, gj) in kt:
                continue
            for s, gen2 in enumerate(gens):
                for (l2, g2) in gen2:
                    if l2 != lj:
                        continue
                    E = G.mul(G.MATS[gj], np.conj(G.MATS[g2]).T)
                    rs = G.rowsupp(E)
                    if not (rs == 2).all():
                        continue
                    k = G.key(E)
                    if k not in out:
                        out[k] = E
                        if len(out) >= maxE:
                            return list(out.values())
    return list(out.values())

def build_lineage(ft, m, rng):
    G = ft.G; reps = ft.geo_rep_letters
    gens = [G.inp(0), G.inp(1)]
    v1 = rng.randrange(1, 1 << m); b1 = rng.choice(reps[v1])
    gens.append(ft.move([(1, G.lmul(b1, gens[0]))]))
    # move 2: route G1 one step AND inject x2; retry choices until a double advance exists afterwards
    for attempt in range(30):
        R = list(ft.routes(ft.lead_terms(gens[2])[0][1]).values())
        r2 = rng.choice(R)[0]
        v2 = rng.randrange(1, 1 << m); b2 = rng.choice(reps[v2]); s = rng.choice([1, -1])
        u = ft.move([(1, G.lmul(r2, gens[2])), (s, G.lmul(b2, gens[1]))])
        if double_advance(ft, u, rng, samples=3000):
            break
    else:
        return None
    gens.append(u)
    for t in range(3, m + 2):
        B = double_advance(ft, gens[-1], rng, samples=4000)
        if not B:
            return None
        gens.append(ft.move([(1, G.lmul(rng.choice(B), gens[-1]))]))
    return gens

def main():
    m, n, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    maxE = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    out = sys.argv[5] if len(sys.argv) > 5 else None
    t0 = time.time(); rng = random.Random(seed)
    ft = setup(m, seed)
    st = {"lineages": 0, "fail_build": 0, "x1_pass": 0, "x2_pass": 0, "both": 0, "cap": 0, "E_tested": 0}
    hits = []
    for it in range(n):
        if len(ft.G.MATS) > 40000:
            ft = setup(m, seed + 1000 + it)
        gens = build_lineage(ft, m, rng)
        if gens is None:
            st["fail_build"] += 1; continue
        st["lineages"] += 1
        ok = []
        for l in range(2):
            passed = False
            for E in junk_E(ft, gens, l, maxE):
                st["E_tested"] += 1
                r, sol = ft.G.span_test(gens, l, 1500, E)
                if r == "cap":
                    st["cap"] += 1
                if r == "yes":
                    passed = True; break
            ok.append(passed)
        st["x1_pass"] += ok[0]; st["x2_pass"] += ok[1]
        if all(ok):
            st["both"] += 1; hits.append(it)
    st["seconds"] = round(time.time() - t0, 1)
    print(f"stag2 m={m} n={n} seed={seed} maxE={maxE}: {st}", flush=True)
    if out:
        json.dump({"args": sys.argv[1:], "stats": st, "hits": hits}, open(out, "w"))

if __name__ == "__main__":
    main()
