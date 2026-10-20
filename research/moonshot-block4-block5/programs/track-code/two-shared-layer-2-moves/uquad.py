"""uquad.py -- w5-twoshared (resumed): for private cores, U_x = L0 cap Lag(D mu_x D) (the 2-space of L0 that defines
the quad through L0 and Lag(D mu_x D)).  An s_j -> v solution survives D2 and D2' only if the neighbour N_c of L0
lies in both quads, which needs dim(U_e + U_j) <= 3, i.e. U_e cap U_j != 0.
Numpy-free replica of sweep.core(rng, 'priv') (same rng consumption: lift draws two randrange(16) per letter), so the
cores are the same as in esearch / rd2batch.  Usage: uquad.py ncores seed ..."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sp
D = sp.D; comp = sp.comp; lag = sp.lag
def rand_full(rng):
    S, P, Bm = sp.letter_full(rng); rng.randrange(16); rng.randrange(16)
    return (S, None)
def route(rng):
    while True:
        b = rand_full(rng); s = sp.mv(b[0], sp.I8)
        for _ in range(50):
            a = rand_full(rng); p = sp.mv(a[0], s)
            if sp.lam(p) == 2: break
        else: continue
        for _ in range(600):
            g = rand_full(rng); q = sp.mv(g[0], p)
            if sp.lam(q) == 3 and sp.dist(sp.lag(q), sp.LK) == 1:
                for _ in range(3000):
                    mu = rand_full(rng)
                    if sp.lag(sp.mv(mu[0], q)) == sp.LK:
                        return dict(b=b, a=a, g=g, mu=mu, s=s, p=p, q=q)
                break
def tags_priv(R0, R1):
    tags = []
    for e, (Re, Rj) in enumerate(((R0, R1), (R1, R0))):
        J = comp(comp(D, Re['mu'][0]), Rj['q']); L = sp.lam(J)
        if L == 4: tags.append('K' if lag(J) == sp.LK else '4')
        elif L == 3: tags.append('3V' if lag(Rj['q']) == lag(J) else '3')
        else: tags.append(str(L))
    return tags
def core_priv(rng, tries=400):
    for _ in range(tries):
        R0 = route(rng); R1 = route(rng)
        tg = tags_priv(R0, R1)
        if any(t in ('4', '3', 'K') for t in tg): continue
        return R0, R1, tg
def dim(S): return len(S).bit_length() - 1
if __name__ == "__main__":
    tot = {}
    ncores = int(sys.argv[1])
    for seed in [int(a) for a in sys.argv[2:]]:
        rng = random.Random(seed)
        for c in range(ncores):
            R0, R1, tg = core_priv(rng)
            U = [sp.L0 & lag(comp(D, comp(r['mu'][0], D))) for r in (R0, R1)]
            key = (dim(U[0]), dim(U[1]), dim(U[0] & U[1]))
            tot[key] = tot.get(key, 0) + 1
            print(f"seed {seed} core {c}: tags {tg}, dim U_0 {key[0]}, dim U_1 {key[1]}, dim(U_0 cap U_1) {key[2]}", flush=True)
    print("TALLY (dim U_0, dim U_1, dim cap):", tot, flush=True)
