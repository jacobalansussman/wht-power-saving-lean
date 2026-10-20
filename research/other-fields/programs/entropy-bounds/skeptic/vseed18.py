# vseed18.py -- replay nc-gain's seed-18 block-4 screen (w3cut.rand_3k2 with random.Random(18)) on my Gamma_4 and
# classify the first eight layouts (the claim: 7 compatible, #8 = [(2051,594),(2031,712),(777,952)]).
import sys, random
import numpy as np
sys.path.insert(0, '.')
import vlemmas
E = np.load('gamma4_edges_dpgidx.npy'); D = np.load('gamma4_D_dpgidx.npy').astype(np.int16)
class G: pass
g = G(); g.m = 4; g.n = D.shape[0]; g.D = D; g.opp = [np.nonzero(D[s] == 4)[0].tolist() for s in range(g.n)]
def rand_3k2(g, rng, W=3):      # copied logic of nc-gain w3cut.rand_3k2 (same RNG calls)
    while True:
        P = []; used = set()
        for _ in range(W):
            for _t in range(100):
                s = rng.randrange(g.n); t = rng.choice(g.opp[s])
                if s not in used and t not in used: break
            P.append((s, t)); used |= {s, t}
        if len(used) == 2 * W: return P
rng = random.Random(18)
for k in range(10):
    P = rand_3k2(g, rng); e = vlemmas.edef(g, P)[0]
    print(k + 1, P, 'e =', e)
