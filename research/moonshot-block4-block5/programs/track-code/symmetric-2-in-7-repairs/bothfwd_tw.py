"""bothfwd_tw.py -- mw2-absorb4: same question as bothfwd.py, for TWISTED absorbed shared moves (every admissible twist
beta, including twists with a non-trivial address part).  Usage: bothfwd_tw.py lim.pkl tsh nsamp"""
import sys, pickle
from twist import *
d = pickle.load(open(sys.argv[1], 'rb')); tsh = int(sys.argv[2]); ns = int(sys.argv[3])
m, W, H = d.get('m', 4), d.get('W', 2), d.get('H', 'Z2swap')
rec = d['recipe']
base = build_T(m, W, H, rec, d['skseed'], {})
G, R = base.G, base.R; n = G.n
kind, summ, Dm = rec[tsh]
rng = random.Random(2)
cands = beta_candidates(G, R, Dm, rng, nrand=20000)
def xpart_trivial(M):
    sig, q = G.mono_of(M)
    s0 = sig[0]
    return all(sig[x] == (x ^ s0) for x in range(n))      # translation (Pauli-type) permutation
print("admissible twists", len(cands), "with non-translation address part:", sum(1 for (b, o, bt) in cands if not xpart_trivial(b)))
letters = []
for it in range(ns):
    if it % 2: letters.append(rand_free_b(G, rng, m))
    else:
        sig = list(range(n)); rng.shuffle(sig); letters.append(G.mono(sig, [rng.randrange(4) for _ in range(n)]))
for k, (bt, o, btot) in enumerate(cands):
    if xpart_trivial(bt) and k > 3:
        continue
    rec2 = list(rec); rec2[tsh] = ('T', summ, Dm)
    sk = build_T(m, W, H, rec2, d['skseed'], {tsh: bt}, R)
    if not sk.moves[tsh]['absorbed']:
        continue
    T = sk.gens[sk.moves[tsh]['rep']]
    h0 = [g for (l, g) in T if l == 0]; h1 = [g for (l, g) in T if l == 1]
    one = 0; dists = {}
    for b in letters:
        g0 = G.intern(G.mul(G.D1, b, G.MATS[h0[0]]))
        if G.dist(g0) != 0: continue
        one += 1
        x = G.dist(G.intern(G.mul(G.D1, b, G.MATS[h1[0]]))); dists[x] = dists.get(x, 0) + 1
    print(f"twist {k} (order(beta_tot) {o}, translation part only: {xpart_trivial(bt)}): shared dists a0 {[G.dist(g) for g in h0]} a1 {[G.dist(g) for g in h1]}; forward on a0: {one}; a1 distances {dists}", flush=True)
