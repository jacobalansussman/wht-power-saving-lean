# lemma_calib.py -- calibration of Lemma W (two-class stretch): on every sampled two-class-stretchable layout the full
# undirected cut LP must equal W*m; also report the cut LP on non-stretchable layouts (where it may fall short).
import sys, random, time
import w3cut, w3image, stretch2
m = int(sys.argv[1]); W = int(sys.argv[2]); N = int(sys.argv[3]); seed = int(sys.argv[4])
g = w3cut.G(m); rng = random.Random(seed); t0 = time.time()
res = {True: [], False: []}
for _ in range(N):
    P = stretch2.rand_layout(g, rng, W, distinct=False)
    if len({s for s, t in P}) + len({t for s, t in P}) < 2 and False: continue
    st = stretch2.two_class(g.D, m, P) is not None
    v = w3image.cut_value(g.n, g.edges, P)[0]
    res[st].append(v)
    if st and v < W * m - 1e-6: print("LEMMA VIOLATION", P, v, flush=True)
for st in (True, False):
    vs = res[st]
    if vs: print(f"m={m} W={W} stretchable={st}: {len(vs)} layouts, cut LP min {min(vs):.4f} max {max(vs):.4f} (W*m = {W*m}); below W*m: {sum(1 for v in vs if v < W*m - 1e-6)}")
print(f"({time.time()-t0:.0f}s)")
