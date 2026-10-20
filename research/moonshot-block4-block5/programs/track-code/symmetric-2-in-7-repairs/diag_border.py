"""diag_border.py -- mw2-absorb4 step 1: anatomy of mw-sym4's border 2-in-7 in the frame-sum relaxation.
Loads a skeleton (initial recipe of (seed,start), or a pickled recipe from rs2.py), prints its structure, runs LM in
rounds and logs the coefficient path (which alpha diverge, at what rate against the residual).
Usage: diag_border.py seed start [recipe.pkl]"""
import sys, pickle, time
from relax_search import *

m, W, H, shape, menu = 4, 2, 'Z2swap', 'P:inp,P:prev,L:prev,P:prev'.split(','), 2
seed, s = int(sys.argv[1]), int(sys.argv[2])
if len(sys.argv) > 3:
    d = pickle.load(open(sys.argv[3], 'rb'))
    sk = build_from_recipe(m, W, H, d['recipe'], d['skseed'])
    rec = d['recipe']
else:
    sk, cost, rec = record_recipe(m, W, H, shape, menu, seed * 100 + s)
G = sk.G
print("orbits", [1 + len(mv["images"]) for mv in sk.moves], "terms per generator", [len(T) for T in sk.gens])
for t, (k, sm, dd) in enumerate(rec):
    print(f"move {t}: kind {k}, summands {[(src) for (src, b) in sm]}, rep generator {sk.moves[t]['rep']}, images {sk.moves[t]['images']}")
print("alpha index -> (move, summand):")
amap = {}
a = 0
for t, (k, sm, dd) in enumerate(rec):
    for j in range(len(sm)):
        amap[a] = (t, j); a += 1
print(amap, "fixed", sk.fixed)
fidx = {}
for i, T in enumerate(sk.gens):
    for (l, g) in T:
        key = (G.frm(g), l)
        if key not in fidx: fidx[key] = len(fidx)
Kkey = (G.frm(G.KID), 0); nrow = len(fidx)
dist = {}
for i, T in enumerate(sk.gens):
    for (l, g) in T:
        dist[fidx[(G.frm(g), l)]] = G.dist(g)
rowmap = [np.array([fidx[(G.frm(g), l)] for (l, g) in T]) for T in sk.gens]
N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
print("rows (frame,array)", nrow, "free alphas", nf, "K row", fidx[Kkey])
def full(x):
    al = np.ones(N, dtype=complex); al[free_idx] = x[:nf] + 1j * x[nf:]; return al
def Amat(al):
    cv = sk.coefs(al)
    A = np.zeros((nrow, len(sk.gens)), dtype=complex)
    for i in range(len(sk.gens)): np.add.at(A[:, i], rowmap[i], cv[i])
    return A
def rr(al):
    A = Amat(al)
    t = np.zeros(nrow, dtype=complex); t[fidx[Kkey]] = 1
    beta = np.linalg.lstsq(A, t, rcond=1e-12)[0]
    return A @ beta - t, beta
rng = random.Random(3)
x = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
path = []
for rnd in range(int(sys.argv[4]) if len(sys.argv) > 4 else 8):
    x, c = lm(lambda y: rr(full(y))[0], x, iters=40)
    al = full(x)
    r, beta = rr(al)
    path.append((c, al.copy(), beta.copy()))
    print(f"round {rnd}: res^2 {c:.3e} |alpha| {np.round(np.abs(al), 3)} |beta| {np.round(np.abs(beta), 3)}", flush=True)
pickle.dump({"recipe": rec, "path": path, "skseed": seed * 100 + s, "amap": amap, "fixed": sk.fixed}, open(f"diag_{seed}_{s}.pkl", "wb"))
# rate fit: log|alpha_a| vs log res
cs = np.array([p[0] for p in path]); A_ = np.array([np.abs(p[1]) for p in path])
ok = cs > 0
for a in range(N):
    if a in sk.fixed: continue
    if A_[-1, a] > 10 * A_[0, a] and A_[-1, a] > 5:
        sl = np.polyfit(np.log(cs[ok][-4:]), np.log(A_[ok, a][-4:]), 1)[0]
        print(f"alpha[{a}] {amap[a]} diverges: |alpha| {A_[-1,a]:.3g}, log|alpha| ~ {sl:.3f} log(res^2)")
