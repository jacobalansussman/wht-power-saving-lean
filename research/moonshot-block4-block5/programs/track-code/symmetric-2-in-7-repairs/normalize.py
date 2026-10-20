"""normalize.py -- mw2-absorb4: projective normalisation of a skeleton's content coefficients.
Each move's content is only defined up to a scalar (rescaling a generator is undone by the later contents), so the true
parameter space is a product of projective spaces.  Given a recipe + alpha (rs2.py / diag pickles), propagate the scale
generator by generator (u_t -> u_t / sigma_t, later coefficients multiplied by sigma_src), put the dominant summand first
(coefficient exactly 1), and DROP summands whose relative coefficient is below tol.  The result is the projective limit
point the diverging LM path is heading for.  Then: relaxation residual at the limit (LM-polished), closure size, full test.
Usage: normalize.py in.pkl tol fullstarts [cap]"""
import sys, pickle, time
from relax_search import *
from relax_sym import relax_matrix

d = pickle.load(open(sys.argv[1], 'rb'))
tol = float(sys.argv[2]); fullstarts = int(sys.argv[3]); cap = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
m, W, H = d.get('m', 4), d.get('W', 2), d.get('H', 'Z2swap')
rec, al = d['recipe'], np.asarray(d['alpha'])
sk = build_from_recipe(m, W, H, rec, d['skseed'])
sigma = {i: 1.0 + 0j for i in range(W)}
newrec = []; newal_by_move = []
a = 0
for t, (kind, summ, Dm) in enumerate(rec):
    w = [al[a + k] * sigma[src] for k, (src, b) in enumerate(summ)]
    a += len(summ)
    kd = int(np.argmax(np.abs(w))); st = w[kd]
    rel = [x / st for x in w]
    keep = [k for k in range(len(summ)) if abs(rel[k]) >= tol]
    keep.sort(key=lambda k: (k != kd))
    newrec.append((kind, [summ[k] for k in keep], Dm))
    newal_by_move.append([rel[k] for k in keep])
    print(f"move {t} ({kind}): scale {abs(st):.3g}, relative |coef| {np.round(np.abs(rel), 5)} -> keep {keep}")
    mv = sk.moves[t]
    for g in [mv['rep']] + mv['images']:
        sigma[g] = st
sk2 = build_from_recipe(m, W, H, newrec, d['skseed'])
N = sk2.nalpha; free_idx = [x for x in range(N) if x not in sk2.fixed]; nf = len(free_idx)
al0 = np.concatenate([np.array(v, dtype=complex) for v in newal_by_move])
assert len(al0) == N
print("limit skeleton: terms per generator", [len(T) for T in sk2.gens], "free alphas", nf)
# relaxation residual at the limit point and after polishing
G = sk2.G
fidx = {}
for i, T in enumerate(sk2.gens):
    for (l, g) in T:
        key = (G.frm(g), l)
        if key not in fidx: fidx[key] = len(fidx)
Kkey = (G.frm(G.KID), 0); nrow = len(fidx)
rowmap = [np.array([fidx[(G.frm(g), l)] for (l, g) in T]) for T in sk2.gens]
def rr(alv):
    cv = sk2.coefs(alv)
    A = np.zeros((nrow, len(sk2.gens)), dtype=complex)
    for i in range(len(sk2.gens)): np.add.at(A[:, i], rowmap[i], cv[i])
    t = np.zeros(nrow, dtype=complex); t[fidx[Kkey]] = 1
    beta = np.linalg.lstsq(A, t, rcond=1e-12)[0]
    return A @ beta - t
r0 = float(np.sum(np.abs(rr(al0)) ** 2))
x0 = np.concatenate([al0[free_idx].real, al0[free_idx].imag])
def full(x):
    v = np.ones(N, dtype=complex); v[free_idx] = x[:nf] + 1j * x[nf:]; return v
x1, c1 = lm(lambda y: rr(full(y)), x0, iters=100) if nf else (x0, r0)
print(f"relaxation at the limit point: res^2 {r0:.3e}; polished {c1:.3e}; max|alpha| {np.abs(full(x1)).max():.3g}")
out = sys.argv[1].replace('.pkl', '_lim.pkl')
pickle.dump({"m": m, "W": W, "H": H, "recipe": newrec, "skseed": d['skseed'], "alpha": full(x1), "relax": c1}, open(out, "wb"))
A, rowinfo, Krow, fi = relax_matrix(sk2)
names = [f"a{x}" for x in fi]
print("relaxation matrix (rows at frame distance 0 and the K row):")
for r, (dd, l, g) in enumerate(rowinfo):
    if dd == 0:
        ents = [(i, A[(r, i)]) for i in range(len(sk2.gens)) if (r, i) in A]
        print(f"  row {r} arr {l}: " + " | ".join(f"g{i}: {p.__repr__(names)}" for i, p in ents))
t0 = time.time()
ok = sk2.setup_final(cap)
print("closure dims", [(nr, nc) for (nr, nc, E) in sk2.final] if ok else f"cap {cap}", f"({time.time()-t0:.0f}s)", flush=True)
if ok and fullstarts:
    prep_fast(sk2)
    rng = random.Random(11)
    best = np.inf
    starts = [x1] + [np.array([rng.gauss(0, 1) for _ in range(2 * nf)]) for _ in range(fullstarts - 1)]
    for st, xs in enumerate(starts):
        x, c = vp_lm(sk2, N, free_idx, xs, iters=150)
        best = min(best, c)
        print(f"full start {st}{' (from the relaxation limit)' if st == 0 else ''}: res^2 {c:.4e}  max|alpha| {np.abs(full(x)).max():.3g}", flush=True)
    print(f"FULL best res^2 {best:.4e}")
