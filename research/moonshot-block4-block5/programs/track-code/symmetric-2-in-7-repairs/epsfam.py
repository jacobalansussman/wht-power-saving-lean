"""epsfam.py -- mw2-absorb4 (resumed run), question (b): does the FULL problem carry the relaxation's border family?
The relaxation path of the seed-66 skeleton tends to the limit skeleton (rs66_1_final_lim: F1,F2,L1,F1) by sending the
relative weight eps of one summand of the shared move L to 0.  Rebuild the one-parameter family: the limit skeleton
with that summand put back (weight eps), every other coefficient free.  Full exact-structure test (closure over frames,
variable-projection LM over all content coefficients including eps): if the full residual could be pushed to 0 only as
eps -> 0 (or as some coefficient diverges) the full problem carries the border; if its minimum stays at the limit's
value (2/3), it does not.  Also scans eps = 1, 1e-1, ..., 1e-4 with the other coefficients optimised.
Usage: epsfam.py rs66_1_final.pkl cap starts"""
import sys, pickle, time
from twist import *

d = pickle.load(open(sys.argv[1], 'rb'))
cap, starts = int(sys.argv[2]), int(sys.argv[3])
m, W, H = d['m'], d['W'], d['H']
rec, al = d['recipe'], np.asarray(d['alpha'])
sk0 = build_from_recipe(m, W, H, rec, d['skseed'])
# projective normalisation as in normalize.py, but keep EVERY summand of the shared move (kind 'L')
sigma = {i: 1.0 + 0j for i in range(W)}
newrec = []; rel_by_move = []
a = 0
for t, (kind, summ, Dm) in enumerate(rec):
    w = [al[a + k] * sigma[src] for k, (src, b) in enumerate(summ)]
    a += len(summ)
    kd = int(np.argmax(np.abs(w))); st = w[kd]
    rel = [x / st for x in w]
    keep = [k for k in range(len(summ)) if abs(rel[k]) >= 1e-3 or kind == 'L']
    keep.sort(key=lambda k: (k != kd))
    newrec.append((kind, [summ[k] for k in keep], Dm))
    rel_by_move.append([rel[k] for k in keep])
    print(f"move {t} ({kind}): kept {len(keep)} of {len(summ)}; relative |coef| {np.round(np.abs([rel[k] for k in keep]), 6)}")
    mv = sk0.moves[t]
    for g in [mv['rep']] + mv['images']:
        sigma[g] = st
sk = build_from_recipe(m, W, H, newrec, d['skseed'])
N = sk.nalpha; free_idx = [x for x in range(N) if x not in sk.fixed]; nf = len(free_idx)
al0 = np.concatenate([np.array(v, dtype=complex) for v in rel_by_move])
print("family skeleton: summands per move", [len(s) for (_, s, _) in newrec], "terms per generator", [len(T) for T in sk.gens],
      "free coefficients", free_idx)
pickle.dump({"m": m, "W": W, "H": H, "recipe": newrec, "skseed": d['skseed'], "alpha": al0}, open("epsfam_66.pkl", "wb"))
t0 = time.time()
ok = sk.setup_final(cap)
print("closure dims", [(nr, nc) for (nr, nc, E) in sk.final] if ok else f"cap {cap} exceeded", f"({time.time()-t0:.0f}s)", flush=True)
if not ok:
    sys.exit()
prep_fast(sk)
# which free coefficient is eps: the non-leading summand of the L move
aL = sum(len(s) for (_, s, _) in newrec[:2]) + 1
print("eps index", aL, "value on the path", al0[aL])
def full(x):
    v = np.ones(N, dtype=complex); v[free_idx] = x[:nf] + 1j * x[nf:]; return v
rng = random.Random(9)
# (1) unconstrained minimum over all coefficients, several starts
for st in range(starts):
    x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)]) if st else np.concatenate([al0[free_idx].real, al0[free_idx].imag])
    x, c = vp_lm(sk, N, free_idx, x0, iters=200)
    v = full(x)
    print(f"free start {st}: full res^2 {c:.6f}  |eps| {abs(v[aL]):.3g}  max|alpha| {np.abs(v).max():.3g}", flush=True)
# (2) eps scan: eps fixed, the rest optimised
fi2 = [x for x in free_idx if x != aL]; nf2 = len(fi2)
for eps in [1.0, 1e-1, 1e-2, 1e-3, 1e-4]:
    best = np.inf
    for st in range(max(2, starts // 2)):
        x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf2)])
        def fullc(y, eps=eps):
            v = np.ones(N, dtype=complex); v[fi2] = y[:nf2] + 1j * y[nf2:]; v[aL] = eps; return v
        # local LM on the reduced variables
        def resid(y):
            r, _ = vp_eval(sk, fullc(y), [], want_J=False); return r
        y, c = lm(resid, x0, iters=60)
        best = min(best, c)
    print(f"eps {eps:g}: full res^2 (others optimised, {max(2, starts // 2)} starts) {best:.6f}", flush=True)
