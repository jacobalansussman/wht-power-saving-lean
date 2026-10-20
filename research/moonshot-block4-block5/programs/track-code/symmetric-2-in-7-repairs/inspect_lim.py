"""inspect_lim.py -- mw2-absorb4: anatomy of a (small) limit skeleton: generators' terms, the tau action on the shared
move (fixed up to which free factor), the full closure for array 0 with symbolic entries, and the least-squares residual
split by rows.  Usage: inspect_lim.py lim.pkl"""
import sys, pickle
from relax_search import *
from relax_sym import sym_coefs

d = pickle.load(open(sys.argv[1], 'rb'))
m, W, H = d['m'], d['W'], d['H']
sk = build_from_recipe(m, W, H, d['recipe'], d['skseed'])
G, R = sk.G, sk.R
cv, free_idx = sym_coefs(sk)
names = [f"a{x}" for x in free_idx]
print("free", names, "alpha", np.round(d['alpha'][free_idx], 4))
ids = {}
def nm(g):
    if g not in ids: ids[g] = f"h{len(ids)}"
    return ids[g]
for i, T in enumerate(sk.gens):
    print(f"g{i} ({sk.coef_of[i][0]}): " + "  ".join(f"[{nm(g)} d{G.dist(g)}]x{l}*({cv[i][j].__repr__(names)})" for j, (l, g) in enumerate(T)))
# tau on the shared move(s): which free factor beta with tau(u) = beta u
for t, mv in enumerate(sk.moves):
    if mv['images']:
        continue
    u = {tt: 1 for tt in sk.gens[mv['rep']]}
    for (c, pi, cp) in R.Hc[1:]:
        img = R.tau((c, pi, cp), u)
        (l0, g0) = next(iter(u))
        found = None
        for (l1, g1) in img:
            if l1 != l0: continue
            b = G.mul(G.MATS[g1], np.conj(G.MATS[g0]).T)
            if G.is_monomial(b):
                bu = G.lmul(b, {k: 1 for k in u})
                if set(bu) == set(img):
                    found = b; break
        if found is None:
            print(f"move {t}: tau(u) is NOT a free multiple of u (terms {len(u)} vs {len(img)})")
        else:
            sig, q = G.mono_of(found)
            o = 1; M = found.copy()
            while not np.allclose(M, np.eye(G.n)) and o < 64:
                M = G.mul(M, found); o += 1
            sc = np.allclose(found, found[0, 0] * np.eye(G.n))
            print(f"move {t}: tau(u) = beta u, beta {'SCALAR' if sc else 'NON-scalar'} free letter of order {o}; perm {sig} phases {q}")
tok = []; gd = []
for i, T in enumerate(sk.gens):
    dd = {}
    for j, tt in enumerate(T):
        dd[tt] = len(tok) + 1; tok.append((i, j))
    gd.append(dd)
for l in range(W):
    cl = G.closure(gd, (l, G.KID), 20000)
    rows, cols, ent = cl
    print(f"array {l}: closure {len(rows)} x {len(cols)}")
    rl = {v: k for k, v in rows.items()}; clm = {v: k for k, v in cols.items()}
    if len(rows) <= 40:
        for r in range(len(rows)):
            lt, gt = rl[r]
            es = [(c, tok[co - 1]) for (rr_, c, co) in ent if rr_ == r]
            print(f"  row {r} [{nm(gt)} d{G.dist(gt)}]x{lt}: " + " ; ".join(f"col{c}(g{tk[0]}): {cv[tk[0]][tk[1]].__repr__(names)}" for c, tk in es))
        for c in range(len(cols)):
            kind, i, b = clm[c]
            print(f"  col {c}: generator g{i} times free letter {nm(b)}")
    alv = d['alpha']; cvn = sk.coefs(alv)
    A = np.zeros((len(rows), len(cols)), dtype=complex)
    for (r, c, co) in ent:
        i, j = tok[co - 1]; A[r, c] += cvn[i][j]
    t = np.zeros(len(rows), dtype=complex); t[0] = 1
    beta = np.linalg.lstsq(A, t, rcond=None)[0]
    res = A @ beta - t
    print(f"  residual^2 {np.sum(np.abs(res)**2):.4f}; per row {np.round(np.abs(res)**2, 4)}")
