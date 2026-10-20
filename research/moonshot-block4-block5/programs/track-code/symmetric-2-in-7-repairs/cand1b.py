import sys
from relax_search import *
m, W, H, shape, menu, seed, s = 4, 2, 'Z2swap', 'P:inp,P:prev,L:prev,P:prev'.split(','), 2, 66, 1
sk, cost, rec = record_recipe(m, W, H, shape, menu, seed * 100 + s)
G = sk.G
fidx = {}
for i, T in enumerate(sk.gens):
    for (l, g) in T:
        key = (G.frm(g), l)
        if key not in fidx: fidx[key] = len(fidx)
Kkey = (G.frm(G.KID), 0); nrow = len(fidx)
rowmap = [np.array([fidx[(G.frm(g), l)] for (l, g) in T]) for T in sk.gens]
N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
def full(x):
    al = np.ones(N, dtype=complex); al[free_idx] = x[:nf] + 1j * x[nf:]; return al
def rr(al):
    cv = sk.coefs(al)
    A = np.zeros((nrow, len(sk.gens)), dtype=complex)
    for i in range(len(sk.gens)): np.add.at(A[:, i], rowmap[i], cv[i])
    t = np.zeros(nrow, dtype=complex); t[fidx[Kkey]] = 1
    beta = np.linalg.lstsq(A, t, rcond=1e-12)[0]
    return A @ beta - t, beta
print("unknowns", nf, "rows", nrow, "generators", len(sk.gens))
rng = random.Random(3)
x = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
for rnd in range(6):
    x, c = lm(lambda y: rr(full(y))[0], x, iters=40)
    r, beta = rr(full(x))
    print(f"round {rnd}: res^2 {c:.3e}  |alpha| max {np.abs(full(x)).max():.3e}  |beta| max {np.abs(beta).max():.3e}", flush=True)
print("alpha:", np.round(full(x), 4))
