"""chidiag.py -- w3-blueprint: character-resolved frame-sum relaxation of a bp.py skeleton.
For chi = omega^r: project every generator onto chi (cell orbits under tau_*), then list each free-orbit
representative's cell orbits by distance, and test e_K in the span (random coefficients)."""
import sys
from bp import *
W, H, variant, seed = int(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
sk, cost = build(W, H, seed, variant)
G, R = sk.G, sk.R
# cell orbit id: frame key of (g c'^-k) for array k -> canonical: frame of the array-0 representative
# tau_h maps (l, g) -> (pi_h l, c g c'); frame of c g c' = frame(g c'). So cell (F, l) with l = pi_h(0) has rep frame F c'_h^-1.
pis = [h[1] for h in R.Hc]
inv = {}
for hj, h in enumerate(R.Hc):
    inv[h[1][0]] = hj          # array l = pi_h(0) -> element h
def cell_rep(l, g):
    h = R.Hc[inv[l]]
    cpi = np.conj(h[2]).T
    return G.frm(G.intern(G.mul(G.MATS[g], cpi))), hj_of(l)
def hj_of(l):
    return inv[l]
# characters: for Z7 h_j = tau^j (closure order may differ) -> use the array index as the group label for regular H
labels = list(range(W))
if H == 'Z7rot':
    chars = [(r, [np.exp(2j * np.pi * r * l / 7) for l in range(7)]) for r in range(7)]
else:
    chars = [(r, [(-1) ** bin(r & l).count('1') for l in range(8)]) for r in range(8)]
rng = np.random.default_rng(1)
alpha = rng.normal(size=sk.nalpha) + 1j * rng.normal(size=sk.nalpha)
for a in sk.fixed:
    alpha[a] = 1
cv = sk.coefs(alpha)
names = ['x%d' % l for l in range(W)]
for mi, mv in enumerate(sk.moves):
    names.append('M%d' % mi)
    names += ['M%d.%d' % (mi, k + 1) for k in range(len(mv['images']))]
reps = list(range(W)) + [mv['rep'] for mv in sk.moves]
KF = G.frm(G.KID)
for (r, chi) in chars[:3]:
    cols = []; allcells = {}
    vecs = []
    for i in reps:
        v = {}
        for (j, (l, g)) in enumerate(sk.gens[i]):
            f, _ = cell_rep(l, g)
            v[f] = v.get(f, 0) + cv[i][j] * np.conj(chi[l])
        v = {f: c for f, c in v.items() if abs(c) > 1e-9}
        vecs.append(v)
        for f in v:
            allcells.setdefault(f, len(allcells))
    if KF not in allcells:
        allcells[KF] = len(allcells)
    A = np.zeros((len(allcells), len(vecs)), dtype=complex)
    for c_, v in enumerate(vecs):
        for f, c in v.items():
            A[allcells[f], c_] = c
    t = np.zeros(len(allcells), dtype=complex); t[allcells[KF]] = 1
    beta = np.linalg.lstsq(A, t, rcond=1e-10)[0]
    res = np.linalg.norm(A @ beta - t) ** 2
    # distances of cells per vector
    dist_of = {}
    for i in reps:
        for (l, g) in sk.gens[i]:
            f, _ = cell_rep(l, g); dist_of[f] = G.dist(g)
    print(f"chi r={r}: residual^2 at random alpha {res:.3f}; vectors (name: #cells by distance 0..4):")
    for nm, i, v in zip([names[i] for i in reps], reps, vecs):
        hist = [sum(1 for f in v if dist_of.get(f) == d) for d in range(5)]
        print(f"   {nm}: {hist}  K-cell coef {abs(v.get(KF, 0)):.3f}")
    # cells shared by >= 2 vectors
    share = {}
    for c_, v in enumerate(vecs):
        for f in v:
            share.setdefault(f, []).append(c_)
    nsh = sum(1 for f, l in share.items() if len(l) >= 2)
    print(f"   cells shared by >= 2 vectors: {nsh} of {len(share)}")
