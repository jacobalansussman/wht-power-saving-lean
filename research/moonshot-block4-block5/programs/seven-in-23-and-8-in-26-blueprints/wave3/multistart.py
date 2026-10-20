"""multistart.py -- w3-blueprint: can a MULTI-term shared start serve the nontrivial characters? (NOTE.md 5b loophole)

Route A through a shared start F = delta . sum_h psi(h) tau_h(sum_i alpha_i b_i x_(l_i)) puts on the input-distance-4
cells of O3's representative, in the chi-component, the vector  omega . E_phi . beta_chi  where
  omega   = sum over paths of (path weight) [A delta]       (a vector on the cells [K r], r in B: H_max\\B, 1024 cells)
  E_phi   = sum_h phi(h) c_h,   phi = psi chi-bar           (right action of the symmetry's c_h)
  beta_chi = sum_i alpha_i chi-bar(l_i) b_i                 (right action of the content letters)
In F O O O F nothing else has input-distance-4 cells in a nontrivial character, so serving chi needs
omega . E_phi . beta_chi = c_chi [K] with c_chi != 0, for ONE omega common to all served chi.
This script builds the 1024 cells [K r] (engine cosets, keyed by bm's frame key, exact for Clifford elements), the right
action of c_h and of random content letters b_i, and computes, for random contents (r terms on random arrays,
random alpha, psi = 1 or not), the largest set of nontrivial characters that ONE omega (ANY omega: no path constraint,
a relaxation of the relaxation) can serve.
Usage: multistart.py W H ntrials seed"""
import sys, time
from symsolve import *
W, H, ntr, seed = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
R = Run(4, W, H, seed); G = R.G; rng = random.Random(seed)
t0 = time.time()
# ---- generators of the engine free group B (right multiplication)
def addr(cols):
    sig = []
    for x in range(16):
        y = 0
        for j in range(4):
            if (x >> j) & 1: y ^= cols[j]
        sig.append(y)
    return G.mono(sig, [0] * 16)
gensB = []
for i in range(4):
    for j in range(4):
        if i != j:
            cols = [1 << k for k in range(4)]; cols[i] = (1 << i) | (1 << j)
            gensB.append(addr(cols))
for j in range(4):
    q = [((x >> j) & 1) for x in range(16)]
    gensB.append(G.mono(list(range(16)), q))                  # i^{x_j}
    gensB.append(G.Xm(1 << j))                                 # shifts
for i in range(4):
    for j in range(i + 1, 4):
        q = [2 * (((x >> i) & 1) * ((x >> j) & 1)) for x in range(16)]
        gensB.append(G.mono(list(range(16)), q))              # CZ
# ---- BFS over cells [K r]
def key(M):
    return G.frame_key(M)
reps = [G.I]; keys = {key(G.K): 0}
frontier = [G.I]
while frontier:
    nf = []
    for r in frontier:
        for g in gensB:
            r2 = G.mul(r, g)
            k = key(G.mul(G.K, r2))
            if k not in keys:
                keys[k] = len(reps); reps.append(r2); nf.append(r2)
    frontier = nf
N = len(reps)
print(f"cells [K r]: {N} ({time.time()-t0:.0f}s)", flush=True)
def ract(b):
    """permutation p with [K r_t] . b = [K r_p[t]]"""
    p = np.zeros(N, dtype=np.int64)
    for t, r in enumerate(reps):
        k = key(G.mul(G.K, r, b))
        p[t] = keys[k]
    return p
def pmat(p, coef=1.0):
    M = np.zeros((N, N), dtype=complex)
    M[np.arange(N), p] = coef
    return M
# symmetry elements c_h and characters (array label of h = pi_h(0))
lab = [h[1][0] for h in R.Hc]
Pc = [ract(h[0]) for h in R.Hc]
if H == 'Z7rot':
    chars = [[np.exp(2j * np.pi * r * l / 7) for l in range(7)] for r in range(7)]
else:
    chars = [[(-1) ** bin(r & l).count('1') for l in range(8)] for r in range(8)]
K0 = 0
# right <c_h>-orbits on cells, for the isotypic bases
orbit_reps = []
seen = set()
for t in range(N):
    if t in seen:
        continue
    orb = {int(p[t]) for p in Pc}
    seen |= orb
    orbit_reps.append(t)
def served(terms, psi_idx):
    """Per character (decoupled, NOTE 5c): chi is servable iff [K] lies in Iso_phi . beta_chi, phi = psi chi-bar,
    Iso_phi spanned by mu_t = sum_h phi(h) e_(t . c_h).  Returns (#servable nontrivial chars, list)."""
    psi = chars[psi_idx]
    Pb = [(a, ract(b), l) for (a, b, l) in terms]
    ok = []
    for r in range(1, W):
        chi = chars[r]
        phi = [psi[lab[hj]] * np.conj(chi[lab[hj]]) for hj in range(len(Pc))]
        rows = []
        for t in orbit_reps:
            mu = np.zeros(N, dtype=complex)
            for hj, p in enumerate(Pc):
                mu[p[t]] += phi[hj]
            if np.abs(mu).max() < 1e-12:
                continue
            out = np.zeros(N, dtype=complex)
            for (a, p, l) in Pb:
                np.add.at(out, p, a * np.conj(chi[l]) * mu)      # (mu . b)[t . b] = mu[t]
            rows.append(out)
        if not rows:
            continue
        M = np.array(rows)                       # rows span Iso_phi . beta_chi
        e = np.zeros(N, dtype=complex); e[K0] = 1
        x = np.linalg.lstsq(M.T, e, rcond=None)[0]
        res = np.linalg.norm(M.T @ x - e)
        if res < 1e-8:
            ok.append(r)
    return len(ok), ok
print("sanity: single term b = 1, psi = trivial ->", served([(1.0, G.I, 0)], 0), flush=True)
print("sanity: single term b = 1, psi = char 1 ->", served([(1.0, G.I, 0)], 1), " (Theorem R: exactly [1])", flush=True)
results = []
for tr in range(ntr):
    nterm = rng.choice([1, 2, 2, 3, 3, 4])
    terms = []
    for i in range(nterm):
        b = rand_free_b(G, rng, 4)
        l = rng.randrange(W) if i else 0
        a = complex(rng.gauss(0, 1), rng.gauss(0, 1))
        terms.append((a, b, l))
    psi_idx = rng.choice([0, 0, 1])
    t1 = time.time()
    k, which = served(terms, psi_idx)
    results.append((k, nterm, psi_idx))
    print(f"trial {tr}: {nterm} terms (arrays {[t[2] for t in terms]}), psi {psi_idx}: nontrivial characters servable by one omega: {k} {which} ({time.time()-t1:.0f}s)", flush=True)
from collections import Counter
print("histogram (max servable, nterms):", sorted(Counter((k, n) for (k, n, p) in results).items()), flush=True)
