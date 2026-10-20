# c5_labels.py -- the hook seen from the singular subgame.
# (1) V_56 restricted to O(q0) = S8: which S8-irreducibles (character inner products with Young characters).
# (2) B'-invariant subspaces of V_56 with labels: Y(L7, L21) = (U7 (x) L7) + (U21 (x) L21), L in P^1 (real and complex
#     scans) or 0 or C^2; ratio f(K)/f(h). Which labels give 2.8, and does any beat it?
import numpy as np, time, itertools
from lib514 import *
from libs8 import *
from libcut import *
t0 = time.time(); rng = np.random.default_rng(11)
R = Irreps()
refl = [TRANS[v] for v in NS]
Gp = gen_closure(refl)
# (1) character of V_56 on O(q0) by cycle type (one representative per S8 class), inner products with Young characters
from math import factorial
from collections import Counter as _C
def ctype(p):
    seen = set(); ls = []
    for i in range(8):
        if i in seen: continue
        L = 0; j = i
        while j not in seen: seen.add(j); j = int(p[j]); L += 1
        ls.append(L)
    return tuple(sorted(ls, reverse=True))
reps = {}; csize = _C()
for g in Gp:
    p = to_perm(g); t = ctype(p); csize[t] += 1
    if t not in reps: reps[t] = (g, p)
print('S8 classes found:', len(reps), ' sizes sum', sum(csize.values()))
tr56 = {t: np.trace(R.r56(g)) for t, (g, p) in reps.items()}
dec = {}
for lam in partitions(8):
    Y = YoungRep(lam)
    m = sum(csize[t] * tr56[t] * np.trace(Y(reps[t][1])) for t in reps) / 40320
    if abs(m) > 1e-6: dec[lam] = round(m, 6)
print('V_56 restricted to S8:', dec, '[%.0fs]' % (time.time() - t0))
# (2) labels
alt = [g for g in gen_closure(A_GENS) if all(q0(int(g[v])) == q0(v) for v in range(NV))]
Bp_gens = [levi([[1, 1, 0], [0, 1, 0], [0, 0, 1]]), levi([[0, 0, 1], [1, 0, 0], [0, 1, 0]])] + [a for a in alt if not np.array_equal(a, ID)][:3]
Bp = gen_closure(Bp_gens); assert len(Bp) == 1344
mats = [R.r56(b) for b in Bp]
comps = isotypic(mats, rng)
gens_m = [R.r56(b) for b in Bp_gens]
def split_pair(comp):
    """for a multiplicity-2 component: two isometric embeddings Q1, Q2X intertwining the B' action"""
    B = comp['basis']; d = comp['d']
    # re-split B into two irreducible pieces with a random commuting operator restricted to the component
    S = rng.standard_normal((B.shape[1], B.shape[1])); S = S + S.T
    Rm = sum((B.T @ M @ B) @ S @ (B.T @ M @ B).T for M in mats) / len(mats)
    w, V = np.linalg.eigh((Rm + Rm.T) / 2)
    Q1 = B @ V[:, :d]; Q2 = B @ V[:, d:]
    eqs = [np.kron(np.eye(d), (Q1.T @ G @ Q1).T) - np.kron(Q2.T @ G @ Q2, np.eye(d)) for G in gens_m]
    _, sv, vt = np.linalg.svd(np.vstack(eqs)); X = vt[-1].reshape(d, d)
    X = X / np.sqrt(np.trace(X @ X.T) / d)
    assert sv[-1] < 1e-8 and np.allclose(X @ X.T, np.eye(d), atol=1e-8)
    Q2X = Q2 @ X
    for G in gens_m: assert np.allclose(Q1.T @ G @ Q1, Q2X.T @ G @ Q2X, atol=1e-8)
    return Q1, Q2X
pairs = {c['d']: split_pair(c) for c in comps}
print("B'-types in V_56:", sorted((c['d'], c['m']) for c in comps))
rh, rK, rh2 = R.r56(h1), R.r56(K), R.r56(comp(h1, h2))
def Ysub(lab):
    cols = []
    for d, L in lab.items():
        Q1, Q2 = pairs[d]
        if L == 'all': cols += [Q1, Q2]
        elif L is None: continue
        else:
            a, b = L; cols.append(a * Q1 + b * Q2)
    return np.hstack(cols) if cols else np.zeros((56, 0))
# which labels give the hook? project the hook Y (P0-piece of the 28 multiplicity-one weights) onto the pairs
dat = np.load('v56_basis.npz'); Vw = dat['V']; Yhook = Vw[:, dat['o28']]
for d, (Q1, Q2) in pairs.items():
    M = np.vstack([Q1.T @ Yhook, Q2.T @ Yhook])     # (2d) x 28
    # the hook's tau-part is U_tau (x) [a:b]: find [a:b] by the rank-d structure
    u, s, vt = np.linalg.svd(np.vstack([np.linalg.norm(Q1.T @ Yhook), np.linalg.norm(Q2.T @ Yhook)]).reshape(1, 2))
    A1 = Q1.T @ Yhook; A2 = Q2.T @ Yhook
    # find (a,b) with a*A2 - b*A1 = 0 (component of hook in tau is image of a Q1 + b Q2)
    _, sv, vt = np.linalg.svd(np.hstack([A2.reshape(-1, 1), -A1.reshape(-1, 1)]))
    print(f'type {d}: hook label (a, b) =', np.round(vt[-1], 6), ' residual', sv[-1])
def ratio(lab):
    Y = Ysub(lab)
    c = f_of(Y, rh)
    return (f_of(Y, rK) / c if c else 0.0, Y.shape[1], c, f_of(Y, rK), f_of(Y, rh2))
# real scan
best = []
angs = np.linspace(0, np.pi, 181)[:-1]
for t7 in angs:
    for t21 in angs:
        r = ratio({7: (np.cos(t7), np.sin(t7)), 21: (np.cos(t21), np.sin(t21))})
        best.append((r[0], r[1], r[2], r[3], r[4], round(t7, 4), round(t21, 4)))
best.sort(key=lambda t: -t[0])
print('real-line scan (180 x 180): top ratios', [(round(b[0], 4), b[1:5], b[5], b[6]) for b in best[:6]])
from collections import Counter
print('   distribution of (cap, demand):', Counter((b[2], b[3]) for b in best).most_common(12))
# mixed: one type all/none, other a line
for d0 in (7, 21):
    d1 = 21 if d0 == 7 else 7
    for L0 in ('all', None):
        rs = []
        for t in angs:
            lab = {d0: L0, d1: (np.cos(t), np.sin(t))}
            rs.append((ratio(lab), round(t, 4)))
        rs.sort(key=lambda x: -x[0][0])
        print(f'type {d0} = {L0}, type {d1} a real line: best', [(round(x[0][0], 4), x[0][1:], x[1]) for x in rs[:3]])
# complex scan (random)
cbest = []
for _ in range(3000):
    z7 = rng.standard_normal(2) + 1j * rng.standard_normal(2); z21 = rng.standard_normal(2) + 1j * rng.standard_normal(2)
    Y = np.hstack([z7[0] * pairs[7][0] + z7[1] * pairs[7][1], z21[0] * pairs[21][0] + z21[1] * pairs[21][1]])
    c = rank(np.hstack([Y, rh @ Y])) - 28; dm = rank(np.hstack([Y, rK @ Y])) - 28
    cbest.append((dm / c, c, dm))
cbest.sort(reverse=True)
print('complex random lines: top', cbest[:3], ' (cap, demand) counts', Counter((c[1], c[2]) for c in cbest).most_common(5))
print('[%.0fs]' % (time.time() - t0))
