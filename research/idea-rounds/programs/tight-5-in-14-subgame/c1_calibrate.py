# c1_calibrate.py -- calibration: build V_56, check it is a representation with the right character,
# find the hook Y (the P0-piece of the 28 radical weights), reproduce cap 10, profile (10, 18, 28), ratio 14/5,
# and the cut profile of one plain 3-move chain.
import numpy as np, json, time
from lib514 import *
t0 = time.time()
R = Irreps(verbose=True)
rng = np.random.default_rng(1)
gens = [TRANS[w] for w in range(1, NV)]
# homomorphism check on random words
err = 0
for _ in range(20):
    a = gens[rng.integers(63)]; b = gens[rng.integers(63)]; c = gens[rng.integers(63)]
    g = comp(a, comp(b, c))
    err = max(err, np.abs(R.r56(g) - R.r56(a) @ R.r56(b) @ R.r56(c)).max())
print('homomorphism error', err)
print('orthogonality error', max(np.abs(R.r56(g) @ R.r56(g).T - np.eye(56)).max() for g in gens[:10]))
print('trace at identity', np.trace(R.r56(ID)), ' at a transvection', np.trace(R.r56(TRANS[9])), '(table: 56, -24)')
print('trace of h1 h2 (expect class 2b value -8 or 2c value 8):', round(np.trace(R.r56(comp(h1, h2))), 6),
      ' trace of K:', round(np.trace(R.r56(K)), 6))
# radical weights
A = [R.r56(a) for a in A_GENS]
# simultaneous diagonalisation: generic combination
Mix = sum(rng.normal() * a for a in A)
w, V = np.linalg.eigh((Mix + Mix.T) / 2)
signs = np.array([[np.round(V[:, i] @ a @ V[:, i]) for a in A] for i in range(56)])
for a in A: assert np.abs(V.T @ a @ V - np.diag(np.diag(V.T @ a @ V))).max() < 1e-8
wts = [tuple(int(s) for s in row) for row in signs]
print('distinct radical weights:', len(set(wts)), 'of 56')
# Levi orbits on weights: l acts on V_56; weight vector v -> r56(l) v, whose weight is the transported one
wt_index = {wt: i for i, wt in enumerate(wts)}
Lg = [levi([[1, 1, 0], [0, 1, 0], [0, 0, 1]]), levi([[0, 0, 1], [1, 0, 0], [0, 1, 0]])]
def wt_of(vec): return tuple(int(np.round(vec @ a @ vec / (vec @ vec)))) if False else tuple(int(np.round((vec @ (a @ vec)) / (vec @ vec))) for a in A)
parent = list(range(56))
def find(i):
    while parent[i] != i: parent[i] = parent[parent[i]]; i = parent[i]
    return i
for l in Lg:
    rl = R.r56(l)
    for i in range(56):
        j = wt_index[wt_of(rl @ V[:, i])]
        parent[find(i)] = find(j)
orbits = {}
for i in range(56): orbits.setdefault(find(i), []).append(i)
sizes = sorted(len(o) for o in orbits.values()); print('Levi orbits on the 56 radical weights:', sizes)
from collections import Counter
mult = Counter(wts)
o28 = [i for i in range(56) if mult[wts[i]] == 1]; o7 = [i for i in range(56) if mult[wts[i]] == 4]; o21 = []
print('weight multiplicities:', sorted(Counter(mult.values()).items()), ' |multiplicity-one weights| =', len(o28))
Y = V[:, o28]; Yp = V[:, o7]
def f(Ybas, g): return rank(np.hstack([Ybas, R.r56(g) @ Ybas])) - Ybas.shape[1]
hh2 = comp(h1, h2)
print('hook Y = 28-piece: f(h) =', f(Y, h1), ' f(h1h2) =', f(Y, hh2), ' f(K) =', f(Y, K), '  -> ratio', f(Y, K) / f(Y, h1))
print('complement (7+21, the 7 weights of multiplicity 4): f(h) =', f(Yp, h1), ' f(h1h2) =', f(Yp, hh2), ' f(K) =', f(Yp, K))
# P0 invariance of Y: Levi gens and A gens
for g in Lg + A_GENS:
    assert rank(np.hstack([Y, R.r56(g) @ Y])) == 28
print('Y is kept by the P0 generators')
# plain chain profile: registers x, h x, h2 h x (= h1 h2 x), K x with the chain h1, then h2, then h3 (all commute)
def Zdim(regs):
    """dim of {y in Y : r y in Y for all r in regs} (regs = list of group elements)"""
    # y = Y c; need (I - P_Y) r Y c = 0
    PY = Y @ Y.T
    M = np.vstack([(np.eye(56) - PY) @ R.r56(g) @ Y for g in regs]) if regs else np.zeros((0, 28))
    return 28 - rank(M)
chain = [ID, h1, comp(h2, h1), comp(h3, comp(h2, h1))]
print('plain chain h1,h2,h3: dims of Z after each move:', [Zdim(chain[1:k + 1]) for k in range(4)])
chain2 = [ID, h1, comp(TRANS[EE(2)], h1)]
print('elapsed %.1fs' % (time.time() - t0))
np.savez('v56_basis.npz', B7=R.B7, B27=R.B27, B56=R.B56, V=V, o28=np.array(o28), o7=np.array(o7), o21=np.array(o21),
         wts=np.array(wts))
